#!/usr/bin/env python3
#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Timed monolithic GTEP driver.

This script builds and solves a single full monolithic/extensive-form GTEP model
using the data/model/solver sections from a Progressive Hedging YAML config.

It reports time spent in major phases:

* config loading,
* data loading,
* cost-data loading,
* model construction,
* GDP transformation,
* solve,
* solution export.

It streams solver output and writes run artifacts under an output directory.

Example
-------

From the repository root:

    python examples/monolithic/run_monolithic_gtep_timed.py \\
        --config examples/progressive_hedging/gtep_ph_123bus_config.yaml \\
        --output-dir /scratch/jkskolf/gtep_monolithic_runs/123bus_4rep \\
        --log-level INFO

Notes
-----
This script does not run PH. It builds one full model containing all
representative periods specified in the config.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
import traceback
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pyomo.environ as pyo

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_solution import ExpansionPlanningSolution

from gtep.algorithms.progressive_hedging.config import load_ph_config
from gtep.algorithms.progressive_hedging.scenario_data import (
    build_cost_data,
    build_full_data,
)
from gtep.algorithms.progressive_hedging.solver import (
    create_solver,
    get_active_objective_value,
    get_result_bound,
    compute_relative_gap,
)

logger = logging.getLogger("gtep.examples.monolithic.run_monolithic_gtep_timed")


@dataclass
class PhaseTiming:
    """Timing record for one driver phase."""

    name: str
    start_time_sec: float
    end_time_sec: float
    elapsed_sec: float


@dataclass
class MonolithicRunSummary:
    """JSON-serializable monolithic run summary."""

    config_path: str
    output_dir: str
    success: bool
    timings: list[dict[str, Any]]
    solver: dict[str, Any]
    objective: float | None
    gtep_results_dir: str | None
    model_pprint_file: str | None
    large_coefficients_file: str | None
    error: str | None = None


class Timer:
    """Simple phase timer."""

    def __init__(self):
        self.records: list[PhaseTiming] = []

    @contextmanager
    def phase(self, name: str):
        logger.info("START phase: %s", name)
        start = time.perf_counter()
        try:
            yield
        finally:
            end = time.perf_counter()
            elapsed = end - start
            self.records.append(
                PhaseTiming(
                    name=name,
                    start_time_sec=start,
                    end_time_sec=end,
                    elapsed_sec=elapsed,
                )
            )
            logger.info("END phase: %s elapsed=%.3f sec", name, elapsed)

    def to_jsonable(self) -> list[dict[str, Any]]:
        return [asdict(record) for record in self.records]


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    _configure_logging(args.log_level)

    try:
        run_monolithic(args)
    except Exception:
        logger.exception("Monolithic GTEP run failed.")
        return 1

    return 0


def run_monolithic(args: argparse.Namespace) -> MonolithicRunSummary:
    """Build, solve, and save a monolithic GTEP model."""
    timer = Timer()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    gtep_results_dir_written = None

    summary_path = output_dir / "run_summary.json"
    solution_json_path = output_dir / "gtep_solution.json"
    solver_log_path = output_dir / "solver.log"
    model_pprint_path = output_dir / "model.pprint.txt"
    large_coefficients_path = output_dir / "large_coefficients.json"

    cfg = None
    model_object = None
    solve_results = None
    objective_value = None
    solution_json_written = None
    model_pprint_written = None
    large_coefficients_written = None

    try:
        with timer.phase("load_config"):
            cfg = load_ph_config(args.config)

        with timer.phase("create_gtep_results_directory"):
            if args.skip_gtep_results:
                sol_object = None
                gtep_results_dir = None
            else:
                sol_object = ExpansionPlanningSolution(str(cfg.data.data_path))
                gtep_results_dir = Path(args.gtep_results_dir).resolve()
                dir_name = sol_object.create_results_directory(str(gtep_results_dir))
                gtep_results_dir_written = str(dir_name)

        with timer.phase("prepare_solver_environment"):
            _apply_solver_environment_from_config(cfg.solver)

        with timer.phase("load_data"):
            data_object = build_full_data(cfg)

        with timer.phase("load_cost_data"):
            cost_data = build_cost_data(cfg)

        with timer.phase("instantiate_model_object"):
            model_object = ExpansionPlanningModel(
                data=data_object,
                cost_data=cost_data,
            )

            for key, value in cfg.model.as_model_options().items():
                model_object.config[key] = value

        if sol_object is not None:
            with timer.phase("write_gtep_model_config_csv"):
                sol_object.save_model_config_to_csv(model_object, dir_name)

        with timer.phase("construct_model"):
            model_object.create_model()

        with timer.phase("gdp_transformation"):
            transformation = pyo.TransformationFactory(cfg.model.gdp_transformation)
            if transformation is None:
                raise RuntimeError(
                    f"Could not create GDP transformation {cfg.model.gdp_transformation!r}."
                )
            transformation.apply_to(model_object.model)

        if args.write_model_pprint:
            with timer.phase("write_model_pprint"):
                model_object.report_model(outfile=str(model_pprint_path))
                model_pprint_written = str(model_pprint_path)

        if args.report_large_coefficients:
            with timer.phase("report_large_coefficients"):
                model_object.report_large_coefficients(
                    outfile=str(large_coefficients_path),
                    magnitude_cutoff=args.large_coefficient_cutoff,
                )
                large_coefficients_written = str(large_coefficients_path)

        with timer.phase("solve"):
            opt = create_solver(cfg.solver)
            solve_kwargs: dict[str, Any] = {
                "tee": cfg.solver.tee,
                "logfile": str(solver_log_path),
            }

            solve_results = opt.solve(model_object.model, **solve_kwargs)
            model_object.results = solve_results

        with timer.phase("extract_objective"):
            objective_value = get_active_objective_value(model_object.model)


        if sol_object is not None:
            with timer.phase("write_gtep_result_json_files"):
                sol_object.save_results_in_json_files(
                    model_object,
                    dir_name,
                    value_threshold=args.value_threshold,
                )

        solver_summary = _solver_summary(
            cfg.solver.name,
            solve_results,
            model_object.model,
            solver_log_path,
        )

        summary = MonolithicRunSummary(
            config_path=str(Path(args.config).resolve()),
            output_dir=str(output_dir),
            success=True,
            timings=timer.to_jsonable(),
            solver=solver_summary,
            objective=objective_value,
            gtep_results_dir=gtep_results_dir_written,
            model_pprint_file=model_pprint_written,
            large_coefficients_file=large_coefficients_written,
            error=None,
        )

        _write_summary(summary_path, summary)
        _print_summary(summary)
        return summary

    except Exception as err:
        error_text = "".join(
            traceback.format_exception(type(err), err, err.__traceback__)
        )

        solver_summary = {}
        if cfg is not None and solve_results is not None and model_object is not None:
            solver_summary = _solver_summary(
                cfg.solver.name,
                solve_results,
                model_object.model,
                solver_log_path,
            )

        summary = MonolithicRunSummary(
            config_path=str(Path(args.config).resolve()),
            output_dir=str(output_dir),
            success=False,
            timings=timer.to_jsonable(),
            solver=solver_summary,
            objective=objective_value,
            gtep_results_dir=gtep_results_dir_written,
            model_pprint_file=model_pprint_written,
            large_coefficients_file=large_coefficients_written,
            error=error_text,
        )

        _write_summary(summary_path, summary)
        raise


def _apply_solver_environment_from_config(solver_config) -> None:
    """Apply solver environment variables.

    This duplicates the relevant behavior from the PH solver path so the
    monolithic driver sees the same XPRESS license and temp settings.
    """
    for key, value in solver_config.environment.items():
        os.environ[str(key)] = str(value)
        logger.info("Set solver environment variable %s from config.", key)

    if solver_config.license_file is not None:
        license_path = Path(solver_config.license_file).expanduser()

        if not license_path.exists():
            raise FileNotFoundError(
                f"Configured solver.license_file does not exist: {license_path}"
            )

        env_var = solver_config.license_env_var
        if env_var is None:
            if solver_config.name.lower() in {
                "xpress",
                "xpress_direct",
                "xpress_persistent",
            }:
                env_var = "XPAUTH_PATH"
            else:
                raise ValueError(
                    "solver.license_file was provided but no license_env_var was set "
                    f"and no default is known for solver {solver_config.name!r}."
                )

        os.environ[env_var] = str(license_path)
        logger.info("Set %s from solver.license_file.", env_var)


def _solver_summary(
    solver_name: str,
    results,
    model,
    solver_log_path: Path,
) -> dict[str, Any]:
    """Create solver summary dictionary."""
    status = None
    termination_condition = None
    message = None

    if results is not None:
        status = str(results.solver.status)
        termination_condition = str(results.solver.termination_condition)
        message = str(getattr(results.solver, "message", None))

    lower_bound = get_result_bound(results, "lower_bound") if results is not None else None
    upper_bound = get_result_bound(results, "upper_bound") if results is not None else None

    return {
        "name": solver_name,
        "status": status,
        "termination_condition": termination_condition,
        "message": message,
        "objective": get_active_objective_value(model),
        "lower_bound": lower_bound,
        "upper_bound": upper_bound,
        "gap": compute_relative_gap(lower_bound, upper_bound),
        "log_file": str(solver_log_path),
    }


def _write_summary(path: Path, summary: MonolithicRunSummary) -> None:
    """Write run summary JSON."""
    with path.open("w", encoding="utf-8") as f:
        json.dump(asdict(summary), f, indent=2)


def _print_summary(summary: MonolithicRunSummary) -> None:
    """Print concise run summary to stdout."""
    print()
    print("=" * 80)
    print("Monolithic GTEP Run Summary")
    print("=" * 80)
    print(f"Success:      {summary.success}")
    print(f"Output dir:   {summary.output_dir}")
    print(f"Objective:    {summary.objective}")
    print(f"GTEP results: {summary.gtep_results_dir}")
    print(f"Solver log:   {summary.solver.get('log_file')}")
    print()
    print("Timings:")
    for record in summary.timings:
        print(f"  {record['name']:35s} {record['elapsed_sec']:12.3f} sec")
    print("=" * 80)
    print()


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and solve a timed monolithic GTEP model."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to PH-style YAML config whose data/model/solver sections are reused.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where monolithic run artifacts will be written.",
    )

    parser.add_argument(
        "--write-model-pprint",
        action="store_true",
        help="Write a full Pyomo model pprint file. Can be very large.",
    )

    parser.add_argument(
        "--report-large-coefficients",
        action="store_true",
        help="Write large coefficient report using ExpansionPlanningModel.report_large_coefficients().",
    )

    parser.add_argument(
        "--large-coefficient-cutoff",
        type=float,
        default=1.0e5,
        help="Magnitude cutoff for large coefficient report.",
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Python logging level.",
    )

    parser.add_argument(
        "--gtep-results-dir",
        default="naerm_one_off_tests",
        help=(
            "Local directory where GTEP-format result JSON files will be written "
            "using ExpansionPlanningSolution. Default: naerm_one_off_tests"
        ),
    )

    parser.add_argument(
        "--value-threshold",
        type=float,
        default=1.0e-3,
        help=(
            "Minimum variable value to save in GTEP result JSON files. "
            "Default: 1e-3"
        ),
    )

    parser.add_argument(
        "--skip-gtep-results",
        action="store_true",
        help="Skip writing GTEP-format result JSON files.",
    )

    return parser.parse_args(argv)


def _configure_logging(log_level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


if __name__ == "__main__":
    sys.exit(main())