#!/usr/bin/env python3
"""Profiled monolithic GTEP driver.

This driver is intended to diagnose why monolithic GTEP model construction is
slow. It separates:

* data loading,
* cost-data loading,
* model construction,
* GDP transformation,
* solve,
* solution extraction,
* output writing,

and can run cProfile around the expensive phases.

It also enables Pyomo construction timing logs so you can see which Pyomo
components take the longest to construct.

Example:

    python -u examples/monolithic/run_monolithic_gtep_profiled.py \\
      --config examples/progressive_hedging/gtep_ph_2030_pcm_config.yaml \\
      --output-root naerm_one_off_tests \\
      --enable-cprofile \\
      --pyomo-construction-log-level INFO \\
      --skip-solve
"""

from __future__ import annotations

import argparse
import cProfile
import csv
import io
import json
import logging
import math
import os
import pstats
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pyomo.environ as pyo
from pyomo.common.timing import report_timing
from pyomo.core.base.componentuid import ComponentUID

from gtep.gtep_model import ExpansionPlanningModel
from gtep.algorithms.progressive_hedging.config import load_ph_config
from gtep.algorithms.progressive_hedging.scenario_data import (
    build_cost_data,
    build_full_data,
)
from gtep.algorithms.progressive_hedging.solver import (
    compute_relative_gap,
    create_solver,
    get_active_objective_value,
    get_result_bound,
)

logger = logging.getLogger("gtep.examples.monolithic.run_monolithic_gtep_profiled")


@dataclass
class PhaseRecord:
    name: str
    start_perf_counter: float
    end_perf_counter: float
    elapsed_sec: float


class PhaseTimer:
    def __init__(self):
        self.records: list[PhaseRecord] = []

    @contextmanager
    def phase(self, name: str):
        logger.info("START PHASE %s", name)
        start = time.perf_counter()
        try:
            yield
        finally:
            end = time.perf_counter()
            elapsed = end - start
            self.records.append(
                PhaseRecord(
                    name=name,
                    start_perf_counter=start,
                    end_perf_counter=end,
                    elapsed_sec=elapsed,
                )
            )
            logger.info("END PHASE %s elapsed=%.3f sec", name, elapsed)

    def to_jsonable(self):
        return [asdict(r) for r in self.records]


class MaybeProfiler:
    """Small wrapper around cProfile."""

    def __init__(self, enabled: bool):
        self.enabled = enabled
        self.profiler = cProfile.Profile() if enabled else None

    def enable(self):
        if self.profiler is not None:
            self.profiler.enable()

    def disable(self):
        if self.profiler is not None:
            self.profiler.disable()

    @contextmanager
    def active(self):
        self.enable()
        try:
            yield
        finally:
            self.disable()

    def dump(self, output_dir: Path, limit: int):
        if self.profiler is None:
            return

        raw_path = output_dir / "cprofile_raw.prof"
        self.profiler.dump_stats(str(raw_path))

        for sort_key, filename in [
            ("cumulative", "cprofile_stats_by_cumulative.txt"),
            ("time", "cprofile_stats_by_time.txt"),
        ]:
            stream = io.StringIO()
            stats = pstats.Stats(self.profiler, stream=stream)
            stats.strip_dirs().sort_stats(sort_key).print_stats(limit)

            with (output_dir / filename).open("w", encoding="utf-8") as f:
                f.write(stream.getvalue())


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    run_dir = make_output_dir(args)
    configure_logging(args.log_level, run_dir / "driver.log")

    logger.info("Writing profiling artifacts to %s", run_dir)

    timer = PhaseTimer()
    profiler = MaybeProfiler(args.enable_cprofile)

    try:
        run_profiled(args, run_dir, timer, profiler)
    except Exception:
        logger.exception("Profiled monolithic run failed.")
        profiler.dump(run_dir, args.cprofile_limit)
        write_phase_timing(run_dir, timer)
        return 1

    profiler.dump(run_dir, args.cprofile_limit)
    write_phase_timing(run_dir, timer)
    return 0


def run_profiled(
    args: argparse.Namespace,
    run_dir: Path,
    timer: PhaseTimer,
    profiler: MaybeProfiler,
) -> None:
    cfg = None
    model_object = None
    solve_results = None

    with timer.phase("load_config"):
        cfg = load_ph_config(args.config)

    with timer.phase("prepare_solver_environment"):
        apply_solver_environment(cfg.solver)

    if args.pyomo_construction_log_level is not None:
        with timer.phase("enable_pyomo_construction_timing"):
            level = getattr(logging, args.pyomo_construction_log_level)
            # report_timing installs Pyomo construction timers/logging.
            report_timing(stream=True)
            logging.getLogger("pyomo.common.timing").setLevel(level)
            logging.getLogger("pyomo.core").setLevel(max(logging.WARNING, level))

    with profiler.active():
        with timer.phase("load_data"):
            data_object = build_full_data(cfg)

        with timer.phase("load_cost_data"):
            cost_data = build_cost_data(cfg)

        with timer.phase("construct_model_total"):
            model_object = ExpansionPlanningModel(
                data=data_object,
                cost_data=cost_data,
            )

            for key, value in cfg.model.as_model_options().items():
                model_object.config[key] = value

            with timer.phase("construct_model_create_model"):
                model_object.create_model()

        # with timer.phase("write_model_size_pre_gdp"):
        #     write_model_size_artifacts(
        #         model_object.model,
        #         run_dir,
        #         prefix="pre_gdp",
        #     )

        if not args.skip_gdp:
            with timer.phase("gdp_transformation"):
                transformation = pyo.TransformationFactory(cfg.model.gdp_transformation)
                if transformation is None:
                    raise RuntimeError(
                        f"Could not create GDP transformation "
                        f"{cfg.model.gdp_transformation!r}."
                    )
                transformation.apply_to(model_object.model)

            # with timer.phase("write_model_size_post_gdp"):
            #     write_model_size_artifacts(
            #         model_object.model,
            #         run_dir,
            #         prefix="post_gdp",
            #     )

    if args.skip_solve:
        logger.info("Skipping solve because --skip-solve was provided.")
        write_run_summary(
            run_dir,
            cfg=cfg,
            model=model_object.model if model_object is not None else None,
            results=None,
            objective=None,
            skipped_solve=True,
        )
        return

    with timer.phase("solve"):
        opt = create_solver(cfg.solver)
        solve_kwargs = {
            "tee": cfg.solver.tee,
            "logfile": str(run_dir / "solver.log"),
        }
        solve_results = opt.solve(model_object.model, **solve_kwargs)
        model_object.results = solve_results

    with timer.phase("extract_objective"):
        objective = get_active_objective_value(model_object.model)

    if args.write_basic_solution:
        with timer.phase("extract_basic_solution"):
            basic_solution = extract_basic_solution(
                model_object.model,
                nonzero_tolerance=args.nonzero_tolerance,
                max_records=args.max_variable_records,
            )

        with timer.phase("write_basic_solution"):
            write_json(run_dir / "basic_solution.json", basic_solution)

    with timer.phase("write_run_summary"):
        write_run_summary(
            run_dir,
            cfg=cfg,
            model=model_object.model,
            results=solve_results,
            objective=objective,
            skipped_solve=False,
        )


def write_model_size_artifacts(model: pyo.ConcreteModel, run_dir: Path, prefix: str):
    """Write model size and component count artifacts."""
    size = {
        "variables_total": 0,
        "variables_binary": 0,
        "variables_integer": 0,
        "variables_continuous": 0,
        "constraints_total": 0,
        "constraints_active": 0,
        "objectives_total": 0,
        "objectives_active": 0,
        "expressions_total": 0,
        "blocks_total": 0,
        "disjuncts_total": 0,
        "disjunctions_total": 0,
    }

    for var in model.component_data_objects(pyo.Var, descend_into=True):
        size["variables_total"] += 1
        if var.is_binary():
            size["variables_binary"] += 1
        elif var.is_integer():
            size["variables_integer"] += 1
        else:
            size["variables_continuous"] += 1

    for con in model.component_data_objects(pyo.Constraint, descend_into=True):
        size["constraints_total"] += 1
        if con.active:
            size["constraints_active"] += 1

    for obj in model.component_data_objects(pyo.Objective, descend_into=True):
        size["objectives_total"] += 1
        if obj.active:
            size["objectives_active"] += 1

    for _ in model.component_data_objects(pyo.Expression, descend_into=True):
        size["expressions_total"] += 1

    for _ in model.component_data_objects(pyo.Block, descend_into=True):
        size["blocks_total"] += 1

    # GDP components may not always be imported through pyo.
    try:
        from pyomo.gdp import Disjunct, Disjunction

        for _ in model.component_data_objects(Disjunct, descend_into=True):
            size["disjuncts_total"] += 1

        for _ in model.component_data_objects(Disjunction, descend_into=True):
            size["disjunctions_total"] += 1
    except Exception:
        pass

    write_json(run_dir / f"{prefix}_model_size.json", size)

    component_counts = component_count_table(model)
    write_json(run_dir / f"{prefix}_component_counts.json", component_counts)

    with (run_dir / f"{prefix}_component_counts.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "component_type",
                "component_name",
                "active",
                "data_count",
            ],
        )
        writer.writeheader()
        writer.writerows(component_counts)


def component_count_table(model: pyo.ConcreteModel) -> list[dict[str, Any]]:
    """Return counts by indexed component."""
    rows = []

    component_types = [
        pyo.Var,
        pyo.Constraint,
        pyo.Expression,
        pyo.Objective,
        pyo.Param,
        pyo.Set,
        pyo.Block,
    ]

    try:
        from pyomo.gdp import Disjunct, Disjunction

        component_types.extend([Disjunct, Disjunction])
    except Exception:
        pass

    for ctype in component_types:
        for comp in model.component_objects(ctype, descend_into=True):
            try:
                data_count = len(comp)
            except TypeError:
                data_count = 1

            rows.append(
                {
                    "component_type": ctype.__name__,
                    "component_name": comp.name,
                    "active": bool(getattr(comp, "active", True)),
                    "data_count": data_count,
                }
            )

    rows.sort(
        key=lambda r: (r["component_type"], -r["data_count"], r["component_name"])
    )
    return rows


def extract_basic_solution(
    model: pyo.ConcreteModel,
    *,
    nonzero_tolerance: float,
    max_records: int | None,
) -> dict[str, Any]:
    records = []

    for var in model.component_data_objects(pyo.Var, active=True, descend_into=True):
        value = pyo.value(var, exception=False)
        if value is None:
            continue

        value = float(value)
        if abs(value) <= nonzero_tolerance:
            continue

        records.append(
            {
                "name": var.name,
                "component_uid": str(ComponentUID(var)),
                "value": value,
                "is_binary": bool(var.is_binary()),
                "is_integer": bool(var.is_integer()),
                "fixed": bool(var.fixed),
            }
        )

        if max_records is not None and len(records) >= max_records:
            break

    return {
        "nonzero_tolerance": nonzero_tolerance,
        "max_records": max_records,
        "nonzero_variables": records,
    }


def write_run_summary(
    run_dir: Path,
    *,
    cfg,
    model,
    results,
    objective,
    skipped_solve: bool,
):
    solver_summary = {}

    if results is not None and model is not None:
        lower_bound = get_result_bound(results, "lower_bound")
        upper_bound = get_result_bound(results, "upper_bound")
        solver_summary = {
            "name": cfg.solver.name,
            "status": str(results.solver.status),
            "termination_condition": str(results.solver.termination_condition),
            "message": str(getattr(results.solver, "message", None)),
            "objective": objective,
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "gap": compute_relative_gap(lower_bound, upper_bound),
            "log_file": str(run_dir / "solver.log"),
        }

    summary = {
        "config_path": (
            str(Path(cfg.config_path).resolve())
            if cfg.config_path is not None
            else None
        ),
        "run_dir": str(run_dir),
        "skipped_solve": skipped_solve,
        "solver": solver_summary,
        "objective": objective,
    }

    write_json(run_dir / "run_summary.json", summary)


def write_phase_timing(run_dir: Path, timer: PhaseTimer):
    data = timer.to_jsonable()
    write_json(run_dir / "phase_timing.json", {"timings": data})

    with (run_dir / "phase_timing.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "name",
                "start_perf_counter",
                "end_perf_counter",
                "elapsed_sec",
            ],
        )
        writer.writeheader()
        writer.writerows(data)


def apply_solver_environment(solver_config):
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
                    "solver.license_file was provided but no license_env_var was set."
                )

        os.environ[env_var] = str(license_path)
        logger.info("Set %s from solver.license_file.", env_var)


def make_output_dir(args: argparse.Namespace) -> Path:
    if args.output_dir is not None:
        output_root = Path(args.output_dir)
    else:
        output_root = Path("naerm_one_off_tests")

    if args.no_timestamp:
        run_dir = output_root
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = output_root / f"profiled_{timestamp}"

    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir.resolve()


def write_json(path: Path, data: Any):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(json_sanitize(data), f, indent=2, allow_nan=False)


def json_sanitize(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_sanitize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_sanitize(v) for v in value]
    if isinstance(value, tuple):
        return [json_sanitize(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def configure_logging(log_level: str, log_file: Path):
    root = logging.getLogger()
    root.setLevel(getattr(logging, log_level))

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")

    console = logging.StreamHandler()
    console.setLevel(getattr(logging, log_level))
    console.setFormatter(formatter)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(getattr(logging, log_level))
    file_handler.setFormatter(formatter)

    root.handlers.clear()
    root.addHandler(console)
    root.addHandler(file_handler)


def parse_args(argv: list[str] | None):
    parser = argparse.ArgumentParser(
        description="Profile monolithic GTEP model construction and solve."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="PH-style YAML config.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output root. Default: naerm_one_off_tests. A timestamped "
            "subdirectory is created unless --no-timestamp is used."
        ),
    )
    parser.add_argument(
        "--no-timestamp",
        action="store_true",
        help="Write directly to output-dir without creating timestamped subdir.",
    )
    parser.add_argument(
        "--enable-cprofile",
        action="store_true",
        help="Enable cProfile around data/model/GDP construction phases.",
    )
    parser.add_argument(
        "--cprofile-limit",
        type=int,
        default=100,
        help="Number of cProfile rows to print.",
    )
    parser.add_argument(
        "--pyomo-construction-log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"],
        help="Enable Pyomo construction timing logs at this level. Use NONE to disable.",
    )
    parser.add_argument(
        "--skip-gdp",
        action="store_true",
        help="Skip GDP transformation.",
    )
    parser.add_argument(
        "--skip-solve",
        action="store_true",
        help="Skip solver call.",
    )
    parser.add_argument(
        "--write-basic-solution",
        action="store_true",
        help="Write basic nonzero variable solution after solve.",
    )
    parser.add_argument(
        "--nonzero-tolerance",
        type=float,
        default=1.0e-8,
    )
    parser.add_argument(
        "--max-variable-records",
        type=int,
        default=100000,
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )

    args = parser.parse_args(argv)

    if args.pyomo_construction_log_level == "NONE":
        args.pyomo_construction_log_level = None

    return args


if __name__ == "__main__":
    sys.exit(main())
