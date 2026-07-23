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
"""Local synchronous driver for GTEP Progressive Hedging.

This script is a development and smoke-test driver. It runs the same PH
scenario-solve and PH-update logic as the Torc dynamic orchestrator, but without
Torc. Scenario solves are executed sequentially in the local Python process.

Purpose
-------
Use this script to validate:

* YAML config loading,
* output-directory creation,
* initial PH state creation,
* single-representative-period scenario construction,
* GDP transformation on small scenario models,
* nonanticipative variable collection using ComponentUID-backed VariableID,
* PH objective augmentation,
* solver execution,
* scenario-result writing,
* consensus and multiplier updates,
* convergence-history writing,
* final solution writing.

This script intentionally does not use model serialization, persistent solvers,
or warm starts.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from gtep.algorithms.progressive_hedging.config import (
    convergence_history_csv_path,
    convergence_history_json_path,
    final_solution_dir,
    iteration_summary_path,
    load_ph_config,
    prepare_output_directory,
    scenario_result_path,
    state_path,
    write_effective_config,
)
from gtep.algorithms.progressive_hedging.convergence import (
    compute_weighted_expected_objective,
    process_completed_iteration,
)
from gtep.algorithms.progressive_hedging.scenario_data import (
    build_full_data,
    build_scenario_infos,
    get_representative_period_weights,
    get_scenario_probabilities,
    scenario_infos_to_jsonable,
)
from gtep.algorithms.progressive_hedging.solution_io import (
    append_convergence_history,
    export_durable_results,
    make_iteration_summary,
    read_scenario_results,
    scenario_values_from_results,
    write_final_solution,
    write_iteration_summary,
)
from gtep.algorithms.progressive_hedging.solve_scenario import (
    solve_scenario_from_files,
)
from gtep.algorithms.progressive_hedging.state import (
    PHState,
    create_initial_state,
    write_state,
)

logger = logging.getLogger("gtep.examples.progressive_hedging.run_ph_local")


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    args = _parse_args(argv)
    _configure_logging(args.log_level)

    try:
        run_ph_local(
            config_path=args.config,
            stop_after_iteration=args.stop_after_iteration,
            force_resolve=args.force_resolve,
        )
    except Exception:
        logger.exception("Local PH run failed.")
        return 1

    return 0


def run_ph_local(
    *,
    config_path: str | Path,
    stop_after_iteration: int | None = None,
    force_resolve: bool = False,
) -> None:
    """Run Progressive Hedging locally and synchronously.

    Parameters
    ----------
    config_path:
        User-facing PH YAML configuration.
    stop_after_iteration:
        Optional debug limit. If provided, stop after processing this iteration,
        even if PH has not converged.
    force_resolve:
        If true, rerun scenario solves even when result files already exist.
    """
    cfg = load_ph_config(config_path)
    prepare_output_directory(cfg, clean=cfg.run.clean_output_dir)
    write_effective_config(cfg)

    state = _initialize_state_zero(cfg)
    current_state_path = state_path(cfg, state.iteration)
    write_state(state, current_state_path)

    while True:
        logger.info("Starting local PH iteration %s", state.iteration)

        _solve_iteration_scenarios(
            cfg=cfg,
            state=state,
            current_state_path=current_state_path,
            force_resolve=force_resolve,
        )

        scenario_result_paths = [
            scenario_result_path(cfg, state.iteration, scenario_id)
            for scenario_id in state.scenario_ids
        ]

        scenario_results = read_scenario_results(scenario_result_paths)

        unacceptable = [
            scenario_id
            for scenario_id, result in sorted(scenario_results.items())
            if not result.acceptable
        ]
        if unacceptable:
            raise RuntimeError(
                "One or more scenario solves were not acceptable for PH. "
                f"Unacceptable scenario id(s): {unacceptable}"
            )

        scenario_values = scenario_values_from_results(
            scenario_results,
            require_acceptable=True,
        )

        update_result = process_completed_iteration(
            state,
            scenario_values,
            convergence_tolerance=cfg.progressive_hedging.convergence_tolerance,
            metadata={
                "driver": "local",
                "processed_state_file": str(current_state_path),
            },
        )

        expected_objective = _compute_expected_objective_if_available(
            scenario_results=scenario_results,
            scenario_probabilities=state.scenario_probabilities,
            scenario_ids=state.scenario_ids,
        )

        effective_converged = (
            update_result.converged
            and state.iteration >= cfg.progressive_hedging.min_iterations
        )

        if update_result.converged and not effective_converged:
            logger.info(
                "PH residual convergence reached at iteration %s, but "
                "min_iterations=%s, so continuing.",
                state.iteration,
                cfg.progressive_hedging.min_iterations,
            )

        summary = make_iteration_summary(
            iteration=state.iteration,
            converged=effective_converged,
            residuals=update_result.residuals,
            scenario_results=scenario_results,
            expected_objective=expected_objective,
        )

        write_iteration_summary(
            iteration_summary_path(cfg, state.iteration),
            summary,
        )

        append_convergence_history(
            summary=summary,
            history_json_path=(
                convergence_history_json_path(cfg)
                if "json" in cfg.output.convergence_history_formats
                else None
            ),
            history_csv_path=(
                convergence_history_csv_path(cfg)
                if "csv" in cfg.output.convergence_history_formats
                else None
            ),
        )

        next_state = update_result.next_state
        next_state_path = state_path(cfg, next_state.iteration)
        write_state(next_state, next_state_path)

        logger.info(
            "Completed PH iteration %s: converged=%s max_abs_residual=%s",
            state.iteration,
            update_result.converged,
            update_result.residuals.max_abs_residual,
        )

        if effective_converged:
            logger.info("PH converged after iteration %s.", state.iteration)

            if cfg.output.save_final_solution:
                write_final_solution(
                    final_solution_dir(cfg),
                    final_state=next_state,
                    scenario_results=scenario_results,
                    summary=summary,
                    output_config=cfg.output,
                )

            return

        if stop_after_iteration is not None and state.iteration >= stop_after_iteration:
            logger.info(
                "Stopping local PH run after requested iteration %s.",
                stop_after_iteration,
            )
            return

        if next_state.iteration > cfg.progressive_hedging.max_iterations:
            logger.warning(
                "PH reached max_iterations=%s without convergence.",
                cfg.progressive_hedging.max_iterations,
            )

            if cfg.output.save_final_solution:
                write_final_solution(
                    final_solution_dir(cfg),
                    final_state=next_state,
                    scenario_results=scenario_results,
                    summary=summary,
                    output_config=cfg.output,
                )

            if cfg.output.durable_output_dir is not None:
                durable_paths = export_durable_results(
                    run_output_dir=cfg.run.output_dir,
                    durable_output_dir=cfg.output.durable_output_dir,
                )
                logger.info("Exported durable PH results: %s", durable_paths)

            return

        state = next_state
        current_state_path = next_state_path


def _initialize_state_zero(cfg) -> PHState:
    """Build full data metadata and create initial unregularized PH state."""
    full_data = build_full_data(cfg)

    scenario_weights = get_representative_period_weights(full_data)
    scenario_probabilities = get_scenario_probabilities(
        full_data,
        normalize=cfg.progressive_hedging.normalize_scenario_probabilities,
    )
    scenario_infos = build_scenario_infos(
        full_data,
        normalize_probabilities=(
            cfg.progressive_hedging.normalize_scenario_probabilities
        ),
    )

    scenario_ids = [info.scenario_id for info in scenario_infos]

    metadata = {
        "driver": "local",
        "scenario_infos": scenario_infos_to_jsonable(scenario_infos),
        "config_path": str(cfg.config_path) if cfg.config_path is not None else None,
    }

    return create_initial_state(
        scenario_ids=scenario_ids,
        scenario_probabilities=scenario_probabilities,
        scenario_weights=scenario_weights,
        rho=cfg.progressive_hedging.rho,
        metadata=metadata,
    )


def _solve_iteration_scenarios(
    *,
    cfg,
    state: PHState,
    current_state_path: Path,
    force_resolve: bool,
) -> None:
    """Solve all scenario subproblems for one PH iteration sequentially."""
    for scenario_id in state.scenario_ids:
        result_path = scenario_result_path(cfg, state.iteration, scenario_id)

        if result_path.exists() and not force_resolve:
            logger.info(
                "Skipping scenario %s iteration %s because result already exists: %s",
                scenario_id,
                state.iteration,
                result_path,
            )
            continue

        logger.info(
            "Solving scenario %s for local PH iteration %s",
            scenario_id,
            state.iteration,
        )

        solve_scenario_from_files(
            config_path=cfg.config_path,
            state_path=current_state_path,
            scenario_id=scenario_id,
            iteration=state.iteration,
        )


def _compute_expected_objective_if_available(
    *,
    scenario_results,
    scenario_probabilities: dict[int, float],
    scenario_ids: list[int],
) -> float | None:
    """Compute probability-weighted objective if all objectives are present."""
    objective_by_scenario = {
        int(scenario_id): result.objective
        for scenario_id, result in scenario_results.items()
    }

    if any(value is None for value in objective_by_scenario.values()):
        return None

    return compute_weighted_expected_objective(
        {
            scenario_id: float(objective)
            for scenario_id, objective in objective_by_scenario.items()
            if objective is not None
        },
        scenario_probabilities,
        scenario_ids,
    )


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run GTEP Progressive Hedging locally without Torc."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to user-facing PH YAML configuration.",
    )
    parser.add_argument(
        "--stop-after-iteration",
        type=int,
        default=None,
        help=(
            "Optional debug limit. Stop after processing this PH iteration, "
            "even if convergence has not been reached."
        ),
    )
    parser.add_argument(
        "--force-resolve",
        action="store_true",
        help="Rerun scenario solves even if scenario result files already exist.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Python logging level.",
    )

    return parser.parse_args(argv)


def _configure_logging(log_level: str) -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


if __name__ == "__main__":
    sys.exit(main())
