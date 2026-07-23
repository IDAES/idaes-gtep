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
"""Torc dynamic orchestrator for GTEP Progressive Hedging.

This module implements the dynamic Torc orchestration loop for representative-
period Progressive Hedging.

Execution model
---------------
Each orchestrator invocation is responsible for one of two actions:

1. If scenario results for the current PH state iteration do not exist, spawn
   one scenario-solve job per representative period plus an orchestrator
   continuation depending on all scenario jobs.

2. If scenario results for the current PH state iteration exist, process those
   results, compute consensus/residuals/multiplier updates, write the next PH
   state, and either:
      * call ``orch.converge(...)`` if PH converged or max iterations were
        reached, or
      * spawn the next iteration's scenario jobs plus another continuation.

Design constraints
------------------
* Scenario solves rebuild and GDP-transform one representative-period model per
  iteration.
* No Pyomo model serialization.
* No persistent solver interface requirement.
* No warm starts.
* Nonanticipative variable identity uses ``VariableID`` / ``ComponentUID``.
* This module does not match variables by names or string paths.
"""

from __future__ import annotations

import argparse
import logging
import shlex
import sys
from pathlib import Path
from typing import Any

from gtep.algorithms.progressive_hedging.config import (
    PHConfig,
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
from gtep.algorithms.progressive_hedging.state import (
    PHState,
    create_initial_state,
    read_state,
    write_state,
)

logger = logging.getLogger("gtep.algorithms.progressive_hedging.torc_orchestrator")


try:
    from torc import Orchestrator, SpawnJobModel
except ImportError as err:  # pragma: no cover
    Orchestrator = None
    SpawnJobModel = None
    _TORC_IMPORT_ERROR = err
else:
    _TORC_IMPORT_ERROR = None


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    args = _parse_args(argv)
    _configure_logging(args.log_level)

    try:
        run_torc_orchestrator(
            config_path=args.config,
            state_file=args.state,
            lineage_fallback=args.lineage,
        )
    except Exception:  # pragma: no cover - command-line failure path
        logger.exception("PH Torc orchestrator failed.")
        return 1

    return 0


def run_torc_orchestrator(
    *,
    config_path: str | Path,
    state_file: str | Path | None = None,
    lineage_fallback: str | None = None,
) -> None:
    """Run one Torc dynamic-orchestrator generation.

    Parameters
    ----------
    config_path:
        User-facing PH YAML configuration.
    state_file:
        Optional PH state JSON file. If omitted, the orchestrator initializes
        iteration 0 state.
    lineage_fallback:
        Optional Torc lineage fallback for seed invocations.
    """
    _require_torc()

    cfg = load_ph_config(config_path)

    # Only the seed orchestrator invocation should clean the output directory.
    # Continuations must not clean because they need prior state and scenario
    # result files.
    is_seed_invocation = state_file is None
    prepare_output_directory(
        cfg, clean=(cfg.run.clean_output_dir and is_seed_invocation)
    )

    write_effective_config(cfg)

    orch = Orchestrator.from_env(lineage_fallback=lineage_fallback)

    logger.info(
        "Starting PH orchestrator generation=%s lineage=%s state_file=%s",
        orch.generation,
        orch.lineage,
        state_file,
    )

    if state_file is None:
        state = _initialize_state_zero(cfg)
        state_file = state_path(cfg, state.iteration)
        write_state(state, state_file)
    else:
        state = read_state(state_file)

    scenario_result_paths = [
        scenario_result_path(cfg, state.iteration, scenario_id)
        for scenario_id in state.scenario_ids
    ]

    result_status = _scenario_result_file_status(scenario_result_paths)

    if result_status == "none":
        _spawn_iteration(
            orch=orch,
            cfg=cfg,
            state=state,
            state_file=Path(state_file),
        )
        return

    if result_status == "partial":
        missing = [
            str(path) for path in scenario_result_paths if not Path(path).exists()
        ]
        raise RuntimeError(
            "Some, but not all, scenario results exist for PH iteration "
            f"{state.iteration}. Missing result file(s): {missing}"
        )

    _process_existing_iteration_results(
        orch=orch,
        cfg=cfg,
        state=state,
        state_file=Path(state_file),
        scenario_result_paths=scenario_result_paths,
    )


def _initialize_state_zero(cfg: PHConfig) -> PHState:
    """Build data metadata and create initial unregularized PH state."""
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


def _process_existing_iteration_results(
    *,
    orch: Any,
    cfg: PHConfig,
    state: PHState,
    state_file: Path,
    scenario_result_paths: list[Path],
) -> None:
    """Process completed scenario solves for the current PH state iteration."""
    scenario_results = read_scenario_results(scenario_result_paths)

    expected_scenario_ids = set(state.scenario_ids)
    actual_scenario_ids = set(scenario_results)

    if actual_scenario_ids != expected_scenario_ids:
        raise RuntimeError(
            "Scenario result set does not match PH state scenario ids. "
            f"Expected {sorted(expected_scenario_ids)}, got {sorted(actual_scenario_ids)}."
        )

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
            "processed_state_file": str(state_file),
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
    next_state_file = state_path(cfg, next_state.iteration)
    write_state(next_state, next_state_file)

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

        orch.converge(
            state={
                "converged": True,
                "final_iteration": state.iteration,
                "state_file": str(next_state_file),
                "summary_file": str(iteration_summary_path(cfg, state.iteration)),
            }
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

        orch.converge(
            state={
                "converged": False,
                "reason": "max_iterations_reached",
                "last_completed_iteration": state.iteration,
                "state_file": str(next_state_file),
                "summary_file": str(iteration_summary_path(cfg, state.iteration)),
            }
        )
        return

    _spawn_iteration(
        orch=orch,
        cfg=cfg,
        state=next_state,
        state_file=next_state_file,
    )


def _spawn_iteration(
    *,
    orch: Any,
    cfg: PHConfig,
    state: PHState,
    state_file: Path,
) -> None:
    """Spawn all scenario solves for one PH iteration plus continuation."""
    iteration = state.iteration

    iteration_directory = cfg.run.output_dir / "iterations" / f"iter_{iteration:03d}"
    iteration_directory.mkdir(parents=True, exist_ok=True)

    scenario_jobs = []
    scenario_job_names = []

    python_executable = shlex.quote(sys.executable)
    config_arg = shlex.quote(str(cfg.config_path))
    state_arg = shlex.quote(str(state_file))

    for scenario_id in state.scenario_ids:
        job_name = _scenario_job_name(
            lineage=orch.lineage,
            iteration=iteration,
            scenario_id=scenario_id,
        )
        scenario_job_names.append(job_name)

        command = (
            f"{python_executable} -m gtep.algorithms.progressive_hedging.solve_scenario "
            f"--config {config_arg} "
            f"--state {state_arg} "
            f"--scenario-id {scenario_id} "
            f"--iteration {iteration}"
        )

        scenario_jobs.append(
            SpawnJobModel(
                name=job_name,
                command=command,
                resource_requirements=cfg.torc.resource_requirements.scenario,
                priority=cfg.torc.scenario_priority,
            )
        )

    continuation_name = _orchestrator_job_name(
        lineage=orch.lineage,
        generation=orch.generation + 1,
        iteration=iteration,
    )

    continuation_command = (
        f"{python_executable} -m gtep.algorithms.progressive_hedging.torc_orchestrator "
        f"--config {config_arg} "
        f"--state {state_arg}"
    )

    continuation_job = SpawnJobModel(
        name=continuation_name,
        command=continuation_command,
        resource_requirements=cfg.torc.resource_requirements.orchestrator,
        priority=cfg.torc.orchestrator_priority,
        depends_on=scenario_job_names,
        cancel_on_blocking_job_failure=cfg.torc.cancel_on_scenario_failure,
    )

    response = orch.spawn(
        jobs=scenario_jobs + [continuation_job],
        state={
            "iteration": iteration,
            "state_file": str(state_file),
            "scenario_jobs": scenario_job_names,
            "continuation_job": continuation_name,
        },
    )

    logger.info(
        "Spawned PH iteration %s: %s scenario job(s), continuation=%s, "
        "torc_iteration=%s",
        iteration,
        len(scenario_jobs),
        continuation_name,
        getattr(response, "iteration", None),
    )


def _scenario_result_file_status(paths: list[Path]) -> str:
    """Return 'none', 'partial', or 'all' for expected result files."""
    exists = [Path(path).exists() for path in paths]

    if not any(exists):
        return "none"

    if all(exists):
        return "all"

    return "partial"


def _compute_expected_objective_if_available(
    *,
    scenario_results: MappingLike,
    scenario_probabilities: dict[int, float],
    scenario_ids: list[int],
) -> float | None:
    """Compute expected objective if all scenario objective values are present."""
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


def _scenario_job_name(
    *,
    lineage: str,
    iteration: int,
    scenario_id: int,
) -> str:
    """Return Torc job name for one scenario solve."""
    return f"ph_solve_{lineage}_i{iteration:03d}_s{scenario_id:03d}"


def _orchestrator_job_name(
    *,
    lineage: str,
    generation: int,
    iteration: int,
) -> str:
    """Return Torc job name for an orchestrator continuation."""
    return f"ph_orch_{lineage}_g{generation:03d}_after_i{iteration:03d}"


def _require_torc() -> None:
    """Raise a clear error if torc-client is unavailable."""
    if Orchestrator is None or SpawnJobModel is None:  # pragma: no cover
        raise ImportError(
            "torc-client is required to run the Progressive Hedging Torc "
            "orchestrator."
        ) from _TORC_IMPORT_ERROR


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run one Torc dynamic-orchestrator generation for GTEP PH."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to user-facing PH YAML configuration.",
    )
    parser.add_argument(
        "--state",
        default=None,
        help=(
            "Path to PH state JSON. Omit for seed invocation, which initializes "
            "iteration 0 state."
        ),
    )
    parser.add_argument(
        "--lineage",
        default=None,
        help=(
            "Optional Torc lineage fallback for seed invocations. Spawned "
            "continuations inherit lineage from Torc environment."
        ),
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


# Avoid importing typing.Protocol just for one internal annotation.
MappingLike = Any


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
