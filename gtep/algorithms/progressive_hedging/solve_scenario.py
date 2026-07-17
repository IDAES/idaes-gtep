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
"""Scenario-solve entry point for GTEP Progressive Hedging.

Each invocation solves one representative-period PH scenario for one PH
iteration. This module is intended to be launched by the Torc dynamic
orchestrator, but it can also be run directly for debugging.

Responsibilities
----------------
1. Load user-facing YAML PH configuration.
2. Load machine-written PH iteration state JSON.
3. Build full GTEP data and cost data.
4. Slice one representative period into a one-scenario data object.
5. Build the GTEP scenario model.
6. Apply the configured GDP transformation to the small scenario model.
7. Collect nonanticipative variables using ``VariableID`` / ``ComponentUID``.
8. Add the PH scenario objective.
9. Solve with Pyomo ``SolverFactory``.
10. Write required per-scenario JSON output.

Design constraints
------------------
* No Pyomo model serialization.
* No persistent solver interface requirement.
* No warm starts.
* No matching variables by names or strings.
* Variable identity is represented by ``VariableID`` backed by Pyomo
  ``ComponentUID``.
"""

from __future__ import annotations

import argparse
import logging
import sys
import traceback
from pathlib import Path
from typing import Any

import pyomo.environ as pyo

from gtep.gtep_model import ExpansionPlanningModel

from gtep.algorithms.progressive_hedging.config import (
    PHConfig,
    load_ph_config,
    prepare_output_directory,
    scenario_log_path,
    scenario_result_path,
)
from gtep.algorithms.progressive_hedging.nonanticipativity import (
    assert_same_variable_id_set,
    collect_nonanticipative_variables,
)
from gtep.algorithms.progressive_hedging.ph_model import (
    PHObjectiveOptions,
    add_progressive_hedging_block,
)
from gtep.algorithms.progressive_hedging.scenario_data import (
    build_cost_data,
    build_full_data,
    get_scenario_info,
    make_single_representative_period_data,
)
from gtep.algorithms.progressive_hedging.solution_io import (
    write_failed_scenario_result,
    write_scenario_result,
)
from gtep.algorithms.progressive_hedging.solver import solve_model
from gtep.algorithms.progressive_hedging.state import PHState, read_state

logger = logging.getLogger("gtep.algorithms.progressive_hedging.solve_scenario")


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    args = _parse_args(argv)
    _configure_logging(args.log_level)

    try:
        solve_scenario_from_files(
            config_path=args.config,
            state_path=args.state,
            scenario_id=args.scenario_id,
            iteration=args.iteration,
        )
    except Exception:  # pragma: no cover - command-line failure path
        logger.exception("Scenario solve failed with an unhandled exception.")
        return 1

    return 0


def solve_scenario_from_files(
    *,
    config_path: str | Path,
    state_path: str | Path,
    scenario_id: int,
    iteration: int,
) -> Path:
    """Solve one PH scenario using configuration and state files.

    Parameters
    ----------
    config_path:
        User-facing YAML PH configuration path.
    state_path:
        Machine-written PH state JSON path for this iteration.
    scenario_id:
        One-based representative-period scenario id.
    iteration:
        PH iteration index to solve.

    Returns
    -------
    pathlib.Path
        Path to the written scenario-result JSON file.
    """
    cfg = load_ph_config(config_path)

    # Do not clean the output directory from a scenario solve. The state file
    # and other scenario results are owned by the driver/orchestrator.
    prepare_output_directory(cfg, clean=False)

    state = read_state(state_path)

    state = read_state(state_path)
    _validate_requested_iteration(state, iteration)
    _validate_requested_scenario(state, scenario_id)

    output_path = scenario_result_path(cfg, iteration, scenario_id)

    scenario_info = None
    try:
        full_data = build_full_data(cfg)
        cost_data = build_cost_data(cfg)

        scenario_info = get_scenario_info(
            full_data,
            scenario_id,
            normalize_probabilities=(
                cfg.progressive_hedging.normalize_scenario_probabilities
            ),
        )

        scenario_data = make_single_representative_period_data(
            full_data,
            scenario_id,
            copy_representative_model_data=False,
        )

        model_object = _build_scenario_model(
            cfg=cfg,
            scenario_data=scenario_data,
            cost_data=cost_data,
        )

        model = model_object.model

        _apply_gdp_transformation(model, cfg)

        nonant_var_map = collect_nonanticipative_variables(
            model,
            cfg.progressive_hedging.nonanticipativity,
        )

        _validate_or_initialize_nonanticipative_state(
            state=state,
            nonant_var_ids=nonant_var_map.keys(),
        )

        _add_ph_objective(
            model=model,
            state=state,
            scenario_id=scenario_id,
            scenario_probability=scenario_info.probability,
            nonant_var_map=nonant_var_map,
        )

        debug_model_path = (
            cfg.run.output_dir
            / "debug_models"
            / f"iter_{iteration:03d}_scenario_{scenario_id:03d}.lp"
        )
        debug_model_path.parent.mkdir(parents=True, exist_ok=True)
        model.write(
            str(debug_model_path),
            io_options={"symbolic_solver_labels": True},
        )
        logger.info("Wrote debug LP model to %s", debug_model_path)

        log_file = _resolve_solver_log_file(cfg, iteration, scenario_id)
        solve_outcome = solve_model(
            model,
            cfg.solver,
            log_file=log_file,
        )

        if not solve_outcome.acceptable:
            message = (
                "Solver termination condition is not acceptable for PH: "
                f"status={solve_outcome.solver_status}, "
                f"termination_condition={solve_outcome.termination_condition}"
            )
            return write_failed_scenario_result(
                output_path,
                scenario_id=scenario_id,
                representative_period=scenario_info.representative_period,
                representative_date=scenario_info.representative_date,
                iteration=iteration,
                solver_name=cfg.solver.name,
                solver_status=solve_outcome.solver_status,
                termination_condition=solve_outcome.termination_condition,
                message=message,
            )

        return write_scenario_result(
            output_path,
            scenario_id=scenario_id,
            representative_period=scenario_info.representative_period,
            representative_date=scenario_info.representative_date,
            iteration=iteration,
            solve_outcome=solve_outcome,
            nonant_var_map=nonant_var_map,
            save_nonanticipative_metadata=cfg.output.save_nonanticipative_metadata,
        )

    except Exception as err:
        logger.error("Scenario %s iteration %s failed: %s", scenario_id, iteration, err)
        logger.debug("Traceback:\n%s", traceback.format_exc())

        representative_period = scenario_id
        representative_date = "unknown"
        if scenario_info is not None:
            representative_period = scenario_info.representative_period
            representative_date = scenario_info.representative_date

        write_failed_scenario_result(
            output_path,
            scenario_id=scenario_id,
            representative_period=representative_period,
            representative_date=representative_date,
            iteration=iteration,
            solver_name=cfg.solver.name if "cfg" in locals() else "unknown",
            message=f"{type(err).__name__}: {err}",
        )
        raise


def _build_scenario_model(
    *,
    cfg: PHConfig,
    scenario_data: Any,
    cost_data: Any,
) -> ExpansionPlanningModel:
    """Build an ``ExpansionPlanningModel`` for one PH scenario."""
    model_object = ExpansionPlanningModel(
        data=scenario_data,
        cost_data=cost_data,
    )

    for key, value in cfg.model.as_model_options().items():
        model_object.config[key] = value

    model_object.create_model()
    return model_object


def _apply_gdp_transformation(
    model: pyo.ConcreteModel,
    cfg: PHConfig,
) -> None:
    """Apply configured GDP transformation to the scenario model."""
    transformation_name = cfg.model.gdp_transformation

    logger.info("Applying GDP transformation %s", transformation_name)

    transformation = pyo.TransformationFactory(transformation_name)
    if transformation is None:
        raise RuntimeError(
            f"Could not create GDP transformation {transformation_name!r}."
        )

    transformation.apply_to(model)


def _validate_or_initialize_nonanticipative_state(
    *,
    state: PHState,
    nonant_var_ids: Any,
) -> None:
    """Validate collected variable IDs against state when PH terms are enabled.

    Iteration 0 has PH terms disabled and may not yet contain variable IDs.
    Later iterations must exactly match the variable IDs recorded in state.
    """
    collected_ids = list(nonant_var_ids)

    if state.ph_terms_enabled:
        assert_same_variable_id_set(
            state.variable_ids,
            collected_ids,
            first_label="PH state",
            second_label="scenario model",
        )


def _add_ph_objective(
    *,
    model: pyo.ConcreteModel,
    state: PHState,
    scenario_id: int,
    scenario_probability: float,
    nonant_var_map: MappingLike,
) -> None:
    """Add PH objective terms to the scenario model."""
    if state.ph_terms_enabled:
        xbar = state.xbar
        multipliers = state.get_scenario_multipliers(scenario_id)
        rho = state.get_rho()
        options = PHObjectiveOptions(
            include_multiplier_terms=True,
            include_regularization_terms=True,
            deactivate_existing_objectives=True,
        )
    else:
        # Initial unregularized iteration. Still use PH objective construction
        # so investment costs are scenario-probability scaled and deficit-penalty
        # terms are excluded.
        xbar = {variable_id: 0.0 for variable_id in nonant_var_map}
        multipliers = {variable_id: 0.0 for variable_id in nonant_var_map}
        rho = state.get_rho()
        options = PHObjectiveOptions(
            include_multiplier_terms=False,
            include_regularization_terms=False,
            deactivate_existing_objectives=True,
        )

    add_progressive_hedging_block(
        model,
        nonant_var_map,
        scenario_probability=scenario_probability,
        xbar=xbar,
        multipliers=multipliers,
        rho=rho,
        options=options,
    )


def _resolve_solver_log_file(
    cfg: PHConfig,
    iteration: int,
    scenario_id: int,
) -> Path | None:
    """Resolve the solver log file for this scenario solve.

    If ``solver.log_file`` is omitted, a per-scenario default log file is used.
    If ``solver.log_file`` is provided and contains ``{iteration}`` or
    ``{scenario_id}``, those fields are formatted. Otherwise the configured path
    is used exactly.
    """
    configured = cfg.solver.log_file

    if configured is None:
        return scenario_log_path(cfg, iteration, scenario_id)

    configured_text = str(configured)
    if "{iteration}" in configured_text or "{scenario_id}" in configured_text:
        return Path(
            configured_text.format(
                iteration=iteration,
                scenario_id=scenario_id,
            )
        )

    return Path(configured)


def _validate_requested_iteration(state: PHState, iteration: int) -> None:
    """Validate requested iteration against loaded state."""
    if iteration != state.iteration:
        raise ValueError(
            f"Requested iteration {iteration} does not match loaded PH state "
            f"iteration {state.iteration}."
        )


def _validate_requested_scenario(state: PHState, scenario_id: int) -> None:
    """Validate requested scenario id against loaded state."""
    if scenario_id not in state.scenario_ids:
        raise ValueError(
            f"Requested scenario id {scenario_id} is not present in PH state "
            f"scenario ids {state.scenario_ids}."
        )


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Solve one representative-period GTEP PH scenario."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to user-facing PH YAML configuration.",
    )
    parser.add_argument(
        "--state",
        required=True,
        help="Path to machine-written PH state JSON for this iteration.",
    )
    parser.add_argument(
        "--scenario-id",
        required=True,
        type=int,
        help="One-based representative-period scenario id.",
    )
    parser.add_argument(
        "--iteration",
        required=True,
        type=int,
        help="PH iteration index to solve.",
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


# Avoid importing typing.Protocol just for one internal annotation; this keeps
# runtime behavior simple while allowing the function signatures above to remain
# readable.
MappingLike = Any


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
