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
"""Solver utilities for GTEP Progressive Hedging.

The initial PH implementation uses Pyomo's standard ``SolverFactory`` interface.
It does not require APPSI, persistent solver interfaces, serialized solver
state, serialized Pyomo models, or warm starts.

Primary target solver
---------------------
XPRESS is the initial target solver. This module supports the required solver
controls:

* relative MIP gap,
* time limit,
* threads,
* solver log file,
* tee output.

No nonanticipative-variable matching is performed in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import os

import logging
import time

import pyomo.environ as pyo
from pyomo.opt import SolverFactory
from pyomo.opt import SolverResults
from pyomo.opt import SolverStatus
from pyomo.opt import TerminationCondition

from gtep.algorithms.progressive_hedging.config import SolverConfig

logger = logging.getLogger("gtep.algorithms.progressive_hedging.solver")


@dataclass(frozen=True)
class SolveOutcome:
    """Summary of one Pyomo solve."""

    results: SolverResults
    solver_name: str
    solver_status: str
    termination_condition: str
    acceptable: bool
    solve_time_sec: float
    objective: float | None
    lower_bound: float | None
    upper_bound: float | None
    gap: float | None
    log_file: str | None

    def to_jsonable(self) -> dict[str, Any]:
        """Return a JSON-serializable solve summary."""
        return {
            "name": self.solver_name,
            "status": self.solver_status,
            "termination_condition": self.termination_condition,
            "acceptable": self.acceptable,
            "objective": self.objective,
            "lower_bound": self.lower_bound,
            "upper_bound": self.upper_bound,
            "gap": self.gap,
            "solve_time_sec": self.solve_time_sec,
            "log_file": self.log_file,
        }


def create_solver(solver_config: SolverConfig):
    """Create and configure a Pyomo solver.

    Parameters
    ----------
    solver_config:
        Solver configuration from the PH YAML config.

    Returns
    -------
    pyomo.opt.base.solvers.OptSolver
        Configured Pyomo solver object.

    Raises
    ------
    RuntimeError
        If the requested solver is unavailable.
    """
    apply_solver_environment(solver_config)

    solver_name = solver_config.name
    opt = SolverFactory(solver_name)

    if opt is None or not opt.available(exception_flag=False):
        raise RuntimeError(
            f"Requested solver {solver_name!r} is not available through "
            "Pyomo SolverFactory."
        )

    configure_solver_options(opt, solver_config)
    return opt


def apply_solver_environment(solver_config: SolverConfig) -> None:
    """Apply solver-related environment variables before solver construction.

    This supports both generic solver environment variables and a convenience
    license-file setting. For XPRESS, the default license environment variable
    is ``XPAUTH_PATH`` unless ``solver.license_env_var`` is explicitly set.
    """
    for key, value in solver_config.environment.items():
        os.environ[str(key)] = str(value)
        logger.info("Set solver environment variable %s from PH config.", key)

    if solver_config.license_file is None:
        return

    license_path = Path(solver_config.license_file).expanduser()

    if not license_path.exists():
        raise FileNotFoundError(
            f"Configured solver.license_file does not exist: {license_path}"
        )

    env_var = solver_config.license_env_var
    if env_var is None:
        env_var = default_license_env_var(solver_config.name)

    if env_var is None:
        raise ValueError(
            "solver.license_file was provided, but no default license "
            f"environment variable is known for solver {solver_config.name!r}. "
            "Please set solver.license_env_var explicitly."
        )

    os.environ[env_var] = str(license_path)
    logger.info(
        "Set solver license environment variable %s from solver.license_file.",
        env_var,
    )


def default_license_env_var(solver_name: str) -> str | None:
    """Return default license environment variable for a solver, if known."""
    normalized = solver_name.lower()

    if normalized in {"xpress", "xpress_direct", "xpress_persistent"}:
        return "XPAUTH_PATH"

    return None


def configure_solver_options(opt: Any, solver_config: SolverConfig) -> None:
    """Apply configured solver options to a Pyomo solver object.

    The initial implementation targets XPRESS. For other solvers, no generic
    cross-solver option mapping is attempted.
    """
    solver_name = solver_config.name.lower()

    if solver_name == "xpress":
        _configure_xpress_options(opt, solver_config)
    else:
        logger.warning(
            "No solver-option mapping is implemented for solver %s. "
            "Proceeding without applying mip_gap, time_limit, or threads.",
            solver_config.name,
        )


def solve_model(
    model: pyo.ConcreteModel,
    solver_config: SolverConfig,
    *,
    log_file: str | Path | None = None,
) -> SolveOutcome:
    """Solve a Pyomo model and return a structured solve outcome.

    Parameters
    ----------
    model:
        Pyomo model to solve.
    solver_config:
        Solver configuration.
    log_file:
        Optional per-solve log file. If omitted, ``solver_config.log_file`` is
        used. If both are omitted, no logfile argument is supplied to Pyomo.

    Returns
    -------
    SolveOutcome
        Structured solve metadata and raw Pyomo results.
    """
    opt = create_solver(solver_config)

    effective_log_file = log_file
    if effective_log_file is None:
        effective_log_file = solver_config.log_file

    if effective_log_file is not None:
        effective_log_file = Path(effective_log_file)
        effective_log_file.parent.mkdir(parents=True, exist_ok=True)

    solve_kwargs: dict[str, Any] = {
        "tee": solver_config.tee,
    }

    if effective_log_file is not None:
        solve_kwargs["logfile"] = str(effective_log_file)

    logger.info(
        "Solving PH scenario model with solver=%s, tee=%s, logfile=%s",
        solver_config.name,
        solver_config.tee,
        effective_log_file,
    )

    start = time.perf_counter()
    results = opt.solve(model, **solve_kwargs)
    solve_time_sec = time.perf_counter() - start

    solver_status = _solver_status_to_string(results)
    termination_condition = _termination_condition_to_string(results)
    acceptable = is_acceptable_termination_condition(
        results,
        solver_config.acceptable_termination_conditions,
    )

    objective = get_active_objective_value(model)
    lower_bound = get_result_bound(results, "lower_bound")
    upper_bound = get_result_bound(results, "upper_bound")
    gap = compute_relative_gap(lower_bound, upper_bound)

    logger.info(
        "Solve complete: status=%s termination=%s acceptable=%s time=%.3f sec "
        "objective=%s lower_bound=%s upper_bound=%s gap=%s",
        solver_status,
        termination_condition,
        acceptable,
        solve_time_sec,
        objective,
        lower_bound,
        upper_bound,
        gap,
    )

    return SolveOutcome(
        results=results,
        solver_name=solver_config.name,
        solver_status=solver_status,
        termination_condition=termination_condition,
        acceptable=acceptable,
        solve_time_sec=solve_time_sec,
        objective=objective,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        gap=gap,
        log_file=None if effective_log_file is None else str(effective_log_file),
    )


def is_acceptable_termination_condition(
    results: SolverResults,
    acceptable_termination_conditions: list[str] | tuple[str, ...] | set[str],
) -> bool:
    """Return true if a Pyomo solve termination condition is PH-acceptable.

    This function checks solver status and termination condition only. It does
    not guarantee that all nonanticipative variable values are available. The
    scenario-result writer should still require successful value extraction.
    """
    status = results.solver.status
    termination_condition = results.solver.termination_condition

    acceptable_conditions = {
        _normalize_condition_name(condition)
        for condition in acceptable_termination_conditions
    }

    normalized_tc = _normalize_condition_name(termination_condition)

    status_ok = status in {
        SolverStatus.ok,
        SolverStatus.warning,
    }

    return bool(status_ok and normalized_tc in acceptable_conditions)


def get_active_objective_value(model: pyo.ConcreteModel) -> float | None:
    """Return the value of the sole active objective, if available.

    If multiple active objectives are present, this raises an error because PH
    scenario models should have exactly one active objective after augmentation.
    """
    active_objectives = list(model.component_data_objects(pyo.Objective, active=True))

    if not active_objectives:
        return None

    if len(active_objectives) > 1:
        raise RuntimeError(
            "Expected at most one active objective after PH augmentation, but "
            f"found {len(active_objectives)} active objectives."
        )

    value = pyo.value(active_objectives[0], exception=False)
    if value is None:
        return None

    return float(value)


def get_result_bound(results: SolverResults, bound_name: str) -> float | None:
    """Extract a numeric bound from Pyomo results if available.

    Parameters
    ----------
    results:
        Pyomo solver results.
    bound_name:
        Usually ``"lower_bound"`` or ``"upper_bound"``.

    Returns
    -------
    float | None
        Extracted bound, or ``None`` if unavailable.
    """
    if not hasattr(results, "problem") or len(results.problem) == 0:
        return None

    problem = results.problem[0]

    if not hasattr(problem, bound_name):
        return None

    bound = getattr(problem, bound_name)

    if bound is None:
        return None

    try:
        return float(bound)
    except (TypeError, ValueError):
        return None


def compute_relative_gap(
    lower_bound: float | None,
    upper_bound: float | None,
) -> float | None:
    """Compute a relative MIP gap from lower and upper bounds.

    Returns ``None`` if either bound is unavailable.
    """
    if lower_bound is None or upper_bound is None:
        return None

    denominator = max(1.0, abs(upper_bound))
    return abs(upper_bound - lower_bound) / denominator


def result_has_acceptable_status_and_values(
    outcome: SolveOutcome,
    model: pyo.ConcreteModel,
) -> bool:
    """Return true if the solve outcome is acceptable and objective-valued.

    This is a conservative helper for callers that want to reject time-limit
    solves without an incumbent objective.
    """
    if not outcome.acceptable:
        return False

    return get_active_objective_value(model) is not None


def _configure_xpress_options(opt: Any, solver_config: SolverConfig) -> None:
    """Apply XPRESS-specific options."""
    if solver_config.mip_gap is not None:
        opt.options["miprelstop"] = solver_config.mip_gap

    if solver_config.time_limit is not None:
        opt.options["maxtime"] = solver_config.time_limit

    if solver_config.threads is not None:
        opt.options["threads"] = solver_config.threads


def _solver_status_to_string(results: SolverResults) -> str:
    """Return normalized solver status string."""
    return str(results.solver.status)


def _termination_condition_to_string(results: SolverResults) -> str:
    """Return normalized termination-condition string."""
    return str(results.solver.termination_condition)


def _normalize_condition_name(condition: Any) -> str:
    """Normalize a termination-condition enum or string for comparison."""
    if isinstance(condition, TerminationCondition):
        return str(condition)

    return str(condition)
