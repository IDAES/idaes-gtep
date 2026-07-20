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
"""Solution and result I/O for GTEP Progressive Hedging.

This module writes and reads the machine-generated PH result files:

* per-scenario result JSON files,
* per-iteration summary JSON files,
* convergence history JSON/CSV files,
* final consensus and scenario-solution JSON files.

Design constraints
------------------
Nonanticipative variables are represented with ``VariableID`` objects backed by
Pyomo ``ComponentUID``. JSON files serialize variable IDs as records containing
serialized ``ComponentUID`` values. This module does not match variables by
names or manually parsed strings.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping
import shutil

import csv
import json
import logging

from gtep.algorithms.progressive_hedging.config import OutputConfig
from gtep.algorithms.progressive_hedging.convergence import ResidualSummary
from gtep.algorithms.progressive_hedging.nonanticipativity import (
    VariableID,
    deserialize_variable_value_records,
    extract_nonanticipative_metadata,
    extract_nonanticipative_values,
    serialize_nonanticipative_metadata,
    serialize_variable_value_records,
)
from gtep.algorithms.progressive_hedging.solver import SolveOutcome
from gtep.algorithms.progressive_hedging.state import PHState

logger = logging.getLogger("gtep.algorithms.progressive_hedging.solution_io")


@dataclass(frozen=True)
class ScenarioResult:
    """Parsed result from one PH scenario solve."""

    scenario_id: int
    representative_period: int
    representative_date: str
    iteration: int
    acceptable: bool
    solver: dict[str, Any]
    nonanticipative_values: dict[VariableID, float]
    nonanticipative_metadata: list[dict[str, Any]] | None = None
    message: str | None = None

    @property
    def objective(self) -> float | None:
        """Return solver objective value, if present."""
        value = self.solver.get("objective")
        return None if value is None else float(value)

    def to_jsonable(self) -> dict[str, Any]:
        """Return a JSON-serializable scenario result."""
        data: dict[str, Any] = {
            "scenario_id": self.scenario_id,
            "representative_period": self.representative_period,
            "representative_date": self.representative_date,
            "iteration": self.iteration,
            "acceptable": self.acceptable,
            "solver": _json_sanitize(self.solver),
            "nonanticipative_values": serialize_variable_value_records(
                self.nonanticipative_values
            ),
        }

        if self.nonanticipative_metadata is not None:
            data["nonanticipative_metadata"] = _json_sanitize(
                self.nonanticipative_metadata
            )

        if self.message is not None:
            data["message"] = self.message

        return data

    @classmethod
    def from_jsonable(cls, data: Mapping[str, Any]) -> "ScenarioResult":
        """Construct from a JSON-loaded mapping."""
        required = {
            "scenario_id",
            "representative_period",
            "representative_date",
            "iteration",
            "acceptable",
            "solver",
            "nonanticipative_values",
        }
        missing = required - set(data)
        if missing:
            raise ValueError(
                "Scenario result JSON is missing required field(s): "
                + ", ".join(sorted(missing))
            )

        return cls(
            scenario_id=int(data["scenario_id"]),
            representative_period=int(data["representative_period"]),
            representative_date=str(data["representative_date"]),
            iteration=int(data["iteration"]),
            acceptable=bool(data["acceptable"]),
            solver=dict(data["solver"]),
            nonanticipative_values=deserialize_variable_value_records(
                data["nonanticipative_values"]
            ),
            nonanticipative_metadata=data.get("nonanticipative_metadata"),
            message=data.get("message"),
        )


@dataclass(frozen=True)
class IterationSummary:
    """Summary for one completed PH iteration."""

    iteration: int
    converged: bool
    all_scenarios_acceptable: bool
    residuals: ResidualSummary
    scenario_ids: list[int]
    objective_by_scenario: dict[int, float | None]
    expected_objective: float | None = None
    message: str | None = None

    def to_jsonable(self) -> dict[str, Any]:
        """Return JSON-serializable iteration summary."""
        data: dict[str, Any] = {
            "iteration": self.iteration,
            "converged": self.converged,
            "all_scenarios_acceptable": self.all_scenarios_acceptable,
            "scenario_ids": list(self.scenario_ids),
            "objective_by_scenario": {
                str(scenario_id): objective
                for scenario_id, objective in self.objective_by_scenario.items()
            },
            "expected_objective": self.expected_objective,
            "residuals": self.residuals.to_jsonable(),
        }

        if self.message is not None:
            data["message"] = self.message

        return _json_sanitize(data)

    @classmethod
    def from_jsonable(cls, data: Mapping[str, Any]) -> "IterationSummary":
        """Construct from a JSON-loaded mapping."""
        required = {
            "iteration",
            "converged",
            "all_scenarios_acceptable",
            "scenario_ids",
            "objective_by_scenario",
            "expected_objective",
            "residuals",
        }
        missing = required - set(data)
        if missing:
            raise ValueError(
                "Iteration summary JSON is missing required field(s): "
                + ", ".join(sorted(missing))
            )

        residual_data = data["residuals"]
        residuals = ResidualSummary(
            max_abs_residual=float(residual_data["max_abs_residual"]),
            weighted_l1_residual=float(residual_data["weighted_l1_residual"]),
            weighted_l2_residual=float(residual_data["weighted_l2_residual"]),
            weighted_linf_residual=float(residual_data["weighted_linf_residual"]),
            num_nonanticipative_variables=int(
                residual_data["num_nonanticipative_variables"]
            ),
            num_scenarios=int(residual_data["num_scenarios"]),
        )

        return cls(
            iteration=int(data["iteration"]),
            converged=bool(data["converged"]),
            all_scenarios_acceptable=bool(data["all_scenarios_acceptable"]),
            residuals=residuals,
            scenario_ids=[int(sid) for sid in data["scenario_ids"]],
            objective_by_scenario={
                int(sid): (None if value is None else float(value))
                for sid, value in data["objective_by_scenario"].items()
            },
            expected_objective=(
                None
                if data["expected_objective"] is None
                else float(data["expected_objective"])
            ),
            message=data.get("message"),
        )


def write_scenario_result(
    output_path: str | Path,
    *,
    scenario_id: int,
    representative_period: int,
    representative_date: str,
    iteration: int,
    solve_outcome: SolveOutcome,
    nonant_var_map: Mapping[VariableID, Any],
    save_nonanticipative_metadata: bool = False,
    message: str | None = None,
) -> Path:
    """Write required per-scenario PH result JSON.

    Parameters
    ----------
    output_path:
        Destination JSON file.
    scenario_id:
        One-based PH scenario id.
    representative_period:
        One-based representative-period id.
    representative_date:
        Representative-period timestamp.
    iteration:
        PH iteration index.
    solve_outcome:
        Structured solver outcome.
    nonant_var_map:
        Mapping from ``VariableID`` to Pyomo variable data object.
    save_nonanticipative_metadata:
        Whether to include optional nonanticipative variable metadata.
    message:
        Optional diagnostic message.

    Returns
    -------
    pathlib.Path
        Path written.
    """
    nonanticipative_values = extract_nonanticipative_values(
        nonant_var_map,
        require_values=True,
    )

    metadata_json = None
    if save_nonanticipative_metadata:
        metadata_json = serialize_nonanticipative_metadata(
            extract_nonanticipative_metadata(nonant_var_map, include_values=True)
        )

    result = ScenarioResult(
        scenario_id=int(scenario_id),
        representative_period=int(representative_period),
        representative_date=str(representative_date),
        iteration=int(iteration),
        acceptable=bool(solve_outcome.acceptable),
        solver=solve_outcome.to_jsonable(),
        nonanticipative_values=nonanticipative_values,
        nonanticipative_metadata=metadata_json,
        message=message,
    )

    return write_scenario_result_object(output_path, result)


def write_failed_scenario_result(
    output_path: str | Path,
    *,
    scenario_id: int,
    representative_period: int,
    representative_date: str,
    iteration: int,
    solver_name: str,
    message: str,
    solver_status: str | None = None,
    termination_condition: str | None = None,
) -> Path:
    """Write a scenario-result JSON file for an unacceptable failed solve."""
    result = ScenarioResult(
        scenario_id=int(scenario_id),
        representative_period=int(representative_period),
        representative_date=str(representative_date),
        iteration=int(iteration),
        acceptable=False,
        solver={
            "name": solver_name,
            "status": solver_status,
            "termination_condition": termination_condition,
            "acceptable": False,
            "objective": None,
            "lower_bound": None,
            "upper_bound": None,
            "gap": None,
            "solve_time_sec": None,
            "log_file": None,
        },
        nonanticipative_values={},
        nonanticipative_metadata=None,
        message=message,
    )

    return write_scenario_result_object(output_path, result)


def write_scenario_result_object(
    output_path: str | Path,
    result: ScenarioResult,
) -> Path:
    """Write a ``ScenarioResult`` object to JSON."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        json.dump(result.to_jsonable(), f, indent=2, allow_nan=False)

    logger.info(
        "Wrote scenario result: iteration=%s scenario=%s path=%s",
        result.iteration,
        result.scenario_id,
        path,
    )
    return path


def read_scenario_result(path: str | Path) -> ScenarioResult:
    """Read one per-scenario PH result JSON file."""
    input_path = Path(path)

    with input_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    return ScenarioResult.from_jsonable(data)


def read_scenario_results(
    paths: Iterable[str | Path],
) -> dict[int, ScenarioResult]:
    """Read multiple scenario-result files indexed by scenario id."""
    results: dict[int, ScenarioResult] = {}

    for path in paths:
        result = read_scenario_result(path)

        if result.scenario_id in results:
            raise ValueError(
                f"Duplicate scenario result for scenario {result.scenario_id}."
            )

        results[result.scenario_id] = result

    return results


def scenario_values_from_results(
    scenario_results: Mapping[int, ScenarioResult],
    *,
    require_acceptable: bool = True,
) -> dict[int, dict[VariableID, float]]:
    """Extract nonanticipative values from parsed scenario results."""
    values: dict[int, dict[VariableID, float]] = {}

    for scenario_id, result in scenario_results.items():
        if require_acceptable and not result.acceptable:
            raise ValueError(
                f"Scenario {scenario_id} result is not acceptable for PH. "
                f"Message: {result.message}"
            )

        if not result.nonanticipative_values:
            raise ValueError(
                f"Scenario {scenario_id} result has no nonanticipative values."
            )

        values[int(scenario_id)] = dict(result.nonanticipative_values)

    return values


def write_iteration_summary(
    output_path: str | Path,
    summary: IterationSummary,
) -> Path:
    """Write an iteration summary JSON file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        json.dump(summary.to_jsonable(), f, indent=2, allow_nan=False)

    logger.info(
        "Wrote iteration summary for iteration %s to %s", summary.iteration, path
    )
    return path


def read_iteration_summary(path: str | Path) -> IterationSummary:
    """Read an iteration summary JSON file."""
    input_path = Path(path)

    with input_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    return IterationSummary.from_jsonable(data)


def append_convergence_history(
    *,
    summary: IterationSummary,
    history_json_path: str | Path | None = None,
    history_csv_path: str | Path | None = None,
) -> None:
    """Append one iteration summary to convergence history files."""
    if history_json_path is not None:
        append_convergence_history_json(history_json_path, summary)

    if history_csv_path is not None:
        append_convergence_history_csv(history_csv_path, summary)


def append_convergence_history_json(
    path: str | Path,
    summary: IterationSummary,
) -> Path:
    """Append one iteration summary to JSON convergence history."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        with output_path.open("r", encoding="utf-8") as f:
            history = json.load(f)
        if not isinstance(history, list):
            raise ValueError(f"Convergence history JSON is not a list: {output_path}")
    else:
        history = []

    history.append(summary.to_jsonable())

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(history, f, indent=2, allow_nan=False)

    return output_path


def append_convergence_history_csv(
    path: str | Path,
    summary: IterationSummary,
) -> Path:
    """Append one iteration summary to CSV convergence history."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    write_header = not output_path.exists()

    fieldnames = [
        "iteration",
        "converged",
        "all_scenarios_acceptable",
        "expected_objective",
        "max_abs_residual",
        "weighted_l1_residual",
        "weighted_l2_residual",
        "weighted_linf_residual",
        "num_nonanticipative_variables",
        "num_scenarios",
    ]

    row = {
        "iteration": summary.iteration,
        "converged": summary.converged,
        "all_scenarios_acceptable": summary.all_scenarios_acceptable,
        "expected_objective": summary.expected_objective,
        "max_abs_residual": summary.residuals.max_abs_residual,
        "weighted_l1_residual": summary.residuals.weighted_l1_residual,
        "weighted_l2_residual": summary.residuals.weighted_l2_residual,
        "weighted_linf_residual": summary.residuals.weighted_linf_residual,
        "num_nonanticipative_variables": (
            summary.residuals.num_nonanticipative_variables
        ),
        "num_scenarios": summary.residuals.num_scenarios,
    }

    with output_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)

        if write_header:
            writer.writeheader()

        writer.writerow(row)

    return output_path


def make_iteration_summary(
    *,
    iteration: int,
    converged: bool,
    residuals: ResidualSummary,
    scenario_results: Mapping[int, ScenarioResult],
    expected_objective: float | None = None,
    message: str | None = None,
) -> IterationSummary:
    """Build an ``IterationSummary`` from scenario results and residuals."""
    scenario_ids = sorted(int(sid) for sid in scenario_results)
    all_acceptable = all(result.acceptable for result in scenario_results.values())

    objective_by_scenario = {
        scenario_id: scenario_results[scenario_id].objective
        for scenario_id in scenario_ids
    }

    return IterationSummary(
        iteration=int(iteration),
        converged=bool(converged),
        all_scenarios_acceptable=bool(all_acceptable),
        residuals=residuals,
        scenario_ids=scenario_ids,
        objective_by_scenario=objective_by_scenario,
        expected_objective=expected_objective,
        message=message,
    )


def write_final_solution(
    final_dir: str | Path,
    *,
    final_state: PHState,
    scenario_results: Mapping[int, ScenarioResult],
    summary: IterationSummary,
    output_config: OutputConfig,
) -> dict[str, Path]:
    """Write final converged PH solution artifacts.

    Files written
    -------------
    ``nonanticipative_values.json``
        Final consensus values from ``final_state.xbar``.

    ``scenario_solutions.json``
        Final per-scenario result summaries.

    ``ph_summary.json``
        High-level PH convergence summary.

    Returns
    -------
    dict[str, pathlib.Path]
        Mapping of artifact label to written path.
    """
    directory = Path(final_dir)
    directory.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}

    paths["nonanticipative_values"] = _write_final_nonanticipative_values(
        directory / "nonanticipative_values.json",
        final_state=final_state,
        output_config=output_config,
    )

    paths["scenario_solutions"] = _write_final_scenario_solutions(
        directory / "scenario_solutions.json",
        scenario_results=scenario_results,
    )

    paths["ph_summary"] = _write_final_ph_summary(
        directory / "ph_summary.json",
        final_state=final_state,
        summary=summary,
    )

    return paths


def _write_final_nonanticipative_values(
    path: Path,
    *,
    final_state: PHState,
    output_config: OutputConfig,
) -> Path:
    """Write final consensus nonanticipative values."""
    records = []

    for variable_id, value in sorted(
        final_state.xbar.items(),
        key=lambda item: item[0].serialize(),
    ):
        record: dict[str, Any] = {
            "variable_id": variable_id.serialize(),
            "xbar": float(value),
            "final_value": float(value),
            "rounded": False,
            "rounding_error": 0.0,
        }

        rounded_value = _round_if_binary_like(
            value,
            tolerance=output_config.binary_rounding_tolerance,
        )
        if rounded_value is not None:
            record["final_value"] = rounded_value
            record["rounded"] = True
            record["rounding_error"] = abs(float(value) - rounded_value)

        records.append(record)

    data = {
        "iteration": final_state.iteration,
        "converged_state": True,
        "nonanticipative_values": records,
    }

    with path.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(data), f, indent=2, allow_nan=False)

    return path


def _write_final_scenario_solutions(
    path: Path,
    *,
    scenario_results: Mapping[int, ScenarioResult],
) -> Path:
    """Write final per-scenario solution summaries."""
    data = {
        "scenario_solutions": [
            {
                "scenario_id": result.scenario_id,
                "representative_period": result.representative_period,
                "representative_date": result.representative_date,
                "iteration": result.iteration,
                "acceptable": result.acceptable,
                "solver": result.solver,
                "nonanticipative_values": serialize_variable_value_records(
                    result.nonanticipative_values
                ),
            }
            for _, result in sorted(scenario_results.items())
        ]
    }

    with path.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(data), f, indent=2, allow_nan=False)

    return path


def _write_final_ph_summary(
    path: Path,
    *,
    final_state: PHState,
    summary: IterationSummary,
) -> Path:
    """Write high-level final PH summary."""
    data = {
        "converged": summary.converged,
        "final_iteration": summary.iteration,
        "state_iteration": final_state.iteration,
        "num_scenarios": len(final_state.scenario_ids),
        "num_nonanticipative_variables": len(final_state.variable_ids),
        "expected_objective": summary.expected_objective,
        "residuals": summary.residuals.to_jsonable(),
        "scenario_ids": list(final_state.scenario_ids),
    }

    with path.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(data), f, indent=2, allow_nan=False)

    return path


def _round_if_binary_like(value: float, *, tolerance: float) -> float | None:
    """Round values within tolerance of 0 or 1.

    This function does not infer variable type. It only provides a conservative
    final-output convenience for consensus values that are numerically very
    close to binary values.
    """
    value = float(value)

    if abs(value) <= tolerance:
        return 0.0

    if abs(value - 1.0) <= tolerance:
        return 1.0

    return None


def _json_sanitize(value: Any) -> Any:
    """Convert nested values into JSON-safe built-in types.

    Raises
    ------
    ValueError
        If a non-finite float is encountered.
    """
    if isinstance(value, dict):
        return {str(key): _json_sanitize(item) for key, item in value.items()}

    if isinstance(value, list):
        return [_json_sanitize(item) for item in value]

    if isinstance(value, tuple):
        return [_json_sanitize(item) for item in value]

    if isinstance(value, Path):
        return str(value)

    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError(f"Cannot serialize non-finite float value: {value}")
        return value

    return value

def export_durable_results(
    *,
    run_output_dir: str | Path,
    durable_output_dir: str | Path,
    include_iteration_summaries: bool = True,
    include_state: bool = True,
) -> dict[str, Path]:
    """Copy final PH results from scratch to durable storage.

    Parameters
    ----------
    run_output_dir:
        Scratch/run output directory.
    durable_output_dir:
        Persistent destination directory, e.g. under /projects or /home.
    include_iteration_summaries:
        Whether to copy iteration summary JSON files.
    include_state:
        Whether to copy PH state JSON files.

    Returns
    -------
    dict[str, pathlib.Path]
        Paths copied/created in durable storage.
    """
    source_root = Path(run_output_dir)
    dest_root = Path(durable_output_dir)
    dest_root.mkdir(parents=True, exist_ok=True)

    copied: dict[str, Path] = {}

    # Always copy effective config if present.
    config_source = source_root / "config_effective.yaml"
    if config_source.exists():
        config_dest = dest_root / "config_effective.yaml"
        shutil.copy2(config_source, config_dest)
        copied["config_effective"] = config_dest

    # Copy final directory.
    final_source = source_root / "final"
    if final_source.exists():
        final_dest = dest_root / "final"
        if final_dest.exists():
            shutil.rmtree(final_dest)
        shutil.copytree(final_source, final_dest)
        copied["final"] = final_dest

    # Copy convergence history.
    convergence_source = source_root / "convergence"
    if convergence_source.exists():
        convergence_dest = dest_root / "convergence"
        if convergence_dest.exists():
            shutil.rmtree(convergence_dest)
        shutil.copytree(convergence_source, convergence_dest)
        copied["convergence"] = convergence_dest

    # Copy iteration summaries only, not full scenario outputs/logs.
    if include_iteration_summaries:
        summaries_dest = dest_root / "iteration_summaries"
        if summaries_dest.exists():
            shutil.rmtree(summaries_dest)
        summaries_dest.mkdir(parents=True, exist_ok=True)

        iterations_source = source_root / "iterations"
        if iterations_source.exists():
            for summary_path in iterations_source.glob("iter_*/summary.json"):
                iter_dest = summaries_dest / summary_path.parent.name
                iter_dest.mkdir(parents=True, exist_ok=True)
                shutil.copy2(summary_path, iter_dest / "summary.json")

        copied["iteration_summaries"] = summaries_dest

    # Copy state files if requested. These are usually small and useful for
    # restart/debugging.
    if include_state:
        state_source = source_root / "state"
        if state_source.exists():
            state_dest = dest_root / "state"
            if state_dest.exists():
                shutil.rmtree(state_dest)
            shutil.copytree(state_source, state_dest)
            copied["state"] = state_dest

    return copied