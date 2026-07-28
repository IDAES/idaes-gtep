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
"""Convergence and PH update utilities for GTEP Progressive Hedging.

This module computes:

* scenario-probability-weighted consensus values,
* nonanticipativity residuals,
* PH multiplier updates,
* convergence summaries,
* next-iteration ``PHState`` objects.

Design constraints
------------------
Nonanticipative variables are represented only by ``VariableID`` objects backed
by Pyomo ``ComponentUID``. This module never matches variables by names, string
paths, or manually parsed strings.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import math

from gtep.algorithms.progressive_hedging.nonanticipativity import (
    VariableID,
    assert_same_variable_id_set,
    serialize_variable_value_records,
)
from gtep.algorithms.progressive_hedging.state import (
    PHState,
    create_enabled_state,
    get_rho_value,
    initialize_zero_multipliers,
)


@dataclass(frozen=True)
class ResidualSummary:
    """Summary of nonanticipativity residuals for one PH iteration."""

    max_abs_residual: float
    weighted_l1_residual: float
    weighted_l2_residual: float
    weighted_linf_residual: float
    num_nonanticipative_variables: int
    num_scenarios: int

    def converged(self, tolerance: float) -> bool:
        """Return true if the max absolute residual satisfies tolerance."""
        return self.max_abs_residual <= tolerance

    def to_jsonable(self) -> dict[str, Any]:
        """Return a JSON-serializable residual summary."""
        return {
            "max_abs_residual": self.max_abs_residual,
            "weighted_l1_residual": self.weighted_l1_residual,
            "weighted_l2_residual": self.weighted_l2_residual,
            "weighted_linf_residual": self.weighted_linf_residual,
            "num_nonanticipative_variables": self.num_nonanticipative_variables,
            "num_scenarios": self.num_scenarios,
        }


@dataclass(frozen=True)
class PHUpdateResult:
    """Result of processing one completed PH iteration."""

    iteration_completed: int
    next_iteration: int
    xbar: dict[VariableID, float]
    multipliers: dict[int, dict[VariableID, float]]
    residuals: ResidualSummary
    converged: bool
    next_state: PHState

    def to_jsonable(self) -> dict[str, Any]:
        """Return a JSON-serializable update summary."""
        return {
            "iteration_completed": self.iteration_completed,
            "next_iteration": self.next_iteration,
            "converged": self.converged,
            "residuals": self.residuals.to_jsonable(),
            "xbar": serialize_variable_value_records(self.xbar),
        }


def compute_consensus(
    scenario_values: Mapping[int, Mapping[VariableID, float]],
    scenario_probabilities: Mapping[int, float],
    scenario_ids: Iterable[int] | None = None,
) -> dict[VariableID, float]:
    """Compute scenario-probability-weighted consensus values.

    Parameters
    ----------
    scenario_values:
        Mapping from scenario id to nonanticipative variable values.
    scenario_probabilities:
        Mapping from scenario id to normalized scenario probability.
    scenario_ids:
        Optional explicit scenario id ordering/subset. If omitted, keys from
        ``scenario_values`` are used.

    Returns
    -------
    dict[VariableID, float]
        Consensus values \( \bar{x}_i = \sum_s \pi_s x_{s,i} \).
    """
    ordered_scenario_ids = _ordered_scenario_ids(scenario_values, scenario_ids)
    _validate_scenario_probabilities(ordered_scenario_ids, scenario_probabilities)
    variable_ids = _validate_scenario_value_sets(scenario_values, ordered_scenario_ids)

    consensus = {variable_id: 0.0 for variable_id in variable_ids}

    for scenario_id in ordered_scenario_ids:
        probability = float(scenario_probabilities[scenario_id])
        for variable_id in variable_ids:
            consensus[variable_id] += probability * float(
                scenario_values[scenario_id][variable_id]
            )

    return consensus


def compute_residuals(
    scenario_values: Mapping[int, Mapping[VariableID, float]],
    xbar: Mapping[VariableID, float],
    scenario_probabilities: Mapping[int, float],
    scenario_ids: Iterable[int] | None = None,
) -> ResidualSummary:
    """Compute nonanticipativity residuals for scenario variable values.

    Residual for scenario \(s\), variable \(i\):

    \[
        r_{s,i} = x_{s,i} - \bar{x}_i.
    \]

    The primary convergence criterion is ``max_abs_residual``.
    """
    ordered_scenario_ids = _ordered_scenario_ids(scenario_values, scenario_ids)
    _validate_scenario_probabilities(ordered_scenario_ids, scenario_probabilities)
    variable_ids = _validate_scenario_value_sets(scenario_values, ordered_scenario_ids)

    assert_same_variable_id_set(
        variable_ids,
        xbar.keys(),
        first_label="scenario_values",
        second_label="xbar",
    )

    max_abs_residual = 0.0
    weighted_l1_residual = 0.0
    weighted_l2_residual_squared = 0.0
    weighted_linf_residual = 0.0

    for scenario_id in ordered_scenario_ids:
        probability = float(scenario_probabilities[scenario_id])
        scenario_linf = 0.0

        for variable_id in variable_ids:
            residual = float(scenario_values[scenario_id][variable_id]) - float(
                xbar[variable_id]
            )
            abs_residual = abs(residual)

            max_abs_residual = max(max_abs_residual, abs_residual)
            weighted_l1_residual += probability * abs_residual
            weighted_l2_residual_squared += probability * residual * residual
            scenario_linf = max(scenario_linf, abs_residual)

        weighted_linf_residual += probability * scenario_linf

    return ResidualSummary(
        max_abs_residual=max_abs_residual,
        weighted_l1_residual=weighted_l1_residual,
        weighted_l2_residual=math.sqrt(weighted_l2_residual_squared),
        weighted_linf_residual=weighted_linf_residual,
        num_nonanticipative_variables=len(variable_ids),
        num_scenarios=len(ordered_scenario_ids),
    )


def update_multipliers(
    previous_state: PHState,
    scenario_values: Mapping[int, Mapping[VariableID, float]],
    xbar: Mapping[VariableID, float],
) -> dict[int, dict[VariableID, float]]:
    """Update PH multipliers.

    For each scenario \(s\) and nonanticipative variable \(i\):

    \[
        w_{s,i}^{k+1}
        =
        w_{s,i}^{k}
        +
        \rho_i (x_{s,i}^{k} - \bar{x}_i^{k}).
    \]

    If the previous state has PH terms disabled, previous multipliers are
    treated as zero.
    """
    scenario_ids = list(previous_state.scenario_ids)
    variable_ids = _validate_scenario_value_sets(scenario_values, scenario_ids)

    assert_same_variable_id_set(
        variable_ids,
        xbar.keys(),
        first_label="scenario_values",
        second_label="xbar",
    )

    if previous_state.ph_terms_enabled:
        previous_multipliers = previous_state.multipliers
    else:
        previous_multipliers = initialize_zero_multipliers(scenario_ids, variable_ids)

    updated: dict[int, dict[VariableID, float]] = {}

    for scenario_id in scenario_ids:
        scenario_updated: dict[VariableID, float] = {}
        for variable_id in variable_ids:
            previous_w = float(previous_multipliers[scenario_id][variable_id])
            rho_value = get_rho_value(previous_state.rho, variable_id)
            deviation = float(scenario_values[scenario_id][variable_id]) - float(
                xbar[variable_id]
            )
            scenario_updated[variable_id] = previous_w + rho_value * deviation

        updated[scenario_id] = scenario_updated

    return updated


def process_completed_iteration(
    previous_state: PHState,
    scenario_values: Mapping[int, Mapping[VariableID, float]],
    *,
    convergence_tolerance: float,
    metadata: Mapping[str, Any] | None = None,
) -> PHUpdateResult:
    """Process scenario solves from one PH iteration.

    Parameters
    ----------
    previous_state:
        State used to solve the completed iteration.
    scenario_values:
        Scenario nonanticipative values from the completed iteration.
    convergence_tolerance:
        Tolerance for max absolute nonanticipativity residual.
    metadata:
        Optional metadata to attach to the next state.

    Returns
    -------
    PHUpdateResult
        Consensus, residuals, updated multipliers, convergence flag, and next
        PH state.
    """
    previous_state.validate()

    if convergence_tolerance < 0:
        raise ValueError(
            "convergence_tolerance must be nonnegative. "
            f"Received {convergence_tolerance}."
        )

    scenario_ids = list(previous_state.scenario_ids)
    variable_ids = _validate_scenario_value_sets(scenario_values, scenario_ids)

    xbar = compute_consensus(
        scenario_values,
        previous_state.scenario_probabilities,
        scenario_ids,
    )

    residuals = compute_residuals(
        scenario_values,
        xbar,
        previous_state.scenario_probabilities,
        scenario_ids,
    )

    converged = residuals.converged(convergence_tolerance)

    multipliers = update_multipliers(
        previous_state,
        scenario_values,
        xbar,
    )

    next_iteration = previous_state.iteration + 1

    next_metadata = dict(previous_state.metadata)
    if metadata:
        next_metadata.update(metadata)
    next_metadata.update(
        {
            "previous_iteration": previous_state.iteration,
            "previous_iteration_converged": converged,
            "previous_iteration_residuals": residuals.to_jsonable(),
        }
    )

    next_state = create_enabled_state(
        iteration=next_iteration,
        scenario_ids=scenario_ids,
        scenario_probabilities=previous_state.scenario_probabilities,
        scenario_weights=previous_state.scenario_weights,
        variable_ids=variable_ids,
        xbar=xbar,
        multipliers=multipliers,
        rho=previous_state.rho,
        metadata=next_metadata,
    )

    return PHUpdateResult(
        iteration_completed=previous_state.iteration,
        next_iteration=next_iteration,
        xbar=xbar,
        multipliers=multipliers,
        residuals=residuals,
        converged=converged,
        next_state=next_state,
    )


def compute_weighted_expected_objective(
    scenario_objectives: Mapping[int, float],
    scenario_probabilities: Mapping[int, float],
    scenario_ids: Iterable[int] | None = None,
) -> float:
    """Compute a probability-weighted expected objective from scenario values."""
    if scenario_ids is None:
        ordered_scenario_ids = sorted(int(sid) for sid in scenario_objectives)
    else:
        ordered_scenario_ids = [int(sid) for sid in scenario_ids]

    missing_objectives = set(ordered_scenario_ids) - set(scenario_objectives)
    if missing_objectives:
        raise ValueError(
            "Missing scenario objective value(s) for scenario id(s): "
            + ", ".join(str(sid) for sid in sorted(missing_objectives))
        )

    _validate_scenario_probabilities(ordered_scenario_ids, scenario_probabilities)

    return sum(
        float(scenario_probabilities[scenario_id])
        * float(scenario_objectives[scenario_id])
        for scenario_id in ordered_scenario_ids
    )


def scenario_values_from_records_by_scenario(
    records_by_scenario: Mapping[int, Iterable[Mapping[str, Any]]],
) -> dict[int, dict[VariableID, float]]:
    """Deserialize scenario value records into PH calculation structure.

    This helper expects records in the format produced by
    ``serialize_variable_value_records``.
    """
    from gtep.algorithms.progressive_hedging.nonanticipativity import (
        deserialize_variable_value_records,
    )

    return {
        int(scenario_id): deserialize_variable_value_records(records)
        for scenario_id, records in records_by_scenario.items()
    }


def _ordered_scenario_ids(
    scenario_values: Mapping[int, Mapping[VariableID, float]],
    scenario_ids: Iterable[int] | None,
) -> list[int]:
    """Return explicit or inferred ordered scenario ids."""
    if scenario_ids is None:
        ordered = sorted(int(sid) for sid in scenario_values)
    else:
        ordered = [int(sid) for sid in scenario_ids]

    if not ordered:
        raise ValueError("At least one scenario id is required.")

    missing = set(ordered) - set(scenario_values)
    if missing:
        raise ValueError(
            "Missing scenario values for scenario id(s): "
            + ", ".join(str(sid) for sid in sorted(missing))
        )

    extra = set(scenario_values) - set(ordered)
    if extra:
        raise ValueError(
            "Scenario values contain scenario id(s) not requested: "
            + ", ".join(str(sid) for sid in sorted(extra))
        )

    return ordered


def _validate_scenario_probabilities(
    scenario_ids: Iterable[int],
    scenario_probabilities: Mapping[int, float],
) -> None:
    """Validate probabilities for the requested scenario ids."""
    scenario_id_list = list(scenario_ids)

    missing = set(scenario_id_list) - set(scenario_probabilities)
    if missing:
        raise ValueError(
            "Missing scenario probabilities for scenario id(s): "
            + ", ".join(str(sid) for sid in sorted(missing))
        )

    probability_sum = 0.0
    for scenario_id in scenario_id_list:
        probability = float(scenario_probabilities[scenario_id])
        if probability < 0:
            raise ValueError(
                f"Scenario probability for scenario {scenario_id} is negative: "
                f"{probability}."
            )
        probability_sum += probability

    if probability_sum <= 0:
        raise ValueError(
            "Scenario probabilities must have positive sum. "
            f"Received sum={probability_sum}."
        )

    if abs(probability_sum - 1.0) > 1.0e-8:
        # This is not fatal because the configuration supports non-normalized
        # probabilities, but the intended PH configuration uses normalized
        # representative-period weights. The calculations use the values exactly
        # as provided.
        pass


def _validate_scenario_value_sets(
    scenario_values: Mapping[int, Mapping[VariableID, float]],
    scenario_ids: Iterable[int],
) -> list[VariableID]:
    """Validate that all scenarios contain the same variable ID set."""
    scenario_id_list = list(scenario_ids)

    if not scenario_id_list:
        raise ValueError("At least one scenario id is required.")

    for scenario_id in scenario_id_list:
        if scenario_id not in scenario_values:
            raise ValueError(f"Missing scenario values for scenario {scenario_id}.")

    reference_scenario_id = scenario_id_list[0]
    reference_ids = list(scenario_values[reference_scenario_id].keys())

    if not reference_ids:
        raise ValueError(
            f"Scenario {reference_scenario_id} has no nonanticipative values."
        )

    if len(set(reference_ids)) != len(reference_ids):
        raise ValueError(
            f"Scenario {reference_scenario_id} contains duplicate variable IDs."
        )

    for scenario_id in scenario_id_list[1:]:
        assert_same_variable_id_set(
            reference_ids,
            scenario_values[scenario_id].keys(),
            first_label=f"scenario {reference_scenario_id}",
            second_label=f"scenario {scenario_id}",
        )

    for variable_id in reference_ids:
        if not isinstance(variable_id, VariableID):
            raise TypeError(
                "Scenario value mappings must be keyed by VariableID objects. "
                f"Received {type(variable_id)}."
            )

    return reference_ids
