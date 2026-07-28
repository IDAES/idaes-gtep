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
"""State management for GTEP Progressive Hedging.

This module owns the machine-written PH iteration state. User-facing
configuration remains YAML; PH state is stored as JSON.

Design constraints
------------------
* Nonanticipative variable identity is represented by ``VariableID`` objects
  backed by Pyomo ``ComponentUID``.
* JSON files serialize variable IDs as records containing serialized
  ``ComponentUID`` values.
* JSON files do not use Pyomo component names as matching keys.
* Variable matching after reloading must be performed through
  ``ComponentUID.find_component_on(model)`` via ``VariableID.find_on(model)``.

PH initialization convention
----------------------------
Iteration 0 is intended to be an unregularized scenario solve:

* multipliers are zero,
* no PH regularization or multiplier terms are required,
* ``ph_terms_enabled`` is false.

After iteration 0 scenario solves complete, the orchestrator computes the first
consensus solution and updated multipliers, writes iteration 1 state, and PH
terms are enabled for subsequent scenario solves.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import json
import logging

from gtep.algorithms.progressive_hedging.nonanticipativity import (
    VariableID,
    deserialize_variable_value_records,
    serialize_variable_value_records,
)

logger = logging.getLogger("gtep.algorithms.progressive_hedging.state")


@dataclass
class PHState:
    """Progressive Hedging iteration state.

    Attributes
    ----------
    iteration:
        PH iteration index associated with this state.
    scenario_ids:
        One-based representative-period scenario identifiers.
    scenario_probabilities:
        Normalized scenario probabilities indexed by scenario id.
    scenario_weights:
        Original representative-period weights indexed by scenario id.
    variable_ids:
        Ordered nonanticipative variable IDs. This may be empty for the initial
        unregularized iteration before variable IDs have been discovered from
        scenario solutions.
    xbar:
        Consensus values indexed by ``VariableID``.
    multipliers:
        PH multipliers indexed first by scenario id and then by ``VariableID``.
    rho:
        Either scalar rho or variable-specific rho indexed by ``VariableID``.
    ph_terms_enabled:
        If false, scenario solves should use only the scenario base objective.
        If true, scenario solves should include PH multiplier and
        regularization terms.
    metadata:
        Optional JSON-serializable metadata.
    """

    iteration: int
    scenario_ids: list[int]
    scenario_probabilities: dict[int, float]
    scenario_weights: dict[int, float]

    variable_ids: list[VariableID] = field(default_factory=list)
    xbar: dict[VariableID, float] = field(default_factory=dict)
    multipliers: dict[int, dict[VariableID, float]] = field(default_factory=dict)

    rho: float | dict[VariableID, float] = 1.0
    ph_terms_enabled: bool = False

    metadata: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        """Validate internal state consistency."""
        if self.iteration < 0:
            raise ValueError(
                f"PH iteration must be nonnegative. Received {self.iteration}."
            )

        if not self.scenario_ids:
            raise ValueError("PH state must contain at least one scenario id.")

        if len(set(self.scenario_ids)) != len(self.scenario_ids):
            raise ValueError(f"Duplicate scenario ids in PH state: {self.scenario_ids}")

        missing_probabilities = set(self.scenario_ids) - set(
            self.scenario_probabilities
        )
        if missing_probabilities:
            raise ValueError(
                "Missing scenario probabilities for scenario id(s): "
                + ", ".join(str(sid) for sid in sorted(missing_probabilities))
            )

        missing_weights = set(self.scenario_ids) - set(self.scenario_weights)
        if missing_weights:
            raise ValueError(
                "Missing scenario weights for scenario id(s): "
                + ", ".join(str(sid) for sid in sorted(missing_weights))
            )

        probability_sum = sum(
            self.scenario_probabilities[sid] for sid in self.scenario_ids
        )
        if probability_sum <= 0:
            raise ValueError(
                "Scenario probabilities must have positive sum. "
                f"Received sum={probability_sum}."
            )

        for sid in self.scenario_ids:
            probability = self.scenario_probabilities[sid]
            weight = self.scenario_weights[sid]

            if probability < 0:
                raise ValueError(
                    f"Scenario probability for scenario {sid} is negative: {probability}."
                )

            if weight < 0:
                raise ValueError(
                    f"Scenario weight for scenario {sid} is negative: {weight}."
                )

        variable_id_set = set(self.variable_ids)
        if len(variable_id_set) != len(self.variable_ids):
            raise ValueError(
                "PH state contains duplicate nonanticipative variable IDs."
            )

        if self.xbar:
            missing_xbar_ids = variable_id_set - set(self.xbar)
            extra_xbar_ids = set(self.xbar) - variable_id_set

            if missing_xbar_ids:
                raise ValueError(
                    "PH state variable_ids contains IDs missing from xbar: "
                    + _format_variable_ids(missing_xbar_ids)
                )

            if extra_xbar_ids:
                raise ValueError(
                    "PH state xbar contains IDs not listed in variable_ids: "
                    + _format_variable_ids(extra_xbar_ids)
                )

        for sid, scenario_multipliers in self.multipliers.items():
            if sid not in self.scenario_ids:
                raise ValueError(
                    f"PH state multipliers contain unknown scenario id {sid}."
                )

            extra_multiplier_ids = set(scenario_multipliers) - variable_id_set
            if extra_multiplier_ids:
                raise ValueError(
                    f"PH state multipliers for scenario {sid} contain IDs not "
                    "listed in variable_ids: "
                    + _format_variable_ids(extra_multiplier_ids)
                )

            if self.ph_terms_enabled:
                missing_multiplier_ids = variable_id_set - set(scenario_multipliers)
                if missing_multiplier_ids:
                    raise ValueError(
                        f"PH terms are enabled, but multipliers for scenario {sid} "
                        "are missing variable IDs: "
                        + _format_variable_ids(missing_multiplier_ids)
                    )

        if self.ph_terms_enabled:
            missing_multiplier_scenarios = set(self.scenario_ids) - set(
                self.multipliers
            )
            if missing_multiplier_scenarios:
                raise ValueError(
                    "PH terms are enabled, but multipliers are missing for "
                    "scenario id(s): "
                    + ", ".join(
                        str(sid) for sid in sorted(missing_multiplier_scenarios)
                    )
                )

            if not self.variable_ids:
                raise ValueError("PH terms are enabled but variable_ids is empty.")

            if not self.xbar:
                raise ValueError("PH terms are enabled but xbar is empty.")

        if isinstance(self.rho, Mapping):
            rho_ids = set(self.rho)
            if rho_ids != variable_id_set:
                missing_rho = variable_id_set - rho_ids
                extra_rho = rho_ids - variable_id_set

                if missing_rho:
                    raise ValueError(
                        "Variable-specific rho is missing variable IDs: "
                        + _format_variable_ids(missing_rho)
                    )

                if extra_rho:
                    raise ValueError(
                        "Variable-specific rho contains IDs not listed in variable_ids: "
                        + _format_variable_ids(extra_rho)
                    )

            for value in self.rho.values():
                if value < 0:
                    raise ValueError(
                        "All variable-specific rho values must be nonnegative."
                    )
        else:
            if self.rho < 0:
                raise ValueError(
                    f"Scalar rho must be nonnegative. Received {self.rho}."
                )

    def get_scenario_multipliers(self, scenario_id: int) -> dict[VariableID, float]:
        """Return multipliers for one scenario, defaulting to zero if disabled."""
        if scenario_id not in self.scenario_ids:
            raise KeyError(f"Unknown scenario id {scenario_id}.")

        if scenario_id in self.multipliers:
            return dict(self.multipliers[scenario_id])

        if self.ph_terms_enabled:
            raise KeyError(
                f"PH terms are enabled but no multipliers exist for scenario {scenario_id}."
            )

        return {variable_id: 0.0 for variable_id in self.variable_ids}

    def get_rho(self) -> float | dict[VariableID, float]:
        """Return scalar or variable-specific rho."""
        if isinstance(self.rho, Mapping):
            return dict(self.rho)
        return float(self.rho)

    def to_jsonable(self) -> dict[str, Any]:
        """Convert state to a JSON-serializable dictionary."""
        self.validate()

        if isinstance(self.rho, Mapping):
            rho_json: dict[str, Any] = {
                "type": "variable_specific",
                "values": serialize_variable_value_records(self.rho),
            }
        else:
            rho_json = {
                "type": "scalar",
                "value": float(self.rho),
            }

        return {
            "iteration": self.iteration,
            "scenario_ids": list(self.scenario_ids),
            "scenario_probabilities": {
                str(sid): float(self.scenario_probabilities[sid])
                for sid in self.scenario_ids
            },
            "scenario_weights": {
                str(sid): float(self.scenario_weights[sid]) for sid in self.scenario_ids
            },
            "variable_ids": [
                {"variable_id": variable_id.serialize()}
                for variable_id in self.variable_ids
            ],
            "xbar": serialize_variable_value_records(self.xbar),
            "multipliers": {
                str(sid): serialize_variable_value_records(
                    self.multipliers.get(sid, {})
                )
                for sid in self.scenario_ids
            },
            "rho": rho_json,
            "ph_terms_enabled": bool(self.ph_terms_enabled),
            "metadata": self.metadata,
        }

    @classmethod
    def from_jsonable(cls, data: Mapping[str, Any]) -> "PHState":
        """Construct state from a JSON-loaded mapping."""
        required = {
            "iteration",
            "scenario_ids",
            "scenario_probabilities",
            "scenario_weights",
            "variable_ids",
            "xbar",
            "multipliers",
            "rho",
            "ph_terms_enabled",
        }
        missing = required - set(data)
        if missing:
            raise ValueError(
                "PH state JSON is missing required field(s): "
                + ", ".join(sorted(missing))
            )

        scenario_ids = [int(sid) for sid in data["scenario_ids"]]

        scenario_probabilities = {
            int(sid): float(value)
            for sid, value in data["scenario_probabilities"].items()
        }
        scenario_weights = {
            int(sid): float(value) for sid, value in data["scenario_weights"].items()
        }

        variable_ids = [
            VariableID.from_serialized(record["variable_id"])
            for record in data["variable_ids"]
        ]

        xbar = deserialize_variable_value_records(data["xbar"])

        multipliers = {
            int(sid): deserialize_variable_value_records(records)
            for sid, records in data["multipliers"].items()
        }

        rho_data = data["rho"]
        rho_type = rho_data.get("type")
        if rho_type == "scalar":
            rho: float | dict[VariableID, float] = float(rho_data["value"])
        elif rho_type == "variable_specific":
            rho = deserialize_variable_value_records(rho_data["values"])
        else:
            raise ValueError(
                "PH state rho field must have type 'scalar' or "
                f"'variable_specific'. Received {rho_type!r}."
            )

        state = cls(
            iteration=int(data["iteration"]),
            scenario_ids=scenario_ids,
            scenario_probabilities=scenario_probabilities,
            scenario_weights=scenario_weights,
            variable_ids=variable_ids,
            xbar=xbar,
            multipliers=multipliers,
            rho=rho,
            ph_terms_enabled=bool(data["ph_terms_enabled"]),
            metadata=dict(data.get("metadata", {})),
        )

        state.validate()
        return state


def create_initial_state(
    *,
    scenario_ids: Iterable[int],
    scenario_probabilities: Mapping[int, float],
    scenario_weights: Mapping[int, float],
    rho: float,
    metadata: Mapping[str, Any] | None = None,
) -> PHState:
    """Create initial unregularized PH state for iteration 0.

    Iteration 0 scenario solves should use the scenario base objective only.
    Nonanticipative variable IDs are discovered from scenario-solve outputs and
    introduced into the next PH state.
    """
    state = PHState(
        iteration=0,
        scenario_ids=list(scenario_ids),
        scenario_probabilities={
            int(sid): float(prob) for sid, prob in scenario_probabilities.items()
        },
        scenario_weights={
            int(sid): float(weight) for sid, weight in scenario_weights.items()
        },
        variable_ids=[],
        xbar={},
        multipliers={},
        rho=float(rho),
        ph_terms_enabled=False,
        metadata=dict(metadata or {}),
    )
    state.validate()
    return state


def create_enabled_state(
    *,
    iteration: int,
    scenario_ids: Iterable[int],
    scenario_probabilities: Mapping[int, float],
    scenario_weights: Mapping[int, float],
    variable_ids: Iterable[VariableID],
    xbar: Mapping[VariableID, float],
    multipliers: Mapping[int, Mapping[VariableID, float]],
    rho: float | Mapping[VariableID, float],
    metadata: Mapping[str, Any] | None = None,
) -> PHState:
    """Create PH state with multiplier and regularization terms enabled."""
    state = PHState(
        iteration=int(iteration),
        scenario_ids=list(scenario_ids),
        scenario_probabilities={
            int(sid): float(prob) for sid, prob in scenario_probabilities.items()
        },
        scenario_weights={
            int(sid): float(weight) for sid, weight in scenario_weights.items()
        },
        variable_ids=list(variable_ids),
        xbar={variable_id: float(value) for variable_id, value in xbar.items()},
        multipliers={
            int(sid): {
                variable_id: float(value)
                for variable_id, value in scenario_values.items()
            }
            for sid, scenario_values in multipliers.items()
        },
        rho=(
            {variable_id: float(value) for variable_id, value in rho.items()}
            if isinstance(rho, Mapping)
            else float(rho)
        ),
        ph_terms_enabled=True,
        metadata=dict(metadata or {}),
    )
    state.validate()
    return state


def write_state(state: PHState, path: str | Path) -> Path:
    """Write PH state to JSON."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(state.to_jsonable(), f, indent=2)

    logger.info("Wrote PH state for iteration %s to %s", state.iteration, output_path)
    return output_path


def read_state(path: str | Path) -> PHState:
    """Read PH state from JSON."""
    input_path = Path(path)

    with input_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    state = PHState.from_jsonable(data)
    logger.info("Read PH state for iteration %s from %s", state.iteration, input_path)
    return state


def initialize_zero_multipliers(
    scenario_ids: Iterable[int],
    variable_ids: Iterable[VariableID],
) -> dict[int, dict[VariableID, float]]:
    """Create zero multipliers for all scenarios and variable IDs."""
    variable_id_list = list(variable_ids)
    return {
        int(sid): {variable_id: 0.0 for variable_id in variable_id_list}
        for sid in scenario_ids
    }


def initialize_zero_xbar(
    variable_ids: Iterable[VariableID],
) -> dict[VariableID, float]:
    """Create zero consensus values for all variable IDs."""
    return {variable_id: 0.0 for variable_id in variable_ids}


def scalar_rho_to_variable_specific(
    variable_ids: Iterable[VariableID],
    rho: float,
) -> dict[VariableID, float]:
    """Expand scalar rho to variable-specific rho values."""
    if rho < 0:
        raise ValueError(f"rho must be nonnegative. Received {rho}.")
    return {variable_id: float(rho) for variable_id in variable_ids}


def get_rho_value(
    rho: float | Mapping[VariableID, float],
    variable_id: VariableID,
) -> float:
    """Return scalar or variable-specific rho for one variable."""
    if isinstance(rho, Mapping):
        return float(rho[variable_id])
    return float(rho)


def copy_state_with_new_iteration(
    state: PHState,
    *,
    iteration: int,
    xbar: Mapping[VariableID, float],
    multipliers: Mapping[int, Mapping[VariableID, float]],
    ph_terms_enabled: bool = True,
    metadata_update: Mapping[str, Any] | None = None,
) -> PHState:
    """Create a new state from an existing state with updated PH quantities."""
    metadata = dict(state.metadata)
    if metadata_update:
        metadata.update(metadata_update)

    new_state = PHState(
        iteration=int(iteration),
        scenario_ids=list(state.scenario_ids),
        scenario_probabilities=dict(state.scenario_probabilities),
        scenario_weights=dict(state.scenario_weights),
        variable_ids=list(state.variable_ids),
        xbar={variable_id: float(value) for variable_id, value in xbar.items()},
        multipliers={
            int(sid): {
                variable_id: float(value)
                for variable_id, value in scenario_values.items()
            }
            for sid, scenario_values in multipliers.items()
        },
        rho=(
            {variable_id: float(value) for variable_id, value in state.rho.items()}
            if isinstance(state.rho, Mapping)
            else float(state.rho)
        ),
        ph_terms_enabled=bool(ph_terms_enabled),
        metadata=metadata,
    )
    new_state.validate()
    return new_state


def _format_variable_ids(variable_ids: Iterable[VariableID]) -> str:
    """Format variable IDs for error messages."""
    return ", ".join(
        variable_id.serialize()
        for variable_id in sorted(variable_ids, key=lambda item: item.serialize())
    )
