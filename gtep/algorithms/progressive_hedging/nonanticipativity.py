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
"""Nonanticipative-variable utilities for GTEP Progressive Hedging.

Design constraint
-----------------
Nonanticipative variables must never be matched by Pyomo component names,
string paths, or manual string parsing. Durable variable identity is represented
with Pyomo ``ComponentUID`` objects.

Because PH scenario solves run in separate Python processes, variable IDs must
be serialized for JSON transport. The serialized value is treated only as the
serialized representation of a ``ComponentUID``. It is not manually parsed or
used for name-based component lookup. Rebinding to a rebuilt Pyomo model is
always performed through ``ComponentUID.find_component_on(model)``.

This module collects only true investment-state/status/capacity variables:

* renewable investment variables,
* thermal investment-status associated binaries,
* transmission investment-status associated binaries,
* storage investment-status associated binaries.

It intentionally does not collect investment-stage accounting or lower-level
summary variables such as expansion cost, operating cost, curtailment cost, or
quota-deficit variables.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import logging

import pyomo.environ as pyo
from pyomo.common.collections import ComponentMap
from pyomo.core.base.componentuid import ComponentUID
from pyomo.core.base.enums import SortComponents

from gtep.algorithms.progressive_hedging.config import NonanticipativityConfig

logger = logging.getLogger("gtep.algorithms.progressive_hedging.nonanticipativity")


@dataclass(frozen=True)
class VariableID:
    """Durable identifier for a Pyomo variable based on ``ComponentUID``.

    Parameters
    ----------
    serialized_component_uid:
        Serialized representation of a Pyomo ``ComponentUID``.

    Notes
    -----
    The serialized value exists for transport and hashing across processes.
    Matching to a rebuilt Pyomo model must be done only by reconstructing the
    ``ComponentUID`` and calling ``find_component_on(model)``.
    """

    serialized_component_uid: str

    @classmethod
    def from_component(cls, component: Any) -> "VariableID":
        """Create a variable ID from a Pyomo component data object."""
        return cls(str(ComponentUID(component)))

    @classmethod
    def from_serialized(cls, value: str) -> "VariableID":
        """Create a variable ID from a serialized ``ComponentUID`` value."""
        if not isinstance(value, str):
            raise TypeError(
                "Serialized ComponentUID value must be a string. "
                f"Received {type(value)}."
            )
        return cls(value)

    @property
    def component_uid(self) -> ComponentUID:
        """Return this ID as a Pyomo ``ComponentUID`` object."""
        return ComponentUID(self.serialized_component_uid)

    def find_on(self, model: pyo.ConcreteModel) -> Any:
        """Find the referenced component on a Pyomo model.

        Parameters
        ----------
        model:
            Rebuilt Pyomo model on which to locate the variable.

        Returns
        -------
        Any
            Pyomo component data object found by ``ComponentUID``.

        Raises
        ------
        KeyError
            If the component cannot be found on the model.
        """
        component = self.component_uid.find_component_on(model)
        if component is None:
            raise KeyError(
                "Could not find nonanticipative variable on model using "
                f"ComponentUID {self.serialized_component_uid!r}."
            )
        return component

    def serialize(self) -> str:
        """Return JSON-serializable representation."""
        return self.serialized_component_uid


@dataclass(frozen=True)
class NonanticipativeVariableMetadata:
    """Lightweight metadata for one nonanticipative variable."""

    variable_id: VariableID
    is_binary: bool
    is_integer: bool
    is_continuous: bool
    fixed: bool
    lower_bound: float | None
    upper_bound: float | None
    value: float | None = None

    def to_jsonable(self) -> dict[str, Any]:
        """Convert metadata to a JSON-serializable dictionary."""
        return {
            "variable_id": self.variable_id.serialize(),
            "is_binary": self.is_binary,
            "is_integer": self.is_integer,
            "is_continuous": self.is_continuous,
            "fixed": self.fixed,
            "lower_bound": self.lower_bound,
            "upper_bound": self.upper_bound,
            "value": self.value,
        }


def collect_nonanticipative_variables(
    model: pyo.ConcreteModel,
    config: NonanticipativityConfig | None = None,
) -> dict[VariableID, Any]:
    """Collect nonanticipative investment variables from a GTEP model.

    Parameters
    ----------
    model:
        Pyomo GTEP model. The model may be GDP-transformed. For GDP investment
        status disjuncts, this function collects the associated binary variables
        through Pyomo's object API.
    config:
        Nonanticipativity selection configuration. If omitted, defaults are used.

    Returns
    -------
    dict[VariableID, VarData]
        Mapping from durable ``VariableID`` to Pyomo variable data object.

    Notes
    -----
    This function uses explicit model structure and Pyomo component objects. It
    does not scan variable names and does not match variables by strings.
    """
    if config is None:
        config = NonanticipativityConfig()

    if not hasattr(model, "investmentStage"):
        raise ValueError("Model does not contain an investmentStage component.")

    variable_map: dict[VariableID, Any] = {}

    for stage in model.stages:
        investment_block = model.investmentStage[stage]

        if config.include_renewable_investment_variables:
            _collect_renewable_investment_variables(variable_map, investment_block)

        if config.include_thermal_status_binaries:
            _collect_thermal_status_binaries(variable_map, investment_block)

        if config.include_transmission_status_binaries:
            _collect_transmission_status_binaries(variable_map, investment_block)

        if config.include_storage_status_binaries:
            _collect_storage_status_binaries(variable_map, investment_block)

    if not variable_map:
        logger.warning("No nonanticipative variables were collected.")

    return variable_map


def collect_nonanticipative_variable_ids(
    model: pyo.ConcreteModel,
    config: NonanticipativityConfig | None = None,
) -> list[VariableID]:
    """Collect nonanticipative variable IDs from a GTEP model."""
    return list(collect_nonanticipative_variables(model, config).keys())


def bind_variable_ids_to_model(
    model: pyo.ConcreteModel,
    variable_ids: Iterable[VariableID],
) -> dict[VariableID, Any]:
    """Bind serialized/reloaded variable IDs to variables on a rebuilt model.

    Parameters
    ----------
    model:
        Rebuilt Pyomo model.
    variable_ids:
        Iterable of durable variable IDs.

    Returns
    -------
    dict[VariableID, VarData]
        Mapping from each ID to the corresponding model variable.

    Raises
    ------
    TypeError
        If any located component is not a Pyomo variable data object.
    KeyError
        If a variable ID cannot be located on the model.
    """
    bound: dict[VariableID, Any] = {}

    for variable_id in variable_ids:
        component = variable_id.find_on(model)

        if not _is_variable_data(component):
            raise TypeError(
                "ComponentUID resolved to a component that is not a Pyomo "
                f"variable data object: {variable_id.serialize()!r}"
            )

        bound[variable_id] = component

    return bound


def extract_nonanticipative_values(
    nonant_var_map: Mapping[VariableID, Any],
    *,
    require_values: bool = True,
) -> dict[VariableID, float]:
    """Extract numeric values of nonanticipative variables.

    Parameters
    ----------
    nonant_var_map:
        Mapping from ``VariableID`` to Pyomo variable data object.
    require_values:
        If true, raise an error if any value is unavailable. If false, variables
        with unavailable values are skipped.

    Returns
    -------
    dict[VariableID, float]
        Numeric variable values indexed by ``VariableID``.
    """
    values: dict[VariableID, float] = {}

    for variable_id, var in nonant_var_map.items():
        value = pyo.value(var, exception=False)

        if value is None:
            if require_values:
                raise ValueError(
                    "No value is available for nonanticipative variable with "
                    f"ComponentUID {variable_id.serialize()!r}."
                )
            continue

        values[variable_id] = float(value)

    return values


def serialize_variable_value_records(
    values: Mapping[VariableID, float],
) -> list[dict[str, Any]]:
    """Serialize nonanticipative values as JSON records.

    The returned structure is intentionally a list of records, not a dictionary
    keyed by variable name.
    """
    return [
        {
            "variable_id": variable_id.serialize(),
            "value": float(value),
        }
        for variable_id, value in sorted(
            values.items(),
            key=lambda item: item[0].serialize(),
        )
    ]


def deserialize_variable_value_records(
    records: Iterable[Mapping[str, Any]],
) -> dict[VariableID, float]:
    """Deserialize JSON records into values indexed by ``VariableID``."""
    values: dict[VariableID, float] = {}

    for record in records:
        if "variable_id" not in record:
            raise ValueError("Variable-value record is missing 'variable_id'.")
        if "value" not in record:
            raise ValueError("Variable-value record is missing 'value'.")

        variable_id = VariableID.from_serialized(record["variable_id"])
        values[variable_id] = float(record["value"])

    return values


def extract_nonanticipative_metadata(
    nonant_var_map: Mapping[VariableID, Any],
    *,
    include_values: bool = False,
) -> list[NonanticipativeVariableMetadata]:
    """Extract optional metadata for nonanticipative variables."""
    metadata: list[NonanticipativeVariableMetadata] = []

    for variable_id, var in nonant_var_map.items():
        value = None
        if include_values:
            raw_value = pyo.value(var, exception=False)
            value = None if raw_value is None else float(raw_value)

        lb = var.lb
        ub = var.ub

        metadata.append(
            NonanticipativeVariableMetadata(
                variable_id=variable_id,
                is_binary=bool(var.is_binary()),
                is_integer=bool(var.is_integer()),
                is_continuous=bool(var.is_continuous()),
                fixed=bool(var.fixed),
                lower_bound=None if lb is None else float(pyo.value(lb)),
                upper_bound=None if ub is None else float(pyo.value(ub)),
                value=value,
            )
        )

    return metadata


def serialize_nonanticipative_metadata(
    metadata: Iterable[NonanticipativeVariableMetadata],
) -> list[dict[str, Any]]:
    """Serialize nonanticipative metadata records."""
    return [
        item.to_jsonable()
        for item in sorted(
            metadata,
            key=lambda item: item.variable_id.serialize(),
        )
    ]


def make_component_map(
    nonant_var_map: Mapping[VariableID, Any],
) -> ComponentMap:
    """Create a Pyomo ``ComponentMap`` from variable object to ``VariableID``.

    This is useful inside model-construction code where object-based lookup is
    preferable.
    """
    component_map = ComponentMap()
    for variable_id, var in nonant_var_map.items():
        component_map[var] = variable_id
    return component_map


def assert_same_variable_id_set(
    first: Iterable[VariableID],
    second: Iterable[VariableID],
    *,
    first_label: str = "first",
    second_label: str = "second",
) -> None:
    """Assert that two collections contain exactly the same variable IDs."""
    first_set = set(first)
    second_set = set(second)

    missing_from_second = first_set - second_set
    missing_from_first = second_set - first_set

    if missing_from_second or missing_from_first:
        message_parts = ["Nonanticipative variable ID sets do not match."]

        if missing_from_second:
            message_parts.append(
                f"Present in {first_label} but missing from {second_label}: "
                + ", ".join(
                    variable_id.serialize()
                    for variable_id in sorted(
                        missing_from_second,
                        key=lambda item: item.serialize(),
                    )
                )
            )

        if missing_from_first:
            message_parts.append(
                f"Present in {second_label} but missing from {first_label}: "
                + ", ".join(
                    variable_id.serialize()
                    for variable_id in sorted(
                        missing_from_first,
                        key=lambda item: item.serialize(),
                    )
                )
            )

        raise ValueError(" ".join(message_parts))


def _collect_renewable_investment_variables(
    variable_map: dict[VariableID, Any],
    investment_block: Any,
) -> None:
    """Collect renewable investment-state/capacity variables."""
    components = []

    try:
        components.append(investment_block.renewableOperational)
        components.append(investment_block.renewableInstalled)
        components.append(investment_block.renewableRetired)
        components.append(investment_block.renewableExtended)
        components.append(investment_block.renewableDisabled)
    except AttributeError:
        return

    for component in components:
        _add_var_component(variable_map, component)


def _collect_thermal_status_binaries(
    variable_map: dict[VariableID, Any],
    investment_block: Any,
) -> None:
    """Collect associated binaries for thermal investment-status disjuncts."""
    disjunct_components = []

    try:
        disjunct_components.append(investment_block.genOperational)
        disjunct_components.append(investment_block.genInstalled)
        disjunct_components.append(investment_block.genRetired)
        disjunct_components.append(investment_block.genDisabled)
        disjunct_components.append(investment_block.genExtended)
    except AttributeError:
        return

    for disjunct_component in disjunct_components:
        _add_associated_binaries_from_disjunct_component(
            variable_map,
            disjunct_component,
        )


def _collect_transmission_status_binaries(
    variable_map: dict[VariableID, Any],
    investment_block: Any,
) -> None:
    """Collect associated binaries for transmission investment-status disjuncts."""
    disjunct_components = []

    try:
        disjunct_components.append(investment_block.branchOperational)
        disjunct_components.append(investment_block.branchInstalled)
        disjunct_components.append(investment_block.branchRetired)
        disjunct_components.append(investment_block.branchDisabled)
        disjunct_components.append(investment_block.branchExtended)
    except AttributeError:
        return

    for disjunct_component in disjunct_components:
        _add_associated_binaries_from_disjunct_component(
            variable_map,
            disjunct_component,
        )


def _collect_storage_status_binaries(
    variable_map: dict[VariableID, Any],
    investment_block: Any,
) -> None:
    """Collect associated binaries for storage investment-status disjuncts."""
    disjunct_components = []

    try:
        disjunct_components.append(investment_block.storOperational)
        disjunct_components.append(investment_block.storInstalled)
        disjunct_components.append(investment_block.storRetired)
        disjunct_components.append(investment_block.storDisabled)
        disjunct_components.append(investment_block.storExtended)
    except AttributeError:
        return

    for disjunct_component in disjunct_components:
        _add_associated_binaries_from_disjunct_component(
            variable_map,
            disjunct_component,
        )


def _add_var_component(
    variable_map: dict[VariableID, Any],
    var_component: Any,
) -> None:
    """Add all variable data objects from a Var component."""
    for var_data in _iter_component_data(var_component):
        if not _is_variable_data(var_data):
            raise TypeError("Expected a Pyomo variable data object.")
        _add_variable(variable_map, var_data)


def _add_associated_binaries_from_disjunct_component(
    variable_map: dict[VariableID, Any],
    disjunct_component: Any,
) -> None:
    """Add associated binary variables from a Disjunct component."""
    for disjunct_data in _iter_component_data(disjunct_component):
        binary_var = _get_associated_binary_from_disjunct(disjunct_data)
        _add_variable(variable_map, binary_var)


def _get_associated_binary_from_disjunct(disjunct_data: Any) -> Any:
    """Return the binary variable associated with a GDP disjunct.

    The primary path uses the indicator Boolean variable's associated binary.
    A direct ``binary_indicator_var`` fallback is included because Pyomo exposes
    that object on DisjunctData in common GDP workflows.
    """
    binary_var = disjunct_data.indicator_var.get_associated_binary()

    if binary_var is None:
        try:
            binary_var = disjunct_data.binary_indicator_var
        except AttributeError as err:
            raise RuntimeError(
                "Could not locate associated binary variable for a GDP disjunct. "
                "Ensure the GDP transformation or logical transformation has "
                "created associated binary variables before PH terms are added."
            ) from err

    if not _is_variable_data(binary_var):
        raise TypeError(
            "GDP disjunct associated binary is not a Pyomo variable data object."
        )

    if not binary_var.is_binary():
        raise TypeError("GDP disjunct associated variable is not binary.")

    return binary_var


def _add_variable(
    variable_map: dict[VariableID, Any],
    var_data: Any,
) -> None:
    """Add one variable data object to the nonanticipative variable map."""
    variable_id = VariableID.from_component(var_data)

    existing = variable_map.get(variable_id)
    if existing is not None and existing is not var_data:
        raise ValueError(
            "Duplicate ComponentUID resolved to different variable objects: "
            f"{variable_id.serialize()!r}"
        )

    variable_map[variable_id] = var_data


def _iter_component_data(component: Any) -> Iterable[Any]:
    """Iterate over scalar or indexed Pyomo component data objects."""
    if component.is_indexed():
        yield from component.values(sort=SortComponents.ORDERED_INDICES)
    else:
        yield component


def _is_variable_data(component: Any) -> bool:
    """Return true if component behaves as a Pyomo variable data object."""
    return component.is_variable_type()
