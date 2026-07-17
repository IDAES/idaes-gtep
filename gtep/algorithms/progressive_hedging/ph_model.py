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
"""Pyomo model augmentation for GTEP Progressive Hedging.

This module adds MILP-compatible Progressive Hedging terms to a single-scenario
GTEP model.

Design constraints
------------------
* Nonanticipative variables are identified externally with ``VariableID``
  objects backed by Pyomo ``ComponentUID``.
* This module never matches variables by Pyomo component names or string paths.
* PH terms are added using the actual Pyomo variable objects supplied in the
  nonanticipative variable map.
* The initial implementation avoids MIQP subproblems:
    - binary variables receive the standard quadratic PH penalty after using
      \(x^2 = x\), resulting in a linear term,
    - non-binary variables receive an L1 absolute-deviation penalty.

Scenario objective
------------------
For representative-period scenario \(s\), this module creates

\[
    \pi_s C^{investment}(x_s)
    + C^{operating}_s(x_s, y_s)
    + C^{curtailment}_s(x_s, y_s)
    + \text{PH terms}.
\]

The model-wide quota-deficit and deficit-penalty machinery is intentionally not
used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import logging

import pyomo.environ as pyo
from pyomo.environ import units as u

from gtep.algorithms.progressive_hedging.nonanticipativity import VariableID

logger = logging.getLogger("gtep.algorithms.progressive_hedging.ph_model")


@dataclass(frozen=True)
class PHObjectiveOptions:
    """Options controlling PH objective augmentation."""

    include_multiplier_terms: bool = True
    include_regularization_terms: bool = True
    deactivate_existing_objectives: bool = True


def add_progressive_hedging_block(
    model: pyo.ConcreteModel,
    nonant_var_map: Mapping[VariableID, Any],
    *,
    scenario_probability: float,
    xbar: Mapping[VariableID, float] | None = None,
    multipliers: Mapping[VariableID, float] | None = None,
    rho: float | Mapping[VariableID, float] = 1.0,
    options: PHObjectiveOptions | None = None,
) -> pyo.Block:
    """Add PH parameters, terms, and objective to a single-scenario GTEP model.

    Parameters
    ----------
    model:
        Pyomo GTEP scenario model.
    nonant_var_map:
        Mapping from ``VariableID`` to the actual Pyomo variable data object on
        ``model``.
    scenario_probability:
        Normalized representative-period probability \(\pi_s\). This scales
        the first-stage investment cost in the scenario objective.
    xbar:
        Consensus values \(\bar{x}\) indexed by ``VariableID``. Required when
        regularization terms are enabled.
    multipliers:
        PH multiplier values \(w_s\) indexed by ``VariableID``. If omitted,
        multipliers default to zero.
    rho:
        Either a scalar PH penalty or a mapping from ``VariableID`` to penalty.
    options:
        Objective augmentation options.

    Returns
    -------
    pyo.Block
        The added ``model.ph`` block.
    """
    if options is None:
        options = PHObjectiveOptions()

    _validate_inputs(
        nonant_var_map=nonant_var_map,
        scenario_probability=scenario_probability,
        xbar=xbar,
        multipliers=multipliers,
        rho=rho,
        include_regularization_terms=options.include_regularization_terms,
    )

    if hasattr(model, "ph"):
        raise RuntimeError("Model already contains a 'ph' block.")

    model.ph = pyo.Block()
    ph = model.ph

    # Store Python-only lookup tables. These are not used for cross-process
    # matching; they only connect integer PH indices to already-supplied Pyomo
    # variable objects and ComponentUID-backed VariableID objects.
    ordered_variable_ids = list(nonant_var_map.keys())

    ph._index_to_variable_id = {
        idx: variable_id
        for idx, variable_id in enumerate(ordered_variable_ids, start=1)
    }
    ph._variable_id_to_index = {
        variable_id: idx for idx, variable_id in ph._index_to_variable_id.items()
    }
    ph._index_to_var = {
        idx: nonant_var_map[variable_id]
        for idx, variable_id in ph._index_to_variable_id.items()
    }

    ph.NONANT_INDEX = pyo.RangeSet(
        len(ordered_variable_ids),
        doc="Integer index over nonanticipative variables",
    )

    ph.scenario_probability = pyo.Param(
        initialize=float(scenario_probability),
        mutable=False,
        within=pyo.NonNegativeReals,
        units=u.dimensionless,
        doc="Normalized PH scenario probability",
    )

    ph.term = pyo.Block(ph.NONANT_INDEX)

    for idx in ph.NONANT_INDEX:
        variable_id = ph._index_to_variable_id[idx]
        var = ph._index_to_var[idx]

        _add_single_variable_ph_term(
            ph.term[idx],
            var=var,
            variable_id=variable_id,
            xbar_value=_lookup_value(
                xbar,
                variable_id,
                default=0.0,
                required=options.include_regularization_terms,
                label="xbar",
            ),
            multiplier_value=_lookup_value(
                multipliers,
                variable_id,
                default=0.0,
                required=False,
                label="multipliers",
            ),
            rho_value=_lookup_rho_value(rho, variable_id),
            include_multiplier_terms=options.include_multiplier_terms,
            include_regularization_terms=options.include_regularization_terms,
        )

    ph.linear_multiplier_term = pyo.Expression(
        expr=sum(ph.term[idx].linear_multiplier_term for idx in ph.NONANT_INDEX),
        doc="Total linear PH multiplier term",
    )

    ph.regularization_term = pyo.Expression(
        expr=sum(ph.term[idx].regularization_term for idx in ph.NONANT_INDEX),
        doc="Total MILP-compatible PH regularization term",
    )

    ph.total_ph_term = pyo.Expression(
        expr=ph.linear_multiplier_term + ph.regularization_term,
        doc="Total PH augmentation term",
    )

    if options.deactivate_existing_objectives:
        deactivate_active_objectives(model)

    model.ph_scenario_objective = pyo.Objective(
        expr=_scenario_base_objective(model) + ph.total_ph_term,
        sense=pyo.minimize,
        doc="Progressive Hedging scenario objective",
    )

    return ph


def update_progressive_hedging_parameters(
    model: pyo.ConcreteModel,
    *,
    xbar: Mapping[VariableID, float],
    multipliers: Mapping[VariableID, float],
    rho: float | Mapping[VariableID, float] | None = None,
) -> None:
    """Update mutable PH parameters on a model that already has ``model.ph``.

    This helper is included for completeness even though the initial Torc design
    rebuilds each scenario model for every PH iteration.
    """
    if not hasattr(model, "ph"):
        raise RuntimeError("Model does not contain a 'ph' block.")

    ph = model.ph

    for idx in ph.NONANT_INDEX:
        variable_id = ph._index_to_variable_id[idx]
        term = ph.term[idx]

        if variable_id not in xbar:
            raise KeyError(
                "Missing xbar value for nonanticipative variable "
                f"{variable_id.serialize()!r}."
            )
        if variable_id not in multipliers:
            raise KeyError(
                "Missing multiplier value for nonanticipative variable "
                f"{variable_id.serialize()!r}."
            )

        term.xbar.set_value(float(xbar[variable_id]))
        term.w.set_value(float(multipliers[variable_id]))

        if rho is not None:
            term.rho.set_value(float(_lookup_rho_value(rho, variable_id)))


def deactivate_active_objectives(model: pyo.ConcreteModel) -> None:
    """Deactivate all active objectives on a model."""
    active_objectives = list(model.component_data_objects(pyo.Objective, active=True))

    if not active_objectives:
        logger.warning("No active objectives were found to deactivate.")

    for objective in active_objectives:
        objective.deactivate()


def _add_single_variable_ph_term(
    block: pyo.Block,
    *,
    var: Any,
    variable_id: VariableID,
    xbar_value: float,
    multiplier_value: float,
    rho_value: float,
    include_multiplier_terms: bool,
    include_regularization_terms: bool,
) -> None:
    """Add PH terms for one nonanticipative variable to an indexed block."""
    var_units = _get_variable_units(var)
    coefficient_units = _coefficient_units_for_variable(var_units)

    # Store only non-Pyomo metadata on the PH term block. Do not assign the
    # Pyomo variable object as an attribute of this block: that would attempt to
    # re-parent the variable component from its original owning block to the PH
    # block, which Pyomo does not allow.
    object.__setattr__(block, "_variable_id", variable_id)

    block.xbar = pyo.Param(
        initialize=float(xbar_value),
        mutable=True,
        units=var_units,
        doc="PH consensus value for this nonanticipative variable",
    )

    block.w = pyo.Param(
        initialize=float(multiplier_value),
        mutable=True,
        units=coefficient_units,
        doc="PH multiplier for this nonanticipative variable",
    )

    block.rho = pyo.Param(
        initialize=float(rho_value),
        mutable=True,
        within=pyo.NonNegativeReals,
        units=coefficient_units,
        doc="PH penalty coefficient for this nonanticipative variable",
    )

    if include_multiplier_terms:
        linear_multiplier_expr = block.w * var
    else:
        linear_multiplier_expr = 0 * u.USD

    block.linear_multiplier_term = pyo.Expression(
        expr=linear_multiplier_expr,
        doc="Linear PH multiplier contribution",
    )

    if not include_regularization_terms:
        regularization_expr = 0 * u.USD

    elif var.is_binary():
        # Standard PH quadratic penalty:
        #
        #   rho / 2 * (x - xbar)^2
        #
        # For binary x, x^2 = x, so ignoring the constant xbar^2 term gives:
        #
        #   rho / 2 * (1 - 2*xbar) * x
        #
        # This keeps the subproblem linear.
        regularization_expr = 0.5 * block.rho * (1 - 2 * block.xbar) * var

    else:
        # L1 penalty for non-binary variables:
        #
        #   rho * abs(x - xbar)
        #
        # modeled through an auxiliary absolute-deviation variable.
        block.abs_dev = pyo.Var(
            within=pyo.NonNegativeReals,
            units=var_units,
            doc="Absolute deviation from PH consensus",
        )

        block.abs_dev_pos = pyo.Constraint(
            expr=block.abs_dev >= var - block.xbar,
            doc="Positive side of absolute deviation",
        )

        block.abs_dev_neg = pyo.Constraint(
            expr=block.abs_dev >= block.xbar - var,
            doc="Negative side of absolute deviation",
        )

        regularization_expr = block.rho * block.abs_dev

    block.regularization_term = pyo.Expression(
        expr=regularization_expr,
        doc="MILP-compatible PH regularization contribution",
    )


def _scenario_base_objective(model: pyo.ConcreteModel) -> Any:
    """Return the non-PH base objective expression for one PH scenario.

    The expression is:

        scenario_probability * expansionCostTotal
        + operatingCostTotal
        + renewableCurtailmentCostTotal

    If ``renewableCurtailmentCostTotal`` is not present, the expression is
    constructed directly from investment-stage curtailment variables.
    """
    required_components = [
        "expansionCostTotal",
        "operatingCostTotal",
        "investmentStage",
        "stages",
    ]

    missing = [name for name in required_components if not hasattr(model, name)]
    if missing:
        raise AttributeError(
            "Cannot construct PH scenario objective because the model is missing "
            f"required component(s): {', '.join(missing)}"
        )

    return (
        model.ph.scenario_probability * model.expansionCostTotal
        + model.operatingCostTotal
        + _renewable_curtailment_cost_total(model)
    )


def _renewable_curtailment_cost_total(model: pyo.ConcreteModel) -> Any:
    """Return renewable curtailment cost total without deficit-penalty terms."""
    if hasattr(model, "renewableCurtailmentCostTotal"):
        return model.renewableCurtailmentCostTotal

    return sum(
        model.investmentStage[stage].renewableCurtailmentInvestment
        for stage in model.stages
    )


def _validate_inputs(
    *,
    nonant_var_map: Mapping[VariableID, Any],
    scenario_probability: float,
    xbar: Mapping[VariableID, float] | None,
    multipliers: Mapping[VariableID, float] | None,
    rho: float | Mapping[VariableID, float],
    include_regularization_terms: bool,
) -> None:
    """Validate inputs to PH block construction."""
    if scenario_probability < 0:
        raise ValueError(
            f"scenario_probability must be nonnegative. Received {scenario_probability}."
        )

    if not nonant_var_map:
        raise ValueError("Cannot add PH block with an empty nonant_var_map.")

    for variable_id, var in nonant_var_map.items():
        if not isinstance(variable_id, VariableID):
            raise TypeError(
                "nonant_var_map keys must be VariableID objects. "
                f"Received {type(variable_id)}."
            )
        if not var.is_variable_type():
            raise TypeError(
                "nonant_var_map values must be Pyomo variable data objects."
            )

    variable_ids = set(nonant_var_map)

    if include_regularization_terms:
        if xbar is None:
            raise ValueError("xbar is required when regularization terms are enabled.")
        missing_xbar = variable_ids - set(xbar)
        if missing_xbar:
            raise KeyError(
                "Missing xbar values for nonanticipative variable(s): "
                + ", ".join(
                    variable_id.serialize()
                    for variable_id in sorted(
                        missing_xbar,
                        key=lambda item: item.serialize(),
                    )
                )
            )

    if multipliers is not None:
        extra_multipliers = set(multipliers) - variable_ids
        if extra_multipliers:
            raise KeyError(
                "Multiplier values were provided for variable IDs that are not "
                "in the nonanticipative variable map: "
                + ", ".join(
                    variable_id.serialize()
                    for variable_id in sorted(
                        extra_multipliers,
                        key=lambda item: item.serialize(),
                    )
                )
            )

    if isinstance(rho, Mapping):
        missing_rho = variable_ids - set(rho)
        extra_rho = set(rho) - variable_ids

        if missing_rho:
            raise KeyError(
                "Missing rho values for nonanticipative variable(s): "
                + ", ".join(
                    variable_id.serialize()
                    for variable_id in sorted(
                        missing_rho,
                        key=lambda item: item.serialize(),
                    )
                )
            )

        if extra_rho:
            raise KeyError(
                "Rho values were provided for variable IDs that are not in the "
                "nonanticipative variable map: "
                + ", ".join(
                    variable_id.serialize()
                    for variable_id in sorted(
                        extra_rho,
                        key=lambda item: item.serialize(),
                    )
                )
            )

        for value in rho.values():
            if value < 0:
                raise ValueError("All rho values must be nonnegative.")
    else:
        if rho < 0:
            raise ValueError(f"rho must be nonnegative. Received {rho}.")


def _lookup_value(
    mapping: Mapping[VariableID, float] | None,
    variable_id: VariableID,
    *,
    default: float,
    required: bool,
    label: str,
) -> float:
    """Look up a value indexed by ``VariableID``."""
    if mapping is None:
        if required:
            raise KeyError(
                f"{label} value is required for nonanticipative variable "
                f"{variable_id.serialize()!r}."
            )
        return float(default)

    if variable_id not in mapping:
        if required:
            raise KeyError(
                f"Missing {label} value for nonanticipative variable "
                f"{variable_id.serialize()!r}."
            )
        return float(default)

    return float(mapping[variable_id])


def _lookup_rho_value(
    rho: float | Mapping[VariableID, float],
    variable_id: VariableID,
) -> float:
    """Look up scalar or variable-specific rho."""
    if isinstance(rho, Mapping):
        return float(rho[variable_id])
    return float(rho)


def _get_variable_units(var: Any) -> Any:
    """Return Pyomo units for a variable, defaulting to dimensionless."""
    var_units = pyo.units.get_units(var)
    if var_units is None:
        return u.dimensionless
    return var_units


def _coefficient_units_for_variable(var_units: Any) -> Any:
    """Return units for a linear objective coefficient multiplying ``var``."""
    return u.USD / var_units
