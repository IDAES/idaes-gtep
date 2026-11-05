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

"""Advanced hydropower modeling for the GTEP model (``advanced_hydro=True``).

Every hydro unit (``unit_type == "HYDRO"``) is modeled separately from
the other renewables. With ``hydro_budget=False`` (the default) each unit
runs at its ``p_max`` series in every commitment period.

With ``hydro_budget=True`` each unit is dispatched between a minimum and
a maximum output in each commitment period, and has an energy budget
over each representative period:

- maximum: the ``p_max`` series (PMax MW);
- minimum: the ``DAY_AHEAD_hydro_min.csv`` series when the unit has one,
  otherwise ``p_max``;
- inflow: the ``DAY_AHEAD_hydro_budget.csv`` series when the unit has
  one, otherwise ``p_max``.

Over each representative period,

    sum_t  generation_t * dt  +  spill  ==  sum_t  inflow_t * dt     [MWh]

with ``spill >= 0``. Generation can therefore be shifted within the
representative period, and water that is not used is spilled. Budgets
are not carried between representative periods.

With neither file, the minimum equals the maximum and the budget equals
their sum, so output is fixed at ``p_max`` (the same as
``hydro_budget=False``). With only a minimum file, the unit can run
anywhere between its minimum and ``p_max`` and spill the rest.
"""

import logging

import pyomo.environ as pyo
from pyomo.environ import units as u

logger = logging.getLogger("gtep.model_library.hydropower_generation")


def _series_value(gen_data, key, commitment_period):
    """Value of a generator attribute for one commitment period.

    Handles EGRET time series (``{"data_type": "time_series", ...}``) and
    plain scalars. Time series are indexed by commitment period, like the
    other GTEP time series (period c uses the c-th entry of the
    representative day).
    """
    value = gen_data[key]
    if isinstance(value, dict):
        return float(value["values"][commitment_period - 1])
    return float(value)


def fix_hydropower_limits(b, commitmentPeriod):
    """Add hydro output limits (and inflows) to a commitment-period block."""

    m = b.model()
    gen_data = m.md.data["elements"]["generator"]

    max_values = {
        g: _series_value(gen_data[g], "p_max", commitmentPeriod)
        for g in m.hydroGenerators
    }

    if m.config["hydro_budget"]:
        min_values = {}
        budget_values = {}
        for g in m.hydroGenerators:
            min_key = "hydro_p_min" if "hydro_p_min" in gen_data[g] else "p_max"
            value = _series_value(gen_data[g], min_key, commitmentPeriod)
            if value > max_values[g]:
                logger.warning(
                    f"Hydro unit {g}: minimum output {value} MW exceeds maximum "
                    f"{max_values[g]} MW in commitment period {commitmentPeriod} "
                    f"({b.name}); using the maximum as the minimum."
                )
                value = max_values[g]
            min_values[g] = value

            budget_key = "hydro_budget" if "hydro_budget" in gen_data[g] else "p_max"
            value = _series_value(gen_data[g], budget_key, commitmentPeriod)
            budget_values[g] = value
    else:
        min_values = max_values
        budget_values = max_values

    b.hydroCapacityExpected = pyo.Param(
        m.hydroGenerators,
        initialize=max_values,
        units=u.MW,
        doc="Maximum hydro output for each hydro generator in this commitment period",
    )

    b.hydroMinimumExpected = pyo.Param(
        m.hydroGenerators,
        initialize=min_values,
        units=u.MW,
        doc="Minimum hydro output for each hydro generator in this commitment period",
    )

    b.hydroInflow = pyo.Param(
        m.hydroGenerators,
        initialize=budget_values,
        units=u.MW,
        doc=(
            "Reservoir inflow in this commitment period, as average MW of "
            "generation-equivalent energy"
        ),
    )


def add_representative_hydropower(b, repPer):
    """Add the hydro energy-budget constraints to a representative-period block.

    Only used when ``hydro_budget`` is True. Without a budget the minimum
    equals the maximum in every period, so there is nothing to shift or
    spill and the constraints would be redundant.
    """
    m = b.model()
    if not m.config["hydro_budget"]:
        return

    # Relative tolerance for the feasibility check on budgets vs minimum output.
    budget_rtol = 1e-9

    gen_data = m.md.data["elements"]["generator"]

    def dispatch_blocks():
        for c in b.commitmentPeriods:
            for d in b.commitmentPeriod[c].dispatchPeriods:
                yield b.commitmentPeriod[c], b.commitmentPeriod[c].dispatchPeriod[d]

    b.hydroSpill = pyo.Var(
        m.hydroGenerators,
        within=pyo.NonNegativeReals,
        initialize=0,
        units=u.MW * u.hr,
        doc="Unused (spilled) hydro energy over the representative period",
    )

    @b.Expression(
        m.hydroGenerators,
        doc="Inflow energy available over the representative period (MWh)",
    )
    def hydroEnergyBudget(b, g):
        return sum(
            b.commitmentPeriod[c].hydroInflow[g]
            * u.convert(b.commitmentPeriod[c].commitmentPeriodLength, to_units=u.hr)
            for c in b.commitmentPeriods
        )

    @b.Expression(
        m.hydroGenerators,
        doc="Hydro energy generated over the representative period (MWh)",
    )
    def hydroEnergyGenerated(b, g):
        return sum(
            d_blk.hydroGeneration[g]
            * u.convert(d_blk.dispatchPeriodLength, to_units=u.hr)
            for _, d_blk in dispatch_blocks()
        )

    # The minimum-output bounds must leave room inside the budget, or the
    # model is infeasible. Catch that here with a readable message.
    for g in m.hydroGenerators:
        budget = pyo.value(b.hydroEnergyBudget[g])
        min_energy = sum(
            pyo.value(c_blk.hydroMinimumExpected[g])
            * pyo.value(u.convert(d_blk.dispatchPeriodLength, to_units=u.hr))
            for c_blk, d_blk in dispatch_blocks()
        )
        if min_energy > budget + budget_rtol * max(1.0, abs(budget)):
            min_src = (
                "DAY_AHEAD_hydro_min.csv"
                if "hydro_p_min" in gen_data[g]
                else "p_max (no DAY_AHEAD_hydro_min.csv column for this unit)"
            )
            budget_src = (
                "DAY_AHEAD_hydro_budget.csv"
                if "hydro_budget" in gen_data[g]
                else "p_max (no DAY_AHEAD_hydro_budget.csv column for this unit)"
            )
            raise ValueError(
                f"Hydro unit {g} in {b.name} ({b.representative_date}): minimum output "
                f"requires {min_energy:.4g} MWh but the inflow budget is only "
                f"{budget:.4g} MWh. Minimum taken from {min_src}; inflow taken "
                f"from {budget_src}."
            )

    @b.Constraint(
        m.hydroGenerators,
        doc="Hydro generation plus spill equals the inflow budget over the period",
    )
    def hydro_energy_budget(b, g):
        return b.hydroEnergyGenerated[g] + b.hydroSpill[g] == b.hydroEnergyBudget[g]
