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

"""Variables and Constraints for the Dispatch Stage in the Generation
and Transmission Expansion Planning (GTEP) Model

"""

import math

import pyomo.environ as pyo
from pyomo.environ import units as u

import gtep.model_library.gen as gens
import gtep.model_library.storage as stor
import gtep.model_library.transmission as transm


def add_dispatch_variables(b, dispatch_period):
    """This method adds dispatch-associated variables to
    representative period block.

    """

    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    # Add variables and bounds for generators and storage, when needed
    gens.add_dispatch_generators_variables(m, b)

    # Add expressions for thermal and renewable generators
    @b.Expression(
        m.renewableGenerators, doc="Surplus generation per renewable generator in MW"
    )
    def renewableGenerationSurplus(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] - b.renewableCurtailment[renewableGen]
        )

    @b.Expression(m.renewableGenerators, doc="Curtailment cost per generator in $")
    def renewableCurtailmentCost(b, renewableGen):
        return (
            b.renewableCurtailment[renewableGen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.curtailmentCost
        )

    @b.Expression(m.thermalGenerators, doc="Cost per thermal generator in $")
    def thermalGeneratorCost(b, gen):
        return (
            b.thermalGeneration[gen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * (m.fixedCost[gen] + m.varCost[gen])
        )

    @b.Expression(m.renewableGenerators, doc="Cost per renewable generator in $")
    def renewableGeneratorCost(b, gen):
        return (
            b.renewableGeneration[gen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.fixedCost[gen]
        )

    if m.config["flow_model"] == "ACR" or m.config["flow_model"] == "ACP":

        @b.Expression(m.thermalGenerators, doc="Reactive power cost per generator")
        def reactiveGeneratorCost(b, gen):
            return b.thermalReactiveGeneration[gen] * m.fuelCost[gen]

    b.loadShed = pyo.Var(
        m.buses,
        domain=pyo.NonNegativeReals,
        initialize=0,
        units=u.MW,
        doc="Load shed per bus",
    )

    @b.Expression(m.buses, doc="Load shed cost per bus in $")
    def loadShedCost(b, bus):
        return (
            b.loadShed[bus]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.loadShedCostperCurtailment  # $/MWh
        )

    # Track total dispatch values and costs
    @b.Expression(doc="Total surplus power for renewable generators in MW")
    def renewableSurplusDispatch(b):
        return sum(b.renewableGenerationSurplus[gen] for gen in m.renewableGenerators)

    @b.Expression()
    def thermalGenerationCostDispatch(b):
        return sum(b.thermalGeneratorCost[gen] for gen in m.thermalGenerators)

    @b.Expression()
    def renewableGenerationCostDispatch(b):
        return sum(b.renewableGeneratorCost[gen] for gen in m.renewableGenerators)

    # Reactive generation cost
    if m.config["flow_model"] == "ACR" or m.config["flow_model"] == "ACP":
        b.reactiveGenerationCostDispatch = sum(b.reactiveGeneratorCost.values())
    else:
        b.reactiveGenerationCostDispatch = 0

    @b.Expression()
    def loadShedCostDispatch(b):
        return sum(b.loadShedCost[bus] for bus in m.buses)

    @b.Expression()
    def curtailmentCostDispatch(b):
        return sum(b.renewableCurtailmentCost[gen] for gen in m.renewableGenerators)

    if m.config["storage"]:
        # Add storage variables and constraints. It also includes its
        # operational costs variables.
        stor.add_dispatch_storage_variables_and_constraints(m, b)

    # [BLN TODO: Check the config check in the Expression rule.]
    @b.Expression(doc="Total cost for dispatch in $")
    def operatingCostDispatch(b):

        # [ESR WIP: If I don't add the 0 value for storage cost
        # dispatch, the optimal solution has a value of 0. Check why
        # this is happening.]
        if m.config["storage"]:
            storage_term = b.storageCostDispatch
        else:
            storage_term = 0

        return (
            b.thermalGenerationCostDispatch
            + b.reactiveGenerationCostDispatch
            + b.renewableGenerationCostDispatch
            + b.loadShedCostDispatch
            + b.curtailmentCostDispatch
            + storage_term
        )

    @b.Expression(doc="Total curtailment dispatch for renewable generators in MW")
    def renewableCurtailmentDispatch(b):
        return sum(b.renewableCurtailment[gen] for gen in m.renewableGenerators)

    # Restrictions on flow over uninvested lines are enforced in a
    # disjuction below
    def power_flow_limits(b, branch, doc="Bounds on transmission line capacity"):
        return (
            -m.transmissionCapacity[branch],
            m.transmissionCapacity[branch],
        )

    # [TODO: Check units since they might need to be fixed for
    # variable temporal resolution.]
    b.powerFlow = pyo.Var(
        m.transmission,
        domain=pyo.Reals,
        bounds=power_flow_limits,
        initialize=0,
        units=u.MW,
        doc="Power flow in MW",
    )

    # Add transmission lines state disjuncts (in use and not in
    # use). The power flow is calculated here using OPF formulations.
    transm.add_transmission_state_disjuncts(m, b, i_p)

    def spinning_reserve_limits(
        b, thermalGen, doc="Bounds on thermal generator spinning reserve supply"
    ):
        return (
            0 * u.MW,
            m.spinningReserveFraction[thermalGen] * m.thermalCapacity[thermalGen],
        )

    b.spinningReserve = pyo.Var(
        m.thermalGenerators,
        domain=pyo.NonNegativeReals,
        bounds=spinning_reserve_limits,
        initialize=0,
        units=u.MW,
    )

    def quickstart_reserve_limits(
        b, thermalGen, doc="Bounds on thermal generator quickstart reserve supply"
    ):
        return (
            0 * u.MW,
            m.quickstartReserveFraction[thermalGen] * m.thermalCapacity[thermalGen],
        )

    b.quickstartReserve = pyo.Var(
        m.thermalGenerators,
        domain=pyo.NonNegativeReals,
        bounds=quickstart_reserve_limits,
        initialize=0,
        units=u.MW,
    )


def add_dispatch_constraints(b, disp_per):
    """This method adds dispatch-associated inequalities to the
    representative period block.

    """

    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    # Add power flow balance
    if m.config["flow_model"] == "CP":

        @b.Constraint(doc="Energy balance for copper-plate formulation")
        def CP_flow_balance(b):
            balance = 0
            buses = [bus for bus in m.buses]
            loads = [l for l in b.loads]
            gens = [gen for gen in m.generators]
            batts = [bat for bat in m.storage]
            balance += sum(
                b.thermalGeneration[g] for g in gens if g in m.thermalGenerators
            )
            balance += sum(
                b.renewableGeneration[g] for g in gens if g in m.renewableGenerators
            )
            """ Battery Storage added to flow balance constraint """
            balance += sum(b.storageDischarged[bt] for bt in batts)
            balance -= sum(b.storageCharged[bt] for bt in batts)

            balance -= sum(m.loads[l] for l in loads)
            balance += sum(b.loadShed[bus] for bus in buses)

            return balance == 0

    else:

        @b.Constraint(m.buses, doc="Energy balance constraint")
        def flow_balance(b, bus):
            balance = 0
            end_points = [
                line
                for line in m.transmission
                if m.transmission[line]["from_bus"] == bus
            ]
            start_points = [
                line for line in m.transmission if m.transmission[line]["to_bus"] == bus
            ]
            gens = [
                gen
                for gen in m.generators
                if m.md.data["elements"]["generator"][gen]["bus"] == bus
            ]
            batts = []
            if m.config["storage"]:
                batts = [
                    bat
                    for bat in m.storage
                    if m.md.data["elements"]["storage"][bat]["bus"] == bus
                ]
            balance -= sum(b.powerFlow[i] for i in end_points)
            balance += sum(b.powerFlow[i] for i in start_points)
            balance += sum(
                b.thermalGeneration[g] for g in gens if g in m.thermalGenerators
            )
            balance += sum(
                b.renewableGeneration[g] for g in gens if g in m.renewableGenerators
            )
            """ Battery Storage added to flow balance constraint """
            balance += sum(b.storageDischarged[bt] for bt in batts)
            balance -= sum(b.storageCharged[bt] for bt in batts)

            # Add the loads as a parameter (already includes units).
            balance -= m.loads[bus]
            balance += b.loadShed[bus]
            return balance == 0 * u.MW

    # NOTE: In comparison to reference [1], this is "per renewable
    # generator". [TODO: Should we include charging costs from
    # non-colocated plants?]
    @b.Constraint(m.renewableGenerators, doc="Capacity factor constraint")
    def capacity_factor(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] + b.renewableCurtailment[renewableGen]
            == c_p.renewableCapacityExpected[renewableGen]
        )

    # [TODO: Add renewableExtended to this and anywhere else.]
    @b.Constraint(m.renewableGenerators)
    def operational_renewables_only(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen]
            <= i_p.renewableInstalled[renewableGen]
            + i_p.renewableOperational[renewableGen]
            + i_p.renewableExtended[renewableGen]
        )

    # Add legacy equations. These are not used in current version of
    # model but keeping here to see how to integrate them in future
    # versions of the model.
    """
    # Calculate reserve considering total operating (spinning +
    # quickstart). Consider the minimum operating reserve is a
    # percentage of load. [TODO: Reserve enforcement causes
    # infeasibility issues.  We should track reserve shortage in some
    # way and find a way to penalize it. Check how it is done in ISOs?
    # Is it an issue to assign this as a regional?.]
    @b.Constraint(m.regions, doc="Total operating reserve constraint")
    def total_operating_reserve(b, region):
        return sum(
            b.spinningReserve[gen] + b.quickstartReserve[gen]
            for gen in m.thermalGenerators & m.gensAtRegion[region]
        ) >= m.minOperatingReserve[region] * (
            sum(
                (m.loads.get(bus) or 0)
                for bus in m.buses
                if m.md.data["elements"]["bus"][bus]["area"] == region
            )
        )

    @b.Constraint(m.regions, doc="Total spinning reserve constraint")
    def total_spinning_reserve(b, region):
        return sum(
            b.spinningReserve[gen]
            for gen in m.thermalGenerators & m.gensAtRegion[region]
        ) >= m.minSpinningReserve[region] * sum(
            (m.loads.get(bus) or 0)
            for bus in m.buses
            if m.md.data["elements"]["bus"][bus]["area"] == region
        )
    """
