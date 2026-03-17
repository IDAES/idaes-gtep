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


def add_dispatch_variables(b, dispatch_period):
    """Add dispatch-associated variables to representative period block."""

    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    def thermal_generation_limits(
        b, thermalGen, doc="Bounds on active generation of thermal generators"
    ):
        return (0, m.thermalCapacity[thermalGen])

    b.thermalGeneration = pyo.Var(
        m.thermalGenerators,
        domain=pyo.NonNegativeReals,
        bounds=thermal_generation_limits,
        initialize=0,
        units=u.MW,
        doc="Thermal generation",
    )

    """ Battery --> Storage Parameters """
    if m.config["storage"]:

        def storage_capacity_limits(b, bat):
            return (
                m.minStorageChargeLevel[bat],
                m.storageCapacity[bat],
            )  # The lower bound should be > 0 - data input

        # TODO: Note that this does not fix initial battery capacity at the first dispatch period - need to adjust constraint
        def init_storage_capacity(b, bat):
            return m.initStorageChargeLevel[bat]

        b.storageChargeLevel = pyo.Var(
            m.storage,
            domain=pyo.NonNegativeReals,
            bounds=storage_capacity_limits,
            initialize=init_storage_capacity,
            units=u.MW,
        )

        # Define bounds on charging/discharging capability. Note that constraints
        # enforce that there are min & max charge/discharge levels if the bat is in
        # the charging or discharging state
        def storage_charge_limits(b, bat):
            return (0, m.chargeMax[bat])

        def storage_discharge_limits(b, bat):
            return (0, m.dischargeMax[bat])

        b.storageCharged = pyo.Var(
            m.storage,
            domain=pyo.NonNegativeReals,
            bounds=storage_charge_limits,
            initialize=0,
            units=u.MW,
        )

        b.storageDischarged = pyo.Var(
            m.storage,
            domain=pyo.NonNegativeReals,
            bounds=storage_discharge_limits,
            initialize=0,
            units=u.MW,
        )

    if m.config["flow_model"] == "ACR" or m.config["flow_model"] == "ACP":
        # Define bounds on thermal generator reactive generation
        def thermal_reactive_generation_limits(b, thermalGen):
            return (0, m.thermalReactiveCapacity[thermalGen])

        b.thermalReactiveGeneration = pyo.Var(
            m.thermalGenerators,
            domain=pyo.Reals,
            bounds=thermal_reactive_generation_limits,
            initialize=0,
            units=u.MVAR,
            doc="Thermal generation",
        )

    # [ESR WIP: Still deciding if this should be Nameplate]
    def renewable_generation_limits(
        b, renewableGen, doc="Bounds on active generation of renewable generator in MW"
    ):
        return (0, m.renewableCapacityNameplate[renewableGen])

    b.renewableGeneration = pyo.Var(
        m.renewableGenerators,
        domain=pyo.NonNegativeReals,
        bounds=renewable_generation_limits,
        initialize=0,
        units=u.MW,
        doc="Renewable generation",
    )

    # Define bounds on renewable generator curtailment
    def curtailment_limits(
        b, renewableGen, doc="Bounds on renewable generator curtailment in MW"
    ):
        return (0, m.renewableCapacityNameplate[renewableGen])

    b.renewableCurtailment = pyo.Var(
        m.renewableGenerators,
        domain=pyo.NonNegativeReals,
        bounds=curtailment_limits,
        initialize=0,
        units=u.MW,
        doc="Curtailment of renewable generators",
    )

    @b.Expression(
        m.renewableGenerators, doc="Surplus generation per renewable generator in MW"
    )
    def renewableGenerationSurplus(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] - b.renewableCurtailment[renewableGen]
        )

    # [ESR WIP: Added the dispatch period length.]
    @b.Expression(m.renewableGenerators, doc="Curtailment cost per generator in $")
    def renewableCurtailmentCost(b, renewableGen):
        return (
            b.renewableCurtailment[renewableGen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.curtailmentCost
        )

    # [ESR WIP: Added the dispatch period length.]
    @b.Expression(m.thermalGenerators, doc="Cost per thermal generator in $")
    def thermalGeneratorCost(b, gen):
        return (
            b.thermalGeneration[gen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * (m.fixedCost[gen] + m.varCost[gen])
        )

    # [ESR WIP: Added the dispatch period length.]
    @b.Expression(m.renewableGenerators, doc="Cost per renewable generator in $")
    def renewableGeneratorCost(b, gen):
        return (
            b.renewableGeneration[gen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.fixedCost[gen]
        )

    if m.config["flow_model"] == "ACR" or m.config["flow_model"] == "ACP":
        # Per generator reactive power cost
        @b.Expression(m.thermalGenerators)
        def reactiveGeneratorCost(b, gen):
            return b.thermalReactiveGeneration[gen] * m.fuelCost[gen]

    b.loadShed = pyo.Var(
        m.buses,
        domain=pyo.NonNegativeReals,
        initialize=0,
        units=u.MW,
        doc="Load shed per bus",
    )

    # [ESR WIP: Add the dispatch period length.]
    @b.Expression(m.buses, doc="Load shed cost per bus in $")
    def loadShedCost(b, bus):
        return (
            b.loadShed[bus]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.loadShedCostperCurtailment  # $/MWh
        )

    if m.config["storage"]:
        """Per-Battery Operational cost variables"""

        @b.Expression(m.storage)
        def storageChargingCost(b, bat):
            return b.storageCharged[bat] * m.chargingCost[bat]

        # JSC addn Per Battery Discharging Cost
        @b.Expression(m.storage)
        def storageDischargingCost(b, bat):
            return b.storageDischarged[bat] * m.dischargingCost[bat]

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
        """Per-Battery Operational costs"""
        b.chargingCostDispatch = sum(b.storageChargingCost.values())

        b.dischargingCostDispatch = sum(b.storageDischargingCost.values())

        b.storageCostDispatch = b.chargingCostDispatch + b.dischargingCostDispatch
    else:
        b.chargingCostDispatch = 0
        b.dischargingCostDispatch = 0
        b.storageCostDispatch = 0

    # BLN TODO: Check the config check in the Expression rule
    @b.Expression(doc="Total cost for dispatch in $")
    def operatingCostDispatch(b):
        return (
            b.thermalGenerationCostDispatch
            + b.reactiveGenerationCostDispatch
            + b.renewableGenerationCostDispatch
            + b.loadShedCostDispatch
            + b.curtailmentCostDispatch
            + b.chargingCostDispatch
            + b.dischargingCostDispatch
            + b.storageCostDispatch
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

    # (Original) NOTE: this is an abuse of units and needs to be fixed for
    # variable temporal resolution
    b.powerFlow = pyo.Var(
        m.transmission,
        domain=pyo.Reals,
        bounds=power_flow_limits,
        initialize=0,
        units=u.MW,
        doc="Power flow in MW",
    )

    @b.Disjunct(m.transmission)
    def branchInUse(disj, branch):
        b = disj.parent_block()

        def bus_angle_bounds(disj, bus, doc="Voltage angle"):
            return (-1000, 1000)
            return (-math.pi / 6, math.pi / 6)

        # Only create bus angle variables for the buses associated with this
        # branch that is in use
        disj.branch_buses = [
            bb
            for bb in m.buses
            if (
                m.transmission[branch]["from_bus"] == bb
                or m.transmission[branch]["to_bus"] == bb
            )
        ]

        disj.busAngle = pyo.Var(
            disj.branch_buses,
            domain=pyo.Reals,
            initialize=0,
            bounds=bus_angle_bounds,
        )

        def delta_bus_angle_bounds(disj, bus, doc="Voltage angle"):
            return (-math.pi / 6, math.pi / 6)

        def delta_bus_angle_rule(disj, doc="Maximum bus angle discrepancy"):
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            return disj.busAngle[tb] - disj.busAngle[fb]

        disj.deltaBusAngle = pyo.Var(
            domain=pyo.Reals,
            bounds=delta_bus_angle_bounds,
            rule=delta_bus_angle_rule,
        )

        if m.config["flow_model"] == "ACP":
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            resistance = m.md.data["elements"]["branch"][branch].get("resistance", 0.0)
            reactance = m.md.data["elements"]["branch"][branch].get("reactance", 1e-6)

            # Transformer tap ratio and phase shift
            if m.md.data["elements"]["branch"][branch]["branch_type"] == "transformer":
                reactance *= m.md.data["elements"]["branch"][branch][
                    "transformer_tap_ratio"
                ]
                phase_shift = m.md.data["elements"]["branch"][branch][
                    "transformer_phase_shift"
                ]
            else:
                phase_shift = 0

            admittance = 1 / complex(resistance, reactance)
            G = admittance.real
            B = admittance.imag

            # Define voltage magnitude variables for from and to buses
            disj.voltage_from = pyo.Var(bounds=(0, 2))
            disj.voltage_to = pyo.Var(bounds=(0, 2))

            # Define active and reactive power flow variables
            disj.P_flow = pyo.Var(bounds=(-1000, 1000))
            disj.Q_flow = pyo.Var(bounds=(-1000, 1000))

            # Polar Active Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_p(disj):
                return disj.P_flow == (
                    disj.voltage_from**2 * G
                    - disj.voltage_from
                    * disj.voltage_to
                    * (
                        G * pyo.cos(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                        + B
                        * pyo.sin(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                    )
                )

            # Polar Reactive Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_q(disj):
                return disj.Q_flow == (
                    -disj.voltage_from**2 * B
                    - disj.voltage_from
                    * disj.voltage_to
                    * (
                        G * pyo.sin(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                        - B
                        * pyo.cos(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                    )
                )

        if m.config["flow_model"] == "ACR":
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            resistance = m.md.data["elements"]["branch"][branch].get("resistance", 0.0)
            reactance = m.md.data["elements"]["branch"][branch].get("reactance", 1e-6)

            # Transformer tap ratio and phase shift
            if m.md.data["elements"]["branch"][branch]["branch_type"] == "transformer":
                reactance *= m.md.data["elements"]["branch"][branch][
                    "transformer_tap_ratio"
                ]
                phase_shift = m.md.data["elements"]["branch"][branch][
                    "transformer_phase_shift"
                ]
            else:
                phase_shift = 0

            admittance = 1 / complex(resistance, reactance)
            G = admittance.real
            B = admittance.imag

            # Define rectangular voltage variables for from and to buses
            disj.real_voltage_from = pyo.Var(bounds=(0, 2))
            disj.real_voltage_to = pyo.Var(bounds=(-2, 2))
            disj.imag_voltage_from = pyo.Var(bounds=(0, 2))
            disj.imag_voltage_to = pyo.Var(bounds=(-2, 2))

            # Define active and reactive power flow variables
            disj.P_flow = pyo.Var(bounds=(-1000, 1000))
            disj.Q_flow = pyo.Var(bounds=(-1000, 1000))

            # Rectangular Active Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_p(disj):
                Vf_r = disj.real_voltage_from
                Vf_i = disj.imag_voltage_from
                Vt_r = disj.real_voltage_to
                Vt_i = disj.imag_voltage_to

                # Active Power Flow Equation
                return disj.P_flow == (
                    G * (Vf_r**2 + Vf_i**2)
                    - G * (Vf_r * Vt_r + Vf_i * Vt_i)
                    - B * (Vf_r * Vt_i - Vf_i * Vt_r)
                )

            # Rectangular Reactive Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_q(disj):
                Vf_r = disj.real_voltage_from
                Vf_i = disj.imag_voltage_from
                Vt_r = disj.real_voltage_to
                Vt_i = disj.imag_voltage_to

                # Reactive Power Flow Equation
                return disj.Q_flow == (
                    B * (Vf_r**2 + Vf_i**2)
                    + B * (Vf_r * Vt_r + Vf_i * Vt_i)
                    - G * (Vf_r * Vt_i - Vf_i * Vt_r)
                )

        if m.config["flow_model"] == "DC":

            @disj.Constraint()
            def dc_power_flow(disj):
                fb = m.transmission[branch]["from_bus"]
                tb = m.transmission[branch]["to_bus"]
                reactance = m.md.data["elements"]["branch"][branch]["reactance"]
                if (
                    m.md.data["elements"]["branch"][branch]["branch_type"]
                    == "transformer"
                ):
                    reactance *= m.md.data["elements"]["branch"][branch][
                        "transformer_tap_ratio"
                    ]
                    shift = m.md.data["elements"]["branch"][branch][
                        "transformer_phase_shift"
                    ]
                else:
                    shift = 0
                # TODO: Fix the units in this constraint
                return b.powerFlow[branch] / u.MW == (-1 / reactance) * (
                    disj.busAngle[tb] - disj.busAngle[fb] + shift
                )

    @b.Disjunct(m.transmission)
    def branchNotInUse(disj, branch):

        # Fixing power flow to 0 and not creating bus angle variables for
        # branches that are not in use.
        @disj.Constraint()
        def dc_power_flow(disj):
            return b.powerFlow[branch] == 0 * u.MW

        return

    # Branches are either in-use or not. This disjunction may provide the
    # basis for transmission switching in the future.
    @b.Disjunction(m.transmission)
    def branchInUseStatus(disj, branch):
        return [disj.branchInUse[branch], disj.branchNotInUse[branch]]

    if m.config["transmission"]:
        # JSC update - If a branch is in use, it must be active
        # Update this when switching is implemented
        @b.LogicalConstraint(m.transmission)
        def must_use_active_branches(b, branch):
            return b.branchInUse[branch].indicator_var.implies(
                pyo.lor(
                    i_p.branchOperational[branch].indicator_var,
                    i_p.branchInstalled[branch].indicator_var,
                    i_p.branchExtended[branch].indicator_var,
                )
            )

        # JSC update - If a branch is not in use, it must be inactive.
        # Update this when switching is implemented
        @b.LogicalConstraint(m.transmission)
        def cannot_use_inactive_branches(b, branch):
            return b.branchNotInUse[branch].indicator_var.implies(
                pyo.lor(
                    i_p.branchDisabled[branch].indicator_var,
                    i_p.branchRetired[branch].indicator_var,
                )
            )

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
    """Add dispatch-associated inequalities to representative period block."""
    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    # [ESR WIP: Commented for now but think about how to implement this in
    # a better way.]
    # for key in m.loads.keys():
    #     m.loads[key] *= max(0, rng.normal(0.5, 0.2))

    # Energy balance constraint
    if m.config["flow_model"] == "CP":

        @b.Constraint()
        def CP_flow_balance(b):
            balance = 0
            # load = m.loads.get(bus) or 0
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

        @b.Constraint(m.buses)
        def flow_balance(b, bus):
            balance = 0
            # [ESR WIP: Comment load and call all the loads as a
            # parameter instead of the original dictionary. Also, note
            # that the loads are now declared for all the buses in
            # m.buses, and set to 0 for the buses that are not in
            # m.load_buses.]
            # load = c_p.loads.get(bus) or 0

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

            balance -= m.loads[bus]  # add new parameter (already includes units)
            balance += b.loadShed[bus]
            return balance == 0 * u.MW

    # Capacity factor constraint
    # NOTE: In comparison to reference work, this is *per renewable generator*
    # JKS - charging costs from non-colocated plants?
    @b.Constraint(m.renewableGenerators)
    def capacity_factor(b, renewableGen):
        # if m.md.data["elements"]["generator"][renewableGen]["fuel"] == "H":
        #     # print(m.md.data["elements"]["generator"][renewableGen])
        #     return pyo.Constraint.Skip
        return (
            b.renewableGeneration[renewableGen] + b.renewableCurtailment[renewableGen]
            == c_p.renewableCapacityExpected[renewableGen]
        )

    ## TODO: (@jkskolf) add renewableExtended to this and anywhere else
    @b.Constraint(m.renewableGenerators)
    def operational_renewables_only(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen]
            <= i_p.renewableInstalled[renewableGen]
            + i_p.renewableOperational[renewableGen]
            + i_p.renewableExtended[renewableGen]
        )

    # RESERVE -- total operating (spinning + quickstart)
    # Total operating reserve constraint
    ## NOTE: min operating reserve is a percentage of load
    ## FIXME: Reserve enforcement causes infeasibility issues.  We should track
    ## reserve shortage in some way and find a way to penalize it -- how is this
    ## done in ISOs?  Is it an issue to assign this as a regional
    # @b.Constraint(m.regions)
    # def total_operating_reserve(b, region):
    #     return sum(
    #         b.spinningReserve[gen] + b.quickstartReserve[gen]
    #         for gen in m.thermalGenerators & m.gensAtRegion[region]
    #     ) >= m.minOperatingReserve[region] * (
    #         sum(
    #             (m.loads.get(bus) or 0)
    #             for bus in m.buses
    #             if m.md.data["elements"]["bus"][bus]["area"] == region
    #         )
    #     )

    # # Total spinning reserve constraint
    # @b.Constraint(m.regions)
    # def total_spinning_reserve(b, region):
    #     return sum(
    #         b.spinningReserve[gen]
    #         for gen in m.thermalGenerators & m.gensAtRegion[region]
    #     ) >= m.minSpinningReserve[region] * sum(
    #         (m.loads.get(bus) or 0)
    #         for bus in m.buses
    #         if m.md.data["elements"]["bus"][bus]["area"] == region
    #     )
