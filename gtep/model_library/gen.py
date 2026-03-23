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

"""Constraints for the Generators in the Generation and Transmission
Expansion Planning (GTEP) Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u


def add_generators_status_disjuncts(b, thermalgens_set, renewablegens_set):
    """This method declares variables and a Disjunction and its
    disjuncts to model the selection of the renewable and thermal
    generators status. The alternatives are:

    Renewable:
    renewableOperational: Renewable is active and producing power.
    renewableInstalled:   Renewable is newly added and active.
    renewableRetired:     Renewable is removed from service.
    renewableDisabled:    Renewable is temporarily out of service.
    renewableExtended:    Renewable is upgraded beyond its original capacity.

    Thermal:
    genOperational: Generator is active and producing power.
    genInstalled:   Generator is newly added and active.
    genRetired:     Generator is removed from service.
    genDisabled:    Generator is temporarily out of service.
    genExtended:    Generator is upgraded beyond its original capacity.

    """

    # Declare renewable generators status as variables, instead of
    # disjuncts. The alternatives are operational, installed, retired,
    # and extended. [TODO: Convert these to discrete decisions.]
    b.renewableOperational = pyo.Var(
        renewablegens_set, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableInstalled = pyo.Var(
        renewablegens_set, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableRetired = pyo.Var(
        renewablegens_set, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableExtended = pyo.Var(
        renewablegens_set, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableDisabled = pyo.Var(
        renewablegens_set, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )

    @b.Disjunct(thermalgens_set)
    def genOperational(disj, gen):
        return

    @b.Disjunct(thermalgens_set)
    def genInstalled(disj, gen):
        return

    @b.Disjunct(thermalgens_set)
    def genRetired(disj, gen):
        return

    @b.Disjunct(thermalgens_set)
    def genDisabled(disj, gen):
        return

    @b.Disjunct(thermalgens_set)
    def genExtended(disj, gen):
        return

    @b.Disjunction(thermalgens_set)
    def genInvestStatus(disj, gen):
        return [
            disj.genOperational[gen],
            disj.genInstalled[gen],
            disj.genRetired[gen],
            disj.genDisabled[gen],
            disj.genExtended[gen],
        ]


def add_investment_generators_constraints(m, b, investment_stage):

    # These constraints take the "in_service" data and fix the
    # indicator variables of the status disjuncts to define the
    # operation of generators.
    for gen in m.thermalGenerators:
        if (
            m.md.data["elements"]["generator"][gen]["in_service"] == False
            and investment_stage == 1
        ):
            b.genOperational[gen].indicator_var.fix(False)
        elif (
            m.md.data["elements"]["generator"][gen]["in_service"] == True
            and investment_stage == 1
        ):
            b.genOperational[gen].indicator_var.fix(True)

    for gen in m.renewableGenerators:
        if (
            m.md.data["elements"]["generator"][gen]["in_service"] == False
            and investment_stage == 1
        ):
            b.renewableOperational[gen].fix(0)
        elif (
            m.md.data["elements"]["generator"][gen]["in_service"] == True
            and investment_stage == 1
        ):
            b.renewableOperational[gen].fix(m.renewableCapacityNameplate[gen])

    @b.Expression(doc="Generators investment costs in $")
    def generators_investment_cost(b):
        return (
            sum(
                # [ESR: When including the disjunction
                # genInvestStatus, think if we should replace this
                # cost with generatorInstallationCost.]
                m.generatorInvestmentCost[gen]
                * m.thermalCapacity[gen]  # in MW
                * m.capitalMultiplier[gen]
                * b.genInstalled[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.capitalMultiplier[gen]
                * b.renewableInstalled[gen]  # in MW
                for gen in m.renewableGenerators
            )
            + sum(
                # [ESR WIP: When including the disjunction
                # genInvestStatus, think if we should replace this cost
                # with generatorInstallationCost.]
                m.generatorInvestmentCost[gen]
                * m.extensionMultiplier[gen]
                * m.thermalCapacity[gen]
                * b.genExtended[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.extensionMultiplier[gen]
                * b.renewableExtended[gen]
                for gen in m.renewableGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.retirementMultiplier[gen]
                * b.renewableRetired[gen]
                for gen in m.renewableGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.retirementMultiplier[gen]
                * m.thermalCapacity[gen]
                * b.genRetired[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
        )

    # Add legacy equations below. These equations are not used in
    # current versions of the model. [TODO: Determine if we need them
    # in future versions of the model.]

    """
    # NOTE: renewableCapacityValue is a percentage of renewableCapacity
    # TODO: renewableCapacityValue ==> renewableCapacityFactor
    # NOTE: reserveMargin is a percentage of peakLoad
    # TODO: check and re-enable with additional bounding transform before bigm
    # TODO: renewableCapacityValue... should this be time iterated? is it tech based?
    # is it site based? who can say?

    @b.Constraint(doc="Planning reserve requirement constraint")
    def planning_reserve_requirement(b):
        return (
            sum(
                m.renewableCapacityNameplate[gen]
                * m.renewableCapacityValue[gen]
                * (b.renewableOperational[gen] + b.renewableInstalled[gen])
                for gen in m.renewableGenerators
            )
            + sum(
                m.thermalCapacity[gen]
                * (
                    b.genOperational[gen].indicator_var.get_associated_binary()
                    + b.genInstalled[gen].indicator_var.get_associated_binary()
                )
                for gen in m.thermalGenerators
            )
            >= (1 + m.reserveMargin[investment_stage]) * m.peakLoad[investment_stage]
        )

    # NOTE: temporarily disabled maximum investment as default option
    # TODO: These capacities shouldn't be enabled by default since they can
    # easily cause absurd results/possibly even infeasibility.  Will need to add
    # user-defined handling for this.

    @b.Constraint(
        m.regions,
        doc="Maximum investment stage installation for thermal generators",
    )
    def maximum_thermal_investment(b, region):
        return (
            sum(
                m.thermalCapacity[gen]
                * b.genInstalled[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators & m.gensAtRegion[region]
            )
            <= b.maxThermalInvestment[region]
        )

    @b.Constraint(
        m.regions,
        doc="Maximum investment stage installation for renewable generators",
    )
    def maximum_renewable_investment(b, region):
        return (
            sum(
                m.renewableCapacityNameplate[gen]
                * b.genInstalled[gen].indicator_var.get_associated_binary()
                for gen in m.renewableGenerators & m.gensAtRegion[region]
            )
            <= b.maxRenewableInvestment[region]
            if m.renewableGenerators & m.gensAtRegion[region]
            else pyo.Constraint.Skip
        )
  
    """


def add_generators_state_disjuncts(m, b, r_p, i_p, commitment_period):
    """This method defines a Disjunction with disjuncts representing
    the alternatives for generators state operation. The alternatives
    are:

    genOn:        Generator is operating and producing power.
    genStartup:   Generator is starting up during this period.
    genShutdown:  Generator is shutting down during this period.
    genOff:       Generator is offline and not producing power.

    """

    @b.Disjunct(m.thermalGenerators)
    def genOn(disj, generator):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods, doc="Minimum operating limits")
        def operating_limit_min(d, dispatchPeriod):
            return (
                m.thermalMin[generator]
                <= b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
            )

        @disj.Constraint(b.dispatchPeriods, doc="Maximum operating limits")
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalCapacity[generator]
            )

        @disj.Constraint(
            b.dispatchPeriods,
            m.thermalGenerators,
            doc="Ramp up limits for fully-on thermal generators",
        )
        def ramp_up_limits(disj, dispatchPeriod, generator):
            if dispatchPeriod != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    - b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                    <= m.rampUpRates[generator]
                    * b.dispatchPeriod[dispatchPeriod].periodLength
                    * m.thermalCapacity[generator]
                )
            elif dispatchPeriod == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    - r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .thermalGeneration[generator]
                    <= m.rampUpRates[generator]
                    * b.dispatchPeriod[dispatchPeriod].periodLength
                    * m.thermalCapacity[generator]
                )
            else:
                return pyo.Constraint.Skip

        @disj.Constraint(
            b.dispatchPeriods,
            m.thermalGenerators,
            doc="Ramp down limits for fully-on thermal generators",
        )
        def ramp_down_limits(disj, dispatchPeriod, generator):
            if dispatchPeriod != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                    - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    <= m.rampDownRates[generator]  # in MW/min
                    * b.dispatchPeriod[dispatchPeriod].periodLength  # in min
                    * m.thermalCapacity[generator]  # in MW
                )
            elif dispatchPeriod == 1 and commitment_period != 1:
                return (
                    r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .thermalGeneration[generator]
                    - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    <= m.rampDownRates[generator]  # in MW/min
                    * b.dispatchPeriod[dispatchPeriod].periodLength  # in min
                    * m.thermalCapacity[generator]  # in MW
                )
            else:
                return pyo.Constraint.Skip

        # NOTE: maxSpinningReserve is a percentage of thermalCapacity
        @disj.Constraint(
            b.dispatchPeriods, m.thermalGenerators, doc="Maximum spinning reserve"
        )
        def max_spinning_reserve(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.maxSpinningReserve[generator] * m.thermalCapacity[generator]
            )

    # [TODO: add quick start reserve = 0]

    @b.Disjunct(m.thermalGenerators)
    def genStartup(disj, generator):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods, doc="Minimum operating limits")
        def operating_limit_min(d, dispatchPeriod):
            return (
                0 * u.MW
                <= b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
            )

        @disj.Constraint(b.dispatchPeriods, doc="Maximum operating limits")
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalMin[generator]
            )

        # [TODO: Define if the max() function is necessary here.]
        @disj.Constraint(
            b.dispatchPeriods,
            m.thermalGenerators,
            doc="Ramp up constraints for generators starting up",
        )
        def ramp_up_limits(disj, dispatchPeriod, generator):
            if dispatchPeriod != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    - b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                    <= max(
                        pyo.value(m.thermalMin[generator]),
                        # [ESR: Make sure the time units are consistent
                        # here since we are only taking the value]
                        pyo.value(m.rampUpRates[generator])
                        * pyo.value(
                            b.dispatchPeriod[dispatchPeriod].periodLength
                        )  # in minutes
                        * pyo.value(m.thermalCapacity[generator]),
                    )
                    * u.MW
                )
            elif dispatchPeriod == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    - r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .thermalGeneration[generator]
                    <= max(
                        pyo.value(m.thermalMin[generator]),
                        # [ESR: Make sure the time units are consistent
                        # here since we are only taking the value]
                        pyo.value(m.rampUpRates[generator])
                        * pyo.value(
                            b.dispatchPeriod[dispatchPeriod].periodLength
                        )  # in minutes
                        * pyo.value(m.thermalCapacity[generator]),
                    )
                    * u.MW
                )
            else:
                return pyo.Constraint.Skip

    @b.Disjunct(m.thermalGenerators)
    def genShutdown(disj, generator):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods, doc="Minimum operating limits")
        def operating_limit_min(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                >= 0 * u.MW
            )

        @disj.Constraint(b.dispatchPeriods, doc="Maximum operating limits")
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalMin[generator]
            )

        # [RMA: We may need to turn off ramp down constraints for
        # feasibility purposes. We will need to think about this for
        # future work, but commenting this out is probably fine for
        # the purposes of this paper.]

        @disj.Constraint(
            b.dispatchPeriods,
            m.thermalGenerators,
            doc="Ramp down constraints for generators shutting down",
        )
        def ramp_down_limits(disj, dispatchPeriod, generator):
            if dispatchPeriod != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                    - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    <= max(
                        pyo.value(m.thermalMin[generator]),
                        # [ESR: Make sure the time units are consistent
                        # here since we are taking the value only]
                        pyo.value(m.rampDownRates[generator])
                        * pyo.value(
                            b.dispatchPeriod[dispatchPeriod].periodLength
                        )  # in minutes
                        * pyo.value(m.thermalCapacity[generator]),
                    )
                    * u.MW
                )
            elif dispatchPeriod == 1 and commitment_period != 1:
                return (
                    r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .thermalGeneration[generator]
                    - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                    <= max(
                        pyo.value(m.thermalMin[generator]),
                        # [ESR: Make sure the time units are consistent
                        # here since we are taking the value only]
                        pyo.value(m.rampDownRates[generator])
                        * pyo.value(
                            b.dispatchPeriod[dispatchPeriod].periodLength
                        )  # in minutes
                        * pyo.value(m.thermalCapacity[generator]),
                    )
                    * u.MW
                )
            else:
                return pyo.Constraint.Skip

    @b.Disjunct(m.thermalGenerators)
    def genOff(disj, generator):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods, doc="Maximum operating limits")
        def operating_limit_max(disj, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                <= 0 * u.MW
            )

        # NOTE: maxQuickstartReserve is a percentage of
        # thermalCapacity. [TODO: This isn't needed so, instead, we
        # need to set spinning reserve to 0.]
        @disj.Constraint(
            b.dispatchPeriods,
            m.thermalGenerators,
            doc="Maximum quickstart reserve constraint",
        )
        def max_quickstart_reserve(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].quickstartReserve[generator]
                <= m.maxQuickstartReserve[generator] * m.thermalCapacity[generator]
            )

    @b.Disjunction(m.thermalGenerators, doc="Disjunction for generator status")
    def genStatus(disj, generator):
        return [
            disj.genOn[generator],
            disj.genStartup[generator],
            disj.genShutdown[generator],
            disj.genOff[generator],
        ]

    @b.LogicalConstraint(
        m.thermalGenerators,
        doc="Enforces that generators cannot be committed unless they are operational or just installed",
    )
    def commit_active_gens_only(b, generator):
        return pyo.lor(
            b.genOn[generator].indicator_var,
            b.genStartup[generator].indicator_var,
            b.genShutdown[generator].indicator_var,
        ).implies(
            pyo.lor(
                i_p.genOperational[generator].indicator_var,
                i_p.genInstalled[generator].indicator_var,
                i_p.genExtended[generator].indicator_var,
            )
        )


def add_generators_logical_constraints(m):
    """This method defines logical constraints to ensure that thermal
    and renewable generators statuses transitions are operationally
    consistent over time, across the investment stages.

    """

    # [TODO: This needs to be tested.]
    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Thermal generator retirement",
    )
    def thermalgen_retirement(m, stage, gen):
        return (
            (
                m.investmentStage[stage - pyo.value(m.lifetimes[gen])]
                .genOperational[gen]
                .indicator_var
                | m.investmentStage[stage - pyo.value(m.lifetimes[gen])].genInstalled[
                    gen
                ]
            ).implies(
                m.investmentStage[stage].genRetired[gen].indicator_var
                | m.investmentStage[stage].genExtended[gen].indicator_var
            )
            if stage > pyo.value(m.lifetimes[gen])
            else pyo.LogicalConstraint.Skip
        )

    # Total renewable generation (in MW) operational at a given stage
    # is equal to what was operational and/or installed in the previous stage
    # less what was retired in the previous stage
    @m.Constraint(m.stages, m.renewableGenerators)
    def renewable_stats_link(m, stage, gen):
        return (
            m.investmentStage[stage].renewableOperational[gen]
            == m.investmentStage[stage - 1].renewableOperational[gen]
            + m.investmentStage[stage - 1].renewableInstalled[gen]
            - m.investmentStage[stage - 1].renewableExtended[gen]
            - m.investmentStage[stage - 1].renewableRetired[gen]
            if stage != 1
            else pyo.Constraint.Skip
        )

    @m.Constraint(m.stages, m.renewableGenerators)
    def renewable_retirement_link(m, stage, gen):
        return (
            m.investmentStage[stage].renewableRetired[gen]
            == m.investmentStage[stage - 1].renewableRetired[gen]
            - m.investmentStage[stage - 1].renewableDisabled[gen]
            if stage != 1
            else pyo.Constraint.Skip
        )

    @m.Constraint(m.stages, m.renewableGenerators)
    def renewable_extension_link(m, stage, gen):
        return (
            m.investmentStage[stage].renewableExtended[gen]
            == m.investmentStage[stage - 1].renewableExtended[gen]
            - m.investmentStage[stage - 1].renewableRetired[gen]
            if stage != 1
            else pyo.Constraint.Skip
        )

    @m.Constraint(m.stages, m.renewableGenerators)
    def renewable_capacity_enforcement(m, stage, gen):
        return (
            m.investmentStage[stage].renewableOperational[gen]
            + m.investmentStage[stage].renewableInstalled[gen]
            + m.investmentStage[stage].renewableExtended[gen]
            + m.investmentStage[stage].renewableRetired[gen]
            + m.investmentStage[stage].renewableRetired[gen]
            <= m.renewableCapacityNameplate[gen]
        )

    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Enforces that, if a thermal gen is online at time t, it must have been online or installed at time t-1",
    )
    def consistent_operation(m, stage, gen):
        return (
            m.investmentStage[stage]
            .genOperational[gen]
            .indicator_var.implies(
                m.investmentStage[stage - 1].genOperational[gen].indicator_var
                | m.investmentStage[stage - 1].genInstalled[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Enforces that, if a thermal gen is online at time t, it must be online, extended, or retired at time t+1",
    )
    def consistent_operation_future(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genOperational[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genOperational[gen].indicator_var
                | m.investmentStage[stage].genExtended[gen].indicator_var
                | m.investmentStage[stage].genRetired[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Enforces that, if a thermal gen is retired in period t-1, it must be disabled in period t",
    )
    def full_retirement(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genRetired[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genDisabled[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Enforces that, if a thermal gen is disabled at time t-1, it must stay disabled  at time t",
    )
    def consistent_disabled(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genDisabled[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genDisabled[gen].indicator_var
                | m.investmentStage[stage].genInstalled[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Enforces that, if a thermal gen is extended at time t-1, it must stay extended or be retired at time t",
    )
    def consistent_extended(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genExtended[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genExtended[gen].indicator_var
                | m.investmentStage[stage].genRetired[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.thermalGenerators,
        doc="Enforces that, if a thermal gen is installed in period t-1, it must be operational in period t",
    )
    def full_investment(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genInstalled[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genOperational[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )


def add_dispatch_generators_variables(m, b):

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

    if m.config["flow_model"] == "ACR" or m.config["flow_model"] == "ACP":

        def thermal_reactive_generation_limits(
            b, thermalGen, doc="Bounds on thermal generator reactive generation"
        ):
            return (0, m.thermalReactiveCapacity[thermalGen])

        b.thermalReactiveGeneration = pyo.Var(
            m.thermalGenerators,
            domain=pyo.Reals,
            bounds=thermal_reactive_generation_limits,
            initialize=0,
            units=u.MVAR,
            doc="Thermal generation",
        )

    # [ESR: Still deciding if this should be Nameplate]
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
