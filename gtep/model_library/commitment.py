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

"""Variables and Constraints for the Commitment Stage in the
Generation and Transmission Expansion Planning (GTEP) Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u

import gtep.model_library.storage as stor


def add_commitment_disjuncts(b, commitment_period):
    """This method adds discrete alternatives and constraints
    associated to the commitment period block. It defines a
    Disjunction with disjuncts representing the alternatives for
    generators status operation. The alternatives are:

    genOn:        Generator is operating and producing power.
    genStartup:   Generator is starting up during this period.
    genShutdown:  Generator is shutting down during this period.
    genOff:       Generator is offline and not producing power.

    """

    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    # Define disjunction on generator status operation:
    # on/startup/shutdown/off
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
                        * b.dispatchPeriod[dispatchPeriod].periodLength
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
                        * b.dispatchPeriod[dispatchPeriod].periodLength
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
                        * b.dispatchPeriod[dispatchPeriod].periodLength
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
                        * b.dispatchPeriod[dispatchPeriod].periodLength
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

    # Add constraints within disjunctions on battery storage
    # commitment (charging/discharging/off) if storage is needed.
    if m.config["storage"]:
        stor.add_storage_disjuncts(m, b, commitment_period)


def add_commitment_constraints(b, comm_per):
    """This method adds the commitment disjunctions and constraints to
    representative period block.

    """

    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    @b.Expression(doc="Total renewable surplus/deficit for commitment block")
    def renewableSurplusCommitment(b):
        return sum(
            b.dispatchPeriod[disp_per].renewableSurplusDispatch  # in MW
            for disp_per in b.dispatchPeriods
        )

    # Define total operating costs for commitment block. [TODO:
    # Replace this constraint with expressions using bounds transform
    # and check if the costs considered need to be re-assessed and
    # account for missing data.]

    @b.Expression()
    def operatingCostCommitment(b):
        if m.config["include_commitment"]:
            return (
                # [ESR WIP: This term includes the op cost for each
                # 15-min dispatch period.]
                sum(
                    b.dispatchPeriod[disp_per].operatingCostDispatch  # in $
                    for disp_per in b.dispatchPeriods
                )
                + sum(
                    m.fixedCost[gen] * b.commitmentPeriodLength
                    # [ESR WIP: Assuming we are paying for the full
                    # capacity of our generator. Note that a capacity
                    # should be included since the associated binaries
                    # are dimensionless. This makes the constraint
                    # unit consistent.]
                    * m.thermalCapacity[gen]
                    * (
                        b.genOn[gen].indicator_var.get_associated_binary()
                        + b.genShutdown[gen].indicator_var.get_associated_binary()
                        + b.genStartup[gen].indicator_var.get_associated_binary()
                    )
                    for gen in m.thermalGenerators
                )
                + sum(
                    m.startupCost[gen]
                    * b.genStartup[gen].indicator_var.get_associated_binary()
                    for gen in m.thermalGenerators
                )
            )
        else:
            return sum(
                b.dispatchPeriod[disp_per].operatingCostDispatch
                for disp_per in b.dispatchPeriods
            )

    # Define total storage costs for commitment block. [TODO: Replace
    # this constraint with expressions using bounds transform and
    # check if costs considered need to be re-assessed and account for
    # missing data.]

    # Compute Battery Storage cost per dispatch period if storage is
    # needed
    if m.config["storage"]:

        @b.Expression()
        def storageCostCommitment(b):
            return sum(
                b.dispatchPeriod[disp_per].storageCostDispatch
                for disp_per in b.dispatchPeriods
            )

    @b.Expression(doc="Total curtailment for commitment block in MW")
    def renewableCurtailmentCommitment(b):
        return sum(
            b.dispatchPeriod[disp_per].renewableCurtailmentDispatch
            for disp_per in b.dispatchPeriods
        )


def add_investment_commitment_variables(b):
    b.operatingCostInvestment = pyo.Var(
        within=pyo.NonNegativeReals, initialize=0, units=u.USD
    )
    b.renewableCurtailmentInvestment = pyo.Var(
        within=pyo.NonNegativeReals, initialize=0, units=u.USD
    )


def add_investment_commitment_constraints(m, b, investment_stage):

    @b.Constraint(doc="Curtailment penalties for investment period")
    def renewable_curtailment_cost(b):
        renewableCurtailmentRep = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                renewableCurtailmentRep += (
                    m.weights[rep_per]
                    * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .renewableCurtailmentCommitment  # in MW
                    # [ESR Question: Do we need to include this term here?]
                    * m.curtailmentCost
                )  # units are in $
        return (
            b.renewableCurtailmentInvestment  # in $
            == m.investmentFactor[investment_stage] * renewableCurtailmentRep
        )

    @b.Expression(doc="Operating costs for investment period")
    def commitmentOperatingCostInvestment(b):
        return m.investmentFactor[investment_stage] * sum(
            sum(
                m.weights[rep_per]
                * b.representativePeriod[rep_per]
                .commitmentPeriod[com_per]
                .operatingCostCommitment
                for com_per in b.representativePeriod[rep_per].commitmentPeriods
            )
            for rep_per in b.representativePeriods
        )
