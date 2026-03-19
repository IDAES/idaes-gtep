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

"""Variables and Constraints for the Investment Stage in the
Generation and Transmission Expansion Planning (GTEP) Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u

import gtep.model_library.storage as stor
import gtep.model_library.transmission as transm
import gtep.model_library.commitment as commit


def add_investment_variables(b, investment_stage):
    """This method adds variables and disjuncts to the investment
    stage block. It includes a Disjunction and its disjuncts to model
    the selection of the thermal generators status. The alternatives
    are:

    genOperational: Generator is active and producing power.
    genInstalled:   Generator is newly added and active.
    genRetired:     Generator is removed from service.
    genDisabled:    Generator is temporarily out of service.
    genExtended:    Generator is upgraded beyond its original capacity.

    :param b: Investment block
    :param investment_stage: Investment stage index or name
    :return: None

    """

    m = b.model()
    b.investmentStage = investment_stage

    # Track and accumulate costs and penalties
    b.quotaDeficit = pyo.Var(within=pyo.NonNegativeReals, initialize=0, units=u.MW)
    b.expansionCost = pyo.Var(within=pyo.NonNegativeReals, initialize=0, units=u.USD)
    if m.config["include_commitment"]:
        commit.add_investment_commitment_variables(b)
    if m.config["storage"]:
        stor.add_investment_storage_variables(b)

    # Declare thermal generator disjuncts (operational, installed,
    # retired, disabled, extended)
    @b.Disjunct(m.thermalGenerators)
    def genOperational(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genInstalled(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genRetired(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genDisabled(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genExtended(disj, gen):
        return

    @b.Disjunction(m.thermalGenerators)
    def genInvestStatus(disj, gen):
        return [
            disj.genOperational[gen],
            disj.genInstalled[gen],
            disj.genRetired[gen],
            disj.genDisabled[gen],
            disj.genExtended[gen],
        ]

    # Declare renewable generators status as variables, instead of
    # disjuncts. The alternatives are operational, installed, retired,
    # and extended. [TODO: Convert these to discrete decisions.]
    b.renewableOperational = pyo.Var(
        m.renewableGenerators, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableInstalled = pyo.Var(
        m.renewableGenerators, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableRetired = pyo.Var(
        m.renewableGenerators, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableExtended = pyo.Var(
        m.renewableGenerators, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableDisabled = pyo.Var(
        m.renewableGenerators, within=pyo.NonNegativeReals, initialize=0, units=u.MW
    )

    # Add storage status alternatives, if needed
    if m.config["storage"]:
        stor.add_storage_disjuncts(b, m.storage)

    # Add transmission lines status alternatives, if needed
    if m.config["transmission"]:
        transm.add_transmission_disjuncts(b, m.transmission)


def add_investment_constraints(
    b,
    investment_stage,
):
    """This method adds standard inequalities (i.e., those not
    involving disjunctions) to the investment stage block.

    """

    m = b.model()

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

    if m.config["transmission"]:
        transm.add_investment_transmission_constraints(m, b, investment_stage)

    # Fixing in-service batteries initial investment state based on
    # input. [TODO: Also initialize storage level (state of charge)]
    if m.config["storage"]:
        stor.add_investment_storage_constraints(m, b, investment_stage)

    # NOTE: The following constraints can be split into rep_per and
    # invest_stage components if desired.

    # [TODO: The definition of the investment cost needs to be
    # revisited AND possibly depends on data format. NOTE: It is
    # _rare_ for these values to be defined at all, let alone
    # consistently.]
    @b.Expression(doc="Investment costs for the investment period in $")
    def investment_cost(b):
        baseline_cost = (
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
            + (
                b.storage_investment_cost
                if m.config["storage"] == True
                else (0 * u.USD)
            )
            + (
                b.transmission_investment_cost
                if m.config["transmission"] == True
                else (0 * u.USD)
            )
        )
        return m.investmentFactor[investment_stage] * baseline_cost

    if m.config["include_commitment"]:
        commit.add_investment_commitment_constraints(m, b, investment_stage)

    # [BLN: Convert this to a constraint using operatingCostInvestment
    # Var. May also need to move it]
    @b.Constraint(doc="Operating costs for investment period")
    def rule_operatingCostInvestment(b):
        return b.operatingCostInvestment == (
            b.commitmentOperatingCostInvestment
            if m.config["include_commitment"] == True
            else 0
        )

    # [ESR Question: Do we need to add the flag for investment? Should
    # we better put this flag when calling this function?]
    if m.config["include_investment"]:

        # NOTE: This is constraint (13) in reference [1]
        @b.Constraint(doc="Minimum per-stage renewable generation requirement")
        def renewable_generation_requirement(b):
            renewableSurplusRepresentative = 0
            ## TODO: preprocess loads for the appropriate sum here
            ed = 0
            for rep_per in b.representativePeriods:
                for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                    renewableSurplusRepresentative += (
                        m.weights[rep_per]
                        * b.representativePeriod[rep_per]
                        .commitmentPeriod[com_per]
                        .renewableSurplusCommitment
                    )
            return (
                renewableSurplusRepresentative + b.quotaDeficit
                >= m.renewableQuota[investment_stage] * ed
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
