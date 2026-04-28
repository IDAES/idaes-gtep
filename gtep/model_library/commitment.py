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

import gtep.model_library.gen as gens
import gtep.model_library.storage as stor
import gtep.model_library.scaling as scaling


def add_commitment_parameters(b, commitment_period, investmentStage):

    m = b.model()

    b.commitmentPeriodLength = pyo.Param(
        within=pyo.PositiveReals, default=1, units=u.hr
    )
    b.carbonTax = pyo.Param(default=0)

    # [ESR: Corrected to be in the commitment block "b", not in main
    # model "m".]
    b.renewableCapacityExpected = {}

    units_renewable_capacity = u.MW
    for renewableGen in m.renewableGenerators:
        if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float:
            b.renewableCapacityExpected[renewableGen] = 0 * units_renewable_capacity
        else:
            b.renewableCapacityExpected[renewableGen] = (
                m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"][
                    commitment_period - 1
                ]
                * units_renewable_capacity
            )

    # [TODO: Redesign load scaling and allow nature of
    # it as argument.]
    if m.config["scale_loads"]:
        scaling.add_load_scaling(
            m,
            b,
            commitment_period,
            investmentStage,
            scaling_value=10,
        )


def add_commitment_disjuncts(b, commitment_period):
    """This method adds discrete alternatives and constraints
    associated to the commitment period block.

    """

    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    if m.config["include_commitment"]:
        # Add disjunction and disjuncts to define the operational
        # state of generators (on/startup/shutdown/off) and storage
        # (charging/discharging/off), when needed.
        gens.add_generators_state_disjuncts(m, b, r_p, i_p, commitment_period)
    else:

        gens.generators_status_always_on(m, b)

    if m.config["storage"]:
        stor.add_storage_state_disjuncts(m, b, commitment_period)


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

    # [ESR TODO: Replace this constraint with expressions using bounds
    # transform and check if the costs considered need to be
    # re-assessed and account for missing data.]
    @b.Expression(doc="Total operating costs for commitment block in $")
    def operatingCostCommitment(b):
        # [ESR Note: This term includes the op cost for each
        # 15-min dispatch period.]
        op_cost_dispatch = sum(
            b.dispatchPeriod[disp_per].operatingCostDispatch  # in $
            for disp_per in b.dispatchPeriods
        )
        # [ESR: Assuming we are paying for the full capacity of our
        # generator and should be included to have consistent units.]
        if m.config["include_commitment"]:
            op_cost_gen_state = sum(
                m.fixedCost[gen]
                * b.commitmentPeriodLength
                * m.thermalCapacity[gen]
                * (
                    b.genOn[gen].indicator_var.get_associated_binary()
                    + b.genShutdown[gen].indicator_var.get_associated_binary()
                    + b.genStartup[gen].indicator_var.get_associated_binary()
                )
                for gen in m.thermalGenerators
            )
            op_cost_gen_startup = sum(
                m.startupCost[gen]
                * b.genStartup[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )

            return op_cost_dispatch + op_cost_gen_state + op_cost_gen_startup
        else:
            op_cost_gen_state = sum(
                m.fixedCost[gen]
                * b.commitmentPeriodLength
                * m.thermalCapacity[gen]
                * b.genOn[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
            return op_cost_dispatch + op_cost_gen_state  # ESR: Added op_cost_gen_cost

    @b.Expression(doc="Total curtailment for commitment block in MW")
    def renewableCurtailmentCommitment(b):
        return sum(
            b.dispatchPeriod[disp_per].renewableCurtailmentDispatch
            for disp_per in b.dispatchPeriods
        )

    # Compute Battery Storage cost per dispatch period if storage is
    # needed
    if m.config["storage"]:
        stor.add_commitment_storage_constraints(b)


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

    # [BLN: Convert this to a constraint using operatingCostInvestment
    # Var. May also need to move it.]

    @b.Constraint(doc="Operating costs for investment period")
    def operatingCostInvestment_constraint(b):
        return b.operatingCostInvestment == (
            m.investmentFactor[investment_stage]
            * sum(
                sum(
                    m.weights[rep_per]
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .operatingCostCommitment
                    for com_per in b.representativePeriod[rep_per].commitmentPeriods
                )
                for rep_per in b.representativePeriods
            )
        )
