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

"""Defines all relevant variables and constraints for the investment
stage in the Generation and Transmission Expansion Planning (GTEP)
model.

"""

import pyomo.environ as pyo
from pyomo.environ import units as u

import gtep.model_library.gen as gens
import gtep.model_library.storage as stor
import gtep.model_library.transmission as transm
import gtep.model_library.commitment as commit


def add_investment_params_and_variables(b, investment_stage):
    """This method defines variables and disjuncts to the investment
    stage block.

    :param b: Investment block
    :param investment_stage: Investment stage index or name

    """

    m = b.model()

    # Add investment parameters.
    b.maxThermalInvestment = pyo.Param(m.regions, default=1000, units=u.MW)
    b.maxRenewableInvestment = pyo.Param(m.regions, default=1000, units=u.MW)

    # Add variables to track and accumulate costs and penalties
    b.quotaDeficit = pyo.Var(within=pyo.NonNegativeReals, initialize=0, units=u.MW)
    b.expansionCost = pyo.Var(within=pyo.NonNegativeReals, initialize=0, units=u.USD)

    # [TODO: Add this only when storage is included. Commented
    # condition for now since it causes an error during testing.]
    # if m.config["storage"]:
    #     stor.add_investment_storage_variables(b)
    b.storageCostInvestment = pyo.Var(
        within=pyo.NonNegativeReals, initialize=0, units=u.USD
    )

    if m.config["include_commitment"]:
        commit.add_investment_commitment_variables(b)


def add_investment_disjuncts(b):
    """This method adds the status disjuncts for thermal and renewable
    generators and transmission lines and storage, when
    needed. Options are: operational, installed, retired, disabled, or
    extended.

    """

    m = b.model()

    #  NOTE: The disjuncts for generators should always be included in
    #  the model.
    gens.add_generators_status_disjuncts(b, m.thermalGenerators, m.renewableGenerators)

    if m.config["storage"]:
        stor.add_storage_status_disjuncts(b, m.storage)

    if m.config["transmission"]:
        transm.add_transmission_status_disjuncts(b, m.transmission)


def add_investment_constraints(
    b,
    investment_stage,
):
    """This method adds standard inequalities (i.e., those not
    involving disjunctions) to the investment stage block.

    """

    m = b.model()

    # Add constraints for generators and transmission lines and
    # storage, when needed. These constraints include fixing indicator
    # variables based on "in_service" data for the units and the
    # calculation of investment cost.
    gens.add_investment_generators_constraints(m, b, investment_stage)

    if m.config["transmission"]:
        transm.add_investment_transmission_constraints(m, b, investment_stage)

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
            b.generators_investment_cost
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
    # Var. May also need to move it.]
    @b.Constraint(doc="Operating costs for investment period")
    def operatingCostInvestment_constraint(b):
        return b.operatingCostInvestment == (
            b.commitmentOperatingCostInvestment
            if m.config["include_commitment"] == True
            else 0
        )

    # [ESR Question: Do we need to add the flag for investment when
    # inside investment? Should we better put this flag when calling
    # this function?]
    if m.config["include_investment"]:

        # NOTE: This is constraint (13) in reference [1]
        @b.Constraint(
            doc="Minimum per-stage renewable generation requirement",
        )
        def renewable_generation_requirement(b):
            renewableSurplusRepresentative = 0
            # [TODO: preprocess loads for the appropriate sum here.]
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
