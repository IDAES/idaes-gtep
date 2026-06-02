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

"""Defines the objective function component in the Generation and
Transmission Expansion Planning (GTEP) Model

"""


def create_objective_function(m):
    """This method defines the objective function for the GTEP model
    as the minimization of the total cost. The total cost is
    calculated as the sum of operating costs, expansion costs, and
    penalty costs (which account for generation deficits, renewable
    quota deficits, and curtailment).

    :param m: Pyomo GTEP model.

    """

    # NOTE: We add battery storage cost only when "storage" is set to
    # True in the configuration input, otherwise its cost value is 0.
    @m.Expression()
    def operatingCostTotal(m):
        return sum(
            m.investmentStage[stage].operatingCostInvestment for stage in m.stages
        )

    # [ESR: Comment this expression since storageCostInvestment is
    # declared but not calculated in the model.]
    # @m.Expression()
    # def storageCostTotal(m):
    #     if m.config["storage"]:
    #         return sum(
    #             m.investmentStage[stage].storageCostInvestment for stage in m.stages
    #         )
    #     else:
    #         return 0

    # [ESR: expansionCostTotal already includes the storage_investment_cost
    # when storage is True]
    @m.Expression()
    def expansionCostTotal(m):
        return sum(m.investmentStage[stage].investment_cost for stage in m.stages)

    @m.Expression()
    def penaltyCostTotal(m):
        return sum(
            m.deficitPenalty[stage]
            * m.investmentFactor[stage]
            * m.investmentStage[stage].quotaDeficit
            + m.investmentStage[stage].renewableCurtailmentInvestment
            for stage in m.stages
        )

    @m.Objective()
    def total_cost_objective_rule(m):
        return (
            m.operatingCostTotal
            + m.expansionCostTotal
            + m.penaltyCostTotal
            # + m.storageCostTotal
        )
