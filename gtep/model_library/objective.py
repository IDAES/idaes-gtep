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
    renewable curtailment costs.

    :param m: Pyomo GTEP model.

    """

    # operatingCostTotal includes generator operating costs from the
    # commitment stage, calculated using fixed and variable cost
    # terms. NOTE: the operatingCostInvestment is weighted by the
    # representative period weights.
    @m.Expression()
    def operatingCostTotal(m):
        return sum(
            m.investmentStage[stage].operatingCostInvestment for stage in m.stages
        )

    # expansionCostTotal includes capital investment costs for
    # generators, storage (if enabled), and transmission (if
    # enabled). Generator investment costs include both thermal and
    # renewable resources.
    @m.Expression()
    def expansionCostTotal(m):
        return sum(m.investmentStage[stage].investment_cost for stage in m.stages)

    # renewableCurtailmentCostTotal includes curtailment costs from
    # the representative-period operational model. These costs are
    # weighted by representative-period weights through
    # renewableCurtailmentInvestment.
    @m.Expression()
    def renewableCurtailmentCostTotal(m):
        return sum(
            m.investmentStage[stage].renewableCurtailmentInvestment
            for stage in m.stages
        )

    @m.Objective()
    def total_cost_objective_rule(m):
        return (
            m.operatingCostTotal
            + m.expansionCostTotal
            + m.renewableCurtailmentCostTotal
        )
