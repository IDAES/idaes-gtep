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
"""
Unit tests for the function to create the objective function component
"""


import pytest
import pyomo.environ as pyo
from gtep.model_library.objective import create_objective_function


# Mock investmentStage
class StageData:
    def __init__(
        self,
        op_cost=0,
        storage_cost=0,
        invest_cost=0,
        quota_deficit=0,
        renewable_curtailment=0,
        expansion_cost=0,
    ):
        self.operatingCostInvestment = op_cost
        self.storageCostInvestment = storage_cost
        self.investment_cost = invest_cost
        self.quotaDeficit = quota_deficit
        self.renewableCurtailmentInvestment = renewable_curtailment
        self.expansionCost = expansion_cost


# helper function
def create_test_model(
    stages=[1, 2],
    op_cost=[0, 0],
    storage_cost=[0, 0],
    invest_cost=None,
    quota_deficit=[0, 0],
    renewable_curtailment=[0, 0],
    expansion_cost=None,
    deficit_penalty={1: 0, 2: 0},
    investment_factor={1: 0, 2: 0},
    storage=True,
):
    m = pyo.ConcreteModel()
    m.stages = stages

    investment_stage = {}
    if invest_cost is not None:
        for i, stage in enumerate(stages):
            investment_stage[stage] = StageData(
                op_cost=op_cost[i],
                storage_cost=storage_cost[i],
                invest_cost=invest_cost[i],
                quota_deficit=quota_deficit[i],
                renewable_curtailment=renewable_curtailment[i],
            )
    elif expansion_cost is not None:
        for i, stage in enumerate(stages):
            investment_stage[stage] = StageData(
                op_cost=op_cost[i],
                storage_cost=storage_cost[i],
                quota_deficit=quota_deficit[i],
                renewable_curtailment=renewable_curtailment[i],
                expansion_cost=expansion_cost[i],
            )

    m.investmentStage = investment_stage

    m.deficitPenalty = deficit_penalty
    m.investmentFactor = investment_factor

    # Config with storage True
    m.config = {"storage": storage}

    return m


# ----------------------------------UNIT TESTS-----------------------------------------#
def test_objective_multiple_stages_storage_true():
    m = create_test_model(
        op_cost=[100, 150],
        storage_cost=[50, 75],
        invest_cost=[200, 300],
        quota_deficit=[10, 20],
        renewable_curtailment=[5, 10],
        deficit_penalty={1: 2, 2: 3},
        investment_factor={1: 1.5, 2: 2.0},
    )

    create_objective_function(m)

    # Check that the objective was added
    assert hasattr(m, "total_cost_objective_rule")

    # expected costs
    operating_cost = 100 + 150
    storage_cost = 50 + 75
    expansion_cost = 200 + 300
    penalty_cost = sum(
        m.deficitPenalty[stage]
        * m.investmentFactor[stage]
        * m.investmentStage[stage].quotaDeficit
        + m.investmentStage[stage].renewableCurtailmentInvestment
        for stage in m.stages
    )

    assert m.operatingCost == operating_cost
    assert m.storageCost == storage_cost
    assert m.expansionCost == expansion_cost
    assert m.penaltyCost == penalty_cost

    expected_total = operating_cost + storage_cost + expansion_cost + penalty_cost

    # Evaluate the objective expression
    expr = pyo.value(m.total_cost_objective_rule.expr)
    assert expr == pytest.approx(expected_total)


def test_objective_multiple_stages_storage_false():
    m = create_test_model(
        op_cost=[100, 150],
        storage_cost=[50, 75],
        invest_cost=[200, 300],
        quota_deficit=[10, 20],
        renewable_curtailment=[5, 10],
        deficit_penalty={1: 2, 2: 3},
        investment_factor={1: 1.5, 2: 2.0},
        storage=False,
    )

    create_objective_function(m)

    # Check that the objective was added
    assert hasattr(m, "total_cost_objective_rule")
    assert m.storageCost == 0

    # expected costs
    operating_cost = 100 + 150
    expansion_cost = 200 + 300
    penalty_cost = sum(
        m.deficitPenalty[stage]
        * m.investmentFactor[stage]
        * m.investmentStage[stage].quotaDeficit
        + m.investmentStage[stage].renewableCurtailmentInvestment
        for stage in m.stages
    )

    assert m.operatingCost == operating_cost
    assert m.expansionCost == expansion_cost
    assert m.penaltyCost == penalty_cost

    expected_total = operating_cost + expansion_cost + penalty_cost

    # Evaluate the objective expression
    expr = pyo.value(m.total_cost_objective_rule.expr)
    assert expr == pytest.approx(expected_total)


def test_objective_single_stage():
    m = create_test_model(
        stages=[1],
        op_cost=[100],
        expansion_cost=[50],
        invest_cost=[200],
        quota_deficit=[10],
        renewable_curtailment=[5],
        deficit_penalty={1: 2},
        investment_factor={1: 1.5},
        storage=False,
    )

    create_objective_function(m)

    assert hasattr(m, "total_cost_objective_rule")

    expected_total = (
        m.investmentStage[1].operatingCostInvestment
        + m.investmentStage[1].storageCostInvestment
        + m.investmentStage[1].expansionCost
        + m.deficitPenalty[1]
        * m.investmentFactor[1]
        * m.investmentStage[1].quotaDeficit
        + m.investmentStage[1].renewableCurtailmentInvestment
        + m.investmentStage[1].storageCostInvestment
    )

    expr = pyo.value(m.total_cost_objective_rule.expr)
    assert expr == pytest.approx(expected_total)


def test_storage_cost_zero_when_storage_false():
    m = create_test_model(
        op_cost=[100, 150],
        storage_cost=[50, 75],
        invest_cost=[200, 300],
        quota_deficit=[10, 20],
        renewable_curtailment=[5, 10],
        deficit_penalty={1: 2, 2: 3},
        investment_factor={1: 1.5, 2: 2.0},
        storage=False,
    )

    create_objective_function(m)

    # storageCost should be zero
    assert m.storageCost == 0

    # Check objective expression does not include storageCost
    expr = pyo.value(m.total_cost_objective_rule.expr)
    operating_cost = sum(
        m.investmentStage[stage].operatingCostInvestment for stage in m.stages
    )
    expansion_cost = sum(m.investmentStage[stage].investment_cost for stage in m.stages)
    penalty_cost = sum(
        m.deficitPenalty[stage]
        * m.investmentFactor[stage]
        * m.investmentStage[stage].quotaDeficit
        + m.investmentStage[stage].renewableCurtailmentInvestment
        for stage in m.stages
    )
    expected_total = operating_cost + expansion_cost + penalty_cost

    assert expr == pytest.approx(expected_total)


def test_multiple_stages_zero_costs():
    m = create_test_model(invest_cost=[0, 0])
    create_objective_function(m)

    assert hasattr(m, "total_cost_objective_rule")

    expr_val = pyo.value(m.total_cost_objective_rule.expr)
    assert expr_val == pytest.approx(0)


def test_single_stage_zero_costs():
    m = create_test_model(stages=[1], expansion_cost=[0])
    create_objective_function(m)

    assert hasattr(m, "total_cost_objective_rule")

    expr_val = pyo.value(m.total_cost_objective_rule.expr)
    assert expr_val == pytest.approx(0)


def test_many_stages_stress():
    num_stages = 20
    stages = list(range(1, num_stages + 1))
    m = create_test_model(
        stages=stages,
        op_cost=[item * 1.0 for item in stages],
        storage_cost=[item * 0.5 for item in stages],
        invest_cost=[item * 2.0 for item in stages],
        quota_deficit=[item * 0.1 for item in stages],
        renewable_curtailment=[item * 0.05 for item in stages],
        deficit_penalty={item: 1.0 for item in stages},
        investment_factor={item: 1.0 for item in stages},
    )

    create_objective_function(m)

    assert hasattr(m, "total_cost_objective_rule")

    # Calculate expected total cost manually
    operating_cost = sum(
        m.investmentStage[stage].operatingCostInvestment for stage in m.stages
    )
    storage_cost = sum(
        m.investmentStage[stage].storageCostInvestment for stage in m.stages
    )
    expansion_cost = sum(m.investmentStage[stage].investment_cost for stage in m.stages)
    penalty_cost = sum(
        m.deficitPenalty[stage]
        * m.investmentFactor[stage]
        * m.investmentStage[stage].quotaDeficit
        + m.investmentStage[stage].renewableCurtailmentInvestment
        for stage in m.stages
    )

    expected_total = operating_cost + storage_cost + expansion_cost + penalty_cost

    expr_val = pyo.value(m.total_cost_objective_rule.expr)
    assert expr_val == pytest.approx(expected_total)
