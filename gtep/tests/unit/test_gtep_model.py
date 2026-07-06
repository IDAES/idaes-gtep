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


import pytest
import pyomo.common.unittest as unittest
from pyomo.environ import ConcreteModel, SolverFactory, value, LogicalConstraint
from pyomo.environ import units as u
from pyomo.environ import TransformationFactory
from pyomo.util.check_units import (
    assert_units_consistent,
    _component_data_handlers,
    assert_units_equivalent,
)
from gtep.gtep_model import ExpansionPlanningModel
from gtep.tests.unit.utils_for_testing import create_model
from egret.data.model_data import ModelData


@pytest.fixture
def patch_unit_handlers():
    # This is a temporary hack to get the unit consistency tests working.
    # This disables units checking on LogicalConstraints. This code can be
    # deleted once a units checker for LogicalConstraints is added to Pyomo
    if LogicalConstraint in _component_data_handlers:
        print(
            "Warning: Found a unit checker for LogicalConstraints. This testing hack "
            "can likely be removed."
        )
        yield
    else:
        _component_data_handlers[LogicalConstraint] = None
        yield
        _component_data_handlers.pop(LogicalConstraint)


@pytest.mark.usefixtures("patch_unit_handlers")
class TestGTEP(unittest.TestCase):
    def test_model_init(self):
        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with default values, including building a Pyomo.ConcreteModel object
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 3,
                "len_reps": 24,
                "num_commit": 24,
                "num_dispatch": 4,
                "duration_dispatch": 60,
            },
            include_cost_data=False,
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 1)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 3)
        self.assertEqual(modObject.len_reps, 24)
        self.assertEqual(modObject.num_commit, 24)
        self.assertEqual(modObject.num_dispatch, 4)
        self.assertEqual(modObject.duration_dispatch, 60)

        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with non-default values
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 4,
                "len_reps": 16,
                "num_commit": 12,
                "num_dispatch": 12,
                "duration_dispatch": 30,
            },
            include_cost_data=False,
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 2)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 4)
        self.assertEqual(modObject.len_reps, 16)
        self.assertEqual(modObject.num_commit, 12)
        self.assertEqual(modObject.num_dispatch, 12)
        self.assertEqual(modObject.duration_dispatch, 30)

        # We have expansion blocks and they are where and what we think they are
        expansion_blocks = modObject.model.component("investmentStage")
        self.assertEqual(len(expansion_blocks), 2)
        self.assertIs(modObject.model.investmentStage, expansion_blocks)

        # Each expansion block has 4 representative period blocks also, and they make sense
        # Well, at least the first expansion block

        representative_blocks_1 = expansion_blocks[1].component("representativePeriod")
        self.assertEqual(len(representative_blocks_1), 4)
        for p in modObject.model.representativePeriods:
            self.assertIs(
                expansion_blocks[1].representativePeriod[p], representative_blocks_1[p]
            )

        # First representative block should have commitment blocks...
        # As many as it has commitment periods assigned...
        # These may be not the most ideal comparisons
        commitment_block_q = representative_blocks_1[1].component("commitmentPeriod")
        self.assertEqual(
            len(commitment_block_q), len(representative_blocks_1[1].commitmentPeriods)
        )

        dispatch_block_1 = commitment_block_q[1].component("dispatchPeriod")
        self.assertEqual(
            len(dispatch_block_1), len(commitment_block_q[1].dispatchPeriods)
        )

    def test_model_unit_consistency(self):
        # Test that the ExpansionPlanningModel has consistent units and spot check that
        # components have their expected units
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "len_reps": 2,
                "num_commit": 2,
                "num_dispatch": 2,
            },
            include_cost_data=False,
        )
        m = modObject.model

        # Check for consistent units
        assert_units_consistent(m)

        # Check that subset of model components have expected units
        assert_units_equivalent(m.renewable_capacity_enforcement[1, "10_PV"].expr, u.MW)
        assert_units_equivalent(
            m.investmentStage[1].renewable_curtailment_cost.expr, u.USD
        )
        assert_units_equivalent(
            m.investmentStage[1]
            .representativePeriod[1]
            .commitmentPeriod[1]
            .dispatchPeriod[1]
            .flow_balance["bus1"]
            .expr,
            u.MW,
        )
        assert_units_equivalent(m.commitmentPeriodLength, u.h)
        assert_units_equivalent(m.dispatchPeriodLength, u.min)
        assert_units_equivalent(m.rampUpRates, u.dimensionless)
        assert_units_equivalent(m.varCost, u.USD / u.h / u.MW)
        assert_units_equivalent(
            m.investmentStage[1]
            .representativePeriod[1]
            .commitmentPeriod[1]
            .dispatchPeriod[1]
            .spinningReserve,
            u.MW,
        )
        assert_units_equivalent(
            m.investmentStage[1]
            .representativePeriod[1]
            .commitmentPeriod[1]
            .genOn["3_CT"]
            .operating_limit_min[1]
            .expr,
            u.MW,
        )
        assert_units_equivalent(
            m.investmentStage[1]
            .representativePeriod[1]
            .commitmentPeriod[1]
            .genOn["4_STEAM"]
            .max_spinning_reserve[1, "4_STEAM"]
            .expr,
            u.MW,
        )

    def test_solve_bigm(self):
        # Solve the debug model as is
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "len_reps": 1,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_dispatch": 15,
            },
            prescient_data_args={
                "representative_dates": ["2020-01-28 00:00"],
            },
            include_cost_data=False,
        )

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        TransformationFactory("gdp.bound_pretransformation").apply_to(modObject.model)
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        # previous successful objective values: 9207.95, 6078.86, 531860.15, 531883.43, 7977055.4, 7977055.4
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 7977150.30, places=1
        )
        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_no_investment(self):
        # Solve the debug model with no investment
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "len_reps": 1,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_dispatch": 15,
            },
            prescient_data_args={
                "representative_dates": ["2020-01-28 00:00"],
            },
            config={
                "include_investment": False,
            },
            include_cost_data=False,
        )

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        TransformationFactory("gdp.bound_pretransformation").apply_to(modObject.model)
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        # previous successful objective values: 531860.15, 531883.43, 7977055.4, 7977055.4
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 7977150.30, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_commitment(self):
        # Test ExpansionPlanningModel with cost data
        # This model originated from driver_esr.py
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "len_reps": 1,
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_dispatch": 15,
            },
            config={
                "include_investment": True,
                "include_commitment": True,
                "include_redispatch": True,
                "scale_loads": True,
                "transmission": True,
                "storage": False,
                "flow_model": "DC",
            },
            candidate_gens=[
                "Natural Gas_CT",
                "Natural Gas_FE",
                "Solar - Utility PV",
                "Land-Based Wind",
            ],
        )

        modObject.config["include_investment"] = True
        modObject.config["include_commitment"] = True
        modObject.config["include_redispatch"] = True
        modObject.config["scale_loads"] = True
        modObject.config["transmission"] = True
        modObject.config["storage"] = False
        modObject.config["flow_model"] = "DC"
        modObject.config["advanced_hydro"] = False

        modObject.create_model()

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        # Apply transformations to logical terms
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        # previous successful objective values: 1524581869.89, 779334165.7, 779344643.1
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 779486340.91, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_no_commitment(self):

        # This test verifies that the expansion planning model can be
        # built and solved using preprocessed cost data when unit
        # commitment is disabled. The test also checks unit
        # consistency and validates the resulting objective value
        # against an expected benchmark.

        # Test ExpansionPlanningModel with cost data and no commitment
        # This model originated from driver_esr.py
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "len_reps": 1,
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_dispatch": 15,
            },
            config={
                "include_investment": True,
                "include_commitment": False,
                "include_redispatch": True,
                "scale_loads": True,
                "transmission": True,
                "storage": False,
                "flow_model": "DC",
            },
            candidate_gens=[
                "Natural Gas_CT",
                "Natural Gas_FE",
                "Solar - Utility PV",
                "Land-Based Wind",
            ],
        )

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        # Apply transformations to logical terms
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        # previous successful objective values: 1524533561.02, 926187704.4
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 926194856.96, places=1
        )

        modObject.config["include_investment"] = True
        modObject.config["include_commitment"] = False
        modObject.config["include_redispatch"] = True
        modObject.config["scale_loads"] = True
        modObject.config["transmission"] = True
        modObject.config["storage"] = False
        modObject.config["flow_model"] = "DC"
        modObject.config["advanced_hydro"] = False
        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_weights(self):

        # Typically, representative weights should reflect how many
        # actual days each representative day represents. For example,
        # with 2 representative days in a 365-day year, the weights
        # would be approximately [182, 183].  For now, we use
        # simplified test values because the model becomes infeasible
        # with HiGHS when using the full representative-day weights.
        # These values are intended only to verify that the model
        # reads and applies representative weights correctly.
        weights = [10, 12]
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "len_reps": 1,
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_dispatch": 15,
            },
            prescient_data_args={
                "representative_dates": ["2020-01-28 00:00", "2020-04-23 00:00"],
                "representative_weights": weights,
            },
            config={
                "include_investment": True,
                "include_commitment": True,
                "include_redispatch": True,
                "scale_loads": True,
                "transmission": True,
                "storage": False,
                "flow_model": "DC",
            },
            candidate_gens=[
                "Natural Gas_CT",
                "Natural Gas_FE",
                "Solar - Utility PV",
                "Land-Based Wind",
            ],
        )

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        # Apply transformations to logical terms
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 8584301655.08, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_hydro(self):

        # This test verifies that the expansion planning model can be
        # built and solved using preprocessed cost data with advanced
        # hydropower enabled. The test also checks unit consistency
        # and validates the resulting objective value against an
        # expected benchmark.
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "len_reps": 1,
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_dispatch": 15,
            }
        )

        modObject.config["include_investment"] = True
        modObject.config["include_commitment"] = True
        modObject.config["include_redispatch"] = True
        modObject.config["scale_loads"] = True
        modObject.config["transmission"] = True
        modObject.config["storage"] = False
        modObject.config["flow_model"] = "DC"
        modObject.config["advanced_hydro"] = True

        modObject.create_model()

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        # Apply transformations to logical terms
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 779418083.72, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)
        # Check that each representative period in the model is
        # assigned the expected representative weight
        for rep_period, expected_weight in zip(
            modObject.model.representativePeriods, weights
        ):
            self.assertEqual(
                value(modObject.model.weights[rep_period]), expected_weight
            )
