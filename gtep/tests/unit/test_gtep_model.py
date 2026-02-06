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


from os.path import abspath, join, dirname
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
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing
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


# Helper functions
def read_debug_model():
    curr_dir = dirname(abspath(__file__))
    debug_data_path = abspath(join(curr_dir, "..", "..", "data", "5bus"))
    dataObject = ExpansionPlanningData()
    dataObject.load_prescient(debug_data_path)
    return dataObject


@pytest.mark.usefixtures("patch_unit_handlers")
class TestGTEP(unittest.TestCase):
    def test_model_init(self):
        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with default values, including building a Pyomo.ConcreteModel object
        data_object = read_debug_model()
        modObject = ExpansionPlanningModel(data=data_object)
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        modObject.create_model()
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 1)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 3)
        self.assertEqual(modObject.len_reps, 24)
        self.assertEqual(modObject.num_commit, 24)
        self.assertEqual(modObject.num_dispatch, 4)

        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with non-default values
        modObject = ExpansionPlanningModel(
            data=data_object,
            stages=2,
            num_reps=4,
            len_reps=16,
            num_commit=12,
            num_dispatch=12,
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        modObject.create_model()
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 2)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 4)
        self.assertEqual(modObject.len_reps, 16)
        self.assertEqual(modObject.num_commit, 12)
        self.assertEqual(modObject.num_dispatch, 12)

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
        data_object = read_debug_model()
        modObject = ExpansionPlanningModel(
            data=data_object,
            stages=2,
            num_reps=2,
            len_reps=2,
            num_commit=2,
            num_dispatch=2,
        )
        modObject.create_model()
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
        data_object = read_debug_model()
        modObject = ExpansionPlanningModel(
            data=data_object, num_reps=1, len_reps=1, num_commit=1, num_dispatch=1
        )
        modObject.create_model()

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        TransformationFactory("gdp.bound_pretransformation").apply_to(modObject.model)
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        # previous successful objective values: 9207.95, 6078.86, 531860.15, 531883.43
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 531883.43, places=1
        )
        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_no_investment(self):
        # Solve the debug model with no investment
        data_object = read_debug_model()
        modObject = ExpansionPlanningModel(
            data=data_object, num_reps=1, len_reps=1, num_commit=1, num_dispatch=1
        )
        modObject.config["include_investment"] = False
        modObject.create_model()

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(modObject.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        TransformationFactory("gdp.bound_pretransformation").apply_to(modObject.model)
        TransformationFactory("gdp.bigm").apply_to(modObject.model)

        modObject.results = opt.solve(modObject.model)

        # previous successful objective values: 531860.15, 531883.43
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 531883.43, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data(self):
        # Test more ExpansionPlanningModel with additional cost data
        # This model originated from driver_esr.py
        data_object = read_debug_model()

        curr_dir = dirname(abspath(__file__))
        data_path = abspath(join(curr_dir, "..", "..", "data", "costs"))
        bus_data_path = abspath(join(data_path, "Bus_data_gen_weights_mappings.csv"))
        cost_data_path = abspath(
            join(
                data_path,
                "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx",
            )
        )
        candidate_gens = [
            "Natural Gas_CT",
            "Natural Gas_FE",
            "Solar - Utility PV",
            "Land-Based Wind",
        ]

        data_processing_object = DataProcessing()
        data_processing_object.load_gen_data(
            bus_data_path=bus_data_path,
            cost_data_path=cost_data_path,
            candidate_gens=candidate_gens,
        )

        # Populate and create GTEP model
        mod_object = ExpansionPlanningModel(
            stages=2,
            data=data_object,
            cost_data=data_processing_object,
            num_reps=2,
            len_reps=1,
            num_commit=6,
            num_dispatch=4,
            duration_dispatch=15,
        )

        mod_object.config["include_investment"] = True
        mod_object.config["include_commitment"] = True
        mod_object.config["include_redispatch"] = True
        mod_object.config["scale_loads"] = True
        mod_object.config["transmission"] = True
        mod_object.config["storage"] = False
        mod_object.config["flow_model"] = "DC"

        mod_object.create_model()

        # Check for consistent units
        # Note: Need to do this check before applying the GDP transformations
        assert_units_consistent(mod_object.model)

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        # Apply transformations to logical terms
        TransformationFactory("gdp.bigm").apply_to(mod_object.model)

        mod_object.results = opt.solve(mod_object.model)

        # previous successful objective values: 1524581869.89
        self.assertAlmostEqual(
            value(mod_object.model.total_cost_objective_rule), 1524581869.89, places=1
        )

        assert_units_equivalent(mod_object.model.total_cost_objective_rule.expr, u.USD)
