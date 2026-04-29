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

import os
import json

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
def read_debug_model(
    stages=1,
    num_reps=3,
    num_commit=24,
    num_dispatch=4,
    duration_dispatch=15,
):
    curr_dir = dirname(abspath(__file__))
    debug_data_path = abspath(join(curr_dir, "..", "..", "data", "5bus"))
    dataObject = ExpansionPlanningData(
        stages=stages,
        num_reps=num_reps,
        num_commit=num_commit,
        num_dispatch=num_dispatch,
        duration_dispatch=duration_dispatch,
    )
    dataObject.load_prescient(debug_data_path)
    return dataObject


def prepare_model_and_cost_data(
    stages=1,
    num_reps=3,
    num_commit=24,
    num_dispatch=4,
    duration_dispatch=15,
):
    # Prepare model and cost data
    dataObject = read_debug_model(
        stages,
        num_reps,
        num_commit,
        num_dispatch,
        duration_dispatch,
    )
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

    dataProcessingObject = DataProcessing()
    dataProcessingObject.load_gen_data(
        bus_data_path=bus_data_path,
        cost_data_path=cost_data_path,
        candidate_gens=candidate_gens,
    )
    return dataObject, dataProcessingObject


@pytest.mark.usefixtures("patch_unit_handlers")
class TestGTEP(unittest.TestCase):
    def test_model_init(self):
        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with default values, including building a Pyomo.ConcreteModel object
        data_object = read_debug_model(
            stages=1,
            num_reps=3,
            num_commit=6,
            num_dispatch=4,
            duration_dispatch=60,
        )
        modObject = ExpansionPlanningModel(data=data_object)
        self.assertIsInstance(modObject, ExpansionPlanningModel)

        modObject.create_model()

        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 1)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 3)
        self.assertEqual(modObject.num_commit, {1: 6, 2: 6, 3: 6})
        self.assertEqual(
            modObject.num_dispatch,
            {
                1: {1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4},
                2: {1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4},
                3: {1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4},
            },
        )
        self.assertEqual(
            modObject.duration_representative_period, {1: 24, 2: 24, 3: 24}
        )

        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with non-default values
        data_object = read_debug_model(
            stages=2,
            num_reps=4,
            num_commit=12,
            num_dispatch=12,
            duration_dispatch=30,
        )
        modObject = ExpansionPlanningModel(
            data=data_object,
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        modObject.create_model()
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 2)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 4)
        self.assertEqual(modObject.num_commit, {1: 12, 2: 12, 3: 12, 4: 12})
        self.assertEqual(
            modObject.num_dispatch,
            {
                1: {
                    1: 12,
                    2: 12,
                    3: 12,
                    4: 12,
                    5: 12,
                    6: 12,
                    7: 12,
                    8: 12,
                    9: 12,
                    10: 12,
                    11: 12,
                    12: 12,
                },
                2: {
                    1: 12,
                    2: 12,
                    3: 12,
                    4: 12,
                    5: 12,
                    6: 12,
                    7: 12,
                    8: 12,
                    9: 12,
                    10: 12,
                    11: 12,
                    12: 12,
                },
                3: {
                    1: 12,
                    2: 12,
                    3: 12,
                    4: 12,
                    5: 12,
                    6: 12,
                    7: 12,
                    8: 12,
                    9: 12,
                    10: 12,
                    11: 12,
                    12: 12,
                },
                4: {
                    1: 12,
                    2: 12,
                    3: 12,
                    4: 12,
                    5: 12,
                    6: 12,
                    7: 12,
                    8: 12,
                    9: 12,
                    10: 12,
                    11: 12,
                    12: 12,
                },
            },
        )

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
        data_object = read_debug_model(
            stages=2,
            num_reps=2,
            num_commit=2,
            num_dispatch=2,
        )
        modObject = ExpansionPlanningModel(
            data=data_object,
        )
        modObject.create_model()
        m = modObject.model

        # Check for consistent units
        assert_units_consistent(m)

        # Check that subset of model components have expected units
        m_inv = m.investmentStage[1]
        m_rep_period = m_inv.representativePeriod[1]
        m_commit = m_rep_period.commitmentPeriod[1]
        m_disp = m_commit.dispatchPeriod[1]

        assert_units_equivalent(m_rep_period.representativePeriodLength, u.h)
        assert_units_equivalent(m_commit.commitmentPeriodLength, u.h)
        assert_units_equivalent(m_disp.dispatchPeriodLength, u.min)
        assert_units_equivalent(m.renewable_capacity_enforcement[1, "10_PV"].expr, u.MW)
        assert_units_equivalent(m_inv.renewable_curtailment_cost.expr, u.USD)
        assert_units_equivalent(m_disp.flow_balance["bus1"].expr, u.MW)
        assert_units_equivalent(m.rampUpRates, u.dimensionless)
        assert_units_equivalent(m.varCost, u.USD / u.h / u.MW)
        assert_units_equivalent(m_disp.spinningReserve, u.MW)
        assert_units_equivalent(
            m_commit.genOn["3_CT"].operating_limit_min[1].expr,
            u.MW,
        )
        assert_units_equivalent(
            m_commit.genOn["4_STEAM"].max_spinning_reserve[1, "4_STEAM"].expr,
            u.MW,
        )

    def test_solve_bigm(self):
        # Solve the debug model as is
        data_object = read_debug_model(num_reps=1, num_commit=1, num_dispatch=1)
        modObject = ExpansionPlanningModel(data=data_object)
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
        data_object = read_debug_model(
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
        )
        modObject = ExpansionPlanningModel(
            data=data_object,
        )

        modObject = ExpansionPlanningModel(data=data_object)
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

    def test_with_cost_data_and_commitment(self):
        # Test ExpansionPlanningModel with cost data
        # This model originated from driver_esr.py
        dataObject, dataProcessingObject = prepare_model_and_cost_data(
            stages=2,
            num_reps=2,
            num_commit=6,
            num_dispatch=4,
            duration_dispatch=15,
        )

        # Populate and create GTEP model
        modObject = ExpansionPlanningModel(
            data=dataObject,
            cost_data=dataProcessingObject,
        )

        modObject.config["include_investment"] = True
        modObject.config["include_commitment"] = True
        modObject.config["include_redispatch"] = True
        modObject.config["scale_loads"] = True
        modObject.config["transmission"] = True
        modObject.config["storage"] = False
        modObject.config["flow_model"] = "DC"

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

        # previous successful objective values: 1524581869.89
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 1524533561.02, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_no_commitment(self):
        # Test ExpansionPlanningModel with cost data and no commitment
        # This model originated from driver_esr.py
        dataObject, dataProcessingObject = prepare_model_and_cost_data(
            stages=2,
            num_reps=2,
            num_commit=6,
            num_dispatch=4,
            duration_dispatch=15,
        )

        # Populate and create GTEP model
        modObject = ExpansionPlanningModel(
            data=dataObject,
            cost_data=dataProcessingObject,
        )

        modObject.config["include_investment"] = True
        modObject.config["include_commitment"] = False
        modObject.config["include_redispatch"] = True
        modObject.config["scale_loads"] = True
        modObject.config["transmission"] = True
        modObject.config["storage"] = False
        modObject.config["flow_model"] = "DC"

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
            value(modObject.model.total_cost_objective_rule), 1524533561.02, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_period_structure_from_scalars_and_json(self):
        # Test with scalar/list arguments (all periods same)
        dataObject, dataProcessingObject = prepare_model_and_cost_data(
            num_reps=2,
            num_commit=3,
            num_dispatch=4,
            duration_dispatch=15,
        )

        modObject = ExpansionPlanningModel(
            data=dataObject,
            cost_data=dataProcessingObject,
            duration_representative_period=24,
            duration_commitment=1,
            save_period_structure_file=False,
            period_structure_json_file=None,
        )

        # Check that all values are as expected (all periods same)
        self.assertEqual(modObject.num_commit[1], 3)
        self.assertEqual(modObject.num_dispatch[2][3], 4)
        self.assertEqual(modObject.duration_commitment[1][2], 1)
        self.assertEqual(modObject.duration_dispatch[2][3][4], 15)

        # Test custom period structure with irregular values. This
        # dictionary is saved as a .json file and then used to
        # initialize the ExpansionPlanningModel class.
        period_dict = {
            "number_representative": 2,
            "number_commitment": {1: 2, 2: 3},
            "number_dispatch": {1: {1: 3, 2: 2}, 2: {1: 2, 2: 3, 3: 2}},
            "duration_representative_period": {1: 24, 2: 18},
            "duration_commitment": {1: {1: 1, 2: 2}, 2: {1: 1, 2: 1.5, 3: 2}},
            "duration_dispatch": {
                1: {1: {1: 10, 2: 20, 3: 30}, 2: {1: 30, 2: 90}},
                2: {1: {1: 30, 2: 30}, 2: {1: 20, 2: 20, 3: 50}, 3: {1: 60, 2: 60}},
            },
        }
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        json_path = os.path.join(curr_dir, "test_custom_period_structure.json")
        with open(json_path, "w") as f:
            json.dump(period_dict, f, indent=2)

        # Test that the model correctly reads and assigns the custom
        # period structure values. Here we instantiate the model using
        # the .json file.
        modObject = ExpansionPlanningModel(
            data=dataObject,
            cost_data=dataProcessingObject,
            period_structure_json_file=json_path,
            save_period_structure_file=False,
        )

        # Assert that we have the correct reading of the structure
        self.assertEqual(modObject.num_reps, 2)
        self.assertEqual(modObject.num_commit[2], 3)
        self.assertEqual(modObject.num_dispatch[2][2], 3)
        self.assertEqual(modObject.duration_commitment[2][2], 1.5)
        self.assertEqual(modObject.duration_dispatch[2][2][3], 50)
        self.assertEqual(modObject.duration_representative_period[2], 18)
        self.assertEqual(modObject.duration_dispatch[1][1][2], 20)
        self.assertEqual(modObject.duration_commitment[1][2], 2)

        # Remove the .json file after the test
        os.remove(json_path)
