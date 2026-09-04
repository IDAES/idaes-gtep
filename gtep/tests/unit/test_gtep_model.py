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

import json
import pytest
from pathlib import Path

import pyomo.common.unittest as unittest
from pyomo.environ import ConcreteModel, SolverFactory, value, LogicalConstraint
from pyomo.environ import units as u
from pyomo.environ import TransformationFactory
from pyomo.util.check_units import (
    assert_units_consistent,
    _component_data_handlers,
    assert_units_equivalent,
)
from pyomo.common.tempfiles import TempfileManager

from gtep.gtep_model import ExpansionPlanningModel
from gtep.tests.unit.utils_for_testing import create_data, create_model
from egret.data.model_data import ModelData

curr_dir = Path(__file__).resolve().parent
Texas123_case_path = (curr_dir / ".." / ".." / "data" / "123_Bus_Resil_Week").resolve()


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
                "num_commit": 24,
                "num_dispatch": 4,
                "duration_representative_period": 24,
            },
            include_cost_data=False,
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 1)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 3)
        self.assertEqual(modObject.num_commit, {1: 24, 2: 24, 3: 24})
        self.assertEqual(
            modObject.num_dispatch,
            {
                1: {i: 4 for i in range(1, 25)},
                2: {i: 4 for i in range(1, 25)},
                3: {i: 4 for i in range(1, 25)},
            },
        )
        self.assertEqual(
            modObject.duration_representative_period, {1: 24, 2: 24, 3: 24}
        )
        self.assertEqual(
            modObject.duration_commitment,
            {
                1: {i: 1 for i in range(1, 25)},
                2: {i: 1 for i in range(1, 25)},
                3: {i: 1 for i in range(1, 25)},
            },
        )
        self.assertEqual(
            modObject.duration_dispatch,
            {
                1: {com: {disp: 15 for disp in range(1, 5)} for com in range(1, 25)},
                2: {com: {disp: 15 for disp in range(1, 5)} for com in range(1, 25)},
                3: {com: {disp: 15 for disp in range(1, 5)} for com in range(1, 25)},
            },
        )

        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with non-default values
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 4,
                "num_commit": 12,
                "num_dispatch": 12,
                "duration_representative_period": 12,
            },
            include_cost_data=False,
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 2)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.model.md, ModelData)
        self.assertEqual(modObject.num_reps, 4)
        self.assertEqual(modObject.num_commit, {1: 12, 2: 12, 3: 12, 4: 12})
        self.assertEqual(
            modObject.num_dispatch,
            {
                1: {i: 12 for i in range(1, 13)},
                2: {i: 12 for i in range(1, 13)},
                3: {i: 12 for i in range(1, 13)},
                4: {i: 12 for i in range(1, 13)},
            },
        )
        self.assertEqual(
            modObject.duration_representative_period, {1: 12, 2: 12, 3: 12, 4: 12}
        )
        self.assertEqual(
            modObject.duration_commitment,
            {
                1: {i: 1 for i in range(1, 13)},
                2: {i: 1 for i in range(1, 13)},
                3: {i: 1 for i in range(1, 13)},
                4: {i: 1 for i in range(1, 13)},
            },
        )
        self.assertEqual(
            modObject.duration_dispatch,
            {
                1: {com: {disp: 5 for disp in range(1, 13)} for com in range(1, 13)},
                2: {com: {disp: 5 for disp in range(1, 13)} for com in range(1, 13)},
                3: {com: {disp: 5 for disp in range(1, 13)} for com in range(1, 13)},
                4: {com: {disp: 5 for disp in range(1, 13)} for com in range(1, 13)},
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
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "num_commit": 2,
                "num_dispatch": 2,
                "duration_representative_period": 2,
            },
            include_cost_data=False,
        )
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
            m_commit.genOn["4_STEAM"].max_spinning_reserve[1].expr,
        )

    def test_solve_bigm(self):
        # Solve the debug model as is
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_representative_period": 1,
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

        # previous successful objective values: 9207.95, 6078.86, 531860.15, 531883.43, 7977055.4,
        # 7977055.4, 7977150.30, 6986122.88, 7118266.88, 27944303.09
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 28076447.10, places=1
        )
        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_no_investment(self):
        # Solve the debug model with no investment
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_representative_period": 1,
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

        # previous successful objective values: 531860.15, 531883.43, 7977055.4, 7977055.4,
        # 7977150.30, 6986122.88, 7977169.84, 27944303.09
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 31908490.95, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_commitment(self):
        # Test ExpansionPlanningModel with cost data
        # This model originated from driver_esr.py
        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_representative_period": 6,
            },
            config={
                "include_investment": True,
                "include_commitment": True,
                "include_redispatch": True,
                "scale_loads": True,
                "transmission": True,
                "storage": False,
                "flow_model": "DC",
                "advanced_hydro": False,
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

        # previous successful objective values: 1524581869.89, 779334165.7, 779344643.1, 779486340.91,
        # 776765960.70
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 779494096.65, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_no_commitment(self):

        # This test verifies that the expansion planning model can be
        # built and solved using preprocessed cost data when unit
        # commitment is disabled. The test also checks unit
        # consistency and validates the resulting objective value
        # against an expected benchmark.

        modObject = create_model(
            planning_data_args={
                "stages": 2,
                "num_reps": 2,
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_representative_period": 6,
            },
            config={
                "include_investment": True,
                "include_commitment": False,
                "include_redispatch": True,
                "scale_loads": True,
                "transmission": True,
                "storage": False,
                "flow_model": "DC",
                "advanced_hydro": False,
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

        # previous successful objective values: 1524533561.02, 926187704.4, 926194856.96, 923474476.75
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 926202612.69, places=1
        )
        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_with_cost_data_and_weights(self):

        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 12,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_representative_period": 1,
            },
            prescient_data_args={
                "representative_dates": [
                    "2020-01-01 00:00",
                    "2020-02-01 00:00",
                    "2020-03-01 00:00",
                    "2020-04-01 00:00",
                    "2020-05-01 00:00",
                    "2020-06-01 00:00",
                    "2020-07-01 00:00",
                    "2020-08-01 00:00",
                    "2020-09-01 00:00",
                    "2020-10-01 00:00",
                    "2020-11-01 00:00",
                    "2020-12-01 00:00",
                ],
                "representative_weights": {
                    "2020-01-01 00:00": 31,
                    "2020-02-01 00:00": 28,
                    "2020-03-01 00:00": 31,
                    "2020-04-01 00:00": 30,
                    "2020-05-01 00:00": 31,
                    "2020-06-01 00:00": 30,
                    "2020-07-01 00:00": 31,
                    "2020-08-01 00:00": 31,
                    "2020-09-01 00:00": 30,
                    "2020-10-01 00:00": 31,
                    "2020-11-01 00:00": 30,
                    "2020-12-01 00:00": 31,
                },
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

        # Previous values: 8584301655.08, 3545283964.22, 3359204716.88, 3487926967.61, 13435364152.91
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 13573783610.18, places=1
        )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_period_structure_from_scalars(self):
        # Test period structure dictionary created using the provided
        # scalars
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 2,
                "num_commit": 3,
                "num_dispatch": 4,
                "duration_representative_period": 3,
            },
            include_cost_data=True,
        )

        # Assert that all values are as expected
        self.assertEqual(modObject.num_commit[1], 3)
        self.assertEqual(modObject.num_dispatch[2][3], 4)
        self.assertEqual(modObject.duration_representative_period[1], 3)
        self.assertEqual(modObject.duration_commitment[1][2], 1)
        self.assertEqual(modObject.duration_dispatch[2][3][4], 15)

        # Test custom period structure dictionary with irregular
        # values.
        period_dict = {
            "number_representative": 2,
            "number_commitment": {1: 2, 2: 3},
            "number_dispatch": {1: {1: 3, 2: 2}, 2: {1: 2, 2: 3, 3: 2}},
            "duration_representative_period": {1: 24, 2: 18},
            "duration_commitment": {1: {1: 12, 2: 12}, 2: {1: 6, 2: 6, 3: 6}},
            "duration_dispatch": {
                1: {1: {1: 360, 2: 180, 3: 180}, 2: {1: 360, 2: 360}},
                2: {
                    1: {1: 180, 2: 180},
                    2: {1: 120, 2: 120, 3: 120},
                    3: {1: 180, 2: 180},
                },
            },
        }

        with TempfileManager.new_context() as tempfile:
            temp_dir = Path(tempfile.mkdtemp())
            json_path = temp_dir / "test_custom_period_structure.json"
            with open(json_path, "w") as f:
                json.dump(period_dict, f, indent=2)

            # Test that the model correctly reads and assigns the
            # custom period structure values. For this, we instantiate
            # the model using the temp .json file.
            modObject = create_model(
                planning_data_args={
                    "period_structure_json_file": str(json_path),
                },
                include_cost_data=True,
            )

            # Assert that we have the correct reading of the structure
            self.assertEqual(modObject.num_reps, 2)
            self.assertEqual(modObject.num_commit[2], 3)
            self.assertEqual(modObject.num_dispatch[2][2], 3)
            self.assertEqual(modObject.duration_commitment[2][2], 6)
            self.assertEqual(modObject.duration_dispatch[2][2][3], 120)
            self.assertEqual(modObject.duration_representative_period[2], 18)
            self.assertEqual(modObject.duration_dispatch[1][1][2], 180)
            self.assertEqual(modObject.duration_commitment[1][2], 12)

    def test_period_structure_consistency_errors(self):

        # Check that a consistency error is raised if
        # commitment/dispatch period durations do not sum to the
        # representative/commitment period duration when loading from
        # a .json file.  The test first creates and writes the .json
        # file, then checks for the error.
        period_dict = {
            "number_representative": 2,
            "number_commitment": {"1": 6, "2": 6},
            "number_dispatch": {
                "1": {"1": 4, "2": 4, "3": 4, "4": 4, "5": 4, "6": 4},
                "2": {"1": 4, "2": 4, "3": 4, "4": 4, "5": 4, "6": 4},
            },
            "duration_representative_period": {"1": 6, "2": 6},
            "duration_commitment": {
                "1": {"1": 1.0, "2": 1.0, "3": 1.0, "4": 1.0, "5": 1.0, "6": 3.0},
                "2": {"1": 1.0, "2": 1.0, "3": 2.0, "4": 1.0, "5": 1.0, "6": 1.0},
            },
            "duration_dispatch": {
                "1": {
                    "1": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "2": {"1": 25.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "3": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "4": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "5": {"1": 15.0, "2": 15.0, "3": 5.0, "4": 15.0},
                    "6": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                },
                "2": {
                    "1": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "2": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "3": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "4": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "5": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                    "6": {"1": 15.0, "2": 15.0, "3": 15.0, "4": 15.0},
                },
            },
        }

        with TempfileManager.new_context() as tempfile:
            temp_dir = Path(tempfile.mkdtemp())
            json_path = temp_dir / "test_consistency_errors_period_structure.json"
            with open(json_path, "w") as f:
                json.dump(period_dict, f, indent=2)

            # Create the data object using the temporary period
            # structure file, but do not instantiate the model.
            dataObject, dataProcessingObject = create_data(
                planning_data_args={
                    "period_structure_json_file": str(json_path),
                },
                include_cost_data=True,
            )

            # Assert that the sum of commitment durations (20hr) does
            # not match the representative period duration (18hr)
            with self.assertRaises(ValueError) as cm:
                ExpansionPlanningModel(data=dataObject, cost_data=dataProcessingObject)

            self.assertIn(
                "Period structure consistency check failed:\n\nERROR: Found (2) mismatches for commitment period duration:",
                str(cm.exception),
            )

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
                "num_commit": 6,
                "num_dispatch": 4,
                "duration_representative_period": 2,
            },
            config={
                "include_investment": True,
                "include_commitment": True,
                "include_redispatch": True,
                "scale_loads": True,
                "transmission": True,
                "storage": False,
                "flow_model": "DC",
                "advanced_hydro": True,
            },
            include_cost_data=True,
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

        # Previous values: 779418083.72, 767038945.08, 767708848.64
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 255902949.54, places=1
        )
        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)

    def test_Texas123_case(self):

        modObject = create_model(
            input_data_path=Texas123_case_path,
            planning_data_args={
                "stages": 1,
                "num_reps": 4,
                "num_commit": 2,
                "num_dispatch": 1,
                "duration_representative_period": 2,
            },
            prescient_data_args={
                "representative_dates": [
                    "2019-01-28 00:00",
                    "2019-04-23 00:00",
                    "2019-07-05 00:00",
                    "2019-10-14 00:00",
                ],
                "representative_weights": {
                    "2019-01-28 00:00": 115,
                    "2019-04-23 00:00": 95,
                    "2019-07-05 00:00": 50,
                    "2019-10-14 00:00": 105,
                },
            },
            config={
                "include_investment": True,
                "include_commitment": True,
                "include_redispatch": True,
                "scale_loads": False,
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

        # Comment out the assertion for the objective value since
        # different values are observed across local, Linux, and
        # Windows GitHub Actions runs. This is tracked in issue #163.
        # This assertion will be uncommented once the observed issue
        # is resolved.
        # self.assertAlmostEqual(
        #     value(modObject.model.total_cost_objective_rule), 20105684865.29, places=1
        # )

        assert_units_equivalent(modObject.model.total_cost_objective_rule.expr, u.USD)
