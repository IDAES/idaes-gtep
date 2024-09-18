import pyomo.common.unittest as unittest

from pyomo.environ import ConcreteModel, Var, SolverFactory, value
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from egret.data.model_data import ModelData
from pyomo.core import TransformationFactory
from prescient.data.providers import gmlc_data_provider
from prescient.simulator.options import Options
from prescient.simulator.config import PrescientConfig
from pyomo.contrib.appsi.solvers.highs import Highs


import logging
from io import StringIO


# Helper functions
def read_debug_model():
    debug_data_path = "./gtep/data/5bus"
    dataObject = ExpansionPlanningData()
    dataObject.load_prescient(debug_data_path)
    return dataObject.md


class TestGTEP(unittest.TestCase):
    def test_model_init(self):
        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with default values, including building a Pyomo.ConcreteModel object
        md = read_debug_model()
        modObject = ExpansionPlanningModel(data=md)
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        modObject.create_model()
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 1)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.data, ModelData)
        self.assertEqual(modObject.num_reps, 3)
        self.assertEqual(modObject.len_reps, 24)
        self.assertEqual(modObject.num_commit, 24)
        self.assertEqual(modObject.num_dispatch, 4)

        # Test that the ExpansionPlanningModel object can read a default dataset and init
        # properly with non-default values
        modObject = ExpansionPlanningModel(
            data=md, stages=2, num_reps=4, len_reps=16, num_commit=12, num_dispatch=12
        )
        self.assertIsInstance(modObject, ExpansionPlanningModel)
        modObject.create_model()
        self.assertIsInstance(modObject.model, ConcreteModel)
        self.assertEqual(modObject.stages, 2)
        self.assertEqual(modObject.formulation, None)
        self.assertIsInstance(modObject.data, ModelData)
        self.assertEqual(modObject.num_reps, 4)
        self.assertEqual(modObject.len_reps, 16)
        self.assertEqual(modObject.num_commit, 12)
        self.assertEqual(modObject.num_dispatch, 12)

        # We have expansion blocks and they are where and what we think they are
        expansion_blocks = modObject.model.component("investmentStage")
        self.assertEqual(len(expansion_blocks), 2)
        self.assertIs(
            modObject.model.investmentStage[modObject.model.stages[1]],
            expansion_blocks[1],
        )

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

    # Solve the debug model as is.  Objective value should be $6078.86
    def test_solve_bigm(self):
        md = read_debug_model()
        modObject = ExpansionPlanningModel(
            data=md, num_reps=1, len_reps=1, num_commit=1, num_dispatch=1
        )
        modObject.create_model()
        opt = Highs()
        if not opt.available():
            print("Ack, no Highs?")
            print(f"{opt.available() = }")
            raise AssertionError
        TransformationFactory("gdp.bound_pretransformation").apply_to(modObject.model)
        TransformationFactory("gdp.bigm").apply_to(modObject.model)
        modObject.results = opt.solve(modObject.model)
        modObject.model.pprint()
        self.assertAlmostEqual(
            value(modObject.model.total_cost_objective_rule), 6078.86, places=1
        )
        # Is it finally time to fix units
        self.assertEqual(
            str(u.get_units(modObject.model.total_cost_objective_rule.expr)), "USD"
        )
