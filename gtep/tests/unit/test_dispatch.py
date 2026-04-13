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


from pathlib import Path

import pyomo.common.unittest as unittest
import pyomo.environ as pyo

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData

curr_dir = Path(__file__).resolve().parent
input_data_source = (curr_dir / ".." / ".." / "data" / "5bus").resolve()


def get_dispatch_block():
    # create model
    data_object = ExpansionPlanningData()
    data_object.load_prescient(
        str(input_data_source)  # load_prescient should accept pathlib paths
    )
    mod_object = ExpansionPlanningModel(
        stages=2,
        data=data_object,
        num_reps=2,
        len_reps=1,
        num_commit=6,
        num_dispatch=4,
    )
    mod_object.create_model()

    # get first dispatch block
    current_block = mod_object.model
    for component_name in [
        "investmentStage",
        "representativePeriod",
        "commitmentPeriod",
        "dispatchPeriod",
    ]:
        block = current_block.component(component_name)
        first_idx = block.index_set().at(1)
        current_block = block[first_idx]

    return current_block


class TestDispatch(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.b = get_dispatch_block()
        cls.m = cls.b.model()

    def test_add_dispatch_variables(self):

        ### CHECK ALL SETS/VARS/ETC WE EXPECT EXIST ###
        # format: (parent, name, type, index)
        to_check = (
            (self.m, "renewableGenerators", pyo.Set, None),
            (self.m, "dispatchPeriodLengthHours", pyo.Param, None),
            (self.m, "curtailmentCost", pyo.Param, None),
            (self.b, "renewableGeneration", pyo.Var, self.m.renewableGenerators),
            (self.b, "renewableCurtailment", pyo.Var, self.m.renewableGenerators),
            (
                self.b,
                "renewableGenerationSurplus",
                pyo.Expression,
                self.m.renewableGenerators,
            ),
            (
                self.b,
                "renewableCurtailmentCost",
                pyo.Expression,
                self.m.renewableGenerators,
            ),
        )
        for parent, name, expected_type, index in to_check:
            self.assertHasAttr(parent, name)
            self.assertIsInstance(parent.component(name), expected_type)
            if index is None:
                self.assertFalse(parent.component(name).is_indexed())
            else:
                self.assertTrue(parent.component(name).is_indexed())
                self.assertIs(parent.component(name).index_set(), index)

        ### CHECK EXPRESSIONS ###
        for renew in self.m.renewableGenerators:

            # renewableGenerationSurplus
            actual = self.b.renewableGenerationSurplus[renew].expr
            expected = (
                self.b.renewableGeneration[renew] - self.b.renewableCurtailment[renew]
            )
            self.assertExpressionsStructurallyEqual(actual, expected)

            # renewableCurtailmentCost
            actual = self.b.renewableCurtailmentCost[renew].expr
            expected = (
                self.b.renewableCurtailment[renew]
                * self.m.dispatchPeriodLengthHours
                * self.m.curtailmentCost
            )
            self.assertExpressionsStructurallyEqual(actual, expected)
