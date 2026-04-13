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
from pyomo.core.base.block import BlockData

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
        to_check = [
            PyomoCheckHelper(
                self,
                self.b,
                "renewableGenerationSurplus",
                pyo.Expression,
                self.m.renewableGenerators,
                (
                    lambda i: self.b.renewableGeneration[i]
                    - self.b.renewableCurtailment[i]
                ),
            ),
            PyomoCheckHelper(
                self,
                self.b,
                "renewableCurtailmentCost",
                pyo.Expression,
                self.m.renewableGenerators,
                (
                    lambda i: self.b.renewableCurtailment[i]
                    * self.m.dispatchPeriodLengthHours
                    * self.m.curtailmentCost
                ),
            ),
            PyomoCheckHelper(
                self,
                self.b,
                "thermalGeneratorCost",
                pyo.Expression,
                self.m.thermalGenerators,
                (
                    lambda i: self.b.thermalGeneration[i]
                    * self.m.dispatchPeriodLengthHours
                    * (self.m.fixedCost[i] + self.m.varCost[i])
                ),
            ),
        ]

        for check_helper in to_check:
            check_helper.check()


class PyomoCheckHelper:
    def __init__(
        self,
        td: TestDispatch,
        block: BlockData,
        name: str,
        obj_type: type,
        index: pyo.Set = None,
        expr: function = None,
    ):
        """
        Class that stores expected properties for a pyomo object. Calling
        the .check() method runs asserts to check for those properties.
        """
        self.td = td
        self.block = block
        self.name = name
        self.obj_type = obj_type
        self.index = index
        self.expr = expr

    def _check_exists(self):
        self.td.assertTrue(hasattr(self.block, self.name))

    def _check_type(self):
        self.td.assertIsInstance(self.block.component(self.name), self.obj_type)

    def _check_index(self):
        if self.index is None:
            self.td.assertFalse(self.block.component(self.name).is_indexed())
        else:
            self.td.assertTrue(self.block.component(self.name).is_indexed())
            self.td.assertIs(self.block.component(self.name).index_set(), self.index)

    def _check_expr(self):
        if self.expr is None:
            self.td.assertFalse(hasattr(self.block.component(self.name), "expr"))
        else:
            for i in self.index:
                self.td.assertExpressionsStructurallyEqual(
                    self.expr(i),
                    self.block.component(self.name)[i].expr,
                )

    def check(self):
        self._check_exists()
        self._check_type()
        self._check_index()
        self._check_expr()