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
from pyomo.environ import units as u

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
                td=self,
                parent=self.b,
                name="renewableGenerationSurplus",
                obj_type=pyo.Expression,
                index=self.m.renewableGenerators,
                expr=(
                    lambda i: self.b.renewableGeneration[i]
                    - self.b.renewableCurtailment[i]
                ),
            ),
            PyomoCheckHelper(
                td=self,
                parent=self.b,
                name="renewableCurtailmentCost",
                obj_type=pyo.Expression,
                index=self.m.renewableGenerators,
                expr=(
                    lambda i: self.b.renewableCurtailment[i]
                    * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                    * self.m.curtailmentCost
                ),
            ),
            PyomoCheckHelper(
                td=self,
                parent=self.b,
                name="thermalGeneratorCost",
                obj_type=pyo.Expression,
                index=self.m.thermalGenerators,
                expr=(
                    lambda i: self.b.thermalGeneration[i]
                    * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                    * (self.m.fixedCost[i] + self.m.varCost[i])
                ),
            ),
            PyomoCheckHelper(
                td=self,
                parent=self.b,
                name="renewableGeneratorCost",
                obj_type=pyo.Expression,
                index=self.m.renewableGenerators,
                expr=(
                    lambda i: self.b.renewableGeneration[i]
                    * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                    * self.m.fixedCost[i]
                ),
            ),
            PyomoCheckHelper(
                td=self,
                parent=self.b,
                name="reactiveGeneratorCost",
                obj_type=pyo.Expression,
                index=self.m.thermalGenerators,
                expr=(
                    lambda i: self.b.thermalReactiveGeneration[i] * self.m.fuelCost[i]
                ),
                condition=(
                    self.m.config["flow_model"] == "ACR"
                    or self.m.config["flow_model"] == "ACP"
                ),
            ),
            PyomoCheckHelper(
                td=self,
                parent=self.b,
                name="loadShed",
                obj_type=pyo.Var,
                index=self.m.buses,
                # domain=pyo.NonNegativeReals,
                # initialize=0,
                # units=u.MW,
            ),
        ]

        for check_helper in to_check:
            check_helper.check()


class PyomoCheckHelper:
    def __init__(
        self,
        td: TestDispatch,
        parent: BlockData,
        name: str,
        obj_type: type,
        index: pyo.Set = None,
        expr=None,
        condition: bool = True,
    ):
        """
        Class that stores expected properties for a pyomo object. Calling
        the .check() method runs asserts to check for those properties.

        :param td:          TestDispatch instance.
        :param parent:      Expected parent object that holds the pyomo object.
        :param name:        Expected name of the pyomo object.
        :param obj_type:    Expected type of the pyomo object.
        :param index:       Expected indexing object of the pyomo object.
        :param expr:        Expected expression of the pyomo object. Defaults to `None`,
                                in which case the object is expected to have no expression.
        :param condition:   Flag that determines whether this pyomo object should be present.
                                If `condition=True`, checks proceed based on the other arguments.
                                If `condition=False`, we check that `block` does __not__ have an
                                attribute `name`, and no other checks are performed.
                                Defaults to `True`.
        :type td:           TestDispatch
        :type parent:       BlockData
        :type name:         str
        :type obj_type:     type
        :type index:        pyomo.environ.Set | None
        :type expr:         function
        :type condition:    bool
        """
        self.td = td
        self.parent = parent
        self.name = name
        self.obj_type = obj_type
        self.index = index
        self.expr = expr
        self.condition = condition

    def _check_exists(self):
        self.td.assertTrue(hasattr(self.parent, self.name))

    def _check_does_not_exist(self):
        self.td.assertFalse(hasattr(self.parent, self.name))

    def _check_type(self):
        self.td.assertIsInstance(self.parent.component(self.name), self.obj_type)

    def _check_index(self):
        if self.index is None:
            self.td.assertFalse(self.parent.component(self.name).is_indexed())
        else:
            self.td.assertTrue(self.parent.component(self.name).is_indexed())
            self.td.assertIs(self.parent.component(self.name).index_set(), self.index)

    def _check_expr(self):
        if self.expr is None:
            self.td.assertFalse(hasattr(self.parent.component(self.name), "expr"))
        else:
            for i in self.index:
                self.td.assertExpressionsStructurallyEqual(
                    self.expr(i),
                    self.parent.component(self.name)[i].expr,
                )

    def check(self):
        if self.condition:
            self._check_exists()
            self._check_type()
            self._check_index()
            self._check_expr()
        else:
            self._check_does_not_exist()
