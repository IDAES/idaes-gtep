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
    data_object = ExpansionPlanningData(
        stages=1,
        num_reps=1,
        num_commit=1,
        num_dispatch=1,
    )
    data_object.load_prescient(str(input_data_source))

    mod_object = ExpansionPlanningModel(
        data=data_object,
    )
    mod_object.create_model()

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
        check_helper = PyomoCheckHelper(self, self.b)
        check_helper.add_object(
            name="renewableGenerationSurplus",
            units=u.MW,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        check_helper.add_object(
            name="renewableCurtailmentCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        check_helper.add_object(
            name="thermalGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
        )
        check_helper.add_object(
            name="renewableGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        check_helper.add_object(
            name="reactiveGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
            cond=(
                self.m.config["flow_model"] == "ACR"
                or self.m.config["flow_model"] == "ACP"
            ),
        )
        check_helper.add_object(
            name="loadShed",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.buses,
            domain=pyo.NonNegativeReals,
            # initialize=0,
        )
        check_helper.add_object(
            name="loadShedCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.buses,
        )

        check_helper.check_all()


class PyomoCheckHelper:
    def __init__(
        self,
        td: TestDispatch,
        parent: BlockData,
    ):
        """
        Class that stores expected properties for pyomo objects. Calling
        the `.check_all()` method runs asserts to check for
        those properties on each object.

        :param td:          TestDispatch instance.
        :param parent:      Expected parent object that holds the pyomo objects.
        """
        self.td = td
        self.parent = parent
        self.object_properties = []

    def add_object(
        self,
        name: str,
        obj_type: type,
        units,
        index: pyo.Set = None,
        domain: pyo.Set = None,
        expr=None,
        cond: bool = True,
    ):
        """
        Add a pyomo object to check.

        :param name:        Expected name of the pyomo object.
        :param obj_type:    Expected type of the pyomo object.
        :param units:       Expected units of the pyomo object. Defaults to `None`, in which
                                case the units are not checked.
        :param index:       Expected indexing object of the pyomo object.
        :param domain:      Expected domain of the pyomo object. Defaults to `None`, in which
                                case the object is expected to have no `domain` attribute.
        :param expr:        Expected expression of the pyomo object. Defaults to `None`,
                                in which case the object is expected to have no `expr` attribute.
        :param cond:        Flag that determines whether this pyomo object should be present.
                                If `cond=True`, checks proceed based on the other arguments.
                                If `cond=False`, we check that `block` does __not__ have an
                                attribute `name`, and no other checks are performed.
                                Defaults to `True`.
        :type td:           TestDispatch
        :type parent:       BlockData
        :type name:         str
        :type obj_type:     type
        :type index:        pyomo.environ.Set | None
        :type domain:       pyomo.environ.Set | None
        :type expr:         function
        :type condition:    bool
        """
        self.object_properties.append(
            {
                "name": name,
                "obj_type": obj_type,
                "index": index,
                "domain": domain,
                "units": units,
                # "expr": expr,
                "cond": cond,
            }
        )

    def _check_exists(self, properties):
        self.td.assertTrue(hasattr(self.parent, properties["name"]))

    def _check_does_not_exist(self, properties):
        self.td.assertFalse(hasattr(self.parent, properties["name"]))

    def _check_type(self, properties):
        self.td.assertIsInstance(properties["obj"], properties["obj_type"])

    def _check_index(self, properties):
        if properties["index"] is None:
            self.td.assertFalse(properties["obj"].is_indexed())
        else:
            self.td.assertTrue(properties["obj"].is_indexed())
            self.td.assertIs(properties["obj"].index_set(), properties["index"])

    def _iter_func_over_index(self, obj, func):
        for i in obj.index_set():
            func(obj[i])

    def _check_domain(self, properties):
        if properties["domain"] is not None:
            iter_func = lambda x: self.td.assertIs(x.domain, properties["domain"])
            self._iter_func_over_index(properties["obj"], iter_func)
        else:
            iter_func = lambda x: self.td.assertFalse(hasattr(x, "domain"))
            self._iter_func_over_index(properties["obj"], iter_func)

    def _check_units(self, properties):
        iter_func = lambda x: self.td.assertEqual(u.get_units(x), properties["units"])
        self._iter_func_over_index(properties["obj"], iter_func)

    # def _check_expr(self, properties):
    #     if properties["expr"] is None:
    #         self.td.assertFalse(
    #             hasattr(self.parent.component(properties["name"]), "expr")
    #         )
    #     else:
    #         for i in self.parent.component(properties["name"]).index_set():
    #             self.td.assertExpressionsStructurallyEqual(
    #                 properties["expr"](i),
    #                 self.parent.component(properties["name"])[i].expr,
    #             )

    def check_all(self):
        for properties in self.object_properties:
            if properties["cond"]:
                self._check_exists(properties)
                properties["obj"] = self.parent.component(properties["name"])
                self._check_type(properties)
                self._check_index(properties)
                self._check_domain(properties)
                self._check_units(properties)
                # self._check_expr(properties)
            else:
                self._check_does_not_exist(properties)
