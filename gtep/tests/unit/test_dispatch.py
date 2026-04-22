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

    # get first dispatch block -- TODO: should we check all the dispatch blocks?
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
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
            expr=(
                lambda i: self.b.renewableGeneration[i] - self.b.renewableCurtailment[i]
            ),
        )
        check_helper.add_object(
            name="renewableCurtailmentCost",
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
            expr=(
                lambda i: self.b.renewableCurtailment[i]
                * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                * self.m.curtailmentCost
            ),
        )
        check_helper.add_object(
            name="thermalGeneratorCost",
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
            expr=(
                lambda i: self.b.thermalGeneration[i]
                * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                * (self.m.fixedCost[i] + self.m.varCost[i])
            ),
        )
        check_helper.add_object(
            name="renewableGeneratorCost",
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
            expr=(
                lambda i: self.b.renewableGeneration[i]
                * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                * self.m.fixedCost[i]
            ),
        )
        check_helper.add_object(
            name="reactiveGeneratorCost",
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
            expr=(lambda i: self.b.thermalReactiveGeneration[i] * self.m.fuelCost[i]),
            condition=(
                self.m.config["flow_model"] == "ACR"
                or self.m.config["flow_model"] == "ACP"
            ),
        )
        check_helper.add_object(
            name="loadShed",
            obj_type=pyo.Var,
            index=self.m.buses,
            domain=pyo.NonNegativeReals,
            # initialize=0,
            units=u.MW,
        )
        check_helper.add_object(
            name="loadShedCost",
            obj_type=pyo.Expression,
            index=self.m.buses,
            expr=(
                lambda i: self.b.loadShed[i]
                * u.convert(self.m.dispatchPeriodLength, to_units=u.hr)
                * self.m.loadShedCostperCurtailment
            ),
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
        index: pyo.Set = None,
        domain: pyo.Set = None,
        units=u.dimensionless,
        expr=None,
        condition: bool = True,
    ):
        """
        Add a pyomo object to check.

        :param name:        Expected name of the pyomo object.
        :param obj_type:    Expected type of the pyomo object.
        :param index:       Expected indexing object of the pyomo object.
        :param domain:      Expected domain of the pyomo object. Defaults to `None`, in which
                                case the object is expected to have no `domain` attribute.
        :param units:       Expected units of the pyomo object. Defaults to
                                `pyomo.environ.units.dimensionless`, corresponding to the case
                                that the object is dimensionless (or does not have a defined unit).
        :param expr:        Expected expression of the pyomo object. Defaults to `None`,
                                in which case the object is expected to have no `expr` attribute.
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
                "expr": expr,
                "condition": condition,
            }
        )

    def _check_exists(self, properties):
        self.td.assertTrue(hasattr(self.parent, properties["name"]))

    def _check_does_not_exist(self, properties):
        self.td.assertFalse(hasattr(self.parent, properties["name"]))

    def _check_type(self, properties):
        self.td.assertIsInstance(
            self.parent.component(properties["name"]), properties["obj_type"]
        )

    def _check_index(self, properties):
        if properties["index"] is None:
            self.td.assertFalse(self.parent.component(properties["name"]).is_indexed())
        else:
            self.td.assertTrue(self.parent.component(properties["name"]).is_indexed())
            self.td.assertIs(
                self.parent.component(properties["name"]).index_set(),
                properties["index"],
            )

    def _check_domain(self, properties):
        if properties["domain"] is not None:
            for i in self.parent.component(properties["name"]).index_set():
                self.td.assertIs(
                    self.parent.component(properties["name"])[i].domain,
                    properties["domain"],
                )
        else:
            self.td.assertFalse(
                hasattr(self.parent.component(properties["name"]), "domain")
            )

    def _check_units(self, properties):
        self.td.assertEqual(
            u.get_units(self.parent.component(properties["name"])), properties["units"]
        )

    def _check_expr(self, properties):
        if properties["expr"] is None:
            self.td.assertFalse(
                hasattr(self.parent.component(properties["name"]), "expr")
            )
        else:
            for i in self.parent.component(properties["name"]).index_set():
                self.td.assertExpressionsStructurallyEqual(
                    properties["expr"](i),
                    self.parent.component(properties["name"])[i].expr,
                )

    def check_all(self):
        for properties in self.object_properties:
            if properties["condition"]:
                self._check_exists(properties)
                self._check_type(properties)
                self._check_index(properties)
                self._check_domain(properties)
                self._check_units(properties)
                self._check_expr(properties)
            else:
                self._check_does_not_exist(properties)
