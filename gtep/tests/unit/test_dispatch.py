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
from io import StringIO
import re

import pyomo.common.unittest as unittest
import pyomo.environ as pyo
from pyomo.core.base.block import BlockData
from pyomo.environ import units as u

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData

curr_dir = Path(__file__).resolve().parent
input_data_source = (curr_dir / ".." / ".." / "data" / "5bus").resolve()


class TestDispatch(unittest.TestCase):

    def _create_model(self, config={}):
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

        for config_option, config_val in config.items():
            mod_object.config[config_option] = config_val

        return mod_object

    def _get_first_dispatch_block(self):
        current_block = self.m
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

    def _add_properties_for_dispatch_variables(self):
        self.check_helper.add_object(
            name="renewableGenerationSurplus",
            units=u.MW,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        self.check_helper.add_object(
            name="renewableCurtailmentCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        self.check_helper.add_object(
            name="thermalGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
        )
        self.check_helper.add_object(
            name="renewableGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        self.check_helper.add_object(
            name="reactiveGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
            cond=(
                self.m.config["flow_model"] == "ACR"
                or self.m.config["flow_model"] == "ACP"
            ),
        )
        self.check_helper.add_object(
            name="loadShed",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.buses,
        )
        self.check_helper.add_object(
            name="loadShedCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.buses,
        )
        self.check_helper.add_object(
            name="renewableSurplusDispatch",
            units=u.MW,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="thermalGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="renewableGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="renewableGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="reactiveGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="loadShedCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="curtailmentCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="operatingCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="renewableCurtailmentDispatch",
            units=u.MW,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="powerFlow",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.transmission,
        )
        self.check_helper.add_object(
            name="spinningReserve",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.thermalGenerators,
        )
        self.check_helper.add_object(
            name="quickstartReserve",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.thermalGenerators,
        )

    def _extract_terms_from_string_expression(self, expr: str) -> list:
        """
        Helper function for `_parse_constraint_pprint` that extracts
        individual terms from a string representation of an expression.
        """
        expr = expr.replace(" ", "")
        stack = [1]  # sign context
        sign = 1
        term = ""
        out = []

        i = 0
        while i < len(expr):
            c = expr[i]
            if c in "+-":
                if term:
                    out.append((sign * stack[-1], term))
                    term = ""
                sign = 1 if c == "+" else -1
            elif c == "(":
                stack.append(stack[-1] * sign)
                sign = 1
            elif c == ")":
                if term:
                    out.append((sign * stack[-1], term))
                    term = ""
                stack.pop()
            else:
                term += c
            i += 1
        if term:
            out.append((sign * stack[-1], term))

        return out

    def _parse_constraint_pprint(self, constraint: pyo.Constraint) -> dict:
        """
        Parses the output of a constraint's `.pprint()` function, returning
        a dictionary representing the contents of the constraints' expressions,
        of the form:
        ```
        {
            i: {
                "expr": [(sign, term), (sign, term), ...],
                "val": val
            }
        }
        ```
        where
            - `i` is an element of the constraint's index (per its `.pprint()` function)
            - `val` is the value of the expression associated with `i`
            - `term` is an individual term of the expression associated with `i`
            - `sign` is the sign of an individual term of teh expression associated with `i`

        :param constraint:  Constraint to be parsed.
        """
        buf = StringIO()
        constraint.pprint(ostream=buf)
        pprinted = buf.getvalue().replace(self.b.name, "b")

        out = {}
        for index_expr_pprinted in pprinted.split("\n")[3:-1]:
            index_expr_split = index_expr_pprinted.split(":")
            i = index_expr_split[0].strip()
            expr = index_expr_split[2].strip()
            val = index_expr_split[3].strip()
            # print(i, ":", expr)

            out[i] = {
                "expr": self._extract_terms_from_string_expression(expr),
                "val": val,
            }
        return out

    def _check_dispatch_constraints(self):
        """Checks dispatch constraints."""
        c = self.b.flow_balance
        constraints_by_index = self._parse_constraint_pprint(c)

        for i in c.index_set():
            load_terms = []
            for sign, term in constraints_by_index[i]["expr"]:
                match = re.fullmatch(r"loads\[([^\]]+)\]", term)
                if match:
                    load_terms.append((sign, term, match.group(1)))

            self.assertTrue(len(load_terms) == 1)  # assert only one load term
            self.assertTrue(load_terms[0][0] == -1)  # assert term is negative
            self.assertTrue(load_terms[0][2] == i)  # assert term is for this index

    def _coordinate_tests(self, config):
        """Creates a model and runs tests for a given set of config options."""
        self.m = self._create_model(config=config).model
        self.b = self._get_first_dispatch_block()

        self.check_helper = PyomoCheckHelper(self, self.b)
        self._add_properties_for_dispatch_variables()
        self.check_helper.check_all()

        self._check_dispatch_constraints()

    def test_default_config_options(self):
        self._coordinate_tests(config={})
        # self._coordinate_tests(config={config_options})


class PyomoCheckHelper:
    def __init__(
        self,
        td: TestDispatch,
        parent: BlockData,
    ):
        """
        Class that stores expected properties for pyomo objects, including
        Expression and Var instances. Call `add_object` to append
        a new set of properties of an object to be checked. Calling
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
        index: pyo.Set | None = None,
        cond: bool = True,
    ):
        """
        Add a pyomo object to check.

        :param name:            Expected name of the pyomo object.
        :param obj_type:        Expected type of the pyomo object.
        :param units:           Expected units of the pyomo object. Defaults to `None`, in which
                                    case the units are not checked.
        :param index:           Expected indexing object of the pyomo object. Defaults to `None`,
                                    in which case we check that there is no index.
        :param cond:            Flag that determines whether this pyomo object should be present.
                                    If `cond=True`, checks proceed based on the other arguments.
                                    If `cond=False`, we check that `block` does NOT have an
                                    attribute `name`, and no other checks are performed.
                                    Defaults to `True`.
        :type td:               TestDispatch
        :type parent:           BlockData
        :type name:             str
        :type obj_type:         type
        :type index:            pyomo.environ.Set | None, optional
        :type bounds:           tuple | None, optional
        :type condition:        bool, optional
        """ 
        self.object_properties.append(
            {
                "name": name,
                "obj_type": obj_type,
                "units": units,
                "index": index,
                "cond": cond,
            }
        )

    def _check_exists(self, properties: dict):
        """Checks that an attribute with the provided name exists on `self.parent`."""
        self.td.assertTrue(hasattr(self.parent, properties["name"]))

    def _check_does_not_exist(self, properties: dict):
        """Checks that an attribute with the provided name does not exist on `self.parent`."""
        self.td.assertFalse(hasattr(self.parent, properties["name"]))

    def _check_type(self, obj: pyo.Component, properties: dict):
        """Checks the type of the provided object."""
        self.td.assertIsInstance(obj, properties["obj_type"])

    def _iter_func_over_index(self, obj: pyo.Component, func):
        """
        Iterates the provided `func` over each element of `obj`. If `obj`
        is not indexed, then simply calls `func(obj)`.
        """
        if obj.is_indexed():
            for i in obj.index_set():
                func(obj[i])
        else:
            func(obj)

    def _check_units(self, obj: pyo.Component, properties: dict):
        """Checks the object's units."""
        iter_func = lambda x: self.td.assertEqual(u.get_units(x), properties["units"])
        self._iter_func_over_index(obj, iter_func)

    def _check_index(self, obj: pyo.Component, properties: dict):
        """Checks the object's index (including whether it has an index)."""
        if properties["index"] is None:
            self.td.assertFalse(obj.is_indexed())
        else:
            self.td.assertTrue(obj.is_indexed())
            if isinstance(properties["index"], dict):
                self.td.assertEqual(len(properties["index"]), len(obj.index_set()))
                for key in properties["index"]:
                    self.td.assertIn(key, obj.index_set())
            else:
                self.td.assertIs(obj.index_set(), properties["index"])

    def _check_bounds(self, obj: pyo.Component):
        """Checks that the object's bounds are consistent with its domain (if it has them)."""
        def iter_func(x):
            if hasattr(x, "bounds"):
                for bound in x.bounds:
                    if bound is not None:
                        self.td.assertIn(bound, x.domain)
        self._iter_func_over_index(obj, iter_func)

    def check_all(self):
        """Performs checks on all members of `self.object_properties`."""
        for properties in self.object_properties:
            if properties["cond"]:
                self._check_exists(properties)
                obj = self.parent.component(properties["name"])

                self._check_type(obj, properties)
                self._check_units(obj, properties)
                self._check_index(obj, properties)
                self._check_bounds(obj)
            else:
                self._check_does_not_exist(properties)
