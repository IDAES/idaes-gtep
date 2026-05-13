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


from io import StringIO
import re

from pyomo.environ import units as u
from pyomo.core.base.block import BlockData
import pyomo.environ as pyo


class PyomoCheckHelper:
    def __init__(
        self,
        test_class,
        parent: BlockData,
    ):
        """
        Class that stores expected properties for pyomo objects, including
        Expression and Var instances. Call `add_object()` to append
        a set of properties of an object to be checked. Calling
        the `.check_all_objects()` method runs asserts to check for
        those properties on each object.

        :param test_class:  Testing class instance.
        :param parent:      Expected parent object that holds the pyomo objects.
        """
        self.test_class = test_class
        self.parent = parent
        self.object_properties = []

    def add_object(
        self,
        name: str,
        obj_type: type,
        units=None,
        index: pyo.Set | None = None,
        check_func=None,
        cond: bool = True,
    ):
        """
        Add a pyomo object to check.

        :param name:            Expected name of the pyomo object.
        :param obj_type:        Expected type of the pyomo object.
        :param units:           Expected units of the pyomo object.
        :param index:           Expected indexing object of the pyomo object. Defaults to `None`,
                                    in which case we check that there is no index.
        :param check_func:      Function with additional checks to run on the object. For instance,
                                    checking the structure of a constraint expression. Defaults to `None`,
                                    in which case no checks are performed. Must take the testing class instance
                                    and object to be checked as its only two arguments.
        :param cond:            Flag that determines whether this pyomo object should be present.
                                    If `cond=True`, checks proceed based on the other arguments.
                                    If `cond=False`, we check that `block` does NOT have an
                                    attribute `name`, and no other checks are performed.
                                    Defaults to `True`.
        :type name:             str
        :type obj_type:         type
        :type units:            tuple | None, optional
        :type index:            pyomo.environ.Set | None, optional
        :type check_func:       function, optional
        :type cond:             bool, optional
        """
        self.object_properties.append(
            {
                "name": name,
                "obj_type": obj_type,
                "units": units,
                "index": index,
                "check_func": check_func,
                "cond": cond,
            }
        )

    def _check_exists(self, properties: dict):
        """Checks that an attribute with the provided name exists on `self.parent`."""
        self.test_class.assertTrue(hasattr(self.parent, properties["name"]))

    def _check_does_not_exist(self, properties: dict):
        """Checks that an attribute with the provided name does not exist on `self.parent`."""
        self.test_class.assertFalse(hasattr(self.parent, properties["name"]))

    def _check_type(self, obj: pyo.Component, properties: dict):
        """Checks the type of the provided object."""
        self.test_class.assertIsInstance(obj, properties["obj_type"])

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
        iter_func = lambda x: self.test_class.assertEqual(
            u.get_units(x), properties["units"]
        )
        self._iter_func_over_index(obj, iter_func)

    def _check_index(self, obj: pyo.Component, properties: dict):
        """Checks the object's index (including whether it has an index)."""
        if properties["index"] is None:
            self.test_class.assertFalse(obj.is_indexed())
        else:
            self.test_class.assertTrue(obj.is_indexed())
            if isinstance(properties["index"], dict):
                self.test_class.assertEqual(
                    len(properties["index"]), len(obj.index_set())
                )
                for key in properties["index"]:
                    self.test_class.assertIn(key, obj.index_set())
            else:
                self.test_class.assertIs(obj.index_set(), properties["index"])

    def _check_bounds(self, obj: pyo.Component):
        """Checks that the object's bounds are consistent with its domain (if it has them)."""

        def iter_func(x):
            if hasattr(x, "bounds"):
                for bound in x.bounds:
                    if bound is not None:
                        self.test_class.assertIn(bound, x.domain)

        self._iter_func_over_index(obj, iter_func)

    def _run_check_func(self, obj: pyo.Component, properties: dict):
        """Runs the additional check function, if a check function is provided."""
        if properties["check_func"] is not None:
            properties["check_func"](self.test_class, obj)

    def check_all_objects(self):
        """Performs checks on all members of `self.object_properties`."""
        for properties in self.object_properties:
            if properties["cond"]:
                self._check_exists(properties)
                obj = self.parent.component(properties["name"])

                self._check_type(obj, properties)
                self._check_units(obj, properties)
                self._check_index(obj, properties)
                self._check_bounds(obj)
                self._run_check_func(obj, properties)
            else:
                self._check_does_not_exist(properties)

    @classmethod
    def _extract_terms_from_string_expression(cls, expr: str) -> list[tuple]:
        """
        Helper function for `parse_constraint_pprint` that extracts
        individual terms from a string representation of an expression,
        including flattening nested parentheses.

        :param expr:    Expression to be parsed.
        :type expr:     str
        :returns:       List of tuples. Each tuple corresponds to a term in the
                            expression and is of the form `(sign, term)` where `sign`
                            is either `1` or `-1` and `term` is the term itself
                            (including index, e.g. `"m.lines[branch_2_1]"`)
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

    @classmethod
    def parse_constraint_pprint(
        cls, block_name: str, constraint: pyo.Constraint
    ) -> dict:
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
            - `sign` is the sign of an individual term of the expression associated with `i`

        :param block_name:  `.name` attribute of the block that `constraint` is on.
        :param constraint:  Constraint to be parsed.
        :type block_name:   str
        :type constraint:   pyomo.environ.Constraint
        """
        buf = StringIO()
        constraint.pprint(ostream=buf)
        pprinted = buf.getvalue().replace(block_name, "b")

        out = {}
        for index_expr_pprinted in pprinted.split("\n")[3:-1]:
            index_expr_split = index_expr_pprinted.split(":")
            i = index_expr_split[0].strip()
            expr = index_expr_split[2].strip()
            val = index_expr_split[3].strip()
            # print(i, ":", expr)

            out[i] = {
                "expr": cls._extract_terms_from_string_expression(expr),
                "val": val,
            }
        return out

    def check_constraint_for_terms(
        self,
        constraint_expr: list,
        term_to_find: str,
        expected_signs: list[int],
        expected_indices: list[str],
    ):
        """
        Checks that the given constraint expression contains expected term(s). Intended as a helper
        function for constraint check functions.

        :param constraint_expr:     Constraint expression (`"expr"` value from an element of `parse_constraint_pprint`).
        :param term_to_find:        Name of term to find (e.g., `"loads"`).
        :param expected_signs:      Expected signs of matching terms.
        :param expected_indices:    Expected indices of matching terms.
        :type constraint_expr:      dict
        :type term_to_find:         str
        :type expected_signs:       list[int]
        :type expected_indices:     list
        """
        if len(expected_signs) != len(expected_indices):
            raise ValueError(
                "Expected signs and indices must be lists of matching length."
            )

        matching_terms = []
        for sign, term in constraint_expr:
            match = re.fullmatch(rf"{term_to_find}\[([^\]]+)\]", term)
            if match:
                matching_terms.append(
                    {
                        "sign": sign,
                        "index": match.group(1),
                    }
                )

        self.test_class.assertEqual(len(matching_terms), len(expected_signs))
        for s, i in zip(expected_signs, expected_indices):
            matching_term = [
                term
                for term in matching_terms
                if term["sign"] == s and term["index"] == i
            ]
            self.test_class.assertEqual(
                len(matching_term),
                1,
                f"There should be only one term matching {'-' if s == -1 else ''}{term_to_find}[{i}]",
            )
