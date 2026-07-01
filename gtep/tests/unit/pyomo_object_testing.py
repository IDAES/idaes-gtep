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


from pyomo.environ import units as u
from pyomo.core.base.block import BlockData
from pyomo.core.base.component import ComponentData
from pyomo.util.check_units import check_units_equivalent
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
        parent: str | None = None,
        skipped_on: list = [],
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
        :param parent:          Optional, used if the object to be checked is not an attribute
                                    of `self.parent` and is instead an attribute of an attribute
                                    of `self.parent`. Identifies the name of the intermediate
                                    object.
        :param skipped_on:      List of indexes where this pyomo object should have been assigned
                                    `pyo.Constraint.Skip`. Defaults to an empty list. Ignored
                                    if `index=None`.
        :type name:             str
        :type obj_type:         type
        :type units:            tuple | None, optional
        :type index:            pyomo.environ.Set | None, optional
        :type check_func:       function, optional
        :type cond:             bool, optional
        :type parent:           str | None, optional
        :type skipped_on:       func, optional
        """
        self.object_properties.append(
            {
                "name": name,
                "obj_type": obj_type,
                "units": units,
                "index": index,
                "check_func": check_func,
                "cond": cond,
                "parent": parent,
                "skipped_on": skipped_on,
            }
        )

    def _get_parent(self, properties: dict):
        if properties["parent"] is None:
            return self.parent
        parent = self.parent.component(properties["parent"])
        if parent.is_indexed():
            return parent[parent.index_set().first()]
        return parent

    def _check_exists(self, properties: dict):
        """Checks that an attribute with the provided name exists."""
        parent = self._get_parent(properties)
        self.test_class.assertTrue(
            hasattr(parent, properties["name"]),
            f"{parent} does not have the attribute {properties['name']}, but should",
        )
        return parent.component(properties["name"])

    def _check_does_not_exist(self, properties: dict):
        """Checks that an attribute with the provided name does not exist."""
        parent = self._get_parent(properties)
        self.test_class.assertFalse(
            hasattr(parent, properties["name"]),
            f"{parent} has the attribute {properties['name']}, but shouldn't",
        )

    def _check_type(self, obj: pyo.Component, properties: dict):
        """Checks the type of the provided object."""
        self.test_class.assertIsInstance(
            obj,
            properties["obj_type"],
            f"{obj} of unexpected type",
        )

    def _iter_func_over_index(self, func, obj: pyo.Component, **additional_args):
        """
        Iterates the provided `func` over each element of `obj`. If `obj`
        is not indexed, then simply calls `func(obj)`.

        If additional keyword arguments are passed, they are also passed into
        the `func` call. If `obj` is indexed, then each `arg[i]` must exist
        for every `i` in `obj.index_set()`.

        Note: this method silently skips any indices in the index set
        which are unused (e.g., via `pyo.Constraint.Skip`). `_check_skipped`
        is used to verify the expected indices are indeed skipped.
        """
        if obj.is_indexed():
            for i in obj.index_set():
                if i not in obj:
                    continue
                func(
                    obj[i],
                    **{k: v[i] for k, v in additional_args.items()},
                )
        else:
            func(obj, **additional_args)

    def _check_units(self, obj: pyo.Component, properties: dict):
        """Checks the object's units."""

        def units(x):
            try:
                return u.get_units(x.expr if hasattr(x, "expr") else x)
            except Exception as e:
                raise RuntimeError(f"Error getting units for {x}") from e

        iter_func = lambda x: self.test_class.assertTrue(
            check_units_equivalent(units(x), properties["units"]),
            f"Expected {x.name} to have units {properties['units']} but got {units(x)}",
        )
        self._iter_func_over_index(iter_func, obj)

    def _check_index(self, obj: pyo.Component, properties: dict):
        """Checks the object's index (including whether it has an index)."""
        if properties["index"] is None:
            self.test_class.assertFalse(
                obj.is_indexed(),
                f"Expected {obj} to not be indexed, but it is",
            )
        else:
            self.test_class.assertTrue(
                obj.is_indexed(), f"Expected {obj} to be indexed, but it isn't"
            )
            if isinstance(
                properties["index"], dict
            ):  # remove once stop indexing by dicts
                exp = len(properties["index"])
                got = len(obj.index_set())
                self.test_class.assertEqual(
                    exp,
                    got,
                    f"Expected index of {obj} to have {exp} elements, but got {got}",
                )
                for key in properties["index"]:
                    self.test_class.assertIn(
                        key,
                        obj.index_set(),
                        f"{key} not found in index of {obj}: {obj.index_set()}",
                    )
            else:
                exp = properties["index"]
                got = obj.index_set()
                self.test_class.assertIs(
                    exp,
                    got,
                    f"Expected {obj} to have index set {exp} but got {got}",
                )

    def _check_skipped(self, obj: pyo.Component, properties: dict):
        """Checks that the object is `pyo.Constraint.Skip`."""
        not_skipped = []
        for i in obj.index_set():
            if i in properties["skipped_on"]:
                self.test_class.assertNotIn(i, obj)
            else:
                self.test_class.assertIn(i, obj)
                not_skipped.append(i)
        return not_skipped

    def _check_bounds(self, obj: pyo.Component):
        """Checks that the object's bounds are consistent with its domain (if it has them)."""

        def iter_func(x):
            if hasattr(x, "bounds"):
                for bound in x.bounds:
                    if bound is not None:
                        self.test_class.assertIn(
                            bound,
                            x.domain,
                            f"{obj}'s bound ({bound}) not in its domain {(x.domain)}",
                        )

        self._iter_func_over_index(iter_func, obj)

    def _run_check_func(self, obj: pyo.Component, properties: dict):
        """Runs the additional check function, if a check function is provided."""
        if properties["check_func"] is not None:
            properties["check_func"](self.test_class, obj)

    def check_all_objects(self):
        """Performs checks on all members of `self.object_properties`."""
        for properties in self.object_properties:
            if properties["cond"]:
                obj = self._check_exists(properties)
                self._check_type(obj, properties)
                self._check_index(obj, properties)
                if properties["index"] is not None:
                    properties["index"] = self._check_skipped(obj, properties)
                self._check_units(obj, properties)
                self._check_bounds(obj)
                self._run_check_func(obj, properties)
            else:
                self._check_does_not_exist(properties)

    def check_expr_contains(self, c: ComponentData, expected: list | dict):
        """
        Checks that the given component expr has exactly the given params and vars.

        :param c:           Component
        :param expected:    A list containing every param and var expected to be
                                in the expr of this component, or a dict mapping
                                index elements to such lists for indexed components.
        :type c:            ComponentData
        :type expected:     list | dict
        """
        self._iter_func_over_index(self._nonindexed_expr_contains, c, expected=expected)

    def _nonindexed_expr_contains(self, c: ComponentData, expected: list):
        actual = [
            a
            for a in (
                list(pyo.visitor.identify_variables(c.expr))
                + list(pyo.visitor.identify_mutable_parameters(c.expr))
            )
            if isinstance(a, ComponentData)
        ]

        ids = [id(a) for a in actual]
        names = [a.name for a in actual]
        expected_names = [e.name for e in expected]

        self.test_class.assertEqual(
            len(actual), len(expected), f"{expected_names} vs {names}"
        )
        for exp, exp_name in zip(expected, expected_names):
            self.test_class.assertIn(id(exp), ids, f"{exp_name} not in {names}")
