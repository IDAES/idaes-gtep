from io import StringIO

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
        Expression and Var instances. Call `add_object` to append
        a new set of properties of an object to be checked. Calling
        the `.check_all()` method runs asserts to check for
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
        :type test_class:               Testest_classispatch
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

    def check_all_added_objects(self):
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


def extract_terms_from_string_expression(expr: str) -> list:
    """
    Helper function for `_parse_constraint_pprint` that extracts
    individual terms from a string representation of an expression.

    :param expr:    Expression to be parsed.
    :type expr:     str
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


def parse_constraint_pprint(block_name: str, constraint: pyo.Constraint) -> dict:
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
            "expr": extract_terms_from_string_expression(expr),
            "val": val,
        }
    return out
