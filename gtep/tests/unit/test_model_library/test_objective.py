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
"""
Unit tests for the function to create the objective function component
"""

import pyomo.common.unittest as unittest
from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper
import pyomo.environ as pyo
from pyomo.environ import units as u
from gtep.tests.unit.utils_for_testing import create_model, path_9_bus


class TestObjective(unittest.TestCase):
    def _make_expressions(self):
        # add necessary expressions
        self.check_helper.add_object(
            name="operatingCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="storageCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="expansionCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="penaltyCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )

        self.check_helper.add_object(
            name="total_cost_objective_rule",
            units=u.USD,
            obj_type=pyo.Objective,
        )

    def _coordinate_tests(self, config):
        """Creates a model and runs tests for a given set of config options."""
        self.m = create_model(input_data_path=path_9_bus, config=config).model

        self.check_helper = PyomoCheckHelper(self, self.m)
        self._make_expressions()
        self.check_helper.check_all_objects()

    def test_check_expressions_with_storage(self):
        self._coordinate_tests(config={"storage": True})

    def test_check_expressions_no_storage(self):
        # check that each of the expressions were created
        self._coordinate_tests(config={"storage": False})
        # check that storage cost was set to 0
        self.assertEqual(pyo.value(self.m.storageCostTotal), 0)
