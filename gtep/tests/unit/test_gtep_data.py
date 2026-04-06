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

import pyomo.common.unittest as unittest
import pytest
from gtep.gtep_data import ExpansionPlanningData


class TestExpansionPlanningData(unittest.TestCase):
    def setup():
        pass

    def test_data_init(self):
        # Test that the ExpansionPlanningData object initializes properly with default values
        testObject = ExpansionPlanningData()
        self.assertIsInstance(testObject, ExpansionPlanningData)
        self.assertEqual(testObject.stages, 2)
        self.assertEqual(testObject.num_reps, 4)
        self.assertEqual(testObject.len_reps, 1)
        self.assertEqual(testObject.num_commit, 24)
        self.assertEqual(testObject.num_dispatch, 1)
        self.assertEqual(testObject.duration_dispatch, 60)

        # Test that the ExpansionPlanningData object initializes properly with input values
        testObject = ExpansionPlanningData(1, 2, 2, 2, 2, 15)
        self.assertIsInstance(testObject, ExpansionPlanningData)
        self.assertEqual(testObject.stages, 1)
        self.assertEqual(testObject.num_reps, 2)
        self.assertEqual(testObject.len_reps, 2)
        self.assertEqual(testObject.num_commit, 2)
        self.assertEqual(testObject.num_dispatch, 2)
        self.assertEqual(testObject.duration_dispatch, 15)

        # Test that the ExpansionPlanningData object initializes properly with partial input values
        testObject = ExpansionPlanningData(duration_dispatch=15)
        self.assertIsInstance(testObject, ExpansionPlanningData)
        self.assertEqual(testObject.stages, 2)
        self.assertEqual(testObject.num_reps, 4)
        self.assertEqual(testObject.len_reps, 1)
        self.assertEqual(testObject.num_commit, 24)
        self.assertEqual(testObject.num_dispatch, 1)
        self.assertEqual(testObject.duration_dispatch, 15)
