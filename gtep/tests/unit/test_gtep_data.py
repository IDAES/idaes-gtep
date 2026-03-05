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
import pytest
import pyomo.common.unittest as unittest
import gtep_data


class TestExpansionPlanningData(unittest.TestCase):
    def test_load_prescient():
        # test that representative data can be loaded from prescient with all default loading
        pass

    def test_load_prescient_with_options_input():
        # test that an options dictionary input gets used in the data loader
        pass

    def test_load_prescient_generator_in_service_false():
        pass

    def test_load_prescient_branch_in_service_false():
        pass

    def test_load_prescient_storage_in_service_false():
        pass

    def test_load_prescient_4_dates():
        # test that 4 representative dates results in the correct weights
        pass
