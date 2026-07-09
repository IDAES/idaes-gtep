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

"""Tests for the commitment.py file of the model.

Variables and Constraints for the Commitment Stage in the
Generation and Transmission Expansion Planning (GTEP) Model

"""

import pytest
import pyomo.common.unittest as unittest
from gtep.model_library.commitment import (
    add_commitment_constraints,
    add_commitment_parameters,
    add_commitment_disjuncts,
    add_investment_commitment_constraints,
    add_investment_commitment_variables,
)
import pyomo.environ as pyo
from pyomo.environ import units as u
from gtep.tests.unit.utils_for_testing import create_model, path_9_bus, input_data_path
from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper
from unittest.mock import patch


class TestObjective(unittest.TestCase):
    def _create_testing_obj(self, config, data_path):

        self.investment_stage = 1
        self.repr_period = 1
        self.commit_period = 1

        self.m = create_model(input_data_path=data_path, config=config).model
        self.b = (
            self.m.investmentStage[self.investment_stage]
            .representativePeriod[self.repr_period]
            .commitmentPeriod[self.commit_period]
        )

    def _make_commitment_param_objects(self):
        self.check_helper.add_object(
            name="commitmentPeriodLength",
            units=u.hr,
            obj_type=pyo.Param,
        )
        self.check_helper.add_object(
            name="carbonTax",
            units=u.dimensionless,
            obj_type=pyo.Param,
        )

    def test_add_commitment_params_no_hydro(self):
        self._create_testing_obj(
            config={"advanced_hydro": False}, data_path=input_data_path
        )

        # Set up helper before function call
        self.check_helper = PyomoCheckHelper(self, self.b)
        self._make_commitment_param_objects()

        with (
            patch(
                "gtep.model_library.hydropower_generation.fix_hydropower_limits"
            ) as mock_hydro,
            patch("gtep.model_library.scaling.add_load_scaling") as mock_scaling,
        ):
            add_commitment_parameters(self.b, self.commit_period, self.investment_stage)

            # scaling should always be called
            mock_scaling.assert_called_once_with(
                self.m,
                self.b,
                self.commit_period,
                self.investment_stage,
                scaling_value=10,
            )

            # advanced_hydro=False, so this should not be called
            mock_hydro.assert_not_called()

        # verify Pyomo params exist / type / units
        self.check_helper.check_all_objects()

        # renewableCapacityExpected is a plain dict, not a Pyomo component
        assert hasattr(self.b, "renewableCapacityExpected")
        assert isinstance(self.b.renewableCapacityExpected, dict)

        # it should contain all renewable generators as keys
        for g in self.m.renewableGenerators:
            assert g in self.b.renewableCapacityExpected

    def test_add_commitment_params_advanced_hydro(self):
        self._create_testing_obj(
            config={"advanced_hydro": True}, data_path=input_data_path
        )

        self.check_helper = PyomoCheckHelper(self, self.b)
        self._make_commitment_param_objects()

        with (
            patch(
                "gtep.model_library.hydropower_generation.fix_hydropower_limits"
            ) as mock_hydro,
            patch("gtep.model_library.scaling.add_load_scaling") as mock_scaling,
        ):
            add_commitment_parameters(self.b, self.commit_period, self.investment_stage)

            mock_hydro.assert_called_once_with(self.b, self.commit_period)
            mock_scaling.assert_called_once_with(
                self.m,
                self.b,
                self.commit_period,
                self.investment_stage,
                scaling_value=10,
            )

        self.check_helper.check_all_objects()
