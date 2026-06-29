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

"""Unit Tests for functions that define variables and constraints
for the investment stage in the Generation and Transmission Expansion Planning (GTEP)
model.

"""

from gtep.model_library.investment import (
    add_investment_params_and_variables,
    add_investment_disjuncts,
    add_investment_constraints,
)
from pyomo.environ import units as u
import pyomo.environ as pyo
import pytest
from unittest.mock import patch
import pyomo.common.unittest as unittest
from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing
from pathlib import Path

curr_dir = Path(__file__).resolve().parent
input_data_source = (
    curr_dir / ".." / ".." / ".." / "data" / "9_bus_GTEP_dir"
).resolve()

bus_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "Bus_data_gen_weights_mappings.csv"
).resolve()
cost_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
).resolve()
ng_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "Total_Energy_Supply_Disposition_and_Price_Summary.csv"
).resolve()


def get_block(model, layer):
    current_block = model
    for component_name in [
        "investmentStage",
        "representativePeriod",
        "commitmentPeriod",
        "dispatchPeriod",
    ]:
        block = current_block.component(component_name)
        first_idx = block.index_set().at(1)
        current_block = block[first_idx]
        if layer == component_name:
            break

    return current_block


class TestInvestment(unittest.TestCase):
    def _create_model(self, config=None):
        # create model
        data_object = ExpansionPlanningData(
            stages=1,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
        )
        data_object.load_prescient(input_data_source)
        data_object.load_storage_csv(str(input_data_source))

        candidate_gens = [
            "Natural Gas_FE",
            "Solar - Utility PV",
            "Land-Based Wind",
        ]

        data_processing_object = DataProcessing()
        data_processing_object.load_gen_data(
            bus_data_path=bus_data_path,
            cost_data_path=cost_data_path,
            ng_cost_path=ng_data_path,
            candidate_gens=candidate_gens,
            save_csv=False,
        )

        mod_object = ExpansionPlanningModel(
            data=data_object,
            cost_data=data_processing_object,
        )

        if config is not None:
            for config_option, config_val in config.items():
                mod_object.config[config_option] = config_val

        mod_object.create_model()

        return mod_object

    def _make_param_and_var_objects(self):
        self.check_helper.add_object(
            name="maxThermalInvestment",
            units=u.MW,
            obj_type=pyo.Param,
        )
        self.check_helper.add_object(
            name="maxRenewableInvestment",
            units=u.MW,
            obj_type=pyo.Param,
        )
        self.check_helper.add_object(
            name="quotaDeficit",
            units=u.MW,
            obj_type=pyo.Var,
        )
        self.check_helper.add_object(
            name="expansionCost",
            units=u.USD,
            obj_type=pyo.Var,
        )
        self.check_helper.add_object(
            name="storageCostInvestment",
            units=u.USD,
            obj_type=pyo.Var,
        )

    def create_testing_obj(self, config):
        """Creates a model and runs tests for a given set of config options."""
        self.investment_stage = 1
        self.m = self._create_model(config=config).model
        self.b = self.m.investmentStage[self.investment_stage]

    # def test_add_investment_param_and_var(self):
    #     self.create_testing_obj(config={})
    #     self.check_helper = PyomoCheckHelper(self, self.b)
    #     self._make_param_and_var_objects()
    #     self.check_helper.check_all_objects()
    #     self.assertEqual(self.b.maxThermalInvestment._default_val, 1000)
    #     self.assertEqual(self.b.maxRenewableInvestment._default_val, 1000)
    #     self.assertEqual(self.b.quotaDeficit.value, 0)
    #     self.assertEqual(self.b.expansionCost.value, 0)
    #     self.assertEqual(self.b.quotaDeficit.domain.name, "NonNegativeReals")
    #     self.assertEqual(self.b.expansionCost.domain.name, "NonNegativeReals")

    def test_add_investment_disjuncts_all_true(self):
        self.create_testing_obj(config={"storage": True, "transmission": True})

        with patch("gtep.model_library.investment.gens.add_generators_status_disjuncts") as mock_gens, \
         patch("gtep.model_library.investment.stor.add_storage_status_disjuncts") as mock_stor, \
         patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts") as mock_transm:

            # call the function you want to test here
            add_investment_disjuncts(self.b)

            mock_gens.assert_called_once_with(
                self.b, self.m.thermalGenerators, self.m.renewableGenerators
            )

            # Check that storage disjuncts are added
            mock_stor.assert_called_once_with(self.b, self.m.storage)

            # Check that transmission disjuncts are added 
            mock_transm.assert_called_once_with(self.b, self.m.lines)

    def test_add_investment_disjuncts_no_storage(self):
        self.create_testing_obj(config={"storage": False, "transmission": True})

        with patch("gtep.model_library.investment.gens.add_generators_status_disjuncts") as mock_gens, \
         patch("gtep.model_library.investment.stor.add_storage_status_disjuncts") as mock_stor, \
         patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts") as mock_transm:

            # call the function you want to test here
            add_investment_disjuncts(self.b)

            mock_gens.assert_called_once_with(
                self.b, self.m.thermalGenerators, self.m.renewableGenerators
            )

            # Check that storage disjuncts are not added to match config
            mock_stor.assert_not_called()

            # Check that transmission disjuncts are added 
            mock_transm.assert_called_once_with(self.b, self.m.lines)

    def test_add_investment_disjuncts_no_transmission(self):
        self.create_testing_obj(config={"storage": True, "transmission": False})

        with patch("gtep.model_library.investment.gens.add_generators_status_disjuncts") as mock_gens, \
         patch("gtep.model_library.investment.stor.add_storage_status_disjuncts") as mock_stor, \
         patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts") as mock_transm:

            # call the function you want to test here
            add_investment_disjuncts(self.b)

            mock_gens.assert_called_once_with(
                self.b, self.m.thermalGenerators, self.m.renewableGenerators
            )

            # Check that storage disjuncts are added to match config
            mock_stor.assert_called_once_with(self.b, self.m.storage)

            # Check that transmission disjuncts are not added 
            mock_transm.assert_not_called()

    def test_add_investment_disjuncts_no_storage_no_transmission(self):
        self.create_testing_obj(config={"storage": False, "transmission": False})

        with patch("gtep.model_library.investment.gens.add_generators_status_disjuncts") as mock_gens, \
         patch("gtep.model_library.investment.stor.add_storage_status_disjuncts") as mock_stor, \
         patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts") as mock_transm:

            # call the function you want to test here
            add_investment_disjuncts(self.b)

            mock_gens.assert_called_once_with(
                self.b, self.m.thermalGenerators, self.m.renewableGenerators
            )

            # check that config dependent functions are not called
            mock_stor.assert_not_called()
            mock_transm.assert_not_called()


### ADD INVESTMENT CONSTRAINTS ###
    def test_add_investment_constraints_all_true(self):

        self.create_testing_obj(config={"storage": True, "transmission": True})

        with patch("gtep.model_library.investment.gens.add_investment_constraints") as mock_gens, \
         patch("gtep.model_library.investment.stor.add_investment_transmission_constraints") as mock_transm, \
         patch("gtep.model_library.investment.transm.add_investment_storage_constraints") as mock_stor:

            # call the function you want to test here
            add_investment_disjuncts(self.b)

            mock_gens.assert_called_once_with(
                self.b, self.m.thermalGenerators, self.m.renewableGenerators
            )

            # Check that storage disjuncts are added
            mock_stor.assert_called_once_with(self.b, self.m.storage)

            # Check that transmission disjuncts are added 
            mock_transm.assert_called_once_with(self.b, self.m.lines)

# @patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
# @patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
# @patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
# @patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
# def test_add_investment_constraints(mock_commitment, mock_transm, mock_stor, mock_gens):

#     add_investment_constraints(inv, 1)

#     # Check external function calls based on config
#     mock_gens.assert_called_once_with(m, inv, 1)
#     mock_transm.assert_called_once_with(m, inv, 1)
#     mock_stor.assert_called_once_with(m, inv, 1)
#     mock_commitment.assert_called_once_with(m, inv, 1)

#     assert hasattr(inv, "investment_cost")
#     expected_baseline = (
#         inv.generators_investment_cost
#         + inv.storage_investment_cost
#         + inv.transmission_investment_cost
#     )
#     expected_value = m.investmentFactor[1] * expected_baseline
#     val = pyo.value(inv.investment_cost)
#     assert val == expected_value

#     # Check renewable_generation_requirement constraint
#     assert hasattr(inv, "renewable_generation_requirement")
#     con = inv.renewable_generation_requirement
#     mock_surplus = 0
#     for rep_per in inv.representativePeriods:
#         for com_per in inv.representativePeriod[rep_per].commitmentPeriods:
#             mock_surplus += (
#                 m.weights[rep_per]
#                 * inv.representativePeriod[rep_per]
#                 .commitmentPeriod[com_per]
#                 .renewableSurplusCommitment
#             )

#     assert mock_surplus + inv.quotaDeficit >= m.renewableQuota[1] * pyo.value(inv.ed)


# @patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
# @patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
# @patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
# @patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
# def test_add_investment_constraints_no_storage(
#     mock_commitment, mock_transm, mock_stor, mock_gens
# ):

#     add_investment_constraints(inv, 1)

#     mock_gens.assert_called_once_with(m, inv, 1)
#     mock_transm.assert_called_once_with(m, inv, 1)
#     mock_stor.assert_not_called()
#     mock_commitment.assert_called_once_with(m, inv, 1)

#     expected_baseline = (
#         inv.generators_investment_cost + 0 + inv.transmission_investment_cost
#     )
#     expected_value = m.investmentFactor[1] * expected_baseline
#     val = pyo.value(inv.investment_cost)
#     assert val == expected_value


# @patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
# @patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
# @patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
# @patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
# def test_no_transmission(mock_commitment, mock_transm, mock_stor, mock_gens):
#     inv, m = create_test_inv(transmission=False)

#     add_investment_constraints(inv, 1)

#     mock_gens.assert_called_once_with(m, inv, 1)
#     mock_transm.assert_not_called()
#     mock_stor.assert_called_once_with(m, inv, 1)
#     mock_commitment.assert_called_once_with(m, inv, 1)

#     expected_baseline = inv.generators_investment_cost + inv.storage_investment_cost + 0
#     expected_value = m.investmentFactor[1] * expected_baseline
#     val = pyo.value(inv.investment_cost)
#     assert val == expected_value


# @patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
# @patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
# @patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
# @patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
# def test_no_include_investment_constraint(
#     mock_commitment, mock_transm, mock_stor, mock_gens
# ):
#     inv, m = create_test_inv(include_investment=False)

#     add_investment_constraints(inv, 1)

#     # The renewable_generation_requirement constraint should NOT be added
#     assert not hasattr(inv, "renewable_generation_requirement")


# @patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
# @patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
# @patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
# @patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
# def test_zero_investment_costs(mock_commitment, mock_transm, mock_stor, mock_gens):
#     inv, m = create_test_inv()
#     inv.generators_investment_cost = 0
#     inv.storage_investment_cost = 0
#     inv.transmission_investment_cost = 0

#     add_investment_constraints(inv, 1)

#     expected_baseline = 0
#     expected_value = m.investmentFactor[1] * expected_baseline
#     val = pyo.value(inv.investment_cost)
#     assert val == expected_value
