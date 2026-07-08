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
from gtep.tests.unit.utils_for_testing import create_model, path_9_bus

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
        self.m = create_model(input_data_path=path_9_bus, config=config).model
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

        with (
            patch(
                "gtep.model_library.investment.gens.add_generators_status_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.stor.add_storage_status_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.transm.add_transmission_status_disjuncts"
            ) as mock_transm,
        ):

            # call the function
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

        with (
            patch(
                "gtep.model_library.investment.gens.add_generators_status_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.stor.add_storage_status_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.transm.add_transmission_status_disjuncts"
            ) as mock_transm,
        ):

            # call the function
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

        with (
            patch(
                "gtep.model_library.investment.gens.add_generators_status_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.stor.add_storage_status_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.transm.add_transmission_status_disjuncts"
            ) as mock_transm,
        ):

            # call the function
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

        with (
            patch(
                "gtep.model_library.investment.gens.add_generators_status_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.stor.add_storage_status_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.transm.add_transmission_status_disjuncts"
            ) as mock_transm,
        ):

            # call the function
            add_investment_disjuncts(self.b)

            mock_gens.assert_called_once_with(
                self.b, self.m.thermalGenerators, self.m.renewableGenerators
            )

            # check that config dependent functions are not called
            mock_stor.assert_not_called()
            mock_transm.assert_not_called()

    ### ADD INVESTMENT CONSTRAINTS ###
    def test_add_investment_constraints_all_true(self):

        self.create_testing_obj(
            config={"storage": True, "transmission": True, "include_investment": True}
        )

        with (
            patch(
                "gtep.model_library.investment.gens.add_investment_generators_constraints"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.transm.add_investment_transmission_constraints"
            ) as mock_transm,
            patch(
                "gtep.model_library.investment.stor.add_investment_storage_constraints"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.commit.add_investment_commitment_constraints"
            ) as mock_comm,
        ):

            # call the function
            add_investment_constraints(self.b, self.investment_stage)

            # check the functions that are always called were called properly
            mock_gens.assert_called_once_with(self.m, self.b, self.investment_stage)
            mock_comm.assert_called_once_with(self.m, self.b, self.investment_stage)

            # Check that transmission constraints are added
            mock_transm.assert_called_once_with(self.m, self.b, self.investment_stage)

            # Check that storage constraints are added
            mock_stor.assert_called_once_with(self.m, self.b, self.investment_stage)

            # check expression and constraint creation
            assert hasattr(self.b, "investment_cost")
            assert hasattr(self.b, "renewable_generation_requirement")
            assert isinstance(self.b.investment_cost, pyo.Expression)
            assert isinstance(self.b.renewable_generation_requirement, pyo.Constraint)

            # check expression config changes
            expr_str = str(self.b.investment_cost.expr)
            assert "storageInvestmentCost" in expr_str
            assert "transmissionCapacity" in expr_str
            assert "investmentStage[1].storInstalled" in expr_str
            assert "investmentStage[1].branchInstalled" in expr_str
            assert "investmentStage[1].storExtended" in expr_str
            assert "investmentStage[1].branchExtended" in expr_str

    def test_add_investment_constraints_all_False(self):

        self.create_testing_obj(
            config={
                "storage": False,
                "transmission": False,
                "include_investment": False,
            }
        )

        with (
            patch(
                "gtep.model_library.investment.gens.add_investment_generators_constraints"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.transm.add_investment_transmission_constraints"
            ) as mock_transm,
            patch(
                "gtep.model_library.investment.stor.add_investment_storage_constraints"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.commit.add_investment_commitment_constraints"
            ) as mock_comm,
        ):

            # call the function
            add_investment_constraints(self.b, self.investment_stage)

            # check the functions that are always called were called properly
            mock_gens.assert_called_once_with(self.m, self.b, self.investment_stage)
            mock_comm.assert_called_once_with(self.m, self.b, self.investment_stage)

            # Check that transmission constraints are not added
            mock_transm.assert_not_called()

            # Check that storage constraints are not added
            mock_stor.assert_not_called()

            # check expression and constraint creation
            assert hasattr(self.b, "investment_cost")
            assert not hasattr(self.b, "renewable_generation_requirement")
            assert isinstance(self.b.investment_cost, pyo.Expression)

            # check expression config changes
            expr_str = str(self.b.investment_cost.expr)
            assert "storageInvestmentCost" not in expr_str
            assert "transmissionCapacity" not in expr_str
            assert "investmentStage[1].storInstalled" not in expr_str
            assert "investmentStage[1].branchInstalled" not in expr_str
            assert "investmentStage[1].storExtended" not in expr_str
            assert "investmentStage[1].branchExtended" not in expr_str

    def test_add_investment_constraints_storage_true(self):

        self.create_testing_obj(
            config={"storage": True, "transmission": False, "include_investment": False}
        )

        with (
            patch(
                "gtep.model_library.investment.gens.add_investment_generators_constraints"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.transm.add_investment_transmission_constraints"
            ) as mock_transm,
            patch(
                "gtep.model_library.investment.stor.add_investment_storage_constraints"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.commit.add_investment_commitment_constraints"
            ) as mock_comm,
        ):

            # call the function
            add_investment_constraints(self.b, self.investment_stage)

            # check the functions that are always called were called properly
            mock_gens.assert_called_once_with(self.m, self.b, self.investment_stage)
            mock_comm.assert_called_once_with(self.m, self.b, self.investment_stage)

            # Check that transmission constraints are added
            mock_transm.assert_not_called()

            # Check that storage constraints are added
            mock_stor.assert_called_once_with(self.m, self.b, self.investment_stage)

            # check expression and constraint creation
            assert hasattr(self.b, "investment_cost")
            assert not hasattr(self.b, "renewable_generation_requirement")
            assert isinstance(self.b.investment_cost, pyo.Expression)

            # check expression config changes
            expr_str = str(self.b.investment_cost.expr)
            assert "storageInvestmentCost" in expr_str
            assert "transmissionCapacity" not in expr_str
            assert "investmentStage[1].storInstalled" in expr_str
            assert "investmentStage[1].branchInstalled" not in expr_str
            assert "investmentStage[1].storExtended" in expr_str
            assert "investmentStage[1].branchExtended" not in expr_str

    def test_add_investment_constraints_transmission_true(self):

        self.create_testing_obj(
            config={"storage": False, "transmission": True, "include_investment": True}
        )

        with (
            patch(
                "gtep.model_library.investment.gens.add_investment_generators_constraints"
            ) as mock_gens,
            patch(
                "gtep.model_library.investment.transm.add_investment_transmission_constraints"
            ) as mock_transm,
            patch(
                "gtep.model_library.investment.stor.add_investment_storage_constraints"
            ) as mock_stor,
            patch(
                "gtep.model_library.investment.commit.add_investment_commitment_constraints"
            ) as mock_comm,
        ):

            # call the function
            add_investment_constraints(self.b, self.investment_stage)

            # check the functions that are always called were called properly
            mock_gens.assert_called_once_with(self.m, self.b, self.investment_stage)
            mock_comm.assert_called_once_with(self.m, self.b, self.investment_stage)

            # Check that transmission constraints are added
            mock_transm.assert_called_once_with(self.m, self.b, self.investment_stage)

            # Check that storage constraints are added
            mock_stor.assert_not_called()

            # check expression and constraint creation
            assert hasattr(self.b, "investment_cost")
            assert hasattr(self.b, "renewable_generation_requirement")
            assert isinstance(self.b.investment_cost, pyo.Expression)
            assert isinstance(self.b.renewable_generation_requirement, pyo.Constraint)

            # check expression config changes
            expr_str = str(self.b.investment_cost.expr)
            assert "storageInvestmentCost" not in expr_str
            assert "transmissionCapacity" in expr_str
            assert "investmentStage[1].storInstalled" not in expr_str
            assert "investmentStage[1].branchInstalled" in expr_str
            assert "investmentStage[1].storExtended" not in expr_str
            assert "investmentStage[1].branchExtended" in expr_str
