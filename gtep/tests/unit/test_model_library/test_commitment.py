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

    ### ADD COMMITMENT PARAMETERS ###
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

    ### ADD COMMITMENT DISJUNCTS ###
    def test_add_commitment_disjuncts_commitment_and_storage_true(self):
        self._create_testing_obj(
            config={"include_commitment": True, "storage": True}, data_path=path_9_bus
        )  # 9 bus used for storage data

        with (
            patch(
                "gtep.model_library.commitment.gens.add_generators_state_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.commitment.stor.add_storage_state_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.commitment.gens.generators_status_always_on"
            ) as mock_always_on,
        ):
            add_commitment_disjuncts(self.b, self.commit_period)
            # check func calls based on config
            mock_gens.assert_called_once_with(
                self.b.model(),
                self.b,
                self.b.parent_block(),
                self.b.parent_block().parent_block(),
                self.commit_period,
            )
            mock_stor.assert_called_once_with(
                self.b.model(), self.b, self.commit_period
            )
            mock_always_on.assert_not_called()

    def test_add_commitment_disjuncts_no_commitment_and_storage_true(self):
        self._create_testing_obj(
            config={"include_commitment": False, "storage": True}, data_path=path_9_bus
        )  # 9 bus used for storage data

        with (
            patch(
                "gtep.model_library.commitment.gens.add_generators_state_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.commitment.stor.add_storage_state_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.commitment.gens.generators_status_always_on"
            ) as mock_always_on,
        ):
            add_commitment_disjuncts(self.b, self.commit_period)
            # check func calls based on config
            mock_gens.assert_not_called()
            mock_stor.assert_called_once_with(
                self.b.model(), self.b, self.commit_period
            )
            mock_always_on.assert_called_once_with(
                self.b.model(),
                self.b,
                self.b.parent_block(),
                self.b.parent_block().parent_block(),
                self.commit_period,
            )

    def test_add_commitment_disjuncts_commitment_true_no_storage(self):
        self._create_testing_obj(
            config={"include_commitment": True, "storage": False},
            data_path=input_data_path,
        )

        with (
            patch(
                "gtep.model_library.commitment.gens.add_generators_state_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.commitment.stor.add_storage_state_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.commitment.gens.generators_status_always_on"
            ) as mock_always_on,
        ):
            add_commitment_disjuncts(self.b, self.commit_period)
            # check func calls based on config
            mock_gens.assert_called_once_with(
                self.b.model(),
                self.b,
                self.b.parent_block(),
                self.b.parent_block().parent_block(),
                self.commit_period,
            )
            mock_stor.assert_not_called()
            mock_always_on.assert_not_called()

    def test_add_commitment_disjuncts_no_commitment_no_storage(self):
        self._create_testing_obj(
            config={"include_commitment": False, "storage": False},
            data_path=input_data_path,
        )

        with (
            patch(
                "gtep.model_library.commitment.gens.add_generators_state_disjuncts"
            ) as mock_gens,
            patch(
                "gtep.model_library.commitment.stor.add_storage_state_disjuncts"
            ) as mock_stor,
            patch(
                "gtep.model_library.commitment.gens.generators_status_always_on"
            ) as mock_always_on,
        ):
            add_commitment_disjuncts(self.b, self.commit_period)
            # check func calls based on config
            mock_gens.assert_not_called()
            mock_stor.assert_not_called()
            mock_always_on.assert_called_once_with(
                self.b.model(),
                self.b,
                self.b.parent_block(),
                self.b.parent_block().parent_block(),
                self.commit_period,
            )

    # ### ADD COMMITMENT CONSTRAINTS ###
    def test_add_commitment_constraints_all_false(self):
        self._create_testing_obj(
            config={
                "storage": False,
                "include_commitment": False,
                "advanced_hydro": False,
            },
            data_path=input_data_path,
        )

        add_commitment_constraints(self.b, self.commit_period)

        assert hasattr(self.b, "renewableSurplusCommitment")
        assert hasattr(self.b, "operatingCostCommitment")
        assert hasattr(self.b, "renewableCurtailmentCommitment")

        assert isinstance(self.b.renewableSurplusCommitment, pyo.Expression)
        assert isinstance(self.b.operatingCostCommitment, pyo.Expression)
        assert isinstance(self.b.renewableCurtailmentCommitment, pyo.Expression)

        op_cost_str = str(self.b.operatingCostCommitment.expr)
        assert "genShutdown" not in op_cost_str
        assert "genStartup" not in op_cost_str
        assert "storagefixedCost" not in op_cost_str

    def test_add_commitment_constraints_no_storage(self):
        self._create_testing_obj(
            config={
                "advanced_hydro": True,
                "include_commitment": True,
                "storage": False,
            },
            data_path=input_data_path,
        )  # 9 bus used for storage data

        # check function call for storage config option
        with patch(
            "gtep.model_library.commitment.stor.add_commitment_storage_constraints"
        ) as mock_stor:
            add_commitment_constraints(self.b, self.commit_period)

            mock_stor.assert_not_called()

        # check expressions exist
        assert hasattr(self.b, "renewableSurplusCommitment")
        assert hasattr(self.b, "operatingCostCommitment")
        assert hasattr(self.b, "renewableCurtailmentCommitment")

        assert isinstance(self.b.renewableSurplusCommitment, pyo.Expression)
        assert isinstance(self.b.operatingCostCommitment, pyo.Expression)
        assert isinstance(self.b.renewableCurtailmentCommitment, pyo.Expression)

        # check expression strings
        surplus_str = str(self.b.renewableSurplusCommitment.expr)
        curtail_str = str(self.b.renewableCurtailmentCommitment.expr)
        op_cost_str = str(self.b.operatingCostCommitment.expr)

        # check what is included in the expressions
        assert "renewableSurplusDispatch" in surplus_str
        assert "renewableCurtailmentDispatch" in curtail_str
        assert "genShutdown" in op_cost_str
        assert "genStartup" in op_cost_str

        # no storage
        assert "storagefixedCost" not in op_cost_str
        assert "storageCapacity" not in op_cost_str

        # if advanced hydro true
        assert "hydroCapacity" in op_cost_str

    def test_add_commitment_constraints_no_storage_no_commitment(self):
        self._create_testing_obj(
            config={
                "advanced_hydro": True,
                "include_commitment": False,
                "storage": False,
            },
            data_path=input_data_path,
        )  # 9 bus used for storage data

        # check function call for storage config option
        with patch(
            "gtep.model_library.commitment.stor.add_commitment_storage_constraints"
        ) as mock_stor:
            add_commitment_constraints(self.b, self.commit_period)

            mock_stor.assert_not_called()

        # check expressions exist
        assert hasattr(self.b, "renewableSurplusCommitment")
        assert hasattr(self.b, "operatingCostCommitment")
        assert hasattr(self.b, "renewableCurtailmentCommitment")

        assert isinstance(self.b.renewableSurplusCommitment, pyo.Expression)
        assert isinstance(self.b.operatingCostCommitment, pyo.Expression)
        assert isinstance(self.b.renewableCurtailmentCommitment, pyo.Expression)

        # check expression strings
        surplus_str = str(self.b.renewableSurplusCommitment.expr)
        curtail_str = str(self.b.renewableCurtailmentCommitment.expr)
        op_cost_str = str(self.b.operatingCostCommitment.expr)

        # check what is included in the expressions
        assert "renewableSurplusDispatch" in surplus_str
        assert "renewableCurtailmentDispatch" in curtail_str
        assert "genShutdown" not in op_cost_str
        assert "genStartup" not in op_cost_str

        # no storage
        assert "storagefixedCost" not in op_cost_str
        assert "storageCapacity" not in op_cost_str

        # if advanced hydro true
        assert "hydroCapacity" in op_cost_str

    def test_add_commitment_constraints_no_hydro(self):
        self._create_testing_obj(
            config={
                "advanced_hydro": False,
                "include_commitment": True,
                "storage": True,
            },
            data_path=path_9_bus,
        )  # 9 bus used for storage data

        # check function call for storage config option
        with patch(
            "gtep.model_library.commitment.stor.add_commitment_storage_constraints"
        ) as mock_stor:
            add_commitment_constraints(self.b, self.commit_period)

            mock_stor.assert_called_once_with(self.b)

        # check expressions exist
        assert hasattr(self.b, "renewableSurplusCommitment")
        assert hasattr(self.b, "operatingCostCommitment")
        assert hasattr(self.b, "renewableCurtailmentCommitment")

        assert isinstance(self.b.renewableSurplusCommitment, pyo.Expression)
        assert isinstance(self.b.operatingCostCommitment, pyo.Expression)
        assert isinstance(self.b.renewableCurtailmentCommitment, pyo.Expression)

        # check expression strings
        surplus_str = str(self.b.renewableSurplusCommitment.expr)
        curtail_str = str(self.b.renewableCurtailmentCommitment.expr)
        op_cost_str = str(self.b.operatingCostCommitment.expr)

        # check what is included in the expressions
        assert "renewableSurplusDispatch" in surplus_str
        assert "renewableCurtailmentDispatch" in curtail_str
        assert "genShutdown" in op_cost_str
        assert "genStartup" in op_cost_str

        # if storage true
        assert "storagefixedCost" in op_cost_str

        # no advanced hydro
        assert "hydroCapacity" not in op_cost_str

    def test_add_commitment_constraints_no_hydro_no_commitment(self):
        self._create_testing_obj(
            config={
                "advanced_hydro": False,
                "include_commitment": False,
                "storage": True,
            },
            data_path=path_9_bus,
        )  # 9 bus used for storage data

        # check function call for storage config option
        with patch(
            "gtep.model_library.commitment.stor.add_commitment_storage_constraints"
        ) as mock_stor:
            add_commitment_constraints(self.b, self.commit_period)

            mock_stor.assert_called_once_with(self.b)

        # check expressions exist
        assert hasattr(self.b, "renewableSurplusCommitment")
        assert hasattr(self.b, "operatingCostCommitment")
        assert hasattr(self.b, "renewableCurtailmentCommitment")

        assert isinstance(self.b.renewableSurplusCommitment, pyo.Expression)
        assert isinstance(self.b.operatingCostCommitment, pyo.Expression)
        assert isinstance(self.b.renewableCurtailmentCommitment, pyo.Expression)

        # check expression strings
        surplus_str = str(self.b.renewableSurplusCommitment.expr)
        curtail_str = str(self.b.renewableCurtailmentCommitment.expr)
        op_cost_str = str(self.b.operatingCostCommitment.expr)

        # check what is included in the expressions
        assert "renewableSurplusDispatch" in surplus_str
        assert "renewableCurtailmentDispatch" in curtail_str
        assert "genShutdown" not in op_cost_str
        assert "genStartup" not in op_cost_str

        # if storage true
        assert "storagefixedCost" in op_cost_str

        # no advanced hydro
        assert "hydroCapacity" not in op_cost_str

    def test_add_commitment_constraints_no_hydro_no_storage(self):
        self._create_testing_obj(
            config={
                "advanced_hydro": False,
                "include_commitment": True,
                "storage": False,
            },
            data_path=input_data_path,
        )

        # check function call for storage config option
        with patch(
            "gtep.model_library.commitment.stor.add_commitment_storage_constraints"
        ) as mock_stor:
            add_commitment_constraints(self.b, self.commit_period)

            mock_stor.assert_not_called()

        # check expressions exist
        assert hasattr(self.b, "renewableSurplusCommitment")
        assert hasattr(self.b, "operatingCostCommitment")
        assert hasattr(self.b, "renewableCurtailmentCommitment")

        assert isinstance(self.b.renewableSurplusCommitment, pyo.Expression)
        assert isinstance(self.b.operatingCostCommitment, pyo.Expression)
        assert isinstance(self.b.renewableCurtailmentCommitment, pyo.Expression)

        # check expression strings
        surplus_str = str(self.b.renewableSurplusCommitment.expr)
        curtail_str = str(self.b.renewableCurtailmentCommitment.expr)
        op_cost_str = str(self.b.operatingCostCommitment.expr)

        # check what is included in the expressions
        assert "renewableSurplusDispatch" in surplus_str
        assert "renewableCurtailmentDispatch" in curtail_str
        assert "genOn" in op_cost_str
        assert "genShutdown" in op_cost_str
        assert "genStartup" in op_cost_str

        # no storage
        assert "storagefixedCost" not in op_cost_str
        assert "storageCapacity" not in op_cost_str

        # no advanced hydro
        assert "hydroCapacity" not in op_cost_str

    ### ADD INVESTMENT COMMITMENT VARIABLES ###
    def _make_investment_commitment_var_objects(self):
        self.check_helper.add_object(
            name="operatingCostInvestment",
            units=u.USD,
            obj_type=pyo.Var,
        )
        self.check_helper.add_object(
            name="renewableCurtailmentInvestment",
            units=u.USD,
            obj_type=pyo.Var,
        )

    def test_add_investment_commitment_variables(self):
        self._create_testing_obj(config={}, data_path=input_data_path)

        self.check_helper = PyomoCheckHelper(self, self.b)
        self._make_investment_commitment_var_objects()

        add_investment_commitment_variables(self.b)

        self.check_helper.check_all_objects()

    ### ADD INVESTMENT COMMITMENT CONSTRAINTS ###
    def _make_investment_commitment_constraint_objects(self):
        self.check_helper.add_object(
            name="renewableCurtailmentInvestment",
            units=u.USD,
            obj_type=pyo.Var,
        )
        self.check_helper.add_object(
            name="operatingCostInvestment",
            units=u.USD,
            obj_type=pyo.Var,
        )

    def check_commitment_constraint_objects(self):
        self._create_testing_obj(
            config={},
            data_path=input_data_path,
        )
        # this function takes in an investment block not a commitment block
        b = self.m.investmentStage[self.investment_stage]
        # Add variables required by the constraints
        add_investment_commitment_variables(b)
        # Set up helper before function call
        self.check_helper = PyomoCheckHelper(self, b)
        self._make_investment_commitment_constraint_objects()

        add_investment_commitment_constraints(
            self.m,
            b,
            self.investment_stage,
        )

        # Check that the variables are still present and correctly typed
        self.check_helper.check_all_objects()

    def test_add_investment_commitment_constraints(self):
        self._create_testing_obj(
            config={},
            data_path=input_data_path,
        )
        # this function takes in an investment block not a commitment block
        b = self.m.investmentStage[self.investment_stage]

        # add variables needed by the constraints
        add_investment_commitment_variables(b)

        add_investment_commitment_constraints(self.m, b, self.investment_stage)

        # check expression and constraint creation
        assert hasattr(b, "renewable_curtailment_cost")
        assert hasattr(b, "operatingCostInvestment_constraint")
        assert isinstance(b.renewable_curtailment_cost, pyo.Constraint)
        assert isinstance(b.operatingCostInvestment_constraint, pyo.Constraint)

        # check expression structure
        curt_expr_str = str(b.renewable_curtailment_cost.expr)
        op_expr_str = str(b.operatingCostInvestment_constraint.expr)

        # renewable curtailment constraint should use these terms
        assert "renewableCurtailmentInvestment" in curt_expr_str
        assert "curtailmentCost" in curt_expr_str
        assert "investmentFactor" in curt_expr_str
        assert "dispatchPeriod" in curt_expr_str
        assert "renewableCurtailment" in curt_expr_str

        # operating cost constraint should use these terms
        assert "operatingCostInvestment" in op_expr_str
        assert "investmentFactor" in op_expr_str
        assert "fixedCost" in op_expr_str
        assert "startupCost" in op_expr_str
        assert "genOn" in op_expr_str
        assert "genShutdown" in op_expr_str
        assert "genStartup" in op_expr_str
        assert "dispatchPeriod" in op_expr_str
        assert "renewableCurtailment" in op_expr_str

        # check that the structure includes the expected indexing pattern
        assert "representativePeriod[1]" in curt_expr_str
        assert "commitmentPeriod[1]" in curt_expr_str
        assert "representativePeriod[1]" in op_expr_str
        assert "commitmentPeriod[1]" in op_expr_str
