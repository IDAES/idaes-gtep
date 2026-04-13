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
import gtep.model_library.investment
from gtep.model_library.investment import (
    add_investment_params_and_variables,
    add_investment_disjuncts,
    add_investment_constraints,
)
from pyomo.environ import units as u
import pyomo.environ as pyo
import pytest
import pyomo.common.unittest as unittest
from unittest.mock import patch
from pint import UnitRegistry


# helper function
def create_test_inv(include_commit=True, storage=True, transmission=True):
    m = pyo.ConcreteModel("GTEP Model")
    m.stages = [1]
    m.config = {
        "include_commitment": include_commit,
        "storage": storage,
        "transmission": transmission,
    }

    m.regions = pyo.Set(initialize=["R1", "R2"])
    m.thermalGenerators = pyo.Set(initialize=["Gen1", "Gen2"])
    m.renewableGenerators = pyo.Set(initialize=["Ren1", "Ren2"])
    m.storage = pyo.Set(initialize=["Stor1"])
    m.transmission = pyo.Set(initialize=["Line1"])

    m.investmentStage = pyo.Block(m.stages)
    inv = m.investmentStage[1]
    return inv, m


# TODO Unit object is throwing an error
# @patch("gtep.model_library.commitment")
# def test_investment_params_and_variables(mock_commit_variable):
#     with patch.object(gtep.model_library.investment, "u", new=custom_units):
#         b, m = create_test_inv()
#         add_investment_params_and_variables(b, 1)

#         mock_commit_variable.assert_called_once()

#         # Check parameters
#         assert hasattr(b, "maxThermalInvestment")
#         assert hasattr(b, "maxRenewableInvestment")
#         # Check default values for parameters for each region
#         for region in m.regions:
#             assert pyo.value(b.maxThermalInvestment[region]) == 1000
#             assert pyo.value(b.maxRenewableInvestment[region]) == 1000

#         # Check variables
#         assert hasattr(b, "quotaDeficit")
#         assert hasattr(b, "expansionCost")
#         assert hasattr(b, "storageCostInvestment")

#         # Check variable domains and initial values
#         assert b.quotaDeficit.domain == pyo.NonNegativeReals
#         assert pyo.value(b.quotaDeficit) == 0
#         assert b.expansionCost.domain == pyo.NonNegativeReals
#         assert pyo.value(b.expansionCost) == 0
#         assert b.storageCostInvestment.domain == pyo.NonNegativeReals
#         assert pyo.value(b.storageCostInvestment) == 0


# ------------------------------------ADD INVESTMENT DISJUNTS------------------------------------ #
@patch("gtep.model_library.investment.gens.add_generators_status_disjuncts")
@patch("gtep.model_library.investment.stor.add_storage_status_disjuncts")
@patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts")
def test_add_investment_disjuncts(mock_transm, mock_stor, mock_gens):
    inv, m = create_test_inv(storage=True, transmission=True)

    # Call the function under test
    add_investment_disjuncts(inv)

    # Check that generator disjuncts are always added
    mock_gens.assert_called_once_with(inv, m.thermalGenerators, m.renewableGenerators)

    # Check that storage disjuncts are added if storage config is True
    mock_stor.assert_called_once_with(inv, m.storage)

    # Check that transmission disjuncts are added if transmission config is True
    mock_transm.assert_called_once_with(inv, m.transmission)


@patch("gtep.model_library.investment.gens.add_generators_status_disjuncts")
@patch("gtep.model_library.investment.stor.add_storage_status_disjuncts")
@patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts")
def test_add_investment_disjuncts_no_storage(mock_transm, mock_stor, mock_gens):
    inv, m = create_test_inv(storage=False, transmission=True)

    add_investment_disjuncts(inv)

    mock_gens.assert_called_once_with(inv, m.thermalGenerators, m.renewableGenerators)
    mock_stor.assert_not_called()
    mock_transm.assert_called_once_with(inv, m.transmission)


@patch("gtep.model_library.investment.gens.add_generators_status_disjuncts")
@patch("gtep.model_library.investment.stor.add_storage_status_disjuncts")
@patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts")
def test_add_investment_disjuncts_no_transmission(mock_transm, mock_stor, mock_gens):
    inv, m = create_test_inv(storage=True, transmission=False)

    add_investment_disjuncts(inv)

    mock_gens.assert_called_once_with(inv, m.thermalGenerators, m.renewableGenerators)
    mock_stor.assert_called_once_with(inv, m.storage)
    mock_transm.assert_not_called()


@patch("gtep.model_library.investment.gens.add_generators_status_disjuncts")
@patch("gtep.model_library.investment.stor.add_storage_status_disjuncts")
@patch("gtep.model_library.investment.transm.add_transmission_status_disjuncts")
def test_add_investment_disjuncts_no_storage_no_transmission(
    mock_transm, mock_stor, mock_gens
):
    inv, m = create_test_inv(storage=False, transmission=False)

    add_investment_disjuncts(inv)

    mock_gens.assert_called_once_with(inv, m.thermalGenerators, m.renewableGenerators)
    mock_stor.assert_not_called()
    mock_transm.assert_not_called()
