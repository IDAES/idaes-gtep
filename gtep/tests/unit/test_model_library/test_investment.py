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


# Create a custom pint UnitRegistry with USD and MW defined
ureg = UnitRegistry()
ureg.define("USD = [currency]")
ureg.define("MW = megawatt")

# Create Pyomo units container with this registry
custom_units = u._PyomoUnitsContainer(ureg)


# helper function
def create_test_inv(include_commit=True):
    m = pyo.ConcreteModel("GTEP Model")
    m.stages = [1]
    m.config["include_commitment"] = include_commit
    m.config["storage"] = True

    m.regions = pyo.Set(initialize=["R1", "R2"])
    m.thermalGenerators = pyo.Set(initialize=["Gen1", "Gen2"])

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


def test_add_investment_disjuncts():
    pass
