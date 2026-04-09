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


# helper function
def create_test_inv():
    m = pyo.ConcreteModel("GTEP Model")
    m.stages = 1
    m.investmentStage = pyo.Block(m.stages)
    inv = m.investmentStage[1]
    return inv, 1
