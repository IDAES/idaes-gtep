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
from unittest.mock import MagicMock

commitment_period = 1
investment_stage = 1


import pyomo.environ as pyo
from unittest.mock import MagicMock


def create_test_inv(
    include_commit=True,
    storage=True,
    scale_loads=True,
    p_max_data=None,
):
    m = pyo.ConcreteModel()
    m.stages = [1]
    m.config = {
        "include_commitment": include_commit,
        "storage": storage,
        "scale_loads": scale_loads,
    }
    m.investmentFactor = {1: 1.0}
    m.weights = {1: 1.0}

    m.regions = pyo.Set(initialize=["R1", "R2"])
    m.thermalGenerators = pyo.Set(initialize=["Gen1", "Gen2"])
    m.renewableGenerators = pyo.Set(initialize=["Ren1", "Ren2"])

    m.investmentStage = pyo.Block(m.stages)
    inv = m.investmentStage[1]

    # representativePeriods is a list of indices
    inv.representativePeriods = [1]
    inv.representativePeriod = pyo.Block(inv.representativePeriods)

    # Declare commitmentPeriods as an indexed Block inside representativePeriod[1]
    inv.representativePeriod[1].commitmentPeriods = pyo.Block([1])

    # Access the block for commitmentPeriod 1 (this initializes it)
    b_comm = inv.representativePeriod[1].commitmentPeriods[1]

    # Add any attributes you want to this block
    b_comm.renewableSurplusCommitment = 10

    if p_max_data is None:
        p_max_data = {
            "Ren1": 100.0,
            "Ren2": {"values": [50.0, 60.0]},
        }

    m.md = MagicMock()
    m.md.data = {
        "elements": {
            "generator": {
                gen: {"p_max": p_max_data.get(gen, 0)} for gen in m.renewableGenerators
            }
        }
    }

    return b_comm, inv, m


def test_add_commitment_parameters_basic():
    b_comm, inv, m = create_test_inv()

    # Call function with commitment_period=2 (index 1)
    add_commitment_parameters(b_comm, commitment_period, investment_stage)

    # Check parameters exist and have correct defaults
    assert isinstance(b_comm.commitmentPeriodLength, pyo.Param)
    assert b_comm.commitmentPeriodLength.value == 1
    assert isinstance(b_comm.carbonTax, pyo.Param)
    assert b_comm.carbonTax.value == 0

    # Check renewableCapacityExpected keys and values
    assert "Ren1" in b_comm.renewableCapacityExpected
    assert "Ren2" in b_comm.renewableCapacityExpected

    # Ren1 p_max is float, so expected capacity should be 0 * units
    assert b_comm.renewableCapacityExpected["Ren1"].magnitude == 0

    # Ren2 p_max is dict, so expected capacity should be 60.0 * units
    assert b_comm.renewableCapacityExpected["Ren2"].magnitude == 60.0
