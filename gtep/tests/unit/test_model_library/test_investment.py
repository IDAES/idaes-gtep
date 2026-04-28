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


# helper function
def create_test_inv(
    include_commit=True, storage=True, transmission=True, include_investment=True
):
    m = pyo.ConcreteModel("GTEP Model")
    m.stages = [1]
    m.config = {
        "include_commitment": include_commit,
        "storage": storage,
        "transmission": transmission,
        "include_investment": include_investment,
    }
    m.investmentFactor = {1: 1.0}
    m.weights = {1: 1.0}
    m.renewableQuota = {1: 0.5}

    class CommitmentPeriod:
        def __init__(self):
            self.renewableSurplusCommitment = 10

    class RepresentativePeriod:
        def __init__(self):
            self.commitmentPeriods = [1]
            self.commitmentPeriod = {1: CommitmentPeriod()}

    m.regions = pyo.Set(initialize=["R1", "R2"])
    m.thermalGenerators = pyo.Set(initialize=["Gen1", "Gen2"])
    m.renewableGenerators = pyo.Set(initialize=["Ren1", "Ren2"])
    m.storage = pyo.Set(initialize=["Stor1"])
    m.transmission = pyo.Set(initialize=["Line1"])

    m.investmentStage = pyo.Block(m.stages)
    inv = m.investmentStage[1]

    # constraint attributes
    inv.representativePeriods = [1]
    inv.representativePeriod = {1: RepresentativePeriod()}

    inv.generators_investment_cost = 100
    inv.storage_investment_cost = 10
    inv.transmission_investment_cost = 10
    inv.commitmentOperatingCostInvestment = 20
    inv.operatingCostInvestment = pyo.Var(initialize=0)
    inv.quotaDeficit = 0

    inv.ed = pyo.Param(initialize=0, mutable=True)

    return inv, m


# ------------------------------------ADD INVESTMENT PARAMS AND VARIABLES------------------------------------ #
def test_investment_params_and_variables():
    # Ensure USD unit is defined
    u._pint_registry.define("USD = []")

    b, m = create_test_inv()
    add_investment_params_and_variables(b, 1)

    # Check parameters
    assert hasattr(b, "maxThermalInvestment")
    assert hasattr(b, "maxRenewableInvestment")
    # Check default values for parameters for each region
    for region in m.regions:
        assert pyo.value(b.maxThermalInvestment[region]) == 1000
        assert pyo.value(b.maxRenewableInvestment[region]) == 1000

    # Check variables
    assert hasattr(b, "quotaDeficit")
    assert hasattr(b, "expansionCost")
    assert hasattr(b, "storageCostInvestment")

    # Check variable domains and initial values
    assert b.quotaDeficit.domain == pyo.NonNegativeReals
    assert pyo.value(b.quotaDeficit) == 0
    assert b.expansionCost.domain == pyo.NonNegativeReals
    assert pyo.value(b.expansionCost) == 0
    assert b.storageCostInvestment.domain == pyo.NonNegativeReals
    assert pyo.value(b.storageCostInvestment) == 0


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


# ------------------------------------ADD INVESTMENT CONSTRAINTS------------------------------------ #
@patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
@patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
@patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
@patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
def test_add_investment_constraints(mock_commitment, mock_transm, mock_stor, mock_gens):
    inv, m = create_test_inv()

    add_investment_constraints(inv, 1)

    # Check external function calls based on config
    mock_gens.assert_called_once_with(m, inv, 1)
    mock_transm.assert_called_once_with(m, inv, 1)
    mock_stor.assert_called_once_with(m, inv, 1)
    mock_commitment.assert_called_once_with(m, inv, 1)

    assert hasattr(inv, "investment_cost")
    expected_baseline = (
        inv.generators_investment_cost
        + inv.storage_investment_cost
        + inv.transmission_investment_cost
    )
    expected_value = m.investmentFactor[1] * expected_baseline
    val = pyo.value(inv.investment_cost)
    assert val == expected_value

    # Check renewable_generation_requirement constraint
    assert hasattr(inv, "renewable_generation_requirement")
    con = inv.renewable_generation_requirement
    mock_surplus = 0
    for rep_per in inv.representativePeriods:
        for com_per in inv.representativePeriod[rep_per].commitmentPeriods:
            mock_surplus += (
                m.weights[rep_per]
                * inv.representativePeriod[rep_per]
                .commitmentPeriod[com_per]
                .renewableSurplusCommitment
            )

    assert mock_surplus + inv.quotaDeficit >= m.renewableQuota[1] * pyo.value(inv.ed)


@patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
@patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
@patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
@patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
def test_add_investment_constraints_no_storage(
    mock_commitment, mock_transm, mock_stor, mock_gens
):
    inv, m = create_test_inv(storage=False)

    add_investment_constraints(inv, 1)

    mock_gens.assert_called_once_with(m, inv, 1)
    mock_transm.assert_called_once_with(m, inv, 1)
    mock_stor.assert_not_called()
    mock_commitment.assert_called_once_with(m, inv, 1)

    expected_baseline = (
        inv.generators_investment_cost + 0 + inv.transmission_investment_cost
    )
    expected_value = m.investmentFactor[1] * expected_baseline
    val = pyo.value(inv.investment_cost)
    assert val == expected_value


@patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
@patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
@patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
@patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
def test_no_transmission(mock_commitment, mock_transm, mock_stor, mock_gens):
    inv, m = create_test_inv(transmission=False)

    add_investment_constraints(inv, 1)

    mock_gens.assert_called_once_with(m, inv, 1)
    mock_transm.assert_not_called()
    mock_stor.assert_called_once_with(m, inv, 1)
    mock_commitment.assert_called_once_with(m, inv, 1)

    expected_baseline = inv.generators_investment_cost + inv.storage_investment_cost + 0
    expected_value = m.investmentFactor[1] * expected_baseline
    val = pyo.value(inv.investment_cost)
    assert val == expected_value


@patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
@patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
@patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
@patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
def test_no_include_investment_constraint(
    mock_commitment, mock_transm, mock_stor, mock_gens
):
    inv, m = create_test_inv(include_investment=False)

    add_investment_constraints(inv, 1)

    # The renewable_generation_requirement constraint should NOT be added
    assert not hasattr(inv, "renewable_generation_requirement")


@patch("gtep.model_library.investment.gens.add_investment_generators_constraints")
@patch("gtep.model_library.investment.stor.add_investment_storage_constraints")
@patch("gtep.model_library.investment.transm.add_investment_transmission_constraints")
@patch("gtep.model_library.investment.commit.add_investment_commitment_constraints")
def test_zero_investment_costs(mock_commitment, mock_transm, mock_stor, mock_gens):
    inv, m = create_test_inv()
    inv.generators_investment_cost = 0
    inv.storage_investment_cost = 0
    inv.transmission_investment_cost = 0

    add_investment_constraints(inv, 1)

    expected_baseline = 0
    expected_value = m.investmentFactor[1] * expected_baseline
    val = pyo.value(inv.investment_cost)
    assert val == expected_value
