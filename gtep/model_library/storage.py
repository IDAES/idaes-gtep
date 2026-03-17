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

"""Constraints for the Storage Component in the Generation and
Transmission Expansion Planning (GTEP) Model

"""

import pyomo.environ as pyo


def add_storage_params(m):
    """Battery Storage properties read-in from data"""

    m.storageCapacity = {
        bat: m.md.data["elements"]["storage"][bat]["energy_capacity"]
        for bat in m.storage
    }  # maximum storage capacity

    m.initStorageChargeLevel = {
        bat: m.md.data["elements"]["storage"][bat]["initial_state_of_charge"]
        for bat in m.storage
    }  # initial storage capacity

    m.minStorageChargeLevel = {
        bat: m.md.data["elements"]["storage"][bat]["minimum_state_of_charge"]
        for bat in m.storage
    }  # minimum storage capacity

    m.chargingCost = {
        bat: m.md.data["elements"]["storage"][bat]["charge_cost"] for bat in m.storage
    }  # cost to charge per unit electricity

    m.dischargingCost = {
        bat: m.md.data["elements"]["storage"][bat]["discharge_cost"]
        for bat in m.storage
    }  # cost to discharge per unit electricity

    m.dischargeMin = {
        bat: m.md.data["elements"]["storage"][bat]["min_discharge_rate"]
        for bat in m.storage
    }  # minimum amount to discharge per dispatch period when discharging

    m.dischargeMax = {
        bat: m.md.data["elements"]["storage"][bat]["max_discharge_rate"]
        for bat in m.storage
    }  # maximum amount to discharge per dispatch period when discharging

    m.chargeMin = {
        bat: m.md.data["elements"]["storage"][bat]["min_charge_rate"]
        for bat in m.storage
    }  # minimum amount to charge per dispatch period when charging

    m.chargeMax = {
        bat: m.md.data["elements"]["storage"][bat]["max_charge_rate"]
        for bat in m.storage
    }  # maximum amount to charge per dispatch period when charging

    m.storageDischargingRampUpRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_up_output_60min"]
        for bat in m.storage
    }  # maximum amount of ramp up between dispatch periods when discharging.
    # Notice that default EGRET naming convention assumes dispatch periods are 60 minutes

    m.storageDischargingRampDownRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_down_output_60min"]
        for bat in m.storage
    }  # maximum amount of ramp down between dispatch periods when discharging.

    m.storageChargingRampUpRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_up_input_60min"]
        for bat in m.storage
    }  # maximum amount of ramp up between dispatch periods when charging.

    m.storageChargingRampDownRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_down_input_60min"]
        for bat in m.storage
    }  # maximum amount of ramp down between dispatch periods when charging.

    m.storageDischargingEfficiency = {
        bat: m.md.data["elements"]["storage"][bat]["discharge_efficiency"]
        for bat in m.storage
    }  # proportion of energy discharged that is not lost to technological
    # inefficiencies with in dispatch periods and which is usable in the flow balance

    m.storageChargingEfficiency = {
        bat: m.md.data["elements"]["storage"][bat]["charge_efficiency"]
        for bat in m.storage
    }  # proportion of energy charged that is not lost to technological
    # inefficiencies within dispatch periods and which is usable in the flow balance

    m.storageRetentionRate = {
        bat: m.md.data["elements"]["storage"][bat]["retention_rate_60min"]
        for bat in m.storage
    }  # proportion of energy discharged that is not lost to technological
    # inefficiencies between dispatch periods and which is usable in the flow balance

    # (Arbitrary) multiplier for new battery investments corresponds to depreciation schedules
    # for individual technologies; higher values are indicative of slow depreciation
    m.storageCapitalMultiplier = {
        bat: m.md.data["elements"]["storage"][bat]["capital_multiplier"]
        for bat in m.storage
    }

    # Cost of life extension for each battery, expressed as a fraction of initial investment cost
    m.storageExtensionMultiplier = {
        bat: m.md.data["elements"]["storage"][bat]["extension_multiplier"]
        for bat in m.storage
    }

    m.storageInvestmentCost = {
        bat: m.md.data["elements"]["storage"][bat]["investment_cost"]
        for bat in m.storage
    }  # Future not real cost: idealized DoE 10-yr targets or something


def add_storage_constraints(m, b, commitment_period):
    """
    Battery Discharging Constraints
    """

    r_p = b.parent_block()
    i_p = r_p.parent_block()

    @b.Disjunct(m.storage)
    def storDischarging(disj, bat):
        # operating limits
        b = disj.parent_block()

        # Minimum operating Limits if storage unit is on
        @disj.Constraint(b.dispatchPeriods)
        def discharge_limit_min(d, disp_per):
            return (
                m.dischargeMin[bat]  # Assuming dischargeMin is an absolute value (MW)
                <= b.dispatchPeriod[disp_per].storageDischarged[bat]
            )

        # Maximum operating limits
        @disj.Constraint(b.dispatchPeriods)
        def discharge_limit_max(d, disp_per):
            return (
                b.dispatchPeriod[disp_per].storageDischarged[bat] <= m.dischargeMax[bat]
            )

        # Ramp up limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def discharge_ramp_up_limits(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageDischarged[bat]
                    - b.dispatchPeriod[disp_per - 1].storageDischarged[bat]
                    <= m.storageDischargingRampUpRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageDischarged[bat]
                    - r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .storageDischarged[bat]
                    <= m.storageDischargingRampUpRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            else:
                return pyo.Constraint.Skip

        # Ramp down limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def discharge_ramp_down_limits(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per - 1].storageDischarged[bat]
                    - b.dispatchPeriod[disp_per].storageDischarged[bat]
                    <= m.storageDischargingRampDownRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods[-1]]
                    .storageDischarged[bat]
                    - b.dispatchPeriod[disp_per].storageDischarged[bat]
                    <= m.storageDischargingRampDownRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            else:
                return pyo.Constraint.Skip

        # Force no charge when discharging
        @disj.Constraint(b.dispatchPeriods)
        def no_charge(disj, disp_per):
            return b.dispatchPeriod[disp_per].storageCharged[bat] <= 0

        # Batteries that are charging both gain and lose energy
        @disj.Constraint(b.dispatchPeriods)
        def discharging_battery_storage_balance(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageChargeLevel[bat]
                    == m.storageRetentionRate[bat]
                    * b.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]
                    - b.dispatchPeriod[disp_per].storageDischarged[bat]
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageChargeLevel[bat]
                    == m.storageRetentionRate[bat]
                    * r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods[-1]]
                    .storageChargeLevel[bat]
                    - m.storageChargingEfficiency[bat]
                    * b.dispatchPeriod[disp_per].storageDischarged[bat]
                )
            else:
                return pyo.Constraint.Skip

    """
    Battery Charging Constraints
    """

    @b.Disjunct(m.storage)
    def storCharging(disj, bat):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods)
        def charge_limit_min(d, disp_per):
            return (
                m.chargeMin[bat]  # Assuming chargeMin is an absolute value (MW)
                <= b.dispatchPeriod[disp_per].storageCharged[bat]
            )

        # Maximum operating limits
        @disj.Constraint(b.dispatchPeriods)
        def charge_limit_max(d, disp_per):
            return b.dispatchPeriod[disp_per].storageCharged[bat] <= m.chargeMax[bat]

        # Ramp up limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def charge_ramp_up_limits(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageCharged[bat]
                    - b.dispatchPeriod[disp_per - 1].storageCharged[bat]
                    <= m.storageChargingRampUpRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageCharged[bat]
                    - r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .storageCharged[bat]
                    <= m.storageChargingRampUpRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            else:
                return pyo.Constraint.Skip

        # Ramp down limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def charge_ramp_down_limits(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per - 1].storageCharged[bat]
                    - b.dispatchPeriod[disp_per].storageCharged[bat]
                    <= m.storageChargingRampDownRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .storageCharged[bat]
                    - b.dispatchPeriod[disp_per].storageCharged[bat]
                    <= m.storageChargingRampDownRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            else:
                return pyo.Constraint.Skip

        @disj.Constraint(b.dispatchPeriods)
        def no_discharge(disj, disp_per):
            return b.dispatchPeriod[disp_per].storageDischarged[bat] <= 0

        # Batteries that are charging both gain and lose energy
        @disj.Constraint(b.dispatchPeriods)
        def charging_battery_storage_balance(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageChargeLevel[bat]
                    == m.storageRetentionRate[bat]
                    * b.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]
                    + m.storageChargingEfficiency[bat]
                    * b.dispatchPeriod[disp_per].storageCharged[bat]
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageChargeLevel[bat]
                    == m.storageRetentionRate[bat]
                    * r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods[-1]]
                    .storageChargeLevel[bat]
                    + m.storageChargingEfficiency[bat]
                    * b.dispatchPeriod[disp_per].storageCharged[bat]
                )
            else:
                return pyo.Constraint.Skip

    """
    Battery Off Constraints
    """

    @b.Disjunct(m.storage)
    def storOff(disj, bat):
        b = disj.parent_block()

        # If battery is off, it is not discharging in terms of sending energy
        # to the grid
        @disj.Constraint(b.dispatchPeriods)
        def no_discharge(disj, disp_per):
            return b.dispatchPeriod[disp_per].storageDischarged[bat] == 0

        # Batteries that are off cannot charge
        @disj.Constraint(b.dispatchPeriods)
        def no_charge(disj, disp_per):
            return b.dispatchPeriod[disp_per].storageCharged[bat] == 0

        # Batteries that are off still lose energy, and none goes to the grid
        @disj.Constraint(b.dispatchPeriods)
        def off_batteries_lose_storage(disj, disp_per):
            if disp_per != 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageChargeLevel[bat]
                    == m.storageRetentionRate[bat]
                    * b.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]
                )
            elif disp_per == 1 and commitment_period != 1:
                return (
                    b.dispatchPeriod[disp_per].storageChargeLevel[bat]
                    == m.storageRetentionRate[bat]
                    * r_p.commitmentPeriod[commitment_period - 1]
                    .dispatchPeriod[b.dispatchPeriods.last()]
                    .storageChargeLevel[bat]
                )
            else:
                return pyo.Constraint.Skip

    # Batteries are exclusively either Charging, Discharging, or Off
    @b.Disjunction(m.storage)
    def storStatus(disj, bat):
        return [
            disj.storCharging[bat],
            disj.storDischarging[bat],
            disj.storOff[bat],
        ]

    # bats cannot be committed unless they are operational or just installed
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    @b.LogicalConstraint(m.storage)
    def commit_active_batts_only(b, bat):
        return pyo.lor(
            b.storCharging[bat].indicator_var, b.storDischarging[bat].indicator_var
        ).implies(
            pyo.lor(
                i_p.storOperational[bat].indicator_var,
                i_p.storInstalled[bat].indicator_var,
                i_p.storExtended[bat].indicator_var,
            )
        )
