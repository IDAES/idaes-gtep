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

import pandas as pd
import pyomo.environ as pyo
from pyomo.environ import units as u


def add_storage_params(m):
    """This method defines all the battery storage properties from
    data

    """

    # Maximum storage capacity
    m.storageCapacity = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["energy_capacity"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW * u.hr,
        doc="Maximum storage capacity in MWh",
    )
    m.initStorageChargeLevel = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat][
                "initial_state_of_charge"
            ]  # 80% energy capacity
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW * u.hr,
        doc="Initial storage capacity",
    )

    m.minStorageChargeLevel = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat][
                "minimum_state_of_charge"
            ]  # 20% energy capacity
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW * u.hr,
        doc="Minimum storage capacity",
    )

    m.dischargeMin = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["min_discharge_rate"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Minimum amount to discharge per dispatch period when discharging",
    )

    m.dischargeMax = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["max_discharge_rate"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum amount to discharge per dispatch period when discharging",
    )

    m.chargeMin = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["min_charge_rate"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Minimum amount to charge per dispatch period when charging",
    )

    m.chargeMax = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["max_charge_rate"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum amount to charge per dispatch period when charging",
    )

    # NOTE: Default EGRET naming convention assumes dispatch periods
    # are 60 minutes.
    m.storageDischargingRampUpRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_up_output_60min"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum amount of ramp up between dispatch periods when discharging",
    )

    m.storageDischargingRampDownRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_down_output_60min"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum amount of ramp down between dispatch periods when discharging.",
    )

    m.storageChargingRampUpRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_up_input_60min"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum amount of ramp up between dispatch periods when charging",
    )

    m.storageChargingRampDownRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_down_input_60min"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum amount of ramp down between dispatch periods when charging",
    )

    m.storageDischargingEfficiency = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["discharge_efficiency"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.dimensionless,
        doc="Proportion of energy discharged that is not lost to technological inefficiencies with in dispatch periods and which is usable in the flow balance",
    )

    m.storageChargingEfficiency = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["charge_efficiency"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.dimensionless,
        doc="Proportion of energy charged that is not lost to technological inefficiencies within dispatch periods and which is usable in the flow balance",
    )

    m.storageRetentionRate = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["retention_rate_60min"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.dimensionless,
        doc="Proportion of energy discharged that is not lost to technological inefficiencies between dispatch periods and which is usable in the flow balance",
    )

    # TODO: Calculate this instead
    m.storageCapitalMultiplier = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["capital_multiplier"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.dimensionless,
        doc="(Arbitrary) multiplier for new battery investments corresponds to depreciation schedules for individual technologies; higher values are indicative of slow depreciation",
    )

    m.storageExtensionMultiplier = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["extension_multiplier"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.dimensionless,
        doc="Cost of life extension for each battery, expressed as a fraction of initial investment cost",
    )

    m.storageLifetimes = pyo.Param(
        m.storage,
        initialize={stor: 3 for stor in m.storage},
        mutable=True,
        units=u.year,
        doc="Lifetime of each storage unit",
    )

    # Add storage charge/discharge caps
    m.storageChargeLimit = pyo.Param(
        initialize=sum(pyo.value(m.chargeMax[bat]) for bat in m.storage),
        mutable=True,
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum storage charging per dispatch period",
    )

    m.storageDischargeLimit = pyo.Param(
        initialize=sum(pyo.value(m.dischargeMax[bat]) for bat in m.storage),
        mutable=True,
        domain=pyo.NonNegativeReals,
        units=u.MW,
        doc="Maximum storage discharging per dispatch period",
    )

    # Initialize charge/discharge costs and fixed and variable
    # costs. The fixed and var costs are updated during the investment
    # stage based on costs given in m.mc data modeling object.
    m.chargingCost = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["charge_cost"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.USD / (u.MW * u.hr),
        doc="Cost to charge per unit electricity",
    )

    m.dischargingCost = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["discharge_cost"]
            for bat in m.storage
        },
        domain=pyo.NonNegativeReals,
        units=u.USD / (u.MW * u.hr),
        doc="Cost to discharge per unit electricity",
    )
    m.storageFixedCost = pyo.Param(
        m.storage,
        initialize={stor: 0 for stor in m.storage},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Fixed operating costs",
    )
    m.storageVariableCost = pyo.Param(
        m.storage,
        initialize={stor: 0 for stor in m.storage},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Variable costs",
    )

    m.storageInvestmentCost = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["investment_cost"]
            for bat in m.storage
        },
        mutable=True,
        domain=pyo.NonNegativeReals,
        units=u.USD / u.MW,
        doc="Investment cost for storage units",
    )


def add_storage_cost_parameters_from_csv(m, year):
    """This method updates investment cost parameters for storage
    units using data from CSV files loaded into the model object
    (m.mc).

    For the specified year, this function updates:
    - Storage lifetime parameters.
    - Investment cost parameters using annualized capex data, converting from $/MW-yr
      to $/MW using the lifetime and a discount rate.
    - Variable and fixed operating costs.

    The final units (to avoid unit consistency issues) should be:
    - fixed cost = $/MWh
    - var cost = $/MWh
    - inv cost = $/Mw

    """

    # Re-populating lifetimes parameters for branches and generators
    # since we have data in the m.mc model object.
    lifetime_col = f"lifetime_{year}"
    new_storage_lifetimes = {
        row["name"]: int(row.get(lifetime_col, 0) if pd.notna(row[lifetime_col]) else 3)
        for _, row in m.mc.storage_data_target.iterrows()
    }

    for storage in m.storage:
        if storage in new_storage_lifetimes:
            m.storageLifetimes[storage] = new_storage_lifetimes[storage]

    # Re-populate the investment cost parameters for storage units
    # since we have available capex data in m.mc modeling
    # object. NOTE: Since the data is annualized ($/MWyr), we
    # de-annualize it using the lifetime parameter and an assumed
    # discounte rate. The final units are in $/MW.
    def annualized_to_total_capex(annualized_cost, years, discount_rate):
        r = discount_rate
        n = years
        crf = (r * (1 + r) ** n) / ((1 + r) ** n - 1)
        total_cost = annualized_cost / crf
        return total_cost

    if m.mc is not None:
        original_units = u.USD / (u.MW * u.year)
        final_units = u.USD / (u.MW * u.hr)
        final_inv_units = u.USD / u.MW

        for index, row in m.mc.storage_data_target.iterrows():
            storage_uid = row["name"]

            # Read costs for the selected year
            capex_yr = float(row[f"capex_{year}"])  # units in $/MW-year
            fixed_ops_yr = float(row[f"fixed_ops_{year}"])  # units in $/MW-year
            var_ops_yr = float(row[f"var_ops_{year}"])  # units in $/MWh

            inv_cost = capex_yr * (u.USD / u.MW)
            # inv_cost = annualized_to_total_capex(
            #     capex_yr,
            #     years=pyo.value(m.storageLifetimes[storage_uid]),
            #     discount_rate=0.07,
            # )

            fixed_cost = pyo.units.convert(
                fixed_ops_yr * original_units, to_units=final_units
            )
            var_cost = var_ops_yr * final_units  # units in $/MWh

            # Assign to Pyomo parameters (strip units for Pyomo Param)
            m.storageFixedCost[storage_uid] = pyo.value(fixed_cost)
            m.storageVariableCost[storage_uid] = pyo.value(var_cost)
            m.storageInvestmentCost[storage_uid] = pyo.value(inv_cost)

    # m.storageInvestmentCost.display()
    # m.storageFixedCost.display()
    # m.storageVariableCost.display()


def add_storage_state_disjuncts(m, b, commitment_period):
    """This method includes the battery storage charging and
    discharging constraints

    """

    r_p = b.parent_block()
    i_p = r_p.parent_block()

    #########################
    # Discharging constraints
    #########################

    @b.Disjunct(m.storage, doc="Storage discharging operating limits")
    def storDischarging(disj, bat):
        b = disj.parent_block()

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage discharging minimum operating Limits if storage unit is on",
        )
        def discharge_limit_min(d, disp_per):
            return (
                m.dischargeMin[bat]  # Assuming dischargeMin is an absolute value (MW)
                <= b.dispatchPeriod[disp_per].storageDischarged[bat]
            )

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage discharging maximum operating limits",
        )
        def discharge_limit_max(d, disp_per):
            return (
                b.dispatchPeriod[disp_per].storageDischarged[bat] <= m.dischargeMax[bat]
            )

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage discharging ramp up limit when fully on",
        )
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

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage discharging ramp down limit when fully on",
        )
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
                    .dispatchPeriod[b.dispatchPeriods.at(-1)]
                    .storageDischarged[bat]
                    - b.dispatchPeriod[disp_per].storageDischarged[bat]
                    <= m.storageDischargingRampDownRates[
                        bat
                    ]  # battery ramp rates are currently absolute values
                )
            else:
                return pyo.Constraint.Skip

        @disj.Constraint(b.dispatchPeriods, doc="Forces no charge when discharging")
        def no_charge(disj, disp_per):
            return b.dispatchPeriod[disp_per].storageCharged[bat] <= 0

        # @disj.Constraint(
        #     b.dispatchPeriods,
        #     doc="Enforces that batteries that are charging both gain and lose energy",
        # )
        # def discharging_battery_storage_balance(disj, disp_per):
        #     if disp_per != 1 and commitment_period != 1:
        #         return (
        #             b.dispatchPeriod[disp_per].storageChargeLevel[bat]
        #             == m.storageRetentionRate[bat]
        #             * b.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]
        #             - b.dispatchPeriod[disp_per].storageDischarged[bat]
        #         )
        #     elif disp_per == 1 and commitment_period != 1:
        #         return (
        #             b.dispatchPeriod[disp_per].storageChargeLevel[bat]
        #             == m.storageRetentionRate[bat]
        #             * r_p.commitmentPeriod[commitment_period - 1]
        #             .dispatchPeriod[b.dispatchPeriods.at(-1)]
        #             .storageChargeLevel[bat]
        #             - m.storageDischargingEfficiency[bat]
        #             * b.dispatchPeriod[disp_per].storageDischarged[bat]
        #         )
        #     else:
        #         return pyo.Constraint.Skip

    #########################
    # Charging constraints
    #########################

    @b.Disjunct(m.storage)
    def storCharging(disj, bat):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods)
        def charge_limit_min(d, disp_per):
            return (
                m.chargeMin[bat]  # Assuming chargeMin is an absolute value (MW)
                <= b.dispatchPeriod[disp_per].storageCharged[bat]
            )

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage charging maximum operating limits",
        )
        def charge_limit_max(d, disp_per):
            return b.dispatchPeriod[disp_per].storageCharged[bat] <= m.chargeMax[bat]

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage charging ramp up limit when fully on",
        )
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

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Storage charging ramp down limit when fully on",
        )
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

        # @disj.Constraint(
        #     b.dispatchPeriods,
        #     doc="Enforces that batteries that are charging both gain and lose energy",
        # )
        # def charging_battery_storage_balance(disj, disp_per):

        #     if disp_per != 1 and commitment_period != 1:
        #         return (
        #             b.dispatchPeriod[disp_per].storageChargeLevel[bat]
        #             == m.storageRetentionRate[bat]
        #             * b.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]
        #             + m.storageChargingEfficiency[bat]
        #             * b.dispatchPeriod[disp_per].storageCharged[bat]
        #         )
        #     elif disp_per == 1 and commitment_period != 1:
        #         return (
        #             b.dispatchPeriod[disp_per].storageChargeLevel[bat]
        #             == m.storageRetentionRate[bat]
        #             * r_p.commitmentPeriod[commitment_period - 1]
        #             .dispatchPeriod[b.dispatchPeriods.at(-1)]
        #             .storageChargeLevel[bat]
        #             + m.storageChargingEfficiency[bat]
        #             * b.dispatchPeriod[disp_per].storageCharged[bat]
        #         )
        #     else:
        #         return pyo.Constraint.Skip

    #########################
    # Storage Off Constraints
    #########################

    @b.Disjunct(m.storage)
    def storOff(disj, bat):
        b = disj.parent_block()

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Enforces that, if battery is off, it is not discharging in terms of sending energy to the grid",
        )
        def no_discharge(disj, disp_per):
            return b.dispatchPeriod[disp_per].storageDischarged[bat] == 0

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Enforces that batteries that are off cannot charge their status",
        )
        def no_charge(
            disj,
            disp_per,
            doc="Enforces that batteries that are off still lose energy, and none goes to the grid",
        ):
            return b.dispatchPeriod[disp_per].storageCharged[bat] == 0

        # @disj.Constraint(b.dispatchPeriods)
        # def off_batteries_lose_storage(disj, disp_per):
        #     if disp_per != 1 and commitment_period != 1:
        #         return (
        #             b.dispatchPeriod[disp_per].storageChargeLevel[bat]
        #             == m.storageRetentionRate[bat]
        #             * b.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]
        #         )
        #     elif disp_per == 1 and commitment_period != 1:
        #         return (
        #             b.dispatchPeriod[disp_per].storageChargeLevel[bat]
        #             == m.storageRetentionRate[bat]
        #             * r_p.commitmentPeriod[commitment_period - 1]
        #             .dispatchPeriod[b.dispatchPeriods.last()]
        #             .storageChargeLevel[bat]
        #         )
        #     else:
        #         return pyo.Constraint.Skip

    @b.Disjunction(
        m.storage,
        doc="Enforces that batteries are exclusively either Charging, Discharging, or Off",
    )
    def storStatus(disj, bat):
        return [
            disj.storCharging[bat],
            disj.storDischarging[bat],
            disj.storOff[bat],
        ]

    r_p = b.parent_block()
    i_p = r_p.parent_block()

    @b.LogicalConstraint(
        m.storage,
        doc="Enforces that batteries cannot be committed unless they are operational or just installed",
    )
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


def add_investment_storage_constraints(m, b, investment_stage):

    # Fix "in-service" initial investment state for storage units
    # based on input. [TODO: Initialize storage level (state of
    # charge)]
    for stor in m.storage:
        in_service = m.md.data["elements"]["storage"][stor]["in_service"]

        if investment_stage == 1:
            if in_service:
                b.storOperational[stor].indicator_var.fix(True)
            else:
                b.storOperational[stor].indicator_var.fix(False)

    if not m.config["include_investment"]:
        for stor in m.storage:
            is_candidate = str(stor).endswith("-c")

            if is_candidate:
                b.storOperational[stor].indicator_var.fix(False)
                b.storDisabled[stor].indicator_var.fix(True)

    @b.Expression(doc="Storage investment costs in $")
    def storage_investment_cost(b):
        # Using chargeMax for unit consistency, assuming it is
        # equivalent to storageCapacity but in MW
        return sum(
            m.storageInvestmentCost[bat]  # in USD/MW
            # * m.storageCapacity[bat]  # in MWh
            * m.chargeMax[bat]  # in MW
            * m.storageCapitalMultiplier[bat]
            * b.storInstalled[bat].indicator_var.get_associated_binary()
            for bat in m.storage
        ) + sum(
            m.storageInvestmentCost[bat]  # in USD/MW
            # * m.storageCapacity[bat]  # in MWh
            * m.chargeMax[bat]  # in MW
            * m.storageExtensionMultiplier[bat]
            * b.storExtended[bat].indicator_var.get_associated_binary()
            for bat in m.storage
        )

    # Add legacy equation below. This is not used in current version
    # of the model. [TODO: Check if we want to add this constraint in
    # future versions of the model.]
    """ 
    # Initial, untested attempt for enforcing identical storage level
    at beginning and end of representative periods. [TODO: Update to
    use init and end batteryChargeLevel?]

    @b.Constraint(
    b.representativePeriods,
    m.batteryStorageSystems,
    doc="Enforces that we have identical storage level at the beginning and end of representative period",
    )
    def consistent_battery_charge_level_commitment(b, rep_per, bat):
        return (

                b.representativePeriod[rep_per]
                .commitmentPeriod[
                    b.representativePeriod[rep_per]
                    .commitmentPeriods.first()
                    ]
                    .dispatchPeriod[
                        b.representativePeriod[rep_per]
                        .commitmentPeriod[
                            b.representativePeriod[rep_per]
                            .commitmentPeriods.first()
                            ]
                            .dispatchPeriods.first()
                        ]
                        .batteryChargeLevel[bat]
                  ==
                  b.representativePeriod[rep_per]
                  .commitmentPeriod[
                      b.representativePeriod[rep_per]
                      .commitmentPeriods.last()
                      ]
                      .dispatchPeriod[
                          b.representativePeriod[rep_per]
                          .commitmentPeriod[
                              b.representativePeriod[rep_per]
                              .commitmentPeriods.last()
                              ]
                              .dispatchPeriods.last()
                          ]
                          .batteryChargeLevel[bat]
        )
    """


def add_storage_status_disjuncts(b, storage_set):
    """This method implements a Disjunction and its disjuncts to model
    the selection of the storage units status. The possible
    alternatives for each storage unit are represented as a disjunct
    expression within the function. The options are:

    storOperational: Storage is active and transmitting power.
    storInstalled:   Storage is newly added and active.
    storRetired:     Storage is removed from service.
    storDisabled:    Storage is temporarily out of service.
    storExtended:    Storage is upgraded beyond its original capacity.

    """

    @b.Disjunct(storage_set)
    def storOperational(disj, bat):
        return

    @b.Disjunct(storage_set)
    def storInstalled(disj, bat):
        return

    @b.Disjunct(storage_set)
    def storRetired(disj, bat):
        return

    @b.Disjunct(storage_set)
    def storDisabled(disj, bat):
        return

    @b.Disjunct(storage_set)
    def storExtended(disj, bat):
        return

    @b.Disjunction(storage_set)
    def storInvestStatus(disj, bat):
        return [
            disj.storOperational[bat],
            disj.storInstalled[bat],
            disj.storRetired[bat],
            disj.storDisabled[bat],
            disj.storExtended[bat],
        ]


def add_storage_logical_constraints(m):
    """This method defines logical constraints to ensure that storage
    statuses transitions are operationally consistent over time,
    across the investment stages.

    """

    for stor in m.storage:
        if m.md.data["elements"]["storage"][stor]["in_service"] == False:
            m.investmentStage[1].storOperational[stor].indicator_var.fix(False)
            m.investmentStage[1].storExtended[stor].indicator_var.fix(False)
        else:
            m.investmentStage[1].storOperational[stor].indicator_var.fix(True)

    @m.LogicalConstraint(
        m.stages,
        m.storage,
        doc="Enforces that, if a storage device is online at time t, it must have been online or installed at time t-1",
    )
    def consistent_operation_stor(m, stage, gen):
        return (
            m.investmentStage[stage]
            .storOperational[gen]
            .indicator_var.implies(
                m.investmentStage[stage - 1].storOperational[gen].indicator_var
                | m.investmentStage[stage - 1].storInstalled[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.storage,
        doc="Enforces that, if a storage unit is online at time t, it must be online, extended, or retired at time t+1",
    )
    def consistent_operation_future_stor(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .storOperational[gen]
            .indicator_var.implies(
                m.investmentStage[stage].storOperational[gen].indicator_var
                | m.investmentStage[stage].storExtended[gen].indicator_var
                | m.investmentStage[stage].storRetired[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.storage,
        doc="Enforces that, if a storage unit is retired in period t-1 it must be disabled in period t",
    )
    def full_retirement_stor(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .storRetired[gen]
            .indicator_var.implies(
                m.investmentStage[stage].storDisabled[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    # [TODO: Disabling is permanent.  Re investment is a "new" unit.
    # Remove the "or".]
    @m.LogicalConstraint(
        m.stages,
        m.storage,
        doc="Enforces that, if a gen is disabled at time t-1, it must stay disabled  at time t",
    )
    def consistent_disabled_stor(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .storDisabled[gen]
            .indicator_var.implies(
                m.investmentStage[stage].storDisabled[gen].indicator_var
                | m.investmentStage[stage].storInstalled[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.storage,
        doc="Enforces that, if a gen is extended at time t-1, it must stay extended or be retired at time t",
    )
    def consistent_extended_stor(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .storExtended[gen]
            .indicator_var.implies(
                m.investmentStage[stage].storExtended[gen].indicator_var
                | m.investmentStage[stage].storRetired[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.storage,
        doc="Enforces that, if a storage unit is installed in period t-1, it must be operational in period t",
    )
    def full_investment_stor(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .storInstalled[gen]
            .indicator_var.implies(
                m.investmentStage[stage].storOperational[gen].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )


def add_dispatch_storage_variables_and_constraints(m, b):

    c_p = b.parent_block()
    r_p = c_p.parent_block()

    # NOTE: The lower bound should be > 0 (from data input)
    def storage_capacity_limits(b, bat):
        return (
            m.minStorageChargeLevel[bat],
            m.storageCapacity[bat],
        )

    # Define bounds on charging/discharging capability. Note that
    # constraints enforce that there are min & max charge/discharge
    # levels if the bat is in the charging or discharging state
    def storage_charge_limits(b, bat):
        return (0, m.chargeMax[bat])

    def storage_discharge_limits(b, bat):
        return (0, m.dischargeMax[bat])

    b.storageCharged = pyo.Var(
        m.storage,
        domain=pyo.NonNegativeReals,
        bounds=storage_charge_limits,
        initialize=0,
        units=u.MW,
    )

    b.storageDischarged = pyo.Var(
        m.storage,
        domain=pyo.NonNegativeReals,
        bounds=storage_discharge_limits,
        initialize=0,
        units=u.MW,
    )

    # Operational cost variables and expressions per storage unit
    @b.Expression(m.storage, doc="Charging cost per battery")
    def storageChargingCost(b, bat):
        return (
            b.storageCharged[bat]
            * pyo.units.convert(b.dispatchPeriodLength, to_units=u.hr)  # in MWh
            * m.chargingCost[bat]  # in $/MWh
        )

    @b.Expression(m.storage, doc="Discharging cost per battery")
    def storageDischargingCost(b, bat):
        return (
            b.storageDischarged[bat]
            * pyo.units.convert(b.dispatchPeriodLength, to_units=u.hr)  # in MWh
            * m.dischargingCost[bat]  # in $/MWh
        )

    @b.Expression()
    def storageCostDispatch(b):
        return sum(b.storageChargingCost[bat] for bat in m.storage) + sum(
            b.storageDischargingCost[bat] for bat in m.storage
        )  # in $

    # Declare initial values for battery capacity at the first
    # dispatch period.
    def init_storage_capacity(b, bat):
        return m.initStorageChargeLevel[bat]  # in MWh

    b.storageChargeLevel = pyo.Var(
        m.storage,
        domain=pyo.NonNegativeReals,
        bounds=storage_capacity_limits,
        initialize=init_storage_capacity,
        units=u.MW * u.hr,
    )

    @b.Constraint(m.storage, doc="Storage state-of-charge balance")
    def storage_state_of_charge_balance(b, bat):
        disp_per = b.dispatchPeriod
        commitment_period = c_p.commitmentPeriod

        if disp_per != 1 and commitment_period != 1:
            # If this is not the first dispatch period, use the previous
            # dispatch period in the same commitment period.
            previous_soc = c_p.dispatchPeriod[disp_per - 1].storageChargeLevel[bat]

        elif disp_per == 1 and commitment_period != 1:
            # If this is the first dispatch period of this commitment
            # period, but not the first commitment period, use the
            # final dispatch period from the previous commitment
            # period.
            previous_commitment = r_p.commitmentPeriod[commitment_period - 1]
            previous_dispatch = previous_commitment.dispatchPeriods.at(-1)
            previous_soc = previous_commitment.dispatchPeriod[
                previous_dispatch
            ].storageChargeLevel[bat]

        else:
            # If this is the first dispatch period of the first
            # commitment period, use the initial state of charge from
            # data.
            # previous_soc = m.initStorageChargeLevel[bat]
            previous_soc = m.storageCapacity[bat] / 2

        return b.storageChargeLevel[bat] == (
            m.storageRetentionRate[bat] * previous_soc
            + m.storageChargingEfficiency[bat]
            * b.storageCharged[bat]
            * pyo.units.convert(b.dispatchPeriodLength, to_units=u.hr)
            - m.storageDischargingEfficiency[bat]
            * b.storageDischarged[bat]
            * pyo.units.convert(b.dispatchPeriodLength, to_units=u.hr)
        )

    @b.Constraint(doc="Storage cap")
    def total_storage_cap(b):
        return (
            sum(b.storageCharged[bat] for bat in m.storage)  # in MW
            + sum(b.storageDischarged[bat] for bat in m.storage)  # in MW
            <= m.storageDischargeLimit  # in MW
        )


def add_commitment_storage_constraints(b):

    # [TODO: Replace this constraint with expressions using bounds
    # transform and check if costs considered need to be
    # re-assessed and account for missing data.]
    @b.Expression(doc="Total storage costs for commitment block in $")
    def storageCostCommitment(b):
        return sum(
            b.dispatchPeriod[disp_per].storageCostDispatch
            for disp_per in b.dispatchPeriods
        )
