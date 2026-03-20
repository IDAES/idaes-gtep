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
from pyomo.environ import units as u


def add_storage_params(m):
    """This method defines all the battery storage properties from
    data

    """

    # Maximum storage capacity
    m.storageCapacity = {
        bat: m.md.data["elements"]["storage"][bat]["energy_capacity"]
        for bat in m.storage
    }

    # Initial storage capacity
    m.initStorageChargeLevel = {
        bat: m.md.data["elements"]["storage"][bat]["initial_state_of_charge"]
        for bat in m.storage
    }

    # Minimum storage capacity
    m.minStorageChargeLevel = {
        bat: m.md.data["elements"]["storage"][bat]["minimum_state_of_charge"]
        for bat in m.storage
    }

    # Cost to charge per unit electricity
    m.chargingCost = {
        bat: m.md.data["elements"]["storage"][bat]["charge_cost"] for bat in m.storage
    }

    # Cost to discharge per unit electricity
    m.dischargingCost = {
        bat: m.md.data["elements"]["storage"][bat]["discharge_cost"]
        for bat in m.storage
    }

    # Minimum amount to discharge per dispatch period when discharging
    m.dischargeMin = {
        bat: m.md.data["elements"]["storage"][bat]["min_discharge_rate"]
        for bat in m.storage
    }

    # Maximum amount to discharge per dispatch period when discharging
    m.dischargeMax = {
        bat: m.md.data["elements"]["storage"][bat]["max_discharge_rate"]
        for bat in m.storage
    }

    # Minimum amount to charge per dispatch period when charging
    m.chargeMin = {
        bat: m.md.data["elements"]["storage"][bat]["min_charge_rate"]
        for bat in m.storage
    }

    # Maximum amount to charge per dispatch period when charging
    m.chargeMax = {
        bat: m.md.data["elements"]["storage"][bat]["max_charge_rate"]
        for bat in m.storage
    }

    # Maximum amount of ramp up between dispatch periods when
    # discharging. NOTE that default EGRET naming convention assumes
    # dispatch periods are 60 minutes.
    m.storageDischargingRampUpRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_up_output_60min"]
        for bat in m.storage
    }

    # Maximum amount of ramp down between dispatch periods when
    # discharging.
    m.storageDischargingRampDownRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_down_output_60min"]
        for bat in m.storage
    }

    # Maximum amount of ramp up between dispatch periods when
    # charging.
    m.storageChargingRampUpRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_up_input_60min"]
        for bat in m.storage
    }

    # Maximum amount of ramp down between dispatch periods when
    # charging.
    m.storageChargingRampDownRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_down_input_60min"]
        for bat in m.storage
    }

    # Proportion of energy discharged that is not lost to
    # technological inefficiencies with in dispatch periods and which
    # is usable in the flow balance
    m.storageDischargingEfficiency = {
        bat: m.md.data["elements"]["storage"][bat]["discharge_efficiency"]
        for bat in m.storage
    }

    # Proportion of energy charged that is not lost to technological
    # inefficiencies within dispatch periods and which is usable in
    # the flow balance
    m.storageChargingEfficiency = {
        bat: m.md.data["elements"]["storage"][bat]["charge_efficiency"]
        for bat in m.storage
    }

    # Proportion of energy discharged that is not lost to
    # technological inefficiencies between dispatch periods and which
    # is usable in the flow balance
    m.storageRetentionRate = {
        bat: m.md.data["elements"]["storage"][bat]["retention_rate_60min"]
        for bat in m.storage
    }

    # (Arbitrary) multiplier for new battery investments corresponds
    # to depreciation schedules for individual technologies; higher
    # values are indicative of slow depreciation
    m.storageCapitalMultiplier = {
        bat: m.md.data["elements"]["storage"][bat]["capital_multiplier"]
        for bat in m.storage
    }

    # Cost of life extension for each battery, expressed as a fraction
    # of initial investment cost
    m.storageExtensionMultiplier = {
        bat: m.md.data["elements"]["storage"][bat]["extension_multiplier"]
        for bat in m.storage
    }

    # Future not real cost: idealized DoE 10-yr targets or something
    m.storageInvestmentCost = {
        bat: m.md.data["elements"]["storage"][bat]["investment_cost"]
        for bat in m.storage
    }


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
                    .dispatchPeriod[b.dispatchPeriods[-1]]
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

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Enforces that batteries that are charging both gain and lose energy",
        )
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

        @disj.Constraint(
            b.dispatchPeriods,
            doc="Enforces that batteries that are charging both gain and lose energy",
        )
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


# def add_investment_storage_variables(b):
#     b.storageCostInvestment = pyo.Var(
#         within=pyo.NonNegativeReals, initialize=0, units=u.USD
#     )


def add_investment_storage_constraints(m, b, investment_stage):

    # Fix "in-service" batteries initial investment state based on
    # input. [TODO: Initialize storage level (state of charge)]
    for bat in m.storage:
        if (
            m.md.data["elements"]["storage"][bat]["in_service"] == False
            and investment_stage == 1
        ):
            b.storOperational[bat].indicator_var.fix(False)
        elif (
            m.md.data["elements"]["storage"][bat]["in_service"] == True
            and investment_stage == 1
        ):
            b.storOperational[bat].indicator_var.fix(True)

    @b.Expression(doc="Storage investment costs in $")
    def storage_investment_cost(b):
        return sum(
            m.storageInvestmentCost[bat]
            * m.storageCapitalMultiplier[bat]
            * b.storInstalled[bat].indicator_var.get_associated_binary()
            for bat in m.storage
        ) + sum(
            m.storageInvestmentCost[bat]
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

    # NOTE: The lower bound should be > 0 (from data input)
    def storage_capacity_limits(b, bat):
        return (
            m.minStorageChargeLevel[bat],
            m.storageCapacity[bat],
        )

    # [TODO: We need to adjust this constraint since this does not fix
    # initial battery capacity at the first dispatch period.]
    def init_storage_capacity(b, bat):
        return m.initStorageChargeLevel[bat]

    b.storageChargeLevel = pyo.Var(
        m.storage,
        domain=pyo.NonNegativeReals,
        bounds=storage_capacity_limits,
        initialize=init_storage_capacity,
        units=u.MW,
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
        return b.storageCharged[bat] * m.chargingCost[bat]

    @b.Expression(m.storage, doc="Discharging cost per battery")
    def storageDischargingCost(b, bat):
        return b.storageDischarged[bat] * m.dischargingCost[bat]

    @b.Expression()
    def storageCostDispatch(b):
        return sum(b.storageChargingCost[bat] for bat in m.storage) + sum(
            b.storageDischargingCost[bat] for bat in m.storage
        )
