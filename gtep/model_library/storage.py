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
from pyomo.core.base.block import BlockData


def add_storage_params(m: pyo.Model):
    """
    This method defines all the battery storage properties from data.

    :param m:       model object
    :type m:        pyomo.environ.Model
    """
    m.storageCapacity = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["energy_capacity"]
            for bat in m.storage
        },
        units=u.MW * u.hr,
        doc="Maximum storage capacity (in MWh)",
    )

    m.initStorageChargeLevel = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["initial_state_of_charge"]
            for bat in m.storage
        },
        units=u.MW * u.hr,
        doc="Initial state of charge (in MWh)",
    )

    m.minStorageChargeLevel = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["minimum_state_of_charge"]
            for bat in m.storage
        },
        units=u.MW * u.hr,
        doc="Minimum state of charge (in MWh)",
    )

    m.chargingCost = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["charge_cost"]
            for bat in m.storage
        },
        units=u.USD / u.MW / u.hr,
        doc="Cost to charge (in USD/MWh)",
    )

    m.dischargingCost = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["discharge_cost"]
            for bat in m.storage
        },
        units=u.USD / u.MW / u.hr,
        doc="Cost to discharge (in USD/MWh)",
    )

    m.dischargeMin = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["min_discharge_rate"]
            for bat in m.storage
        },
        units=u.MW,
        doc="Minimum discharge rate (in MW)",
    )

    m.dischargeMax = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["max_discharge_rate"]
            for bat in m.storage
        },
        units=u.MW,
        doc="Maximum discharge rate (in MW)",
    )

    m.chargeMin = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["min_charge_rate"]
            for bat in m.storage
        },
        units=u.MW,
        doc="Minimum charging rate (in MW)",
    )

    m.chargeMax = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["max_charge_rate"]
            for bat in m.storage
        },
        units=u.MW,
        doc="Maximum charging rate (in MW)",
    )

    # NOTE that default EGRET naming convention assumes
    # dispatch periods are 60 minutes.
    m.storageDischargingRampUpRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_up_output_60min"]
            for bat in m.storage
        },
        units=u.MW / u.hr,  # units?
        doc="Maximum amount of ramp up between dispatch periods when discharging",
    )

    m.storageDischargingRampDownRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_down_output_60min"]
            for bat in m.storage
        },
        units=u.MW / u.hr,  # units?
        doc="Maximum amount of ramp down between dispatch periods when discharging",
    )

    m.storageChargingRampUpRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_up_input_60min"]
            for bat in m.storage
        },
        units=u.MW / u.hr,  # units?
        doc="Maximum amount of ramp up between dispatch periods when charging",
    )

    m.storageChargingRampDownRates = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["ramp_down_input_60min"]
            for bat in m.storage
        },
        units=u.MW / u.hr,  # units?
        doc="Maximum amount of ramp down between dispatch periods when charging",
    )

    m.storageDischargingEfficiency = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["discharge_efficiency"]
            for bat in m.storage
        },
        units=u.dimensionless,
        doc="Proportion of energy discharged not lost to technological inefficencies within dispatch periods and which is usable in the flow balance",
    )

    m.storageChargingEfficiency = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["charge_efficiency"]
            for bat in m.storage
        },
        units=u.dimensionless,
        doc="Proportion of energy charged not lost to technological inefficiencies within dispatch periods and which is usable in the flow balance",
    )

    m.storageRetentionRate = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["retention_rate_60min"]
            for bat in m.storage
        },
        units=1 / u.hr,
        doc="Proportion of stored energy that is preserved per hour (in 1/hr)",
    )

    m.storageCapitalMultiplier = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["capital_multiplier"]
            for bat in m.storage
        },
        units=u.dimensionless,
        doc="(Arbitrary) multiplier for new battery investments, corresponding to depreciation schedules for individual technologies; higher values indicate slower depreciation",
    )

    m.storageExtensionMultiplier = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["extension_multiplier"]
            for bat in m.storage
        },
        units=u.dimensionless,
        doc="Cost of life extension, as a fraction of initial investment cost",
    )

    m.storageInvestmentCost = pyo.Param(
        m.storage,
        initialize={
            bat: m.md.data["elements"]["storage"][bat]["investment_cost"]
            for bat in m.storage
        },
        units=u.USD / u.MW / u.hr,  # TODO: check that this is indeed what's in the data
        doc="Future not real cost; based on idealized targets, in $/MWh",
    )

    m.storageFixedCost = pyo.Param(
        m.storage,
        initialize={stor: 0 for stor in m.storage},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Storage fixed operating costs, in USD/MWh",
    )

    m.storageVarCost = pyo.Param(
        m.storage,
        initialize={stor: 0 for stor in m.storage},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Storage variable costs, in USD/MWh",
    )


### HELPER FUNCTIONS
def _charge_discharge_limit(b, bat, charging: bool, upper: bool):
    """
    Returns a constraint on battery charging/discharging.

    :param b:           Commitment block
    :param bat:         Battery
    :param charging:    Whether to apply constraint on charge (rather than discharge)
    :param upper:       Whether to apply the upper limit (rather than lower)
    """
    m = b.model()
    descrip = "charge" if charging else "discharge"
    max_or_min = "max" if upper else "min"
    limit = m.component(f"{descrip}{max_or_min.capitalize()}")[bat]

    return pyo.Constraint(
        b.dispatchPeriods,
        expr={
            d: (
                limit * (-1 if upper else 1)
                <= b.dispatchPeriod[d].component(f"storage{descrip.capitalize()}d")[bat]
                * (-1 if upper else 1)
            )
            for d in b.dispatchPeriods
        },
        doc=f"Storage {descrip} {max_or_min} operating limits",
    )


def _ramp_limit_rule(b, bat, disp_per, charging, up):
    m = b.model()
    r_p = b.parent_block()
    comm_per = b.index()

    if (comm_per, disp_per) in r_p.commitDispatchPairsNotFirst:
        cur = (
            b.dispatchPeriod[disp_per].storageCharged[bat]
            if charging
            else b.dispatchPeriod[disp_per].storageDischarged[bat]
        )
        prev = (
            r_p.prevStorageCharged[(comm_per, disp_per), bat]
            if charging
            else r_p.prevStorageDischarged[(comm_per, disp_per), bat]
        )
        limit = m.component(
            f"storage{'C' if charging else 'Disc'}hargingRamp{'Up' if up else 'Down'}Rates"
        )[bat]
        return (cur - prev if up else prev - cur) / u.convert(
            m.dispatchPeriodLength, u.hr
        ) <= limit
    return pyo.Constraint.Skip


def _ramp_limit(b, bat, charging: bool, up: bool):
    """
    Returns a constraint on battery charge/discharge ramping.

    :param b:           Commitment block
    :param bat:         Battery
    :param charging:    Whether to apply constraint on charge (rather than discharge)
    :param up:          Whether to apply ramp up limit (rather than ramp down)
    """
    descrip = "charge" if charging else "discharge"
    up_or_down = "up" if up else "down"

    return pyo.Constraint(
        b.dispatchPeriods,
        expr={
            disp_per: _ramp_limit_rule(b, bat, disp_per, charging, up)
            for disp_per in b.dispatchPeriods
        },
        doc=f"Storage {descrip} ramp {up_or_down} limit",
    )


def _no_charge_discharge(b, bat, charge: bool):
    """
    Returns a constraint to prevent charging/discharging.

    :param b:           Commitment block
    :param bat:         Battery
    :param charging:    Whether to apply constraint on charge (rather than discharge)
    """
    return pyo.Constraint(
        b.dispatchPeriods,
        expr={
            d: (
                b.dispatchPeriod[d].storageCharged[bat]
                if charge
                else b.dispatchPeriod[d].storageDischarged[bat]
            )
            <= 0 * u.MW
            for d in b.dispatchPeriods
        },
        doc=f"Forces no {'charging' if charge else 'discharging'}",
    )


def add_storage_state_disjuncts(b: BlockData):
    """
    This function adds battery storage charging and discharging
    disjuncts and constraints.

    :param b:           Commitment block to add disjuncts to
    :type b:            pyomo.core.base.block.BlockData
    """

    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()
    comm_per = b.index()

    ####################
    # Common constraints
    ####################

    @b.Constraint(m.storage, b.dispatchPeriods)
    def battery_storage_balance(b, bat, disp_per):
        if (comm_per, disp_per) in r_p.commitDispatchPairsNotFirst:
            return b.dispatchPeriod[disp_per].storageChargeLevel[
                bat
            ] == r_p.retainedStorageChargeLevelFromPrev[(comm_per, disp_per), bat] + (
                b.dispatchPeriod[disp_per].storageCharged[bat]
                - b.dispatchPeriod[disp_per].storageDischarged[bat]
            ) * u.convert(
                m.dispatchPeriodLength, u.hr
            )
        return pyo.Constraint.Skip

    #########################
    # Discharging constraints
    #########################

    @b.Disjunct(m.storage, doc="Storage discharging operating limits")
    def storDischarging(disj, bat):
        disj.add_component(
            "discharge_limit_min",
            _charge_discharge_limit(b, bat, charging=False, upper=False),
        )
        disj.add_component(
            "discharge_limit_max",
            _charge_discharge_limit(b, bat, charging=False, upper=True),
        )
        disj.add_component(
            "discharge_ramp_down_limits",
            _ramp_limit(b, bat, charging=False, up=False),
        )
        disj.add_component(
            "discharge_ramp_up_limits",
            _ramp_limit(b, bat, charging=False, up=True),
        )
        disj.add_component("no_charge", _no_charge_discharge(b, bat, charge=True))

    #########################
    # Charging constraints
    #########################

    @b.Disjunct(m.storage)
    def storCharging(disj, bat):
        disj.add_component(
            "charge_limit_min",
            _charge_discharge_limit(b, bat, charging=True, upper=False),
        )
        disj.add_component(
            "charge_limit_max",
            _charge_discharge_limit(b, bat, charging=True, upper=True),
        )
        disj.add_component(
            "charge_ramp_down_limits",
            _ramp_limit(b, bat, charging=True, up=False),
        )
        disj.add_component(
            "charge_ramp_up_limits",
            _ramp_limit(b, bat, charging=True, up=True),
        )
        disj.add_component("no_discharge", _no_charge_discharge(b, bat, charge=False))

    #########################
    # Storage Off Constraints
    #########################

    @b.Disjunct(m.storage)
    def storOff(disj, bat):
        disj.add_component("no_charge", _no_charge_discharge(b, bat, charge=True))
        disj.add_component("no_discharge", _no_charge_discharge(b, bat, charge=False))

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
        m = b.model()

        return sum(
            m.storageInvestmentCost[bat]
            * m.storageCapacity[bat]  # TODO: verify this is correct; was chargeMax
            * m.storageCapitalMultiplier[bat]
            * b.storInstalled[bat].indicator_var.get_associated_binary()
            for bat in m.storage
        ) + sum(
            m.storageInvestmentCost[bat]
            * m.storageCapacity[bat]
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
        units=u.MW * u.hr,
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

    # Operational cost variables and expressions per storage
    # unit. Here we assume the costs are in $/MW. If instead the costs
    # are in $/MWh, the storageCharge/Discharged should be multiplied
    # by b.dispatchPeriodLength, in hours.
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
