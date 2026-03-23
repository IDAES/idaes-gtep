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

"""Constraints for Transmission in the Generation and Transmission
Expansion Planning (GTEP) Model

"""

import math

import pyomo.environ as pyo
from pyomo.environ import units as u


def add_investment_transmission_constraints(m, b, investment_stage):
    for branch in m.transmission:
        if (
            m.md.data["elements"]["branch"][branch]["in_service"] == False
            and investment_stage == 1
        ):
            b.branchOperational[branch].indicator_var.fix(False)
        elif (
            m.md.data["elements"]["branch"][branch]["in_service"] == True
            and investment_stage == 1
        ):
            b.branchOperational[branch].indicator_var.fix(True)

    @b.Expression(doc="Transmission investment costs in $")
    def transmission_investment_cost(b):
        return sum(
            m.branchInvestmentCost[branch]
            * m.branchCapitalMultiplier[branch]
            * b.branchInstalled[branch].indicator_var.get_associated_binary()
            for branch in m.transmission
        ) + sum(
            m.branchInvestmentCost[branch]
            * m.branchExtensionMultiplier[branch]
            * b.branchExtended[branch].indicator_var.get_associated_binary()
            for branch in m.transmission
        )


def add_transmission_status_disjuncts(b, transmission_set):
    """This method implements a Disjunction and its disjuncts to model
    the selection of the transmission lines status. The possible
    alternatives for each transmission line are represented as a
    disjunct expression within the function. The options are:

    branchOperational: Line is active and transmitting power.
    branchInstalled:   Line is newly added and active.
    branchRetired:     Line is removed from service.
    branchDisabled:    Line is temporarily out of service.
    branchExtended:    Line is upgraded beyond its original capacity.

    """

    # For now, mimicking thermal generator disjuncts. [TODO: Check if
    # we need to define different states for transmission lines or if
    # we differentiate between line and transformer investments?]
    @b.Disjunct(transmission_set)
    def branchOperational(disj, branch):
        return

    @b.Disjunct(transmission_set)
    def branchInstalled(disj, branch):
        return

    @b.Disjunct(transmission_set)
    def branchRetired(disj, branch):
        return

    @b.Disjunct(transmission_set)
    def branchDisabled(disj, branch):
        return

    @b.Disjunct(transmission_set)
    def branchExtended(disj, branch):
        return

    # [TODO: Do we differentiate between line and transformer
    # investments?]
    @b.Disjunction(transmission_set)
    def branchInvestStatus(disj, branch):
        return [
            disj.branchOperational[branch],
            disj.branchInstalled[branch],
            disj.branchRetired[branch],
            disj.branchDisabled[branch],
            disj.branchExtended[branch],
        ]


def add_transmission_logical_constraints(m):
    """This method defines logical constraints to ensure that
    transmission lines statuses transitions are operationally
    consistent over time, across the investment stages.

    """

    @m.LogicalConstraint(
        m.stages,
        m.transmission,
        doc="Enforces that, if a branch is online at time t, it must have been online or installed at time t-1",
    )
    def consistent_branch_operation(m, stage, branch):
        return (
            m.investmentStage[stage]
            .branchOperational[branch]
            .indicator_var.implies(
                m.investmentStage[stage - 1].branchOperational[branch].indicator_var
                | m.investmentStage[stage - 1].branchInstalled[branch].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.transmission,
        doc="Enforces that, if a branch is online at time t, it must be online, extended, or retired at time t+1",
    )
    def consistent_branch_operation_future(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchOperational[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchOperational[branch].indicator_var
                | m.investmentStage[stage].branchExtended[branch].indicator_var
                | m.investmentStage[stage].branchRetired[branch].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.transmission,
        doc="Enforces that, if a branch is retired in period t-1, it should be disabled in period t",
    )
    def full_branch_retirement(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchRetired[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchDisabled[branch].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.transmission,
        doc="Enforces that, if a branch is disabled at time t-1, it must stay disabled or be installed at time t",
    )
    def consistent_branch_disabled(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchDisabled[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchDisabled[branch].indicator_var
                | m.investmentStage[stage].branchInstalled[branch].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.transmission,
        doc="Enforces that, if a branch is extended at time t-1, it must stay extended or be retired at time t",
    )
    def consistent_branch_extended(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchExtended[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchExtended[branch].indicator_var
                | m.investmentStage[stage].branchRetired[branch].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )

    @m.LogicalConstraint(
        m.stages,
        m.transmission,
        doc="Enforces that, if a branch is installed in period t-1, it must be operational in period t",
    )
    def full_branch_investment(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchInstalled[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchOperational[branch].indicator_var
            )
            if stage != 1
            else pyo.LogicalConstraint.Skip
        )


def add_transmission_state_disjuncts(m, b, i_p):
    """This method defines constraints and a Disjunction with
    disjuncts representing the alternatives for transmission lines
    state operation. The alternatives are:

    branchInUse:     Line is used.
    branchNotInUse:  Line is not used.

    """

    @b.Disjunct(m.transmission)
    def branchInUse(disj, branch):
        b = disj.parent_block()

        def bus_angle_bounds(disj, bus, doc="Voltage angle"):
            return (-1000, 1000)
            return (-math.pi / 6, math.pi / 6)

        # Create bus angle variables for the buses associated with
        # the branch that is in use
        disj.branch_buses = [
            bb
            for bb in m.buses
            if (
                m.transmission[branch]["from_bus"] == bb
                or m.transmission[branch]["to_bus"] == bb
            )
        ]

        disj.busAngle = pyo.Var(
            disj.branch_buses,
            domain=pyo.Reals,
            initialize=0,
            bounds=bus_angle_bounds,
        )

        def delta_bus_angle_bounds(disj, bus, doc="Voltage angle"):
            return (-math.pi / 6, math.pi / 6)

        def delta_bus_angle_rule(disj, doc="Maximum bus angle discrepancy"):
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            return disj.busAngle[tb] - disj.busAngle[fb]

        disj.deltaBusAngle = pyo.Var(
            domain=pyo.Reals,
            bounds=delta_bus_angle_bounds,
            rule=delta_bus_angle_rule,
        )

        if m.config["flow_model"] == "ACP":
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            resistance = m.md.data["elements"]["branch"][branch].get("resistance", 0.0)
            reactance = m.md.data["elements"]["branch"][branch].get("reactance", 1e-6)

            # Transformer tap ratio and phase shift
            if m.md.data["elements"]["branch"][branch]["branch_type"] == "transformer":
                reactance *= m.md.data["elements"]["branch"][branch][
                    "transformer_tap_ratio"
                ]
                phase_shift = m.md.data["elements"]["branch"][branch][
                    "transformer_phase_shift"
                ]
            else:
                phase_shift = 0

            admittance = 1 / complex(resistance, reactance)
            G = admittance.real
            B = admittance.imag

            # Define voltage magnitude variables for from and to buses
            disj.voltage_from = pyo.Var(bounds=(0, 2))
            disj.voltage_to = pyo.Var(bounds=(0, 2))

            # Define active and reactive power flow variables
            disj.P_flow = pyo.Var(bounds=(-1000, 1000))
            disj.Q_flow = pyo.Var(bounds=(-1000, 1000))

            # Polar Active Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_p(disj):
                return disj.P_flow == (
                    disj.voltage_from**2 * G
                    - disj.voltage_from
                    * disj.voltage_to
                    * (
                        G * pyo.cos(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                        + B
                        * pyo.sin(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                    )
                )

            # Polar Reactive Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_q(disj):
                return disj.Q_flow == (
                    -disj.voltage_from**2 * B
                    - disj.voltage_from
                    * disj.voltage_to
                    * (
                        G * pyo.sin(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                        - B
                        * pyo.cos(disj.busAngle[fb] - disj.busAngle[tb] + phase_shift)
                    )
                )

        if m.config["flow_model"] == "ACR":
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            resistance = m.md.data["elements"]["branch"][branch].get("resistance", 0.0)
            reactance = m.md.data["elements"]["branch"][branch].get("reactance", 1e-6)

            # Transformer tap ratio and phase shift
            if m.md.data["elements"]["branch"][branch]["branch_type"] == "transformer":
                reactance *= m.md.data["elements"]["branch"][branch][
                    "transformer_tap_ratio"
                ]
                phase_shift = m.md.data["elements"]["branch"][branch][
                    "transformer_phase_shift"
                ]
            else:
                phase_shift = 0

            admittance = 1 / complex(resistance, reactance)
            G = admittance.real
            B = admittance.imag

            # Define rectangular voltage variables for from and to buses
            disj.real_voltage_from = pyo.Var(bounds=(0, 2))
            disj.real_voltage_to = pyo.Var(bounds=(-2, 2))
            disj.imag_voltage_from = pyo.Var(bounds=(0, 2))
            disj.imag_voltage_to = pyo.Var(bounds=(-2, 2))

            # Define active and reactive power flow variables
            disj.P_flow = pyo.Var(bounds=(-1000, 1000))
            disj.Q_flow = pyo.Var(bounds=(-1000, 1000))

            # Rectangular Active Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_p(disj):
                Vf_r = disj.real_voltage_from
                Vf_i = disj.imag_voltage_from
                Vt_r = disj.real_voltage_to
                Vt_i = disj.imag_voltage_to

                # Active Power Flow Equation
                return disj.P_flow == (
                    G * (Vf_r**2 + Vf_i**2)
                    - G * (Vf_r * Vt_r + Vf_i * Vt_i)
                    - B * (Vf_r * Vt_i - Vf_i * Vt_r)
                )

            # Rectangular Reactive Power Flow Constraint
            @disj.Constraint()
            def ac_power_flow_q(disj):
                Vf_r = disj.real_voltage_from
                Vf_i = disj.imag_voltage_from
                Vt_r = disj.real_voltage_to
                Vt_i = disj.imag_voltage_to

                # Reactive Power Flow Equation
                return disj.Q_flow == (
                    B * (Vf_r**2 + Vf_i**2)
                    + B * (Vf_r * Vt_r + Vf_i * Vt_i)
                    - G * (Vf_r * Vt_i - Vf_i * Vt_r)
                )

        if m.config["flow_model"] == "DC":

            @disj.Constraint()
            def dc_power_flow(disj):
                fb = m.transmission[branch]["from_bus"]
                tb = m.transmission[branch]["to_bus"]
                reactance = m.md.data["elements"]["branch"][branch]["reactance"]
                if (
                    m.md.data["elements"]["branch"][branch]["branch_type"]
                    == "transformer"
                ):
                    reactance *= m.md.data["elements"]["branch"][branch][
                        "transformer_tap_ratio"
                    ]
                    shift = m.md.data["elements"]["branch"][branch][
                        "transformer_phase_shift"
                    ]
                else:
                    shift = 0

                # [TODO: Fix the units in this constraint.]
                return b.powerFlow[branch] / u.MW == (-1 / reactance) * (
                    disj.busAngle[tb] - disj.busAngle[fb] + shift
                )

    @b.Disjunct(m.transmission)
    def branchNotInUse(disj, branch):

        # Fixing power flow to 0 and not creating bus angle variables
        # for branches that are not in use.
        @disj.Constraint()
        def dc_power_flow(disj):
            return b.powerFlow[branch] == 0 * u.MW

        return

    # Branches are either in-use or not. This disjunction may provide the
    # basis for transmission switching in the future.
    @b.Disjunction(m.transmission)
    def branchInUseStatus(disj, branch):
        return [disj.branchInUse[branch], disj.branchNotInUse[branch]]

    # [ESR: Do we really need this flag here if we are already adding
    # the branches?]
    if m.config["transmission"]:

        # [TODO: Update this when switching is implemented.]
        @b.LogicalConstraint(
            m.transmission,
            doc="Enforces that, if a branch is in use, it must be active",
        )
        def must_use_active_branches(b, branch):
            return b.branchInUse[branch].indicator_var.implies(
                pyo.lor(
                    i_p.branchOperational[branch].indicator_var,
                    i_p.branchInstalled[branch].indicator_var,
                    i_p.branchExtended[branch].indicator_var,
                )
            )

        # JSC update - If a branch is not in use, it must be inactive.
        # Update this when switching is implemented
        @b.LogicalConstraint(m.transmission)
        def cannot_use_inactive_branches(b, branch):
            return b.branchNotInUse[branch].indicator_var.implies(
                pyo.lor(
                    i_p.branchDisabled[branch].indicator_var,
                    i_p.branchRetired[branch].indicator_var,
                )
            )
