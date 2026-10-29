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
Expansion Planning (GTEP) Model.

This module uses the following transmission-line sets:

    m.lines
        All transmission lines.

    m.staticLines
        Existing/non-candidate lines that cannot be retired, disabled, installed,
        or extended. These lines do not receive GDP disjuncts.

    m.nonStaticLines
        Existing/non-candidate lines that may be operational, retired, disabled,
        extended, or installed/reinstalled.

    m.candidateLines
        Candidate lines that may be operational, retired, disabled, extended,
        or installed.

    m.investmentLines
        union of m.nonStaticLines and m.candidateLines

GDP disjuncts are created only for ``m.nonStaticLines`` and
``m.candidateLines``. Static lines receive deterministic transmission-flow
constraints directly on the relevant block.
"""

import math

import pyomo.environ as pyo
from pyomo.environ import units as u

DELTA_BUS_ANGLE_LOWER = -math.pi / 6
DELTA_BUS_ANGLE_UPPER = math.pi / 6


def _delta_bus_angle_inequality(expr):
    return pyo.inequality(
        DELTA_BUS_ANGLE_LOWER,
        expr,
        DELTA_BUS_ANGLE_UPPER,
    )


def _add_bus_angle_variables(m, b):
    """Create bus-angle variables outside all disjuncts when needed."""
    flow_model = m.config["flow_model"]

    if flow_model in ("CP", "transport"):
        return

    if hasattr(b, "busAngle"):
        return

    def bus_angle_bounds(b, bus):
        # Preserve the broad angle bounds used by the previous implementation.
        # The original code also had an unreachable tighter bound after this
        # return statement.
        return (-1000, 1000)

    b.busAngle = pyo.Var(
        m.buses,
        domain=pyo.Reals,
        initialize=0,
        bounds=bus_angle_bounds,
        doc="Voltage angle",
    )

    @b.Expression(m.lines, doc="Voltage angle difference across each branch")
    def deltaBusAngle(b, branch):
        fb = m.from_bus[branch]
        tb = m.to_bus[branch]
        return b.busAngle[tb] - b.busAngle[fb]


def _branch_dc_data(m, branch):
    """Return DC flow parameters for a branch."""
    fb = m.from_bus[branch]
    tb = m.to_bus[branch]
    branch_data = m.md.data["elements"]["branch"][branch]

    reactance = branch_data["reactance"]

    if branch_data["branch_type"] == "transformer":
        reactance *= branch_data["transformer_tap_ratio"]
        shift = branch_data["transformer_phase_shift"]
    else:
        shift = 0

    return fb, tb, reactance, shift


def _branch_ac_data(m, branch):
    """Return common AC flow parameters for a branch."""
    fb = m.from_bus[branch]
    tb = m.to_bus[branch]
    branch_data = m.md.data["elements"]["branch"][branch]

    resistance = branch_data.get("resistance", 0.0)
    reactance = branch_data.get("reactance", 1e-6)

    if branch_data["branch_type"] == "transformer":
        reactance *= branch_data["transformer_tap_ratio"]
        phase_shift = branch_data["transformer_phase_shift"]
    else:
        phase_shift = 0

    admittance = 1 / complex(resistance, reactance)
    G = admittance.real
    B = admittance.imag

    return fb, tb, G, B, phase_shift


def add_investment_transmission_constraints(m, b, investment_stage):
    investment_lines = m.investmentLines
    candidate_lines = m.candidateLines

    for branch in investment_lines:
        if investment_stage == m.stages.first():
            b.branchOperational[branch].indicator_var.fix(
                m.md.data["elements"]["branch"][branch]["in_service"]
            )

        if branch in candidate_lines:
            if not m.config["include_investment"]:
                b.branchDisabled[branch].indicator_var.fix(True)

    @b.Expression(doc="Transmission investment costs in $")
    def transmission_investment_cost(b):
        m = b.model()

        return sum(
            m.transmissionCapacity[branch]
            * m.branchInvestmentCost[branch]
            * m.branchCapitalMultiplier[branch]
            * b.branchInstalled[branch].indicator_var.get_associated_binary()
            for branch in investment_lines
        ) + sum(
            m.transmissionCapacity[branch]
            * m.branchInvestmentCost[branch]
            * m.branchExtensionMultiplier[branch]
            * b.branchExtended[branch].indicator_var.get_associated_binary()
            for branch in investment_lines
        )


def add_transmission_status_disjuncts(b, transmission_set=None):
    """Add investment/status disjuncts for non-static transmission lines.

    Static lines intentionally receive no investment/status disjuncts.

    The possible investment/status alternatives are:

    branchOperational
        Existing line is active and transmitting power.

    branchInstalled
        Line is newly added or restored and active.

    branchRetired
        Existing line is permanently removed from service.

    branchDisabled
        Line is inactive/out of service.

    branchExtended
        Line is upgraded beyond its original capacity.
    """
    m = b.model()

    # Ignore a legacy caller passing m.lines here; static lines must not get
    # investment disjuncts.
    investment_lines = m.investmentLines

    @b.Disjunct(investment_lines)
    def branchOperational(disj, branch):
        return

    @b.Disjunct(investment_lines)
    def branchInstalled(disj, branch):
        return

    @b.Disjunct(investment_lines)
    def branchRetired(disj, branch):
        return

    @b.Disjunct(investment_lines)
    def branchDisabled(disj, branch):
        return

    @b.Disjunct(investment_lines)
    def branchExtended(disj, branch):
        return

    @b.Disjunction(investment_lines)
    def branchInvestStatus(disj, branch):
        return [
            disj.branchOperational[branch],
            disj.branchInstalled[branch],
            disj.branchRetired[branch],
            disj.branchDisabled[branch],
            disj.branchExtended[branch],
        ]


def add_transmission_logical_constraints(m):
    """Add logical constraints for transmission-line status transitions.

    These constraints are applied only to non-static lines, i.e.,
    ``m.nonStaticLines`` and ``m.candidateLines``. Static lines have no status
    disjuncts and therefore are excluded.
    """
    investment_lines = m.investmentLines

    @m.LogicalConstraint(
        m.stages,
        investment_lines,
        doc=(
            "Enforces that, if a branch is online at time t, it must have "
            "been online or installed at time t-1"
        ),
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
        investment_lines,
        doc=(
            "Enforces that, if a branch is online at time t, it must be "
            "online, extended, or retired at time t+1"
        ),
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
        investment_lines,
        doc=(
            "Enforces that, if a branch is retired in period t-1, it should "
            "be disabled in period t"
        ),
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
        investment_lines,
        doc=(
            "Enforces that, if a branch is disabled at time t-1, it must "
            "stay disabled or be installed at time t"
        ),
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
        investment_lines,
        doc=(
            "Enforces that, if a branch is extended at time t-1, it must "
            "stay extended or be retired at time t"
        ),
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
        investment_lines,
        doc=(
            "Enforces that, if a branch is installed in period t-1, it must "
            "be operational in period t"
        ),
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


def _add_static_transmission_constraints(m, b, static_lines):
    """Add deterministic flow constraints for static transmission lines.

    Static lines do not receive GDP disjuncts. Their flow equations are placed
    directly on the operating block.
    """
    flow_model = m.config["flow_model"]

    if flow_model not in ("CP", "transport"):

        @b.Constraint(
            static_lines,
            doc="Voltage angle difference bounds for static transmission lines",
        )
        def static_delta_bus_angle_bounds(b, branch):
            return _delta_bus_angle_inequality(b.deltaBusAngle[branch])

    if flow_model == "DC":

        @b.Constraint(static_lines, doc="DC power-flow constraints for static lines")
        def static_dc_power_flow(b, branch):
            fb, tb, reactance, shift = _branch_dc_data(m, branch)

            return b.powerFlow[branch] / u.MW == (-1 / reactance) * (
                b.busAngle[tb] - b.busAngle[fb] + shift
            )

    elif flow_model == "ACP":
        # Define voltage magnitude variables for from and to buses
        b.static_voltage_from = pyo.Var(static_lines, bounds=(0, 2))
        b.static_voltage_to = pyo.Var(static_lines, bounds=(0, 2))

        # Define active and reactive power flow variables
        b.static_P_flow = pyo.Var(static_lines, bounds=(-1000, 1000))
        b.static_Q_flow = pyo.Var(static_lines, bounds=(-1000, 1000))

        @b.Constraint(
            static_lines, doc="ACP active power-flow constraints for static lines"
        )
        def static_ac_power_flow_p(b, branch):
            fb, tb, G, B, phase_shift = _branch_ac_data(m, branch)

            return b.static_P_flow[branch] == (
                b.static_voltage_from[branch] ** 2 * G
                - b.static_voltage_from[branch]
                * b.static_voltage_to[branch]
                * (
                    G * pyo.cos(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                    + B * pyo.sin(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                )
            )

        @b.Constraint(
            static_lines, doc="ACP reactive power-flow constraints for static lines"
        )
        def static_ac_power_flow_q(b, branch):
            fb, tb, G, B, phase_shift = _branch_ac_data(m, branch)

            return b.static_Q_flow[branch] == (
                -b.static_voltage_from[branch] ** 2 * B
                - b.static_voltage_from[branch]
                * b.static_voltage_to[branch]
                * (
                    G * pyo.sin(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                    - B * pyo.cos(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                )
            )

    elif flow_model == "ACR":
        # Define rectangular voltage variables for from and to buses
        b.static_real_voltage_from = pyo.Var(static_lines, bounds=(0, 2))
        b.static_real_voltage_to = pyo.Var(static_lines, bounds=(-2, 2))
        b.static_imag_voltage_from = pyo.Var(static_lines, bounds=(0, 2))
        b.static_imag_voltage_to = pyo.Var(static_lines, bounds=(-2, 2))

        # Define active and reactive power flow variables
        b.static_P_flow = pyo.Var(static_lines, bounds=(-1000, 1000))
        b.static_Q_flow = pyo.Var(static_lines, bounds=(-1000, 1000))

        @b.Constraint(
            static_lines, doc="ACR active power-flow constraints for static lines"
        )
        def static_ac_power_flow_p(b, branch):
            fb, tb, G, B, phase_shift = _branch_ac_data(m, branch)

            Vf_r = b.static_real_voltage_from[branch]
            Vf_i = b.static_imag_voltage_from[branch]
            Vt_r = b.static_real_voltage_to[branch]
            Vt_i = b.static_imag_voltage_to[branch]

            return b.static_P_flow[branch] == (
                G * (Vf_r**2 + Vf_i**2)
                - G * (Vf_r * Vt_r + Vf_i * Vt_i)
                - B * (Vf_r * Vt_i - Vf_i * Vt_r)
            )

        @b.Constraint(
            static_lines, doc="ACR reactive power-flow constraints for static lines"
        )
        def static_ac_power_flow_q(b, branch):
            fb, tb, G, B, phase_shift = _branch_ac_data(m, branch)

            Vf_r = b.static_real_voltage_from[branch]
            Vf_i = b.static_imag_voltage_from[branch]
            Vt_r = b.static_real_voltage_to[branch]
            Vt_i = b.static_imag_voltage_to[branch]

            return b.static_Q_flow[branch] == (
                B * (Vf_r**2 + Vf_i**2)
                + B * (Vf_r * Vt_r + Vf_i * Vt_i)
                - G * (Vf_r * Vt_i - Vf_i * Vt_r)
            )


def add_transmission_state_disjuncts(m, b, i_p):
    """Add operating-state disjuncts for non-static transmission lines.

    The operating-state alternatives are:

    branchInUse
        Line is used.

    branchNotInUse
        Line is not used.

    Static lines receive no operating-state disjuncts; they instead receive
    deterministic flow constraints directly on block ``b``.
    """
    investment_lines = m.investmentLines
    static_lines = m.staticLines

    _add_bus_angle_variables(m, b)
    _add_static_transmission_constraints(m, b, static_lines)

    @b.Disjunct(investment_lines)
    def branchInUse(disj, branch):

        if m.config["flow_model"] == "transport":
            return

        b = disj.parent_block()
        flow_model = m.config["flow_model"]

        fb = m.from_bus[branch]
        tb = m.to_bus[branch]

        if flow_model not in ("CP", "transport"):

            @disj.Constraint(
                doc="Voltage angle difference bounds when branch is in use",
            )
            def delta_bus_angle_bounds(disj):
                return _delta_bus_angle_inequality(b.deltaBusAngle[branch])

        if flow_model == "ACP":
            fb, tb, G, B, phase_shift = _branch_ac_data(m, branch)

            # Define voltage magnitude variables for from and to buses
            disj.voltage_from = pyo.Var(bounds=(0, 2))
            disj.voltage_to = pyo.Var(bounds=(0, 2))

            # Define active and reactive power flow variables
            disj.P_flow = pyo.Var(bounds=(-1000, 1000))
            disj.Q_flow = pyo.Var(bounds=(-1000, 1000))

            @disj.Constraint()
            def ac_power_flow_p(disj):
                return disj.P_flow == (
                    disj.voltage_from**2 * G
                    - disj.voltage_from
                    * disj.voltage_to
                    * (
                        G * pyo.cos(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                        + B * pyo.sin(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                    )
                )

            @disj.Constraint()
            def ac_power_flow_q(disj):
                return disj.Q_flow == (
                    -disj.voltage_from**2 * B
                    - disj.voltage_from
                    * disj.voltage_to
                    * (
                        G * pyo.sin(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                        - B * pyo.cos(b.busAngle[fb] - b.busAngle[tb] + phase_shift)
                    )
                )

        elif flow_model == "ACR":
            fb, tb, G, B, phase_shift = _branch_ac_data(m, branch)

            # Define rectangular voltage variables for from and to buses
            disj.real_voltage_from = pyo.Var(bounds=(0, 2))
            disj.real_voltage_to = pyo.Var(bounds=(-2, 2))
            disj.imag_voltage_from = pyo.Var(bounds=(0, 2))
            disj.imag_voltage_to = pyo.Var(bounds=(-2, 2))

            # Define active and reactive power flow variables
            disj.P_flow = pyo.Var(bounds=(-1000, 1000))
            disj.Q_flow = pyo.Var(bounds=(-1000, 1000))

            @disj.Constraint()
            def ac_power_flow_p(disj):
                Vf_r = disj.real_voltage_from
                Vf_i = disj.imag_voltage_from
                Vt_r = disj.real_voltage_to
                Vt_i = disj.imag_voltage_to

                return disj.P_flow == (
                    G * (Vf_r**2 + Vf_i**2)
                    - G * (Vf_r * Vt_r + Vf_i * Vt_i)
                    - B * (Vf_r * Vt_i - Vf_i * Vt_r)
                )

            @disj.Constraint()
            def ac_power_flow_q(disj):
                Vf_r = disj.real_voltage_from
                Vf_i = disj.imag_voltage_from
                Vt_r = disj.real_voltage_to
                Vt_i = disj.imag_voltage_to

                return disj.Q_flow == (
                    B * (Vf_r**2 + Vf_i**2)
                    + B * (Vf_r * Vt_r + Vf_i * Vt_i)
                    - G * (Vf_r * Vt_i - Vf_i * Vt_r)
                )

        elif flow_model == "DC":

            @disj.Constraint()
            def dc_power_flow(disj):
                fb, tb, reactance, shift = _branch_dc_data(m, branch)

                return b.powerFlow[branch] / u.MW == (-1 / reactance) * (
                    b.busAngle[tb] - b.busAngle[fb] + shift
                )

        elif flow_model == "transport":
            # In the transport model, being "in use" means the branch is
            # available to carry flow subject to capacity constraints elsewhere.
            return

    @b.Disjunct(investment_lines)
    def branchNotInUse(disj, branch):
        @disj.Constraint()
        def no_power_flow(disj):
            return b.powerFlow[branch] == 0 * u.MW

    @b.Disjunction(investment_lines)
    def branchInUseStatus(disj, branch):
        return [disj.branchInUse[branch], disj.branchNotInUse[branch]]

    if m.config["transmission"]:

        @b.LogicalConstraint(
            investment_lines,
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

        @b.LogicalConstraint(
            investment_lines,
            doc="Enforces that, if a branch is not in use, it must be inactive",
        )
        def cannot_use_inactive_branches(b, branch):
            return b.branchNotInUse[branch].indicator_var.implies(
                pyo.lor(
                    i_p.branchDisabled[branch].indicator_var,
                    i_p.branchRetired[branch].indicator_var,
                )
            )
