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

import pyomo.environ as pyo


def add_transmission_constraints(m, b, investment_stage):
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


def add_transmission_disjuncts(b, transmission_set):
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
    # we need to define different states for transmission lines.]
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
