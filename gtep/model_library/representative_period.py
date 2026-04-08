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

"""Variables and Constraints for the Representative Periods in the
Generation and Transmission Expansion Planning (GTEP) Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u
from math import ceil


def add_representative_period_variables(b, rep_per):
    m = b.model()
    i_p = b.parent_block()

    # [ESR WIP: This variable is never used. Should we remove it?]
    b.renewableSurplusRepresentative = pyo.Var(
        within=pyo.Reals, initialize=0, units=u.USD
    )


def add_representative_period_logical_constraints(b, rep_per):
    m = b.model()
    i_p = b.parent_block()

    # [TODO: This needs to be updated for variable length
    # commitment periods. Do this by (pre) processing the set of
    # commitment periods for req_shutdown_periods.]
    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_shutdown(b, commitmentPeriod, thermalGen):
        req_shutdown_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"])
        )
        return (
            pyo.atmost(
                req_shutdown_periods - 1,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genShutdown[thermalGen]
                    .indicator_var
                    for j in range(min(req_shutdown_periods, commitmentPeriod - 1))
                ],
            ).land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genShutdown[thermalGen]
                .indicator_var
            )
            # | b.commitmentPeriod[commitmentPeriod-1].genOn.indicator_var)
            .implies(
                b.commitmentPeriod[commitmentPeriod]
                .genShutdown[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else pyo.LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_off_after_shutdown(b, commitmentPeriod, thermalGen):
        req_shutdown_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"])
        )
        return (
            pyo.atleast(
                req_shutdown_periods,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genShutdown[thermalGen]
                    .indicator_var
                    for j in range(min(req_shutdown_periods, commitmentPeriod - 1))
                ],
            )
            .land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genShutdown[thermalGen]
                .indicator_var
            )
            .implies(
                b.commitmentPeriod[commitmentPeriod].genOff[thermalGen].indicator_var
            )
            if commitmentPeriod != 1
            else pyo.LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_startup(b, commitmentPeriod, thermalGen):
        req_startup_periods = ceil(
            1
            # / float(m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"])
            / pyo.value(m.rampUpRates[thermalGen])
        )
        return (
            pyo.atmost(
                req_startup_periods - 1,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genStartup[thermalGen]
                    .indicator_var
                    for j in range(min(req_startup_periods, commitmentPeriod - 1))
                ],
            ).land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genStartup[thermalGen]
                .indicator_var
            )
            # | b.commitmentPeriod[commitmentPeriod-1].genOn.indicator_var)
            .implies(
                b.commitmentPeriod[commitmentPeriod]
                .genStartup[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else pyo.LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_on_after_startup(b, commitmentPeriod, thermalGen):
        req_startup_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"])
        )
        return (
            pyo.atleast(
                req_startup_periods,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genStartup[thermalGen]
                    .indicator_var
                    for j in range(min(req_startup_periods, commitmentPeriod - 1))
                ],
            )
            .land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genStartup[thermalGen]
                .indicator_var
            )
            .implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
            )
            if commitmentPeriod != 1
            else pyo.LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_uptime(b, commitmentPeriod, thermalGen):
        return (
            pyo.atmost(
                int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"]) - 1,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genOn[thermalGen]
                    .indicator_var
                    for j in range(
                        min(
                            int(
                                m.md.data["elements"]["generator"][thermalGen][
                                    "min_up_time"
                                ]
                            ),
                            commitmentPeriod - 1,
                        )
                    )
                ],
            )
            .land(
                b.commitmentPeriod[commitmentPeriod - 1].genOn[thermalGen].indicator_var
            )
            .implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
            )
            if commitmentPeriod
            != 1  # int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"])+1
            else pyo.LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_downtime(b, commitmentPeriod, thermalGen):
        return (
            (
                pyo.atmost(
                    int(m.md.data["elements"]["generator"][thermalGen]["min_down_time"])
                    - 1,
                    [
                        b.commitmentPeriod[commitmentPeriod - j - 1]
                        .genOff[thermalGen]
                        .indicator_var
                        for j in range(
                            min(
                                int(
                                    m.md.data["elements"]["generator"][thermalGen][
                                        "min_down_time"
                                    ]
                                ),
                                commitmentPeriod - 1,
                            )
                        )
                    ],
                ).land(
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .genOff[thermalGen]
                    .indicator_var
                )
            ).implies(
                b.commitmentPeriod[commitmentPeriod].genOff[thermalGen].indicator_var
            )
            if commitmentPeriod
            != 1  # >= int(m.md.data["elements"]["generator"][thermalGen]["min_down_time"])+1
            else pyo.LogicalConstraint.Skip
        )

    # Add legacy equation. Check with team if we need to keep this
    # equation.
    """
    # [TODO: Is this constraint necessary?]
    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_shutdown_after_uptime(
        b, commitmentPeriod, thermalGen
    ):
        return (
            (
                atleast(
                    int(
                        m.md.data["elements"]["generator"][thermalGen][
                            "min_up_time"
                        ]
                    ),
                    [
                        b.commitmentPeriod[commitmentPeriod - j - 1]
                        .genOn[thermalGen]
                        .indicator_var
                        for j in range(
                            min(
                                int(
                                    m.md.data["elements"]["generator"][thermalGen][
                                        "min_up_time"
                                    ]
                                ),
                                commitmentPeriod - 1,
                            )
                        )
                    ],
                ).land(
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .genOn[thermalGen]
                    .indicator_var
                )
            ).implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
                | b.commitmentPeriod[commitmentPeriod]
                .genShutdown[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )
    """
