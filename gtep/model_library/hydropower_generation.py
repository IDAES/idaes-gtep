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

"""Advanced implementations of hydropower modeling for GTEP model"""

import pyomo.environ as pyo
from pyomo.environ import units as u


def fix_hydropower_limits(b, commitmentPeriod):

    m = b.model()

    units_renewable_capacity = u.MW

    b.hydroCapacityExpected = pyo.Param(
        m.hydroGenerators,
        initialize={
            hydroGen: (
                m.md.data["elements"]["generator"][hydroGen]["p_max"]["values"][
                    commitmentPeriod - 1
                ]
            )
            for hydroGen in m.hydroGenerators
        },
        # mutable=True,
        units=units_renewable_capacity,
        doc="Expected maximum hydro capacity for each hydro generator",
    )

    b.hydroMinimumExpected = pyo.Param(
        m.hydroGenerators,
        initialize={
            hydroGen: (
                m.md.data["elements"]["generator"][hydroGen]["p_min"]["values"][
                    commitmentPeriod - 1
                ]
            )
            for hydroGen in m.hydroGenerators
        },
        # mutable=True,
        units=units_renewable_capacity,
        doc="Expected minimum hydro generation for each hydro generator",
    )


def add_representative_hydropower_average(b, repPer):
    m = b.model()

    @b.Constraint(
        m.hydroGenerators,
        doc="Enforce total hydropower generation requirements for given representative period",
    )
    def average_hydro_generation(b, hydroGen):
        return sum(
            b.commitmentPeriod[c_p].dispatchPeriod[d_p].hydroGeneration[hydroGen]
            for c_p in b.commitmentPeriods
            for d_p in b.commitmentPeriod[c_p].dispatchPeriods
        ) <= sum(
            b.commitmentPeriod[c_p].hydroAverageExpected[hydroGen]
            for c_p in b.commitmentPeriods
        )
