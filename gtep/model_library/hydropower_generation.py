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

# def add_daily_hydropower_limits(b, commitmentPeriod):
#     m = b.model()

#     # TODO move capacity and minimum constraints to here and let hydro gen variables have overall bounds
#     @b.Constraint(
#         m.hydroGenerators,
#         doc="Enforce commitment period maximum hydropower generation in addition nameplate hydropower capacity",
#     )
#     def strict_hydro_maximum(b, commitmentPeriod, gen):
#         return b.hydroGeneration <= 0

#     @b.Constraint(
#         m.hydroGenerators,
#         doc="Enforce commitment period minimum hydropower generation",
#     )
#     def strict_hydro_minimum(b, commitmentPeriod, gen):
#         return b.hydroGeneration >= 0


def fix_hydropower_limits(b, commitmentPeriod):

    m = b.model()

    b.hydroCapacityExpected = {}
    b.hydroMinimumExpected = {}
    b.hydroAverageExpected = {}

    units_renewable_capacity = u.MW

    for hydroGen in m.hydroGenerators:
        # print(m.md.data["elements"]["generator"][hydroGen])
        b.hydroCapacityExpected[hydroGen] = (
            m.md.data["elements"]["generator"][hydroGen]["p_max"]["values"][
                commitmentPeriod - 1
            ]
            * units_renewable_capacity
        )
        b.hydroMinimumExpected[hydroGen] = (
            m.md.data["elements"]["generator"][hydroGen]["p_min"]["values"][
                commitmentPeriod - 1
            ]
            * units_renewable_capacity
        )
        b.hydroAverageExpected[hydroGen] = (
            m.md.data["elements"]["generator"][hydroGen]["p_average"]["values"][
                commitmentPeriod - 1
            ]
            * units_renewable_capacity
        )


def add_representative_hydropower_average(b, commitmentPeriod):
    m = b.model()

    @b.Constraint(
        m.hydroGenerators,
        doc="Enforce total hydropower generation requirements for given representative period",
    )
    def average_hydro_generation(b, hydroGen):
        return sum(
            b.c_p.d_p.hydroGeneration[hydroGen]
            for c_p in b.commitmentPeriods
            for d_p in c_p.dispatchPeriods
        ) == sum(b.c_p.hydroAverageExpected[hydroGen] for c_p in b.commitmentPeriods)
