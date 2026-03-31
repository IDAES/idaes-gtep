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

"""Scaling for Generation and Transmission Expansion Planning (GTEP)
Model

"""


def add_load_scaling(m, b, commitment_period, investment_stage, scaling_value):

    # Loop over the demand at each bus and scale based on the
    # attributes given in b.
    for load_n in m.load_buses:
        p_load = m.md.data["elements"]["load"][load_n]["p_load"]["values"][
            commitment_period - 1
        ]

        if hasattr(b, "load_scaling"):
            print("**Load scaling exists and will be applied to the load parameters")

            m.loads[load_n] = (
                p_load
                * b.load_scaling[m.md.data["elements"]["load"][load_n]["zone"]].iloc[0]
            )

        else:

            m.loads[load_n] = (
                p_load
                * scaling_value
                * (
                    1
                    + (scaling_value + investment_stage)
                    / (scaling_value + len(m.stages))
                )
            )
