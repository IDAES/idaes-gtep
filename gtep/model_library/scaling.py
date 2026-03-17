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


def add_load_scaling(m, b, commitment_period, investment_stage):

    # Create a false load directory. [TODO: Commented for now since it
    # is not used. Check with team if we want to preserve this for
    # Texas case.]
    # false_loads = []
    # for load in m.md.data["elements"]["load"]:
    #     if type(m.md.data["elements"]["load"][load]) == float:
    #         false_loads.append(load)
    # for load in false_loads:
    #     del m.md.data["elements"]["load"][load]
    #     # del m.loads[load]
    # # print(m.loads)

    # Loop over the demand at each bus and scale based on the case.
    for load_n in m.load_buses:
        p_load = m.md.data["elements"]["load"][load_n]["p_load"]["values"][
            commitment_period - 1
        ]

        if m.config["scale_loads"]:
            temp_scale = 10
            m.loads[load_n] = (
                p_load
                * temp_scale
                * (1 + (temp_scale + investment_stage) / (temp_scale + len(m.stages)))
            )

        elif m.config["scale_texas_loads"]:
            m.loads[load_n] = (
                p_load
                * b.load_scaling[m.md.data["elements"]["load"][load_n]["zone"]].iloc[0]
            )

        else:
            m.loads[load_n] = p_load

        # for key, val in b.loads.items():
        #     # print(f"{key=}")
        #     # print(f"{val=}")
        #     b.loads[key] *= 1/3
        # print(f"total load at time period = {sum(b.loads.values())}")
