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
import pyomo.environ as pyo
from pyomo.environ import units as u

def add_load_scaling(m, b, commitment_period, investment_stage, scaling_value):

    b.loads = pyo.Param(
        m.buses,
        initialize={load_n: 0 for load_n in m.buses},
        mutable=True,
        units=u.MW,
        doc="Demand at each bus",
    )

    b.loads = pyo.Param(
        m.buses,
        initialize={load_n: 0 for load_n in m.buses},
        mutable=True,
        units=u.MW,
        doc="Demand at each bus",
    )

    # Loop over the demand at each bus and scale based on the case.
    for load_n in m.load_buses:
        # print(f"{load_n = }")
        # print(f"{m.md.data["elements"]["load"] = }")
        # print(f"{m.md.data["elements"]["load"][load_n] = }")
        p_load = m.md.data["elements"]["load"][load_n]["p_load"]["values"][
            commitment_period - 1
        ]
        # print(f"{p_load = }")

        if m.config["scale_loads"]:
            temp_scale = 10
            b.loads[load_n] = (
                p_load
                * temp_scale
                * (1 + (temp_scale + investment_stage) / (temp_scale + len(m.stages)))
            )

        elif m.config["scale_texas_loads"]:
            b.loads[load_n] = (
                p_load
                * b.load_scaling[m.md.data["elements"]["load"][load_n]["zone"]].iloc[0]
            )

        else:
            b.loads[load_n] = p_load
    
    #print(f"{b = }")
    #print(f"{pyo.value(sum(b.loads[n] for n in b.loads)) = }")

    # for key, val in b.loads.items():
    #     # print(f"{key=}")
    #     # print(f"{val=}")
    #     b.loads[key] *= 1/3
    # print(f"total load at time period = {sum(b.loads.values())}")
