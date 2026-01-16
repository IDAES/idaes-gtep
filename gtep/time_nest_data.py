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

import awkward as ak
import logging
import pathlib

logger = logging.getLogger("time_nest_data")

# build helper wrapper thingy for awkward arrays
# how to interface with config?
# one thing it should do is send itself to config, most notably the time dict shenanigans


class timeNestData(object):
    def __init__(self):
        pass
