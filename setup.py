#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import os
from setuptools import setup, find_packages

# requires = [
#     "pyomo",
#     "scipy==1.11",
#     "gridx-egret",
#     "gridx-prescient",
#     "anyio==3.1",
#     "pint",
# ]

requires = [
    "pyomo",
    "scipy",
    "gridx-egret",
    "gridx-prescient",
    "anyio",
    "pint",
    "icecream",
]
setup(
    name="gtep",
    version="0.1.dev0",
    python_requires=">=3.7",
    description="a python package",
    packages=find_packages(),
    install_requires=requires,
)
