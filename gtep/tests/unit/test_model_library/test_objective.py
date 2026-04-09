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
"""
Unit tests for the function to create the objective function component
"""


from os.path import abspath, join, dirname
import pytest
import pyomo.environ as pyo
import gtep.model_library.components as comps
from gtep.gtep_model import ExpansionPlanningModel, create_stages
from gtep.gtep_data import ExpansionPlanningData
from gtep.model_library.objective import create_objective_function


# Helper Functions
def read_debug_model():
    curr_dir = dirname(abspath(__file__))
    debug_data_path = abspath(join(curr_dir, "..", "..", "data", "5bus"))
    dataObject = ExpansionPlanningData()
    dataObject.load_prescient(debug_data_path)
    return dataObject


def create_test_model():
    # create the model to use in the tests
    dataObject = read_debug_model()
    modObject = ExpansionPlanningModel(data=dataObject)
    m = pyo.ConcreteModel("GTEP Model")
    comps.add_model_sets(
        m, modObject.stages, rep_per=[i for i in range(1, modObject.num_reps + 1)]
    )

    comps.add_model_parameters(
        m, modObject.num_commit, modObject.num_dispatch, modObject.duration_dispatch
    )

    create_stages(m, modObject.stages)

    return m


def test_create_objective_function():
    pass
