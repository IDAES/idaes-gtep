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

from os.path import abspath, join, dirname, exists
import pyomo.common.unittest as unittest
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs

import pandas as pd
import pandas.testing as pdt
from os import listdir

from gtep.validation import (
    clone_timeseries,
    filter_pointers,
    populate_generators,
    populate_transmission,
)

curr_dir = dirname(abspath(__file__))
input_data_source = abspath(join(curr_dir, "..", "..", "data", "5bus"))
output_data_source = abspath(join(curr_dir, "..", "..", "data", "5bus_out"))


def get_solution_object():
    data_object = ExpansionPlanningData()
    data_object.load_prescient(input_data_source)

    mod_object = ExpansionPlanningModel(
        stages=2,
        data=data_object,
        num_reps=2,
        len_reps=1,
        num_commit=6,
        num_dispatch=4,
    )
    mod_object.create_model()
    TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
    TransformationFactory("gdp.bigm").apply_to(mod_object.model)
    opt = Highs()
    mod_object.results = opt.solve(mod_object.model)

    sol_object = ExpansionPlanningSolution()
    sol_object.load_from_model(mod_object)
    return sol_object


def read_pandas_frames_and_assert_equal(path_1, path_2):
    """
    Read in two csvs as `pandas.DataFrame` objects and assert they
    are equal using `pandas.testing.assert_frame_equal`.

    All the asserts happen in this function.

    :param path_1:
    :param path_2:
    :type path_1: str
    :type path_2: str
    """

    assert exists(path_1)
    assert exists(path_2)
    out_gen_df = pd.read_csv(path_1)
    ref_gen_df = pd.read_csv(path_2)

    # TODO: check if any options should be used here
    pdt.assert_frame_equal(out_gen_df, ref_gen_df)


class TestValidation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.solution = get_solution_object()

    def test_populate_generators_filter_pointers(self):
        populate_generators(input_data_source, self.solution, output_data_source)
        # filter_pointers needs to access the gen.csv file created in populate_generators
        # so these functions need to be tested together
        filter_pointers(input_data_source, output_data_source)

    def test_populate_transmission(self):
        populate_transmission(input_data_source, self.solution, output_data_source)

    def test_clone_timeseries(self):
        clone_timeseries(input_data_source, output_data_source)
