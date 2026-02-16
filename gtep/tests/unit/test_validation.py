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

#############################################################################
# These tests are designed to run over a collection of models, across both
# different data sources and model parameters. `data_paths_dict` and
# `gtep_model_args_list` define the set of collection of models.
#############################################################################

# maps input path to output and reference paths; potentially include other data sources
# ...just 5 bus for now.
data_paths_dict = {
    abspath(join(curr_dir, "..", "..", "data", "5bus")): {
        "out": abspath(join(curr_dir, "..", "..", "data", "testing", "5bus", "out")),
        "ref": abspath(join(curr_dir, "..", "..", "data", "testing", "5bus", "ref")),
    },
}

gtep_model_args_list = [
    {
        "stages": 2,
        "num_reps": 2,
        "len_reps": 1,
        "num_commit": 6,
        "num_dispatch": 4,
    },
    {
        "stages": 3,
        "num_reps": 3,
        "len_reps": 2,
        "num_commit": 7,
        "num_dispatch": 5,
    },
]


# helper functions
def get_solution_objects() -> dict[str, list[ExpansionPlanningSolution]]:
    """
    :return: Dictionary that maps each data source path to a list of solution objects
    (each corresponding to a set of arguments to the model constructor)
    """

    solution_objects_dict = {}

    for data_input_path in data_paths_dict:
        solution_objects_dict[data_input_path] = []

        for gtep_model_args in gtep_model_args_list:

            data_object = ExpansionPlanningData()
            data_object.load_prescient(data_input_path)

            mod_object = ExpansionPlanningModel(**gtep_model_args, data=data_object)
            mod_object.create_model()

            TransformationFactory("gdp.bound_pretransformation").apply_to(
                mod_object.model
            )
            TransformationFactory("gdp.bigm").apply_to(mod_object.model)

            opt = Highs()
            mod_object.results = opt.solve(mod_object.model)

            sol_object = ExpansionPlanningSolution()
            sol_object.load_from_model(mod_object)

            solution_objects_dict[data_input_path].append(sol_object)

    return solution_objects_dict


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
        cls.solution_dict = get_solution_objects()

    def iterate_func_over_data_and_options(self, func: function):
        """
        Loops over every combination of data source and options,
        calling `func` each time.

        :param func: Must accept four positional arguments:
            data input path, solution object, data output path, and data output reference path
        :type func: function
        """

        for data_input_path, solutions in self.solution_dict.items():

            out_path = data_paths_dict[data_input_path]["out"]
            ref_path = data_paths_dict[data_input_path]["ref"]
            for i, solution in enumerate(solutions):

                model_tag = f"model_setup_{i}"
                func(
                    data_input_path,
                    solution,
                    join(out_path, model_tag),
                    join(ref_path, model_tag),
                )

    def test_populate_generators_filter_pointers(self):
        """
        Tests the function `populate_generators`, then the function
        `filter_pointers`. The latter runs on the output of the
        former, so one test is needed to cover both functions.
        """

        def single_check_populate_generators_and_filter_pointers(
            input, soln, output, ref
        ):
            """
            Runs tests on `populate_generators` and `filter_pointers` for a single model setup.
            """

            populate_generators(input, soln, output)
            read_pandas_frames_and_assert_equal(
                join(output, "gen.csv"), join(ref, "gen.csv")
            )

            filter_pointers(input, output)
            read_pandas_frames_and_assert_equal(
                join(output, "timeseries_pointers.csv"),
                join(ref, "timeseries_pointers.csv"),
            )

        self.iterate_func_over_data_and_options(
            single_check_populate_generators_and_filter_pointers
        )

    def test_populate_transmission(self):
        """
        Tests the function `populate_transmission`.
        """

        def single_check_populate_transmission(input, soln, output, ref):
            """
            Runs tests on `populate_transmission` for a single model setup.
            """

            populate_transmission(input, soln, output)
            read_pandas_frames_and_assert_equal(
                join(output, "branch.csv"), join(ref, "branch.csv")
            )

        self.iterate_func_over_data_and_options(single_check_populate_transmission)

    def test_clone_timeseries(self):
        """
        Tests the function `clone_timeseries`.
        """

        def single_check_clone_timeseries(input, _, output, ref):
            """
            Runs tests on `clone_timeseries` for a single model setup.
            """

            clone_timeseries(input, output)

            for file in listdir(input):
                if file in ["branch.csv", "gen.csv", "timeseries_pointers.csv"]:
                    continue  # skip the files handled by the other functions

                read_pandas_frames_and_assert_equal(join(output, file), join(ref, file))

        self.iterate_func_over_data_and_options(single_check_clone_timeseries)
