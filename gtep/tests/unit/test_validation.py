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

from os.path import abspath, join, dirname
import pyomo.common.unittest as unittest
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs

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
        'output': abspath(join(curr_dir, "..", "..", "data", "5bus_out")),
        'reference': abspath(join(curr_dir, "..", "..", "data", "5bus_ref"))
    },
}

gtep_model_args_list = [
    {
        'stages': 2,
        'num_reps': 2,
        'len_reps': 1,
        'num_commit': 6,
        'num_dispatch': 4,
     },
    {
        'stages': 1,
        'num_reps': 2,
        'len_reps': 1,
        'num_commit': 6,
        'num_dispatch': 4,
     }
]

# helper function
def get_solution_objects():
    solution_objects_dict = {}

    for data_input_path in data_paths_dict:
        solution_objects_dict[data_input_path] = []

        for gtep_model_args in gtep_model_args_list:

            data_object = ExpansionPlanningData()
            data_object.load_prescient(data_input_path)

            mod_object = ExpansionPlanningModel(**gtep_model_args, data=data_object)
            mod_object.create_model()

            TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
            TransformationFactory("gdp.bigm").apply_to(mod_object.model)
            
            opt = Highs()
            mod_object.results = opt.solve(mod_object.model)

            sol_object = ExpansionPlanningSolution()
            sol_object.load_from_model(mod_object)

            solution_objects_dict[data_input_path].append(sol_object)

    return solution_objects_dict


class TestValidation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.solution_dict = get_solution_objects()

    def iterate_func_over_data_and_options(self, func: function):
        """
        Loops over every combination of data source and options,
        calling `func` each time.
        
        :param func: Must accept three positional arguments:
            data input path, solution object, and data output path
        :type func: function
        """
        for data_input_path, solutions in self.solution_dict.items():

            out_path = data_paths_dict[data_input_path]['output']
            for solution in solutions:

                func(data_input_path, solution, out_path)

    def test_populate_generators_filter_pointers(self):
        """
        Tests the function `populate_generators`, then the function
        `filter_pointers`. The latter runs on the output of the
        former, so one test is needed to cover both functions.
        """
        def thingy(input, sol, output):
            populate_generators(input, sol, output)
            filter_pointers(input, output)
            assert True

        self.iterate_func_over_data_and_options(thingy)

    def test_populate_transmission(self):
        self.iterate_func_over_data_and_options(populate_transmission)

    def test_clone_timeseries(self):
        self.iterate_func_over_data_and_options(
            lambda input, _, output: clone_timeseries(input, output)
        )