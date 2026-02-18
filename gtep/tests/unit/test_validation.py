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
import pandas as pd
from math import isclose

from gtep.validation import (
    extract_end_variable_values,
    sum_variable_values_by_index,
    safe_extract_variable_index,
    safe_write_dataframe_to_csv,
    copy_prescient_inputs,
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


class TestValidation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.solution = get_solution_object()

    def test_safe_extract_variable_index(self):
        assert safe_extract_variable_index("var[test]") == "test"
        assert safe_extract_variable_index("var") == "var"
        assert safe_extract_variable_index("var]") == "var]"

    def test_extract_end_variable_values(self):
        extract_end_variable_values(self.solution, "gen")

    def test_sum_variable_values_by_index(self):
        input = {"var1[i]": 0.1, "var2[i]": 0.2, "var1[j]": 0.6}
        output = sum_variable_values_by_index(input)
        expected = {"i": 0.3, "j": 0.6}

        assert output.keys() == expected.keys()
        assert all([isclose(output[idx], expected[idx]) for idx in output.keys()])

    def test_safe_write_dataframe_to_csv(self):
        safe_write_dataframe_to_csv(
            pd.DataFrame([[0, 0], [0, 0]]), output_data_source, "test.csv"
        )

    def test_populate_generators_filter_pointers(self):
        populate_generators(input_data_source, self.solution, output_data_source)
        # filter_pointers needs to access the gen.csv file created in populate_generators
        # so these functions need to be tested together
        filter_pointers(input_data_source, output_data_source)

    def test_populate_transmission(self):
        populate_transmission(input_data_source, self.solution, output_data_source)

    def test_copy_prescient_inputs(self):
        copy_prescient_inputs(input_data_source, output_data_source)
