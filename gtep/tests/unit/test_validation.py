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

from os.path import abspath, join, dirname, isdir
from os import remove, rmdir, listdir, makedirs
import pyomo.common.unittest as unittest
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.common.tempfiles import TempfileManager
import pandas as pd
from numbers import Number

from gtep.validation import (
    extract_primals_last_investment_stage,
    extract_variable_values,
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
        input_output_pairs = [("var[i]", "i"), ("var", "var"), ("var]", "var]")]
        for input, output in input_output_pairs:
            self.assertEqual(safe_extract_variable_index(input), output)

    def test_extract_primals_last_investment_stage_and_variable_values(self):
        primals = extract_primals_last_investment_stage(self.solution)
        self.assertIsInstance(primals, dict)
        for key, val in primals.items():
            self.assertIsInstance(key, str)
            self.assertIsInstance(val, dict)

        var_vals = extract_variable_values(primals, "gen")
        self.assertIsInstance(var_vals, dict)
        for key, val in var_vals.items():
            self.assertIsInstance(key, str)
            self.assertIn("gen", key)
            self.assertRegex(key, r"Extended|Operational|Installed")
            self.assertIsInstance(
                val, Number
            )  # or would it ever be something else, e.g. bool?

    def test_sum_variable_values_by_index(self):
        input = {"var1[i]": 0.1, "var2[i]": 0.2, "var1[j]": 0.6}
        output = sum_variable_values_by_index(input)
        expected = {"i": 0.3, "j": 0.6}

        self.assertEqual(output.keys(), expected.keys())
        for idx in output:
            self.assertAlmostEqual(output[idx], expected[idx])

    def helper_test_one_dataframe_write(self, fname, dir):
        self.assertIn(fname, listdir(dir))
        test_csv = pd.read_csv(join(dir, fname))
        self.assertTupleEqual(test_csv.shape, (2, 2))
        for item in test_csv.to_numpy().flatten():
            self.assertAlmostEqual(item, 0)

    def test_safe_write_dataframe_to_csv(self):
        with TempfileManager.new_context() as tempfile:
            temp_dir = tempfile.mkdtemp()

            # test writing dataframe
            safe_write_dataframe_to_csv(
                pd.DataFrame([[0, 0], [0, 0]]), temp_dir, "test.csv"
            )
            self.helper_test_one_dataframe_write("test.csv", temp_dir)

            # test writing dataframe to a directory that doesn't yet exist
            output_subdir = join(temp_dir, "test_dir")
            safe_write_dataframe_to_csv(
                pd.DataFrame([[0, 0], [0, 0]]), output_subdir, "test.csv"
            )
            self.helper_test_one_dataframe_write("test.csv", output_subdir)

    def test_populate_generators_filter_pointers(self):
        # filter_pointers needs to access the gen.csv file created in populate_generators
        # so these functions need to be tested together
        with TempfileManager.new_context() as tempfile:
            temp_dir = tempfile.mkdtemp()
            populate_generators(input_data_source, self.solution, temp_dir)
            self.assertIn("gen.csv", listdir(temp_dir))
            filter_pointers(input_data_source, temp_dir)
            self.assertIn("timeseries_pointers.csv", listdir(temp_dir))

    def test_populate_transmission(self):
        with TempfileManager.new_context() as tempfile:
            temp_dir = tempfile.mkdtemp()
            populate_transmission(input_data_source, self.solution, temp_dir)
            self.assertIn("branch.csv", listdir(temp_dir))

    def test_copy_prescient_inputs(self):
        with TempfileManager.new_context() as tempfile:
            temp_dir = tempfile.mkdtemp()
            copy_prescient_inputs(input_data_source, temp_dir)
            for f in listdir(input_data_source):
                if f not in ["gen.csv", "timeseries_pointers.csv", "branch.csv"]:
                    self.assertIn(f, listdir(temp_dir))
