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

import pytest
import pyomo.common.unittest as unittest
from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper
import pyomo.environ as pyo
from pyomo.environ import units as u
from gtep.model_library.objective import create_objective_function
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing
from pathlib import Path

curr_dir = Path(__file__).resolve().parent
input_data_source = (
    curr_dir / ".." / ".." / ".." / "data" / "9_bus_GTEP_dir"
).resolve()

bus_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "Bus_data_gen_weights_mappings.csv"
).resolve()
cost_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
).resolve()
ng_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "Total_Energy_Supply_Disposition_and_Price_Summary.csv"
).resolve()


class TestObjective(unittest.TestCase):
    def _create_model(self, config={}):
        # create model
        data_object = ExpansionPlanningData(
            stages=1,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
        )
        data_object.load_prescient(input_data_source)
        data_object.load_storage_csv(str(input_data_source))

        candidate_gens = [
            "Natural Gas_FE",
            "Solar - Utility PV",
            "Land-Based Wind",
        ]

        data_processing_object = DataProcessing()
        data_processing_object.load_gen_data(
            bus_data_path=bus_data_path,
            cost_data_path=cost_data_path,
            ng_cost_path=ng_data_path,
            candidate_gens=candidate_gens,
            save_csv=False,
        )

        mod_object = ExpansionPlanningModel(
            data=data_object,
            cost_data=data_processing_object,
        )
        mod_object.create_model()

        for config_option, config_val in config.items():
            mod_object.config[config_option] = config_val

        return mod_object.model

    def _make_expressions(self):
        # add necessary expressions
        self.check_helper.add_object(
            name="operatingCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="storageCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="expansionCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="penaltyCostTotal",
            units=u.USD,
            obj_type=pyo.Expression,
        )

        self.check_helper.add_object(
            name="total_cost_objective_rule",
            obj_type=pyo.Objective,
        )

    def _coordinate_tests(self, config):
        """Creates a model and runs tests for a given set of config options."""
        self.m = self._create_model(config=config)

        self.check_helper = PyomoCheckHelper(self, self.m)
        self._make_expressions()
        # self.check_helper.check_all_objects()

        for properties in self.check_helper.object_properties:
            if properties["cond"]:
                self.check_helper._check_exists(properties)
                obj = self.check_helper.parent.component(properties["name"])
                if properties["name"] != "total_cost_objective_rule":
                    self.check_helper._check_units(obj, properties)

                self.check_helper._check_type(obj, properties)
                self.check_helper._check_index(obj, properties)
                self.check_helper._check_bounds(obj)
                self.check_helper._run_check_func(obj, properties)

    def test_check_expressions_with_storage(self):
        self._coordinate_tests(config={"storage": True})

    def test_check_expressions_no_storage(self):
        # check that each of the expressions were created
        self._coordinate_tests(config={"storage": False})
        # check that storage cost was set to 0
        self.assertEqual(pyo.value(self.m.storageCostTotal), 0)
