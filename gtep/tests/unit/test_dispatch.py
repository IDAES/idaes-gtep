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


from pathlib import Path

import pyomo.common.unittest as unittest
import pyomo.environ as pyo
from pyomo.environ import units as u

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper

curr_dir = Path(__file__).resolve().parent
input_data_source = (curr_dir / ".." / ".." / "data" / "5bus").resolve()


def check_CP_flow_balance(td, c: pyo.Constraint):
    """Checks the copper-plate flow balance constraint."""
    expected = [td.b.thermalGeneration[g] for g in td.m.thermalGenerators]
    expected += [td.b.renewableGeneration[g] for g in td.m.renewableGenerators]
    if td.m.config["storage"]:
        expected += [td.b.storageDischarged[bt] for bt in td.m.storage]
        expected += [td.b.storageCharged[bt] for bt in td.m.storage]
    expected += [td.b.parent_block().loads[bus] for bus in td.m.buses]
    expected += [td.b.loadShed[bus] for bus in td.m.buses]

    td.check_helper.check_expr_contains(c, expected)


def check_flow_balance(td, c: pyo.Constraint):
    """Checks the flow balance constraint."""
    expected = {}
    for i in c.index_set():
        expected[i] = [td.b.parent_block().loads[i], td.b.loadShed[i]]
        expected[i] += [
            td.b.powerFlow[line]
            for line in td.m.lines
            if (td.m.to_bus[line] == i) or (td.m.from_bus[line] == i)
        ]
        expected[i] += [
            td.b.renewableGeneration[g]
            for g in td.m.renewableGenerators
            if td.m.md.data["elements"]["generator"][g]["bus"] == i
        ]
        expected[i] += [
            td.b.thermalGeneration[g]
            for g in td.m.thermalGenerators
            if td.m.md.data["elements"]["generator"][g]["bus"] == i
        ]
        if td.m.config["storage"]:
            expected[i] += [
                td.b.storageDischarged[bat]
                for bat in td.m.storage
                if td.m.md.data["elements"]["storage"][bat]["bus"] == i
            ]
            expected[i] += [
                td.b.storageCharged[bat]
                for bat in td.m.storage
                if td.m.md.data["elements"]["storage"][bat]["bus"] == i
            ]

    td.check_helper.check_expr_contains(c, expected)


def check_capacity_factor(td, c: pyo.Constraint):
    expected = {
        i: [
            td.b.renewableGeneration[i],
            td.b.renewableCurtailment[i],
            td.b.parent_block().renewableCapacityExpected[i],
        ]
        for i in c.index_set()
    }

    td.check_helper.check_expr_contains(c, expected)


class TestDispatch(unittest.TestCase):

    def _create_model(self, config={}):
        # create model
        data_object = ExpansionPlanningData(
            stages=1,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
        )
        data_object.load_prescient(input_data_source)

        mod_object = ExpansionPlanningModel(
            data=data_object,
        )
        for config_option, config_val in config.items():
            mod_object.config[config_option] = config_val

        mod_object.create_model()
        return mod_object

    def _get_first_dispatch_block(self):
        current_block = self.m
        for component_name in [
            "investmentStage",
            "representativePeriod",
            "commitmentPeriod",
            "dispatchPeriod",
        ]:
            block = current_block.component(component_name)
            first_idx = block.index_set().at(1)
            current_block = block[first_idx]

        return current_block

    def _add_expected_properties_for_objects(self):
        self.check_helper.add_object(
            name="renewableGenerationSurplus",
            units=u.MW,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        self.check_helper.add_object(
            name="renewableCurtailmentCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        self.check_helper.add_object(
            name="thermalGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
        )
        self.check_helper.add_object(
            name="renewableGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.renewableGenerators,
        )
        self.check_helper.add_object(
            name="reactiveGeneratorCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.thermalGenerators,
            cond=(
                self.m.config["flow_model"] == "ACR"
                or self.m.config["flow_model"] == "ACP"
            ),
        )
        self.check_helper.add_object(
            name="loadShed",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.buses,
        )
        self.check_helper.add_object(
            name="loadShedCost",
            units=u.USD,
            obj_type=pyo.Expression,
            index=self.m.buses,
        )
        self.check_helper.add_object(
            name="renewableSurplusDispatch",
            units=u.MW,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="thermalGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="renewableGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="renewableGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="reactiveGenerationCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="loadShedCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="curtailmentCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="operatingCostDispatch",
            units=u.USD,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="renewableCurtailmentDispatch",
            units=u.MW,
            obj_type=pyo.Expression,
        )
        self.check_helper.add_object(
            name="powerFlow",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.lines,
        )
        self.check_helper.add_object(
            name="spinningReserve",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.thermalGenerators,
        )
        self.check_helper.add_object(
            name="quickstartReserve",
            units=u.MW,
            obj_type=pyo.Var,
            index=self.m.thermalGenerators,
        )
        self.check_helper.add_object(
            name="CP_flow_balance",
            units=u.MW,
            obj_type=pyo.Constraint,
            index=None,
            check_func=check_CP_flow_balance,
            cond=self.m.config["flow_model"] == "CP",
        )
        self.check_helper.add_object(
            name="flow_balance",
            units=u.MW,
            obj_type=pyo.Constraint,
            index=self.m.buses,
            check_func=check_flow_balance,
            cond=self.m.config["flow_model"] != "CP",
        )
        self.check_helper.add_object(
            name="capacity_factor",
            units=u.MW,
            obj_type=pyo.Constraint,
            index=self.m.renewableGenerators,
            check_func=check_capacity_factor,
        )

    def _coordinate_tests(self, config):
        """Creates a model and runs tests for a given set of config options."""
        self.m = self._create_model(config=config).model
        self.b = self._get_first_dispatch_block()

        self.check_helper = PyomoCheckHelper(self, self.b)
        self._add_expected_properties_for_objects()
        self.check_helper.check_all_objects()

    def test_default_config_options(self):
        self._coordinate_tests(config={})

    def test_copper_plate(self):
        self._coordinate_tests(config={"flow_model": "CP"})

    def _test_other_config(self):
        pass
        # self._coordinate_tests(config={config_options})
