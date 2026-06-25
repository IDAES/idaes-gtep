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

import pyomo.environ as pyo
from pyomo.environ import units as u
from pyomo.common.unittest import TestCase
import pyomo.gdp as gdp

from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing
from gtep.gtep_model import ExpansionPlanningModel
from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper

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


class TestObjective(TestCase):
    def _create_model(self, planning_data_args={}, config={}):
        # create model
        default_data_planning_args = dict(
            stages=1,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
        )
        for arg, val in default_data_planning_args.items():
            if arg not in planning_data_args:
                planning_data_args[arg] = val

        data_object = ExpansionPlanningData(**planning_data_args)
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

        for config_option, config_val in config.items():
            mod_object.config[config_option] = config_val
        mod_object.create_model()

        return mod_object.model

    def _add_storage_params(self):
        self.check_helper.add_object(
            name="storageCapacity",
            obj_type=pyo.Param,
            units=u.MW * u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="initStorageChargeLevel",
            obj_type=pyo.Param,
            units=u.MW * u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="minStorageChargeLevel",
            obj_type=pyo.Param,
            units=u.MW * u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="chargingCost",
            obj_type=pyo.Param,
            units=u.USD / u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="dischargingCost",
            obj_type=pyo.Param,
            units=u.USD / u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="dischargeMin",
            obj_type=pyo.Param,
            units=u.MW,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="dischargeMax",
            obj_type=pyo.Param,
            units=u.MW,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="chargeMin",
            obj_type=pyo.Param,
            units=u.MW,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="chargeMax",
            obj_type=pyo.Param,
            units=u.MW,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageDischargingRampUpRates",
            obj_type=pyo.Param,
            units=u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageDischargingRampDownRates",
            obj_type=pyo.Param,
            units=u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageChargingRampUpRates",
            obj_type=pyo.Param,
            units=u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageChargingRampDownRates",
            obj_type=pyo.Param,
            units=u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageDischargingEfficiency",
            obj_type=pyo.Param,
            units=u.dimensionless,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageChargingEfficiency",
            obj_type=pyo.Param,
            units=u.dimensionless,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageRetentionRate",
            obj_type=pyo.Param,
            units=1 / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageCapitalMultiplier",
            obj_type=pyo.Param,
            units=u.dimensionless,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageExtensionMultiplier",
            obj_type=pyo.Param,
            units=u.dimensionless,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageInvestmentCost",
            obj_type=pyo.Param,
            units=u.USD / u.MW / u.hr,
            index=self.m.storage,
        )

    def _add_storage_state_disjuncts(self):
        continuity_constraint_skipped = [
            d
            for d in self.b.dispatchPeriods
            if (self.b.parent_block().index(), d)
            == self.b.parent_block().commitDispatchPairs.first()
        ]
        self.check_helper.add_object(
            name="storDischarging",
            obj_type=gdp.Disjunct,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="discharge_limit_min",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
        )
        self.check_helper.add_object(
            name="discharge_limit_max",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
        )
        self.check_helper.add_object(
            name="discharge_ramp_up_limits",
            obj_type=pyo.Constraint,
            units=u.MW / u.hr,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
            skipped_on=continuity_constraint_skipped,
        )
        self.check_helper.add_object(
            name="discharge_ramp_down_limits",
            obj_type=pyo.Constraint,
            units=u.MW / u.hr,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
            skipped_on=continuity_constraint_skipped,
        )
        self.check_helper.add_object(
            name="no_charge",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
        )
        self.check_helper.add_object(
            name="discharging_battery_storage_balance",
            obj_type=pyo.Constraint,
            units=u.MW * u.hr,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
            skipped_on=continuity_constraint_skipped,
        )

    def _coordinate_tests(self, planning_data_args, config):
        """Creates a model and runs tests for a given set of config options."""
        self.m = self._create_model(
            planning_data_args=planning_data_args, config=config
        )

        # checks on model attributes
        self.check_helper = PyomoCheckHelper(self, self.m)
        self._add_storage_params()
        self.check_helper.check_all_objects()

        # checks on commitment block attributes
        self.b = self.m.investmentStage[1].representativePeriod[1].commitmentPeriod[1]
        self.check_helper = PyomoCheckHelper(self, self.b)
        self._add_storage_state_disjuncts()
        self.check_helper.check_all_objects()

    def test_(self):
        self._coordinate_tests(
            planning_data_args={"num_commit": 2, "num_dispatch": 2},
            config={"storage": True},
        )
