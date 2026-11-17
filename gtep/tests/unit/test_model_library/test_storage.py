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


import pyomo.environ as pyo
from pyomo.environ import units as u
from pyomo.common.unittest import TestCase
from pyomo.environ import (
    SolverFactory,
    value,
    TransformationFactory,
    check_optimal_termination,
)
import pyomo.gdp as gdp

from gtep.tests.unit.pyomo_object_testing import PyomoCheckHelper
from gtep.tests.unit.utils_for_testing import create_model, path_9_bus


class TestObjective(TestCase):
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
            name="storageLifetimes",
            obj_type=pyo.Param,
            units=u.year,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageChargeLimit",
            obj_type=pyo.Param,
            units=u.MW,
        )
        self.check_helper.add_object(
            name="storageDischargeLimit",
            obj_type=pyo.Param,
            units=u.MW,
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
            name="storageInvestmentCost",
            obj_type=pyo.Param,
            units=u.USD / u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageFixedCost",
            obj_type=pyo.Param,
            units=u.USD / u.MW / u.hr,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="storageVarCost",
            obj_type=pyo.Param,
            units=u.USD / u.MW / u.hr,
            index=self.m.storage,
        )

    def _add_storage_state_disjuncts(self):
        self.check_helper.add_object(
            name="battery_storage_balance",
            obj_type=pyo.Constraint,
            units=u.MW * u.hr,
            index=(self.m.storage, self.b.dispatchPeriods),
        )
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
            name="no_charge",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storDischarging",
        )
        self.check_helper.add_object(
            name="storCharging",
            obj_type=gdp.Disjunct,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="charge_limit_min",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storCharging",
        )
        self.check_helper.add_object(
            name="no_discharge",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storCharging",
        )
        self.check_helper.add_object(
            name="storOff",
            obj_type=gdp.Disjunct,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="no_discharge",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storOff",
        )
        self.check_helper.add_object(
            name="no_charge",
            obj_type=pyo.Constraint,
            units=u.MW,
            index=self.b.dispatchPeriods,
            parent="storOff",
        )
        self.check_helper.add_object(
            name="storStatus",
            obj_type=gdp.Disjunction,
            index=self.m.storage,
        )
        self.check_helper.add_object(
            name="commit_active_batts_only",
            obj_type=pyo.LogicalConstraint,
            index=self.m.storage,
        )

    def _coordinate_tests(self, planning_data_args, config):
        """Creates a model and runs tests for a given set of config options."""
        self.m = create_model(
            input_data_path=path_9_bus,
            planning_data_args=planning_data_args,
            config=config,
        ).model

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

    def test_storage_retention_and_periodicity(self):
        # Test the storage charge-level links in a solved model:
        #   - Between consecutive dispatch periods, the retention rate
        #     is applied to the previous charge level. battery_1a and
        #     battery_1b have a retention rate of 0.95 and 1, respectively
        #   - The first dispatch period of a representative period is
        #     linked to the last one, with the retention rate applied
        #     to the last period's charge level.
        # fix_init_storage is set to False (not the default), so the
        # model chooses the initial charge level.
        # Dispatch periods are 1 hour long, so the retained charge
        # level should be exactly the retention rate times the
        # previous charge level.
        modObject = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "num_commit": 6,
                "num_dispatch": 1,
                "duration_representative_period": 6,
            },
            config={
                "include_investment": True,
                "include_commitment": True,
                "include_redispatch": True,
                "transmission": False,
                "storage": True,
                "fix_init_storage": False,
                "flow_model": "DC",
            },
            include_cost_data=False,
        )
        m = modObject.model

        inv_block = m.investmentStage[m.stages.first()]
        rep_block = inv_block.representativePeriod[inv_block.representativePeriods[0]]

        # The initial charge level should not be fixed
        self.assertFalse(hasattr(inv_block, "fix_init_storage_charge_level"))

        opt = SolverFactory("highs")
        if not opt.available():
            raise unittest.SkipTest("Solver not available")

        # Apply transformations to logical terms
        TransformationFactory("gdp.bigm").apply_to(m)

        modObject.results = opt.solve(m)
        self.assertTrue(check_optimal_termination(modObject.results))

        def dispatch_block(pair):
            comm_per, disp_per = pair
            return rep_block.commitmentPeriod[comm_per].dispatchPeriod[disp_per]

        def net_charge(pair, bat):
            # Net energy added in a dispatch period, in MWh
            disp = dispatch_block(pair)
            hours = value(u.convert(disp.dispatchPeriodLength, u.hr))
            return hours * (
                value(disp.storageCharged[bat])
                * value(m.storageChargingEfficiency[bat])
                - value(disp.storageDischarged[bat])
                / value(m.storageDischargingEfficiency[bat])
            )

        pairs = list(rep_block.commitDispatchPairs)
        first_pair, second_pair, last_pair = pairs[0], pairs[1], pairs[-1]

        for pair in pairs:
            self.assertAlmostEqual(
                value(u.convert(dispatch_block(pair).dispatchPeriodLength, u.hr)),
                1.0,
            )

        expected_retention = {
            "battery_1a": 0.95,
            "battery_1b": 1.00,
        }

        for bat, retention in expected_retention.items():
            with self.subTest(battery=bat):
                self.assertAlmostEqual(value(m.storageRetentionRate[bat]), retention)

                before = value(dispatch_block(first_pair).storageChargeLevel[bat])
                after = value(dispatch_block(second_pair).storageChargeLevel[bat])

                # Charge level retained from the first period into the
                # second period
                self.assertAlmostEqual(
                    value(
                        rep_block.retainedStorageChargeLevelFromPrev[second_pair, bat]
                    ),
                    retention * before,
                    places=5,
                )

                # Removing the energy charged or discharged in the second
                # period, the change in charge level is the retention
                # rate
                self.assertAlmostEqual(
                    after - net_charge(second_pair, bat),
                    retention * before,
                    places=5,
                )

                # The first period is linked to the last period, with
                # the retention rate applied to the last period's charge
                # level
                first_level = value(dispatch_block(first_pair).storageChargeLevel[bat])
                last_level = value(dispatch_block(last_pair).storageChargeLevel[bat])
                self.assertAlmostEqual(
                    first_level - net_charge(first_pair, bat),
                    retention * last_level,
                    places=5,
                )
