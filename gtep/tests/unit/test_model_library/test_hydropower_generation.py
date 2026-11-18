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

"""Unit Tests for the advanced hydropower data and model components in
the Generation and Transmission Expansion Planning (GTEP) model.

These tests use the 5-bus case, which has an existing hydro unit
(1_HYDRO) and a candidate hydro unit (1_HYDRO-c), and the
DAY_AHEAD_hydro_min.csv and DAY_AHEAD_hydro_budget.csv files in that
case.
"""

import pandas as pd

import pyomo.common.unittest as unittest
import pyomo.environ as pyo
from pyomo.environ import units as u
from pyomo.util.check_units import assert_units_equivalent

from gtep.gtep_model import ExpansionPlanningModel
from gtep.tests.unit.utils_for_testing import create_data, input_data_path


class TestHydropowerGeneration(unittest.TestCase):

    def create_testing_obj(self, config, modify=None):
        """Creates a 5-bus model for a given set of config options.

        The model has one representative day (2020-04-23) with six
        1-hour commitment periods, so commitment period c uses Period c
        of the hydro time-series files on that day.

        :param config: Non-default model config options
        :param modify: Optional function that takes a generator data dict
                       and modifies it before the model is built. The
                       model reads generator data from the
                       representative-period copies, so it is applied
                       to md and to every copy.
        """
        self.data_object, _ = create_data(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "num_commit": 6,
                "num_dispatch": 1,
                "duration_representative_period": 6,
            },
            prescient_data_args={"representative_dates": ["2020-04-23 00:00"]},
            include_cost_data=False,
        )

        if modify is not None:
            for md in [self.data_object.md, *self.data_object.representative_data]:
                modify(md.data["elements"]["generator"])

        mod_object = ExpansionPlanningModel(
            config=config,
            data=self.data_object,
            cost_data=None,
        )
        mod_object.create_model()
        self.m = mod_object.model

        # Representative period block and its commitment period blocks
        b_inv = self.m.investmentStage[self.m.stages.first()]
        self.b_rep = b_inv.representativePeriod[b_inv.representativePeriods[0]]
        self.b_comm = [
            self.b_rep.commitmentPeriod[c] for c in self.b_rep.commitmentPeriods
        ]

    ### DATA LOADING ###

    def test_load_hydro_timeseries(self):
        # The hydro minimum and budget files are read into the
        # generator data as time series for both hydro units, covering
        # every model time step, with values matching the files
        data_object, _ = create_data(include_cost_data=False)
        generators = data_object.md.data["elements"]["generator"]
        time_keys = data_object.md.data["system"]["time_keys"]

        self.assertEqual(
            set(data_object.hydro_timeseries_files), {"hydro_p_min", "hydro_budgeted"}
        )

        files = {
            "hydro_p_min": input_data_path / "DAY_AHEAD_hydro_min.csv",
            "hydro_budgeted": input_data_path / "DAY_AHEAD_hydro_budget.csv",
        }
        for key, path in files.items():
            df = pd.read_csv(path)
            for gen in ["1_HYDRO", "1_HYDRO-c"]:
                with self.subTest(key=key, gen=gen):
                    series = generators[gen][key]
                    self.assertEqual(series["data_type"], "time_series")
                    self.assertEqual(len(series["values"]), len(time_keys))

                    # The files and the model both start on 2020-01-01
                    # with 1-hour periods
                    expected = df[gen].astype(float).tolist()[: len(time_keys)]
                    self.assertEqual(series["values"], expected)

        # Non-hydro units get no hydro time series
        for key in files:
            self.assertNotIn(key, generators["4_WIND"])
            self.assertNotIn(key, generators["3_CT"])

    def test_representative_data_hydro_timeseries(self):
        # Each representative day carries the hydro time series for that
        # day only
        self.create_testing_obj(config={})
        gen_data = self.data_object.representative_data[0].data["elements"][
            "generator"
        ]["1_HYDRO"]

        self.assertEqual(len(gen_data["hydro_p_min"]["values"]), 24)
        self.assertEqual(len(gen_data["hydro_budgeted"]["values"]), 24)

        # 2020-04-23, Periods 1-6
        self.assertEqual(gen_data["hydro_p_min"]["values"][:6], [0.0] * 6)
        self.assertEqual(
            gen_data["hydro_budgeted"]["values"][:6],
            [6.35, 6.85, 12.75, 13.2, 12.75, 13.7],
        )

    ### GENERATOR SETS ###

    def test_hydro_sets_advanced(self):
        # With hydro_advanced on, the existing hydro unit is modeled as
        # a hydro generator. Candidate hydro units are not supported, so
        # the candidate is removed from all generator sets and a warning
        # is logged. THe log is then checked to make sure the line was
        # printed to the log
        with self.assertLogs("gtep.model_library.components", level="WARNING") as log:
            self.create_testing_obj(config={"hydro_advanced": True})

        self.assertTrue(any("does not support candidates" in msg for msg in log.output))

        self.assertEqual(set(self.m.hydroGenerators), {"1_HYDRO"})
        for gen in ["1_HYDRO", "1_HYDRO-c"]:
            self.assertNotIn(gen, self.m.renewableGenerators)
            self.assertNotIn(gen, self.m.thermalGenerators)

        self.assertEqual(set(self.m.hydroGeneratorsByBus["bus1"]), {"1_HYDRO"})

    def test_hydro_sets_not_advanced(self):
        # With hydro_advanced off, both hydro units are ordinary
        # renewable generators
        self.create_testing_obj(config={})

        self.assertFalse(hasattr(self.m, "hydroGenerators"))
        self.assertIn("1_HYDRO", self.m.renewableGenerators)
        self.assertIn("1_HYDRO-c", self.m.renewableGenerators)

    def test_hydro_budgeted_without_advanced(self):
        # hydro_budgeted is only used together with hydro_advanced
        self.create_testing_obj(config={"hydro_budgeted": True})

        self.assertFalse(hasattr(self.m, "hydroGenerators"))
        self.assertFalse(hasattr(self.b_comm[0], "hydroMinimumExpected"))
        self.assertFalse(hasattr(self.b_rep, "hydroSpill"))

    ### HYDRO PARAMETERS ###

    def test_hydro_params_budgeted(self):
        # With hydro_advanced and hydro_budgeted on, the minimum output
        # and inflow in each commitment period come from the hydro
        # minimum and budget files, and the maximum output comes from
        # p_max
        self.create_testing_obj(config={"hydro_advanced": True, "hydro_budgeted": True})
        gen_data = self.m.data_list[0].data["elements"]["generator"]["1_HYDRO"]

        # 2020-04-23, Periods 1-6
        expected_min = [0.0] * 6
        expected_budget = [6.35, 6.85, 12.75, 13.2, 12.75, 13.7]

        for c, b_comm in enumerate(self.b_comm):
            with self.subTest(commitment_period=c + 1):
                self.assertAlmostEqual(
                    pyo.value(b_comm.hydroMinimumExpected["1_HYDRO"]), expected_min[c]
                )
                self.assertAlmostEqual(
                    pyo.value(b_comm.hydroInflow["1_HYDRO"]), expected_budget[c]
                )
                self.assertAlmostEqual(
                    pyo.value(b_comm.hydroCapacityExpected["1_HYDRO"]),
                    gen_data["p_max"]["values"][c],
                )

                # Generation is bounded by the minimum and maximum
                b_disp = b_comm.dispatchPeriod[b_comm.dispatchPeriods.first()]
                lb, ub = b_disp.hydroGeneration["1_HYDRO"].bounds
                self.assertAlmostEqual(lb, expected_min[c])
                self.assertAlmostEqual(ub, gen_data["p_max"]["values"][c])

        # The energy budget over the representative day is the sum of the
        # inflows over the six 1-hour commitment periods
        self.assertTrue(hasattr(self.b_rep, "hydroSpill"))
        self.assertTrue(hasattr(self.b_rep, "hydro_energy_budget"))
        self.assertAlmostEqual(
            pyo.value(self.b_rep.hydroEnergyBudget["1_HYDRO"]), sum(expected_budget)
        )
        self.assertEqual(set(self.b_rep.hydro_energy_budget), {"1_HYDRO"})

    def test_hydro_params_not_budgeted(self):
        # With hydro_advanced on and hydro_budgeted off, the minimum
        # output and inflow equal the maximum (p_max), so each unit runs
        # at p_max, and there are no energy-budget components
        self.create_testing_obj(config={"hydro_advanced": True})
        gen_data = self.m.data_list[0].data["elements"]["generator"]["1_HYDRO"]

        for c, b_comm in enumerate(self.b_comm):
            with self.subTest(commitment_period=c + 1):
                p_max = gen_data["p_max"]["values"][c]
                self.assertAlmostEqual(
                    pyo.value(b_comm.hydroCapacityExpected["1_HYDRO"]), p_max
                )
                self.assertAlmostEqual(
                    pyo.value(b_comm.hydroMinimumExpected["1_HYDRO"]), p_max
                )
                self.assertAlmostEqual(pyo.value(b_comm.hydroInflow["1_HYDRO"]), p_max)

        self.assertFalse(hasattr(self.b_rep, "hydroSpill"))
        self.assertFalse(hasattr(self.b_rep, "hydro_energy_budget"))

    def test_hydro_min_exceeds_budget_error(self):
        # If the minimum output needs more energy than the inflow budget
        # over a representative period, building the model raises an
        # error. A minimum above p_max is reduced to p_max, and p_max is
        # larger than the budget on 2020-04-23.
        def raise_min(generators):
            series = generators["1_HYDRO"]["hydro_p_min"]
            series["values"] = [100.0] * len(series["values"])

        with self.assertRaisesRegex(ValueError, "minimum output requires"):
            self.create_testing_obj(
                config={"hydro_advanced": True, "hydro_budgeted": True},
                modify=raise_min,
            )

    def test_hydro_units(self):
        # The hydro parameters are in MW and the energy budget in MWh
        self.create_testing_obj(config={"hydro_advanced": True, "hydro_budgeted": True})

        b_comm = self.b_comm[0]
        assert_units_equivalent(b_comm.hydroMinimumExpected["1_HYDRO"], u.MW)
        assert_units_equivalent(b_comm.hydroCapacityExpected["1_HYDRO"], u.MW)
        assert_units_equivalent(b_comm.hydroInflow["1_HYDRO"], u.MW)
        assert_units_equivalent(self.b_rep.hydroSpill["1_HYDRO"], u.MW * u.hr)
        assert_units_equivalent(self.b_rep.hydroEnergyBudget["1_HYDRO"], u.MW * u.hr)
