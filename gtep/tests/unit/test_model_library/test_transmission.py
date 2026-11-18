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

"""Unit Tests for the transmission sets, disjuncts, and constraints in
the Generation and Transmission Expansion Planning (GTEP) model.

These tests use the 123-bus case (123_Bus_Resil_Week), which has 221
existing lines and 34 candidate lines (UID ending in "-c"), and check
that the transmission sets are built correctly under different
configuration options.
"""

import math
from pathlib import Path

import pyomo.common.unittest as unittest
import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables

from gtep.gtep_model import ExpansionPlanningModel
from gtep.model_library.transmission import _branch_dc_data
from gtep.tests.unit.utils_for_testing import create_data


class TestTransmission(unittest.TestCase):

    def create_testing_obj(self, config, modify=None):
        """Creates a 123-bus model for a given set of config options.

        :param config: Non-default model config options
        :param modify: Optional function that takes a branch data dict
                       and modifies it before the model is built. The
                       model reads branch data from the
                       representative-period copies, so it is applied
                       to md and to every copy.
        """
        model_config = {
            "include_investment": True,
            "include_commitment": False,
            "include_redispatch": True,
            "transmission": True,
            "flow_model": "DC",
        }
        model_config.update(config)

        data_path = Path(__file__).resolve().parents[3] / "data" / "123_Bus_Resil_Week"
        data_object, _ = create_data(
            input_data_path=data_path,
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_representative_period": 1,
            },
            prescient_data_args={"representative_dates": ["2019-01-28 00:00"]},
            include_cost_data=False,
        )

        # Individual tests can modify the data for testing
        if modify is not None:
            for md in [data_object.md, *data_object.representative_data]:
                modify(md.data["elements"]["branch"])

        mod_object = ExpansionPlanningModel(
            config=model_config,
            data=data_object,
            cost_data=None,
        )
        mod_object.create_model()
        self.m = mod_object.model

        # Investment stage block and its first dispatch block
        self.b = self.m.investmentStage[self.m.stages.first()]
        b_rep = self.b.representativePeriod[self.b.representativePeriods[0]]
        b_comm = b_rep.commitmentPeriod[b_rep.commitmentPeriods.first()]
        self.b_disp = b_comm.dispatchPeriod[b_comm.dispatchPeriods.first()]

        # Existing and candidate lines in the input data
        branches = self.m.md.data["elements"]["branch"]
        self.all_lines = set(branches)
        self.candidate_lines = {b for b in branches if b.endswith("-c")}
        self.existing_lines = self.all_lines - self.candidate_lines

    ### TRANSMISSION SETS ###

    def test_sets_transmission_on(self):
        # With transmission on, all lines are in the model. The 123-bus
        # branch.csv has no "Static Line" column, so no line is static
        # and every line can change status.
        self.create_testing_obj(config={})

        self.assertEqual(len(self.existing_lines), 221)
        self.assertEqual(len(self.candidate_lines), 34)

        self.assertEqual(set(self.m.lines), self.all_lines)
        self.assertEqual(set(self.m.candidateLines), self.candidate_lines)
        self.assertEqual(len(self.m.staticLines), 0)
        self.assertEqual(set(self.m.nonStaticLines), self.all_lines)

        # Candidate lines are connected to their buses
        for branch in self.candidate_lines:
            self.assertIn(branch, self.m.branchByFromBus[self.m.from_bus[branch]])
            self.assertIn(branch, self.m.branchByToBus[self.m.to_bus[branch]])

    def test_sets_transmission_off(self):
        # With transmission off, candidate lines are left out and all
        # existing in-service lines are static
        self.create_testing_obj(config={"transmission": False})

        self.assertEqual(set(self.m.lines), self.existing_lines)
        self.assertEqual(set(self.m.staticLines), self.existing_lines)
        self.assertEqual(len(self.m.candidateLines), 0)
        self.assertEqual(len(self.m.nonStaticLines), 0)

        for branch in self.m.lines:
            self.assertFalse(branch.endswith("-c"))

    def test_sets_transmission_off_out_of_service_line(self):
        # With transmission off, an existing line that is not in service
        # is left out of the model
        def take_out_of_service(branches):
            branches["1"]["in_service"] = False

        self.create_testing_obj(
            config={"transmission": False}, modify=take_out_of_service
        )

        self.assertNotIn("1", self.m.lines)
        self.assertEqual(set(self.m.lines), self.existing_lines - {"1"})

    def test_sets_static_lines(self):
        # Lines marked static are in staticLines and not in
        # nonStaticLines. A candidate line marked static is still
        # non-static.
        def mark_static(branches):
            for branch in ["1", "2", "3", "222-c"]:
                branches[branch]["static"] = True

        self.create_testing_obj(config={}, modify=mark_static)

        self.assertEqual(set(self.m.staticLines), {"1", "2", "3"})
        self.assertEqual(set(self.m.nonStaticLines), self.all_lines - {"1", "2", "3"})
        self.assertIn("222-c", self.m.nonStaticLines)
        self.assertEqual(set(self.m.lines), self.all_lines)

    ### INVESTMENT DISJUNCTS ###

    def test_investment_disjuncts_investment_on(self):
        # Investment status disjuncts are built for non-static lines.
        # In the first stage, existing lines are fixed operational and
        # candidate lines are fixed not operational. Candidate lines are
        # not disabled when investment is on.
        self.create_testing_obj(config={})

        self.assertEqual(set(self.b.branchInvestStatus), self.all_lines)

        for branch in self.existing_lines:
            ind = self.b.branchOperational[branch].indicator_var
            self.assertTrue(ind.fixed)
            self.assertTrue(ind.value)

        for branch in self.candidate_lines:
            ind = self.b.branchOperational[branch].indicator_var
            self.assertTrue(ind.fixed)
            self.assertFalse(ind.value)
            self.assertFalse(self.b.branchDisabled[branch].indicator_var.fixed)

        self.assertTrue(hasattr(self.b, "transmission_investment_cost"))

    def test_investment_disjuncts_investment_off(self):
        # With investment off, candidate lines are fixed disabled
        self.create_testing_obj(config={"include_investment": False})

        self.assertEqual(set(self.m.candidateLines), self.candidate_lines)

        for branch in self.candidate_lines:
            ind = self.b.branchDisabled[branch].indicator_var
            self.assertTrue(ind.fixed)
            self.assertTrue(ind.value)

        for branch in self.existing_lines:
            self.assertFalse(self.b.branchDisabled[branch].indicator_var.fixed)

    def test_investment_disjuncts_static_lines(self):
        # Static lines get no investment status disjuncts
        def mark_static(branches):
            for branch in ["1", "2", "3"]:
                branches[branch]["static"] = True

        self.create_testing_obj(config={}, modify=mark_static)

        self.assertEqual(
            set(self.b.branchInvestStatus), self.all_lines - {"1", "2", "3"}
        )
        for branch in {"1", "2", "3"}:
            self.assertNotIn(branch, self.b.branchOperational)

    def test_investment_disjuncts_transmission_off(self):
        # With transmission off, no transmission investment components
        # are built
        self.create_testing_obj(config={"transmission": False})

        self.assertFalse(hasattr(self.b, "branchInvestStatus"))
        self.assertFalse(hasattr(self.b, "transmission_investment_cost"))
        self.assertFalse(hasattr(self.m, "consistent_branch_operation"))

    ### DISPATCH (OPERATING STATE) DISJUNCTS AND FLOW CONSTRAINTS ###

    def test_dispatch_components_transmission_on(self):
        # Non-static lines get in-use/not-in-use disjuncts. Bus angles
        # are defined once on the dispatch block and shared by all
        # lines (not created inside each disjunct).
        self.create_testing_obj(config={})

        self.assertEqual(set(self.b_disp.powerFlow), self.all_lines)
        self.assertEqual(set(self.b_disp.branchInUseStatus), self.all_lines)
        self.assertEqual(len(self.b_disp.static_dc_power_flow), 0)

        self.assertEqual(set(self.b_disp.busAngle), set(self.m.buses))
        self.assertEqual(set(self.b_disp.deltaBusAngle), self.all_lines)

        branch = sorted(self.candidate_lines)[0]
        disj = self.b_disp.branchInUse[branch]

        fb, tb, _, _ = _branch_dc_data(self.m, branch)
        dc_vars = ComponentSet(identify_variables(disj.dc_power_flow.body))
        self.assertIn(self.b_disp.busAngle[fb], dc_vars)
        self.assertIn(self.b_disp.busAngle[tb], dc_vars)

        self.assertEqual(set(self.b_disp.must_use_active_branches), self.all_lines)
        self.assertEqual(set(self.b_disp.cannot_use_inactive_branches), self.all_lines)

    def test_dispatch_components_transmission_off(self):
        # With transmission off, all lines are static: they get DC flow
        # and angle-difference constraints directly on the dispatch
        # block and no in-use/not-in-use disjuncts
        self.create_testing_obj(config={"transmission": False})

        self.assertEqual(set(self.b_disp.powerFlow), self.existing_lines)
        self.assertEqual(len(self.b_disp.branchInUseStatus), 0)
        self.assertEqual(set(self.b_disp.static_dc_power_flow), self.existing_lines)
        self.assertEqual(
            set(self.b_disp.static_delta_bus_angle_bounds), self.existing_lines
        )
        self.assertFalse(hasattr(self.b_disp, "must_use_active_branches"))

        branch = sorted(self.existing_lines)[0]
        c = self.b_disp.static_delta_bus_angle_bounds[branch]
        self.assertAlmostEqual(pyo.value(c.upper), math.pi / 6)
        self.assertAlmostEqual(pyo.value(c.lower), -math.pi / 6)

    def test_dispatch_components_static_lines(self):
        # With transmission on, static lines get flow constraints on the
        # dispatch block and other lines get disjuncts
        def mark_static(branches):
            for branch in ["1", "2", "3"]:
                branches[branch]["static"] = True

        self.create_testing_obj(config={}, modify=mark_static)

        self.assertEqual(set(self.b_disp.static_dc_power_flow), {"1", "2", "3"})
        self.assertEqual(
            set(self.b_disp.branchInUseStatus),
            self.all_lines - {"1", "2", "3"},
        )

    def test_transport_flow_model(self):
        # The transport model has no bus angles or DC flow constraints
        for transmission in (True, False):
            with self.subTest(transmission=transmission):
                self.create_testing_obj(
                    config={"transmission": transmission, "flow_model": "transport"}
                )

                self.assertTrue(hasattr(self.b_disp, "powerFlow"))
                self.assertFalse(hasattr(self.b_disp, "busAngle"))
                self.assertFalse(hasattr(self.b_disp, "static_dc_power_flow"))

    def test_copper_plate_flow_model(self):
        # The copper-plate model has no line flows or line disjuncts
        self.create_testing_obj(config={"flow_model": "CP"})

        self.assertFalse(hasattr(self.b_disp, "powerFlow"))
        self.assertFalse(hasattr(self.b_disp, "branchInUseStatus"))
        self.assertTrue(hasattr(self.b_disp, "CP_flow_balance"))