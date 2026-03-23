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

"""Generation and Transmission Expansion Planning (GTEP) Model

Model equations available at ref[1]


References:

[1] http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

"""

__author__ = "Kyle Skolfield"

import json
import numpy as np
import re

import pyomo.environ as pyo
from pyomo.environ import units as u
from pyomo.common.timing import TicTocTimer
from pyomo.repn.linear import LinearRepnVisitor

from egret.data.model_data import ModelData
from egret.model_library.transmission.tx_utils import scale_ModelData_to_pu

from gtep.config_options import (
    _get_model_config,
    _add_common_configs,
    _add_investment_configs,
)
import gtep.model_library.investment as inv
import gtep.model_library.dispatch as disp
import gtep.model_library.commitment as commit
import gtep.model_library.objective as obj_comp
import gtep.model_library.components as comps
import gtep.model_library.representative_period as rep_period
import gtep.model_library.scaling as scaling
import gtep.model_library.gen as gens
import gtep.model_library.storage as stor
import gtep.model_library.transmission as transm

# Define what a USD is for pyomo units purposes. This will be set to a
# base year and we will do NPV calculations based on automatic Pyomo
# unit transformations.
u.load_definitions_from_strings(["USD = [currency]", "MVAR = [power]"])


## TODO: Egret features


def data_update(investment_stage, storage_object, target_storage_object):
    pass


# This is only used for reporting potentially bad (i.e., large
# magnitude) coefficients and thus only when that argument is passed.
class VisitorConfig(object):
    def __init__(self):
        self.subexpr = {}
        self.var_map = {}
        self.var_order = {}

    def __iter__(self):
        return iter((self.subexpr, self.var_map, self.var_order))


class ExpansionPlanningModel:
    """A generalized generation and transmission expansion planning
    model.

    """

    def __init__(
        self,
        config=None,
        stages=1,
        formulation=None,
        data=None,
        cost_data=None,
        num_reps=3,
        len_reps=24,
        num_commit=24,
        num_dispatch=4,
        duration_dispatch=15,
    ):
        """Initialize generation & expansion planning model object.

        :param stages: integer number of investment periods
        :param formulation: Egret stuff, to be filled
        :param data: full set of model data
        :param cost_data: full set of cost data for all generators
        :param num_reps: integer number of representative periods per investment period
        :param len_reps: (for now integer) length of each representative period (in hours)
        :param num_commit: integer number of commitment periods per representative period
        :param num_dispatch: integer number of dispatch periods per commitment period
        :param duration_dispatch: (for now integer) duration of each dispatch period (in minutes)
        :return: Pyomo model for full GTEP
        """

        self.stages = stages
        self.formulation = formulation
        self.data = data
        self.cost_data = cost_data
        self.num_reps = num_reps
        self.len_reps = len_reps
        self.num_commit = num_commit
        self.num_dispatch = num_dispatch
        self.duration_dispatch = duration_dispatch
        self.config = _get_model_config()
        self.timer = TicTocTimer()

        _add_common_configs(self.config)
        _add_investment_configs(self.config)

    def create_model(self):
        """Create concrete Pyomo model object associated with the
        ExpansionPlanningModel class

        """

        self.timer.tic("Creating GTEP Model")
        m = pyo.ConcreteModel("GTEP Model")

        m.config = self.config

        # NOTE: scale_ModelData_to_pu doesn't account for expansion
        # data -- does it need to?
        # [TODO: checks for active/built/inactive/unbuilt/etc. gen]
        if self.data is None:
            raise
        elif type(self.data.representative_data) is list:
            # If self.data is a list, it is a list of data for
            # representative periods
            m.data_list = self.data.representative_data
            m.md = m.data_list[0]
            m.data = self.data
        else:
            # If self.data is an Egret model data object,
            # representative periods will just copy it unchanged
            m.data_list = None
            m.md = scale_ModelData_to_pu(self.data)
            m.formulation = self.formulation

        # Add cost_data from the DataProcessing class. [TODO: Think
        # about how to do some scaling in this data.]
        m.mc = self.cost_data

        comps.model_set_declaration(
            m, self.stages, rep_per=[i for i in range(1, self.num_reps + 1)]
        )

        comps.model_data_references(
            m, self.num_commit, self.num_dispatch, self.duration_dispatch
        )

        # model_create_investment_stages(m, self.stages)
        create_stages(m, self.stages)

        obj_comp.create_objective_function(m)

        self.model = m

    # [TODO: This method should handle string or i/o object for
    # outfile.]
    def report_model(self, outfile="pretty_model_output.txt"):
        """Pretty prints Pyomo model to outfile.

        :outfile: (str, optional) _description_. Defaults to "pretty_model_output.txt".
        """
        with open(outfile, "w") as outf:
            self.model.pprint(ostream=outf)

    def report_large_coefficients(self, outfile, magnitude_cutoff=1e5):
        """Dump very large magnitude (>= 1e5) coefficients to a json
        file.

        :outfile: should accept filename or open file and write there;
                  (see how we do this in pyomo elsewhere)
        :magnitude_cutoff: magnitude above which to report coefficients

        """
        var_coef_dict = {}
        for e in self.model.component_data_objects(pyo.Constraint):
            cfg = VisitorConfig()
            repn = LinearRepnVisitor(*cfg).walk_expression(e.body)
            repn_dict = repn.linear
            varname_dict = {cfg.var_map[v].name: repn_dict[v] for v in repn_dict.keys()}
            var_coef_dict = dict(var_coef_dict | varname_dict)

        really_bad_var_coef_dict = {
            key: value
            for (key, value) in var_coef_dict.items()
            if abs(value) >= magnitude_cutoff
        }
        really_bad_var_coef_list = sorted(
            really_bad_var_coef_dict.items(), key=lambda x: x[1]
        )
        with open(outfile, "w") as fil:
            json.dump(really_bad_var_coef_list, fil)


####################################
## Model Building Functions Below ##
####################################


def commitment_period_rule(b, commitment_period):
    """Create commitment period block.

    :param b: commitment period block
    :param commitment_period: corresponding commitment period label
    """
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    b.commitmentPeriod = commitment_period
    b.commitmentPeriodLength = pyo.Param(
        within=pyo.PositiveReals, default=1, units=u.hr
    )
    b.dispatchPeriods = pyo.RangeSet(m.numDispatchPeriods[r_p.currentPeriod])
    b.carbonTax = pyo.Param(default=0)
    b.dispatchPeriod = pyo.Block(b.dispatchPeriods)

    # update properties for this time period!!!!
    if m.data_list:
        m.md = m.data_list[i_p.representativePeriods.index(r_p.currentPeriod)]

    # Making an exception for cases where gens were candidates
    # bc their time series reduced to single values. Will probably need to fix
    # this and look at where that reduction is taking place because we need more
    # than a single value if the generator is built. (Probably? Maybe there's a
    # different way to handle candidate renewable data because this assumes
    # knowledge of the future outputs of a candidate... could be captured by scenarios?)
    # Maximum output of each renewable generator

    # [ESR WIP: Corrected to be in the block "b", not in "m". Also,
    # changed its original name "renewableCapacity" to include the
    # word "Expected" since there are two different
    # "renewableCapacity" parameters (the second one is included in
    # model_data_references and includes the word "Nameplate"). The
    # model now distingues between these two.]
    b.renewableCapacityExpected = {}
    units_renewable_capacity = u.MW
    for renewableGen in m.renewableGenerators:
        if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float:
            b.renewableCapacityExpected[renewableGen] = 0 * units_renewable_capacity
        else:
            b.renewableCapacityExpected[renewableGen] = (
                m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"][
                    commitment_period - 1
                ]
                * units_renewable_capacity
            )

    ## TODO: Redesign load scaling and allow nature of it as argument
    scaling.add_load_scaling(m, b, commitment_period, i_p.investmentStage)

    ## TODO: This feels REALLY inelegant and bad.
    ## TODO: Something weird happens if I say periodLength has a unit
    for period in b.dispatchPeriods:
        b.dispatchPeriod[period].periodLength = pyo.Param(
            within=pyo.PositiveReals, default=1
        )
        disp.add_dispatch_variables(b.dispatchPeriod[period], period)

    ## TODO: if commitment is neglected but dispatch is still desired, pull something different here? or simply don't enforce linked commitment constraints?
    if m.config["include_commitment"]:
        commit.add_commitment_disjuncts(b, commitment_period)

    commit.add_commitment_constraints(b, commitment_period)

    for period in b.dispatchPeriods:
        disp.add_dispatch_constraints(b.dispatchPeriod[period], period)


def representative_period_rule(b, representative_period):
    """Create representative period block.

    :b: Representative period block
    :representative_period: corresponding representative period label
    """

    m = b.model()
    i_s = b.parent_block()

    b.representative_date = m.data.representative_dates[representative_period - 1]
    broken_date = list(re.split(r"[-: ]", b.representative_date))
    b.month = int(broken_date[1])
    b.day = int(broken_date[2])
    # b.load_scaling = i_s.load_scaling[
    #    (i_s.load_scaling["month"] == b.month) & (i_s.load_scaling["day"] == b.day)
    # ]

    b.currentPeriod = representative_period

    if m.config["include_commitment"] or m.config["include_redispatch"]:
        b.commitmentPeriods = pyo.RangeSet(
            m.numCommitmentPeriods[representative_period]
        )
        b.commitmentPeriod = pyo.Block(b.commitmentPeriods, rule=commitment_period_rule)

        rep_period.add_representative_period_variables(b, representative_period)
        rep_period.add_representative_period_constraints(b, representative_period)


def investment_stage_rule(b, investment_stage):
    """Creates investment stage block.

    :b: Investment block
    :investment_stage: ID for current investment stage
    """

    m = b.parent_block()

    b.year = m.years[investment_stage - 1]

    print(f"b.year = {b.year}")

    # Declare costs parameters here since they depend on the
    # investment year
    comps.model_data_costs(m, b.year)

    inv.add_investment_params_and_variables(b, investment_stage)

    b.representativePeriods = [
        p
        for p in m.representativePeriods
        # if m.representativePeriodStage[p] == investment_stage
    ]

    b.representativePeriod = pyo.Block(
        b.representativePeriods, rule=representative_period_rule
    )
    inv.add_investment_constraints(b, investment_stage)


def model_create_investment_stages(m, stages):
    """Creates investment blocks and linking constraints for GTEP model.
    Largely manages retirements and links operational units in a given investment stage
    to operational + installed - retired in the previous investment stage.

    :m: Pyomo model object
    :stages: Number of investment stages in planning horizon
    """

    m.investmentStage = pyo.Block(m.stages, rule=investment_stage_rule)

    # Add logical constraints for generators and transmission lines
    # and storage, when needed. These logical constraints ensure that
    # their status are operationally consistent over time
    gens.add_generators_logical_constraints(m)

    if m.config["storage"]:
        stor.add_storage_logical_constraints(m)

    if m.config["transmission"]:
        transm.add_transmission_logical_constraints(m)


def create_stages(m, stages):
    """Creates investment, commitment, and dispatch stages using Blocks

    :m: Pyomo model object
    :stages: Number of investment stages in planning horizon

    """

    # Add investment stage Block and all its equations and
    # variables. this block is the main block where all the stages
    # will be included in a nested way.
    m.investmentStage = pyo.Block(m.stages)

    for investment_stage in m.stages:
        b_inv = m.investmentStage[investment_stage]
        b_inv.year = m.years[investment_stage - 1]
        print(f"{b_inv}.year = {b_inv.year}")

        # Declare costs parameters here since they depend on the
        # investment year
        comps.model_data_costs(m, b_inv.year)

        # Declare investment parameters and variables. This includes the
        # status disjunction for generators and transmission lines and
        # storage, when needed
        inv.add_investment_params_and_variables(b_inv, investment_stage)

        b_inv.representativePeriods = [
            p for p in m.representativePeriods
        ]
        b_inv.representativePeriod = pyo.Block(b_inv.representativePeriods)

        #--------------------------------------------------------------
        # Add representative_period Block and all its variables and
        # equations.
        for representative_period in b_inv.representativePeriods:
            b_rep = b_inv.representativePeriod[representative_period]
            b_rep.representative_date = m.data.representative_dates[
                representative_period - 1
            ]

            # [ESR WIP: Comment out for now since it is not
            # used. Check if we need this in future versions.]
            # broken_date = list(re.split(r"[-: ]", b_rep[per].representative_date))
            # b_rep.month = int(broken_date[1])
            # b_rep.day = int(broken_date[2])
            b_rep.currentPeriod = representative_period

            if m.config["include_commitment"] or m.config["include_redispatch"]:
                b_rep.commitmentPeriods = pyo.RangeSet(
                    m.numCommitmentPeriods[representative_period]
                )
                b_rep.commitmentPeriod = pyo.Block(b_rep.commitmentPeriods)

                # --.--.--.--.--.--.----.--.--.--.--.--.----.--.--.--.--.--.--
                # Add commitment Block and all its equations and
                # constraints
                for commitment_period in b_rep.commitmentPeriods:
                    b_comm = b_rep.commitmentPeriod[commitment_period]
                    b_comm.commitmentPeriod = commitment_period
                    b_comm.dispatchPeriods = pyo.RangeSet(
                        m.numDispatchPeriods[b_rep.currentPeriod]
                    )
                    b_comm.dispatchPeriod = pyo.Block(b_comm.dispatchPeriods)

                    # [TODO: update properties for this time period!]
                    if m.data_list:
                        m.md = m.data_list[
                            b_inv.representativePeriods.index(
                                b_rep.currentPeriod
                            )
                        ]

                    commit.add_params(
                        m, b_comm, commitment_period, b_inv.investmentStage,
                    )

                    #=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
                    # Add dispatch equations
                    
                    # [TODO: This feels REALLY inelegant and
                    # bad. Check a better way of declaring these.]
                    for period in b_comm.dispatchPeriods:
                        b_comm.dispatchPeriod[period].periodLength = pyo.Param(
                            initialize=1,
                            within=pyo.PositiveReals,
                            units=u.minutes,
                        )
                        disp.add_dispatch_variables(
                            b_comm.dispatchPeriod[period], period
                        )
                        disp.add_dispatch_constraints(
                            b_comm.dispatchPeriod[period], period
                        )

                    #=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.

                    # [TODO: If commitment is neglected but dispatch
                    # is still desired, pull something different here?
                    # or simply don't enforce linked commitment
                    # constraints?]
                    if m.config["include_commitment"]:
                        commit.add_commitment_disjuncts(b_comm, commitment_period)

                    commit.add_commitment_constraints(b_comm, commitment_period)

                # --.--.--.--.--.--.----.--.--.--.--.--.----.--.--.--.--.--.--

                rep_period.add_representative_period_variables(
                    b_rep, representative_period
                )
                rep_period.add_representative_period_constraints(
                    b_rep, representative_period
                )
            #--------------------------------------------------------------

    for investment_stage in m.stages:
        inv.add_investment_constraints(m.investmentStage[investment_stage], investment_stage)
    #############

    # Add logical constraints for generators and transmission lines
    # and storage, when needed. These logical constraints ensure that
    # their status are operationally consistent over time
    gens.add_generators_logical_constraints(m)

    if m.config["storage"]:
        stor.add_storage_logical_constraints(m)

    if m.config["transmission"]:
        transm.add_transmission_logical_constraints(m)
