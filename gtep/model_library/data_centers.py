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

"""Constraints for Data Centers in the Generation and Transmission
Expansion Planning (GTEP) Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u

def add_data_center_parameters(m):
    """This method defines parameters related to data centers in the main
    model.

    """

    m.dataCenterOwners = pyo.Set(initialize=(m.md.data["elements"]["data_center"][dc]["owner"] for dc in m.md.data["elements"]["data_center"]), doc="Data center owners")

    m.dataCenterCapacity = pyo.Param(
        m.dataCenters, default=0, units=u.MW, doc="Data center capacity"
    )

    m.dataCenterGenerationCapacity = pyo.Param(
        m.dataCenters, default=0, units=u.MW, doc="Data center co-located generation capacity"
    )

    m.dataCenterOperationalCost = pyo.Param(
        m.dataCenters, default=0, units=u.USD / u.MWh, doc="Data center operational cost"
    )

    m.dataCenterGenerationCost = pyo.Param(
        m.dataCenters,
        default=0,
        units=u.USD / u.MWh,
        doc="Cost of generation from data center co-located generation",
    )

    m.dataCenterInvestmentCost = pyo.Param(
        m.dataCenters,
        default=0,
        units=u.USD / u.MW,
        doc="Capital cost for investing in new data center capacity or upgrading existing capacity",
    )


def add_data_center_status_disjuncts(b, data_centers_set):
    """This method declares variables and a Disjunction and its
    disjuncts to model the selection of data center status. The alternatives are:

    dataCenterOperational: Data center is active and consuming power.
    dataCenterInstalled:   Data center is newly added and active.
    dataCenterRetired:     Data center is removed from service.
    dataCenterDisabled:    Data center is temporarily out of service.
    dataCenterExtended:    Data center is upgraded beyond its original capacity.

    """

    @b.Disjunct(data_centers_set)
    def dataCenterOperational(disj, dc):
        return

    @b.Disjunct(data_centers_set)
    def dataCenterInstalled(disj, dc):
        return

    @b.Disjunct(data_centers_set)
    def dataCenterRetired(disj, dc):
        return

    @b.Disjunct(data_centers_set)
    def dataCenterDisabled(disj, dc):
        return

    @b.Disjunct(data_centers_set)
    def dataCenterExtended(disj, dc):
        return

    @b.Disjunction(data_centers_set)
    def dataCenterInvestStatus(disj, dc):
        return [
            disj.dataCenterOperational[dc],
            disj.dataCenterInstalled[dc],
            disj.dataCenterRetired[dc],
            disj.dataCenterDisabled[dc],
            disj.dataCenterExtended[dc],
        ]


def add_investment_data_centers_constraints(m, b, investment_stage):

    # These constraints take the "in_service" data and fix the
    # indicator variables of the status disjuncts to define the
    # operation of data centers.
    for dc in m.dataCenters:
        if (
            m.md.data["elements"]["data_center"][dc]["in_service"] == False
            and investment_stage == 1
        ):
            b.dataCenterOperational[dc].indicator_var.fix(False)
        elif (
            m.md.data["elements"]["data_center"][dc]["in_service"] == True
            and investment_stage == 1
        ):
            b.dataCenterOperational[dc].indicator_var.fix(True)

    for dc in m.dataCenters:
        if (
            m.md.data["elements"]["data_center"][dc]["in_service"] == False
            and investment_stage == 1
        ):
            b.dataCenterOperational[dc].fix(0)
        elif (
            m.md.data["elements"]["data_center"][dc]["in_service"] == True
            and investment_stage == 1
        ):
            b.dataCenterOperational[dc].fix(m.dataCenterCapacity[dc])

    @b.Expression(doc="Data center investment costs in $")
    def data_centers_investment_cost(b):
        return (
            sum(
                m.dataCenterInvestmentCost[dc]
                * m.dataCenterCapacity[dc]  # in MW
                * m.capitalMultiplier.get(dc, 1)  # Assuming similar multiplier
                * b.dataCenterInstalled[dc].indicator_var.get_associated_binary()
                for dc in m.dataCenters
            )
            + sum(
                m.dataCenterInvestmentCost[dc]
                * m.extensionMultiplier.get(dc, 1)
                * m.dataCenterCapacity[dc]
                * b.dataCenterExtended[dc].indicator_var.get_associated_binary()
                for dc in m.dataCenters
            )
            + sum(
                m.dataCenterInvestmentCost[dc]
                * m.retirementMultiplier.get(dc, 1)
                * b.dataCenterRetired[dc].indicator_var.get_associated_binary()
                for dc in m.dataCenters
            )
        )

def add_representative_period_data_centers_constraints(m, b, rep_per, i_p, commitment_period):

    @b.Constraint(m.dataCenters, doc = "Data center something.")
    def data_center_something_constraint(b, dc):
        return 



def add_data_centers_state_disjuncts(m, b, r_p, i_p, commitment_period):
    """This method defines a Disjunction with disjuncts representing
    the alternatives for data center state operation. The alternatives
    are:

    dataCenterTraining:      Data center is operating and running training tasks.
    dataCenterInference:       Data center is operating and running infeence tasks.
    dataCenterOff:     Data center is offline and not consuming power.

    """
    @b.Disjunct(m.dataCenters)
    def dataCenterTraining(disj, dc):
        b = disj.parent_block()

        @disj.Constraint(doc="Data center requires maximum load while training.")
        def data_center_training_load_constraint(d, dispatchPeriod):
            return b.dataCenterLoad[dc] == m.dataCenterCapacity[dc]

    
    @b.Disjunct(m.dataCenters)
    def dataCenterInference(disj, dc):
        b = disj.parent_block()

        @disj.Constraint(doc = "Data center load is between 50\% and 100\% of capacity while running inference.")
        def data_center_inference_load_constraint(d, dispatchPeriod):
            return b.dataCenterLoad[dc] >= 0.5 * m.dataCenterCapacity[dc]
    
    @b.Disjunct(m.dataCenters)
    def dataCenterOff(disj, dc):
        b = disj.parent_block()

        @disj.Constraint(doc = "Data center load is zero when offline.")
        def data_center_off_load_constraint(d, dispatchPeriod):
            return b.dataCenterLoad[dc] == 0
    
    @b.Disjunction(m.dataCenters)
    def dataCenterOperationalStatus(disj, dc):
        return [
            disj.dataCenterTraining[dc],
            disj.dataCenterInference[dc],
            disj.dataCenterOff[dc],
        ]


def add_commitment_data_centers_constraints(m, b, r_p, i_p, comm_per):
    pass


def add_dispatch_data_centers_variables(m, b):
    """Add data center variables to the dispatch block."""

    def data_center_load_limits(b, dc, doc="Bounds on data center load"):
        return (0, m.dataCenterCapacity[dc])

    b.dataCenterLoad = pyo.Var(
        m.dataCenters,
        domain=pyo.NonNegativeReals,
        bounds=data_center_load_limits,
        initialize=0,
        units=u.MW,
        doc="Data center load consumption",
    )

    def data_center_generation_limits(b, dc, doc="Bounds on data center generation"):
        return (0, m.dataCenterGenerationCapacity[dc])

    b.dataCenterGeneration = pyo.Var(
        m.dataCenters,
        domain=pyo.NonNegativeReals,
        bounds=data_center_generation_limits,
        initialize=0,
        units=u.MW,
        doc="Generation from data center co-located generation",
    )

    def data_center_curtailment_limits(b, dc, doc="Bounds on data center curtailment"):
        return (0, m.dataCenterCapacity[dc])

    b.dataCenterCurtailment = pyo.Var(
        m.dataCenters,
        domain=pyo.NonNegativeReals,
        bounds=data_center_curtailment_limits,
        initialize=0,
        units=u.MW,
        doc="Data center curtailment",
    )

    @b.Expression(m.dataCenters, doc="Data center operational cost in $")
    def dataCenterCost(b, dc):
        return (
            b.dataCenterLoad[dc]
            * pyo.units.convert(b.periodLength, to_units=u.hr)
            * m.dataCenterOperationalCost[dc]
        )

    @b.Expression(m.dataCenters, doc="Data center generation cost in $")
    def dataCenterGenerationCost(b, dc):
        return (
            b.dataCenterGeneration[dc]
            * pyo.units.convert(b.periodLength, to_units=u.hr)
            * m.dataCenterGenerationCost[dc]
        )
