# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

import pyomo.environ as pyo
from pyomo.environ import * 
from pyomo.environ import units as u

# from pyomo.gdp import *

from egret.data.model_data import ModelData
from egret.model_library.transmission.tx_utils import scale_ModelData_to_pu
from pyomo.common.timing import TicTocTimer
from pyomo.repn.linear import LinearRepnVisitor
import json
import numpy as np

import math


from math import ceil
from config_options import _get_model_config

# Define what a USD is for pyomo units purposes
# This will be set to a base year and we will do NPV calculations
# based on automatic pyomo unit transformations
u.load_definitions_from_strings(["USD = [currency]"])

rng = np.random.default_rng(seed=123186)

####################################
########## New Work Here ###########
####################################

## TODO: Egret features


# This is only used for reporting potentially bad (i.e., large magnitude) coefficients
# and thus only when that argument is passed
class VisitorConfig(object):
    def __init__(self):
        self.subexpr = {}
        self.var_map = {}
        self.var_order = {}

    def __iter__(self):
        return iter((self.subexpr, self.var_map, self.var_order))


class ExpansionPlanningModel:
    """A generalized generation and transmission expansion planning model."""

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

    def create_model(self):
        """Create concrete Pyomo model object associated with the ExpansionPlanningModel"""

        self.timer.tic("Creating GTEP Model")
        m = ConcreteModel()

        ## TODO: checks for active/built/inactive/unbuilt/etc. gen
        ## NOTE: scale_ModelData_to_pu doesn't account for expansion data -- does it need to?
        if self.data is None:
            raise
        elif type(self.data) is list:
            # If self.data is a list, it is a list of data for
            # representative periods
            m.data_list = self.data
            m.md = scale_ModelData_to_pu(self.data[0])
        else:
            # If self.data is an Egret model data object, representative periods will just copy it unchanged
            m.data_list = None
            m.md = scale_ModelData_to_pu(self.data)
            m.formulation = self.formulation
        
        # [ESR WIP: Add cost_data]
        # TODO: Think about how to do some scaling in cost data
        if self.cost_data is None:
            raise
        else:
            m.mc = self.cost_data

        model_set_declaration(
            m, self.stages, rep_per=[i for i in range(1, self.num_reps + 1)]
        )
        m.representativePeriodLength = Param(
            m.representativePeriods,
            within=PositiveReals,
            default=24,
            units=u.hr
        )
        m.numCommitmentPeriods = Param(
            m.representativePeriods,
            within=PositiveIntegers,
            default=2,
            initialize=self.num_commit,
        )
        m.numDispatchPeriods = Param(
            m.representativePeriods,
            within=PositiveIntegers,
            default=2,
            initialize=self.num_dispatch,
        )
        m.commitmentPeriodLength = Param(
            within=PositiveReals,
            default=1,
            units=u.hr
        )

        # TODO: index by dispatch period? Certainly index by
        # commitment period
        m.dispatchPeriodLength = Param(
            within=PositiveReals,
            initialize=self.duration_dispatch,
            units=u.minutes
        )

        model_data_references(m)
        model_create_investment_stages(m, self.stages)
        create_objective_function(m)

        self.model = m

    ## TODO: this should handle string or i/o object for outfile
    def report_model(self, outfile="pretty_model_output.txt"):
        """Pretty prints Pyomo model to outfile.

        :outfile: (str, optional) _description_. Defaults to "pretty_model_output.txt".
        """
        with open(outfile, "w") as outf:
            self.model.pprint(ostream=outf)

    def report_large_coefficients(self, outfile, magnitude_cutoff=1e5):
        """Dump very large magnitude (>= 1e5) coefficients to a json file.

        :outfile: should accept filename or open file and write there; see how we do this in pyomo elsewhere
        :magnitude_cutoff: magnitude above which to report coefficients
        """
        var_coef_dict = {}
        for e in self.model.component_data_objects(Constraint):
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


def add_investment_variables(
    b,
    investment_stage,
):
    """Add variables to investment stage block.

    :param b: Investment block
    :param investment_stage: Investment stage index or name
    :return: None
    """

    m = b.model()
    b.investmentStage = investment_stage
    
    # Thermal generator disjuncts (operational, installed, retired,
    # disabled, extended)
    @b.Disjunct(m.thermalGenerators)
    def genOperational(disj, gen):
        return

    # ESR TODO: start adding constraints to installed generators. Not
    # used for now.
    @b.Disjunct(m.thermalGenerators)
    def genInstalled(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genRetired(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genDisabled(disj, gen):
        return

    @b.Disjunct(m.thermalGenerators)
    def genExtended(disj, gen):
        return

    @b.Disjunction(m.thermalGenerators)
    def investStatus(disj, gen):
        return [
            disj.genOperational[gen],
            disj.genInstalled[gen],
            disj.genRetired[gen],
            disj.genDisabled[gen],
            disj.genExtended[gen],
        ]

    # Line disjuncts. For now mimicking thermal generator disjuncts,
    # though different states may need to be defined
    @b.Disjunct(m.transmission)
    def branchOperational(disj, branch):
        return

    @b.Disjunct(m.transmission)
    def branchInstalled(disj, branch):
        return

    @b.Disjunct(m.transmission)
    def branchRetired(disj, branch):
        return

    @b.Disjunct(m.transmission)
    def branchDisabled(disj, branch):
        return

    @b.Disjunct(m.transmission)
    def branchExtended(disj, branch):
        return

    # JSC update (done?)
    # @KyleSkolfield: do we differentiate between line and transformer investments?
    @b.Disjunction(m.transmission)
    def branchInvestStatus(disj, branch):
        return [
            disj.branchOperational[branch],
            disj.branchInstalled[branch],
            disj.branchRetired[branch],
            disj.branchDisabled[branch],
            disj.branchExtended[branch],
        ]

    # Renewable generator MW values (operational, installed, retired, extended)
    b.renewableOperational = Var(
        m.renewableGenerators,
        within=NonNegativeReals,
        initialize=0,
        units=u.MW
    )
    b.renewableInstalled = Var(
        m.renewableGenerators,
        within=NonNegativeReals,
        initialize=0,
        units=u.MW
    )
    b.renewableRetired = Var(
        m.renewableGenerators,
        within=NonNegativeReals,
        initialize=0,
        units=u.MW
    )
    b.renewableExtended = Var(
        m.renewableGenerators,
        within=NonNegativeReals,
        initialize=0,
        units=u.MW

    )

    # Track and accumulate costs and penalties
    b.quotaDeficit = Var(
        within=NonNegativeReals,
        initialize=0,
        units=u.MW
    )
    b.operatingCostInvestment = Var(
        within=Reals,
        initialize=0,
        units=u.USD
    )
    b.expansionCost = Var(
        within=Reals,
        initialize=0,
        units=u.USD
    )
    b.renewableCurtailmentInvestment = Var(
        within=NonNegativeReals,
        initialize=0,
        units=u.USD
    )


def add_investment_constraints(
    b,
    investment_stage,
):
    """Add standard inequalities (i.e., those not involving disjunctions) to investment stage block."""

    m = b.model()

    for gen in m.thermalGenerators:
        if (
            m.md.data["elements"]["generator"][gen]["in_service"] == False
            and investment_stage == 1
        ):
            b.genDisabled[gen].indicator_var.fix(True)
            # b.genDisabled[gen].binary_indicator_var.fix(1)
        elif (
            m.md.data["elements"]["generator"][gen]["in_service"] == True
            and investment_stage == 1
        ):
            b.genOperational[gen].indicator_var.fix(True)
            # b.genInstalled[gen].binary_indicator_var.fix(1)

    for branch in m.transmission:
        if (
            m.md.data["elements"]["branch"][branch]["in_service"] == False
            and investment_stage == 1
        ):
            b.branchDisabled[branch].indicator_var.fix(True)
            # b.branchDisabled[branch].binary_indicator_var.fix(1)
        elif (
            m.md.data["elements"]["branch"][branch]["in_service"] == True
            and investment_stage == 1
        ):
            b.branchOperational[branch].indicator_var.fix(True)
            # b.branchInstalled[branch].binary_indicator_var.fix(1)

    # Planning reserve requirement constraint
    ## NOTE: renewableCapacityValue is a percentage of renewableCapacity
    ## TODO: renewableCapacityValue ==> renewableCapacityFactor
    ## NOTE: reserveMargin is a percentage of peakLoad
    ## TODO: check and re-enable with additional bounding transform before bigm
    ## TODO: renewableCapacityValue... should this be time iterated? is it tech based?
    ## is it site based? who can say?
    """@b.Constraint()
    def planning_reserve_requirement(b):
        return (
            sum(
                m.renewableCapacityNameplate[gen]
                * m.renewableCapacityValue[gen]
                * (b.renewableOperational[gen] + b.renewableInstalled[gen])
                for gen in m.renewableGenerators
            )
            + sum(
                m.thermalCapacity[gen]
                * (
                    b.genOperational[gen].indicator_var.get_associated_binary()
                    + b.genInstalled[gen].indicator_var.get_associated_binary()
                )
                for gen in m.thermalGenerators
            )
            >= (1 + m.reserveMargin[investment_stage]) * m.peakLoad[investment_stage]
        )"""

    # maximum investment stage installation
    ## NOTE: temporarily disabled maximum investment as default option
    ## TODO: These capacities shouldn't be enabled by default since they can
    ## easily cause absurd results/possibly even infeasibility.  Will need to add
    ## user-defined handling for this.
    # @b.Constraint(m.regions)
    # def maximum_thermal_investment(b, region):
    #     return (
    #         sum(
    #             m.thermalCapacity[gen]
    #             * b.genInstalled[gen].indicator_var.get_associated_binary()
    #             for gen in m.thermalGenerators & m.gensAtRegion[region]
    #         )
    #         <= b.maxThermalInvestment[region]
    #     )

    # @b.Constraint(m.regions)
    # def maximum_renewable_investment(b, region):
    #     return (
    #         sum(
    #             m.renewableCapacityNameplate[gen]
    #             * b.genInstalled[gen].indicator_var.get_associated_binary()
    #             for gen in m.renewableGenerators & m.gensAtRegion[region]
    #         )
    #         <= b.maxRenewableInvestment[region]
    #         if m.renewableGenerators & m.gensAtRegion[region]
    #         else Constraint.Skip
    #     )

    ## NOTE: The following constraints can be split into rep_per and invest_stage components if desired

    ## NOTE: Constraint (13) in the reference paper
    @b.Constraint(doc="Minimum per-stage renewable generation requirement")
    def renewable_generation_requirement(b):
        renewableSurplusRepresentative = 0
        ## TODO: preprocess loads for the appropriate sum here
        ed = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                renewableSurplusRepresentative += (
                    m.weights[rep_per]
                    # [ESR WIP: Commented since we are including the
                    # period during dispatch stage.]
                    # * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per].commitmentPeriod[com_per].renewableSurplusCommitment
                )
        return (
            # [ESR WIP: Change quotas to be in MWh to be consistent
            # with the lrs of the equation.]
            renewableSurplusRepresentative + b.quotaDeficit
            >= m.renewableQuota[investment_stage] * ed # [ESR WIP: Q: ed is 0 here, is that correct?] 
        )

    @b.Expression(doc="Operating costs for investment period")
    def operatingCostInvestment(b):
        operatingCostRepresentative = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                operatingCostRepresentative += (
                    m.weights[rep_per]
                    # [ESR WIP: Commented since we are including the
                    # period during dispatch stage.]
                    # * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per].commitmentPeriod[com_per].operatingCostCommitment
                )
        return m.investmentFactor[investment_stage] * operatingCostRepresentative

    ## FIXME: investment cost definition needs to be revisited AND possibly depends on
    ## data format.  It is _rare_ for these values to be defined at all, let alone consistently.
    @b.Constraint(doc="Investment costs for investment period in $")
    def investment_cost(b):
        return b.expansionCost == m.investmentFactor[investment_stage] *(
            sum(
                # [ESR WIP: When including the disjunction
                # investStatus, think if we should replace this cost
                # with generatorInstallationCost.]
                m.generatorInvestmentCost[gen]
                * m.thermalCapacity[gen] # in MW
                # m.generatorInstallationCost[gen]
                * m.capitalMultiplier[gen]
                * b.genInstalled[gen].indicator_var.get_associated_binary() 
                for gen in m.thermalGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.capitalMultiplier[gen]
                * b.renewableInstalled[gen] # in MW
                for gen in m.renewableGenerators
            )
            + sum(
                # [ESR WIP: When including the disjunction
                # investStatus, think if we should replace this cost
                # with generatorInstallationCost.]
                m.generatorInvestmentCost[gen]
                # m.generatorInstallationCost[gen]
                * m.extensionMultiplier[gen]
                # [ESR WIP: Term added for unit consistency]
                * m.thermalCapacity[gen]
                * b.genExtended[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.extensionMultiplier[gen]
                * b.renewableExtended[gen]
                for gen in m.renewableGenerators
            )
            # JSC inprog (done?) - added branch investment costs here
            + sum(
                m.branchInvestmentCost[branch]
                * m.branchCapitalMultiplier[branch]
                * b.branchInstalled[branch].indicator_var.get_associated_binary()
                for branch in m.transmission
            )
            + sum(
                m.branchInvestmentCost[branch]
                * m.branchExtensionMultiplier[branch]
                * b.branchExtended[branch].indicator_var.get_associated_binary()
                for branch in m.transmission
            )
        )

    @b.Constraint(doc="Curtailment penalties for investment period")
    def renewable_curtailment_cost(b):
        renewableCurtailmentRep = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                renewableCurtailmentRep += (
                    m.weights[rep_per] # dimensionless
                    * m.commitmentPeriodLength # in hr
                    * b.representativePeriod[rep_per].commitmentPeriod[com_per].renewableCurtailmentCommitment # in MW
                    # [ESR WIP: Q: Do we need to include this term here?] 
                    * m.curtailmentCost
                ) # units are in $
        return (
            b.renewableCurtailmentInvestment # in $
            == m.investmentFactor[investment_stage] * renewableCurtailmentRep
        )


def add_dispatch_variables(
    b,
    dispatch_period,
):
    """Add dispatch-associated variables to representative period block."""

    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    def thermal_generation_limits(b, thermalGen,
                                  doc="Bounds on active generation of thermal generators"):
        return (0, m.thermalCapacity[thermalGen])

    b.thermalGeneration = Var(
        m.thermalGenerators,
        domain=NonNegativeReals,
        bounds=thermal_generation_limits,
        initialize=0,
        units=u.MW,
        doc="Thermal generation"
    )

    # [ESR WIP: Still deciding if this should be Nameplate]
    def renewable_generation_limits(b, renewableGen,
                                    doc="Bounds on active generation of renewable generator in MW"):
        return (0, m.renewableCapacityNameplate[renewableGen])

    b.renewableGeneration = Var(
        m.renewableGenerators,
        domain=NonNegativeReals,
        bounds=renewable_generation_limits,
        initialize=0,
        units=u.MW,
        doc="Renewable generation"
    )

    def curtailment_limits(b, renewableGen,
                           doc="Bounds on renewable generator curtailment in MW"):
        return (0, m.renewableCapacityNameplate[renewableGen])

    b.renewableCurtailment = Var(
        m.renewableGenerators,
        domain=NonNegativeReals,
        bounds=curtailment_limits,
        initialize=0,
        units=u.MW,
        doc="Curtailment of renewable generators"
    )

    @b.Expression(m.renewableGenerators,
                  doc="Surplus generation per renewable generator in MW")
    def renewableGenerationSurplus(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] - b.renewableCurtailment[renewableGen]
        )

    # [ESR WIP: Add the dispatch period length.]
    @b.Expression(m.renewableGenerators, doc="Curtailment cost per generator in $")
    def renewableCurtailmentCost(b, renewableGen):
        return (
            b.renewableCurtailment[renewableGen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.curtailmentCost
        )            

    # [ESR WIP: Add the dispatch period length.]
    @b.Expression(m.thermalGenerators, doc="Cost per thermal generator in $")
    def thermalGeneratorCost(b, gen):
        return (
            b.thermalGeneration[gen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * (m.fixedCost[gen] + m.varCost[gen])
        )

    @b.Expression(m.renewableGenerators, doc="Cost per renewable generator in $")
    def renewableGeneratorCost(b, gen):
        return (
            b.renewableGeneration[gen]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.fixedCost[gen]
        )

    b.loadShed = Var(
        m.buses,
        domain=NonNegativeReals,
        initialize=0,
        units=u.MW,
        doc="Load shed per bus"
    )

    # [ESR WIP: Add the dispatch period length.]
    @b.Expression(m.buses, doc="Load shed cost per bus in $")
    def loadShedCost(b, bus):
        return (
            b.loadShed[bus]
            * pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            * m.loadShedCostperCurtailment # $/MWh
        )

    # Track total dispatch values and costs
    @b.Expression(doc="Total surplus power for renewable generators in MW")
    def renewableSurplusDispatch(b):
        return sum(b.renewableGenerationSurplus[gen] for gen in m.renewableGenerators)
    @b.Expression()
    def thermalGenerationCostDispatch(b):
        return sum(b.thermalGeneratorCost[gen] for gen in m.thermalGenerators)
    @b.Expression()
    def renewableGenerationCostDispatch(b):
        return sum(b.renewableGeneratorCost[gen] for gen in m.renewableGenerators)
    @b.Expression()
    def loadShedCostDispatch(b):
        return sum(b.loadShedCost[bus] for bus in m.buses)
    @b.Expression()
    def curtailmentCostDispatch(b):
        return sum(b.renewableCurtailmentCost[gen] for gen in m.renewableGenerators)
    @b.Expression(doc="Total cost for dispatch in $")
    def operatingCostDispatch(b):
        return (
            b.thermalGenerationCostDispatch
            + b.renewableGenerationCostDispatch
            + b.loadShedCostDispatch
            + b.curtailmentCostDispatch
        )

    @b.Expression(doc="Total curtailment dispatch for renewable generators in MW")
    def renewableCurtailmentDispatch(b):
        return sum(
            b.renewableCurtailment[gen] for gen in m.renewableGenerators
        )

    # Restrictions on flow over uninvested lines are enforced in a
    # disjuction below
    def power_flow_limits(b, branch, doc="Bounds on transmission line capacity"):
        return (
            -m.transmissionCapacity[branch],
            m.transmissionCapacity[branch],
        )

    # (Original) NOTE: this is an abuse of units and needs to be fixed for
    # variable temporal resolution
    b.powerFlow = Var(
        m.transmission,
        domain=Reals,
        bounds=power_flow_limits,
        initialize=0,
        units=u.MW,
        doc="Power flow in MW"
    )

    @b.Disjunct(m.transmission)
    def branchInUse(disj, branch):
        b = disj.parent_block()

        def bus_angle_bounds(disj, bus, doc="Voltage angle"):
            return (-math.pi / 6, math.pi / 6)

        # Only create bus angle variables for the buses associated with this
        # branch that is in use
        disj.branch_buses = [
            bb
            for bb in m.buses
            if (
                m.transmission[branch]["from_bus"] == bb
                or m.transmission[branch]["to_bus"] == bb
            )
        ]

        disj.busAngle = Var(
            disj.branch_buses,
            domain=Reals,
            initialize=0,
            bounds=bus_angle_bounds
        )

        def delta_bus_angle_bounds(disj, bus, doc="Voltage angle"):
            return (-math.pi / 6, math.pi / 6)

        def delta_bus_angle_rule(disj, doc="Maximum bus angle discrepancy"):
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            return disj.busAngle[tb] - disj.busAngle[fb]

        # @KyleSkolfield - I think this var is unused and commented it
        # out, can we delete?
        disj.deltaBusAngle = Var(
            domain=Reals,
            bounds=delta_bus_angle_bounds,
            rule=delta_bus_angle_rule
        )

        ## FIXME
        # @disj.Constraint()
        # def max_delta_bus_angle(disj):
        #     return abs(disj.deltaBusAngle) <= math.pi/6

        @disj.Constraint()
        def dc_power_flow(disj):
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            reactance = m.md.data["elements"]["branch"][branch]["reactance"]
            if m.md.data["elements"]["branch"][branch]["branch_type"] == "transformer":
                reactance *= m.md.data["elements"]["branch"][branch][
                    "transformer_tap_ratio"
                ]
                shift = m.md.data["elements"]["branch"][branch][
                    "transformer_phase_shift"
                ]
            else:
                shift = 0
            return b.powerFlow[branch] == (-1 / reactance) * (
                disj.busAngle[tb] - disj.busAngle[fb] + shift
            ) * u.MW

    @b.Disjunct(m.transmission)
    def branchNotInUse(disj, branch):

        # JSC update (done?) Fixing power flow to 0 and not creating bus angle
        # variables for branches that are not in use
        @disj.Constraint()
        def dc_power_flow(disj):
            return b.powerFlow[branch] == 0 * u.MW

        return

    # JSC update - Branches are either in-use or not. This disjunction may
    # provide the basis for transmission switching in the future
    @b.Disjunction(m.transmission)
    def branchInUseStatus(disj, branch):
        return [
            disj.branchInUse[branch],
            disj.branchNotInUse[branch],
        ]

    # JSC update - If a branch is in use, it must be active
    # Update this when switching is implemented
    @b.LogicalConstraint(m.transmission)
    def must_use_active_branches(b, branch):
        return b.branchInUse[branch].indicator_var.implies(
            lor(
                i_p.branchOperational[branch].indicator_var,
                i_p.branchInstalled[branch].indicator_var,
                i_p.branchExtended[branch].indicator_var,
            )
        )

    # JSC update - If a branch is not in use, it must be inactive.
    # Update this when switching is implemented
    @b.LogicalConstraint(m.transmission)
    def cannot_use_inactive_branches(b, branch):
        return b.branchNotInUse[branch].indicator_var.implies(
            lor(
                i_p.branchDisabled[branch].indicator_var,
                i_p.branchRetired[branch].indicator_var,
            )
        )

    def spinning_reserve_limits(b, thermalGen,
                                doc="Bounds on thermal generator spinning reserve supply"):
        return (
            0 * u.MW,
            m.spinningReserveFraction[thermalGen] * m.thermalCapacity[thermalGen],
        )

    b.spinningReserve = Var(
        m.thermalGenerators,
        domain=NonNegativeReals,
        bounds=spinning_reserve_limits,
        initialize=0,
        units=u.MW,
    )

    def quickstart_reserve_limits(
            b,
            thermalGen,
            doc="Bounds on thermal generator quickstart reserve supply"
    ):
        return (
            0 * u.MW,
            m.quickstartReserveFraction[thermalGen] * m.thermalCapacity[thermalGen],
        )

    b.quickstartReserve = Var(
        m.thermalGenerators,
        domain=NonNegativeReals,
        bounds=quickstart_reserve_limits,
        initialize=0,
        units=u.MW
    )


def add_dispatch_constraints(b, disp_per):
    """Add dispatch-associated inequalities to representative period block."""
    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    
    # [ESR WIP: Commented for now but think about how to implement this in
    # a better way.]
    # for key in m.loads.keys():
    #     m.loads[key] *= max(0, rng.normal(0.5, 0.2))

    # units_load = u.MW
    @b.Constraint(m.buses, doc="Energy balance")
    def flow_balance(b, bus):
        balance = 0

        # [ESR WIP: Comment load and call all the loads as a parameter
        # instead of the original dictionary. Also, note that the
        # loads are now declared for all the buses in m.buses, and set
        # to 0 for the buses that are not in m.load_buses.]
        # load = value(m.loads.get(bus)) or 0
        
        end_points = [
            line for line in m.transmission if m.transmission[line]["from_bus"] == bus
        ]
        start_points = [
            line for line in m.transmission if m.transmission[line]["to_bus"] == bus
        ]
        gens = [
            gen
            for gen in m.generators
            if m.md.data["elements"]["generator"][gen]["bus"] == bus
        ]
        balance -= sum(b.powerFlow[i] for i in end_points)
        balance += sum(b.powerFlow[i] for i in start_points)
        balance += sum(b.thermalGeneration[g] for g in gens if g in m.thermalGenerators)
        balance += sum(
            b.renewableGeneration[g] for g in gens if g in m.renewableGenerators
        )
        balance -= m.loads[bus] # remove units_load since the parameter already has units
        balance += b.loadShed[bus]
        return balance == 0 * u.MW

    # NOTE: In comparison to reference work, this is *per renewable
    # generator*
    @b.Constraint(m.renewableGenerators, doc="Capacity factor")
    def capacity_factor(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen]
            + b.renewableCurtailment[renewableGen]
            == c_p.renewableCapacityExpected[renewableGen]
        )

    @b.Constraint(m.renewableGenerators)
    def operational_renewables_only(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen]
            <= i_p.renewableInstalled[renewableGen]
            + i_p.renewableOperational[renewableGen]
        )

    # RESERVE -- total operating (spinning + quickstart)
    # Total operating reserve constraint
    ## NOTE: min operating reserve is a percentage of load
    ## FIXME: Reserve enforcement causes infeasibility issues.  We should track
    ## reserve shortage in some way and find a way to penalize it -- how is this
    ## done in ISOs?  Is it an issue to assign this as a regional
    # @b.Constraint(m.regions)
    # def total_operating_reserve(b, region):
    #     return sum(
    #         b.spinningReserve[gen] + b.quickstartReserve[gen]
    #         for gen in m.thermalGenerators & m.gensAtRegion[region]
    #     ) >= m.minOperatingReserve[region] * (
    #         sum(
    #             (m.loads.get(bus) or 0)
    #             for bus in m.buses
    #             if m.md.data["elements"]["bus"][bus]["area"] == region
    #         )
    #     )

    # # Total spinning reserve constraint
    # @b.Constraint(m.regions)
    # def total_spinning_reserve(b, region):
    #     return sum(
    #         b.spinningReserve[gen]
    #         for gen in m.thermalGenerators & m.gensAtRegion[region]
    #     ) >= m.minSpinningReserve[region] * sum(
    #         (m.loads.get(bus) or 0)
    #         for bus in m.buses
    #         if m.md.data["elements"]["bus"][bus]["area"] == region
    #     )


def add_commitment_variables(b, commitment_period):
    """Add variables and disjuncts to commitment period block."""
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    # Define disjunction on generator status: on/startup/shutdown/off
    @b.Disjunct(m.thermalGenerators)
    def genOn(disj, generator):
        # operating limits
        ## NOTE: Reminder: thermalMin is a percentage of thermalCapacity
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods, doc="Minimum operating limits")
        def operating_limit_min(d, dispatchPeriod):
            return (
                m.thermalMin[generator]
                <= b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
            )

        @disj.Constraint(b.dispatchPeriods, doc="Maximum operating limits")
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalCapacity[generator]
            )

        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators,
                         doc="Ramp up limits for fully-on thermal generators")
        def ramp_up_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                <= m.rampUpRates[generator]
                * b.dispatchPeriod[dispatchPeriod].periodLength
                * m.thermalCapacity[generator]
                if dispatchPeriod != 1
                else Constraint.Skip
            )


        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators,
                         doc="Ramp down limits for fully-on thermal generators")
        def ramp_down_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                <= m.rampDownRates[generator] # in MW/min
                * b.dispatchPeriod[dispatchPeriod].periodLength # in min
                * m.thermalCapacity[generator] # in MW
                if dispatchPeriod != 1
                else Constraint.Skip
            )
                    
        ##NOTE: maxSpiningReserve is a percentage of thermalCapacity
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators,
                         doc="Maximum spinning reserve")
        def max_spinning_reserve(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.maxSpinningReserve[generator] * m.thermalCapacity[generator]
            )

    @b.Disjunct(m.thermalGenerators)
    def genStartup(disj, generator):
        b = disj.parent_block()

        # (Original) NOTE: Reminder: thermalMin is a percentage of
        # thermalCapacity
        @disj.Constraint(b.dispatchPeriods, doc="Operating limits")
        def operating_limit_min(d, dispatchPeriod):
            return 0 <= b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]

        @disj.Constraint(b.dispatchPeriods, doc="Maximum operating limits")
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalMin[generator]
            )

        # (Original) TODO: is this max necessary? I would like to
        # remove
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators,
                         doc="Ramp up constraints for generators starting up")
        def ramp_up_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                <= max(value(m.thermalMin[generator]),
                       value(m.rampUpRates[generator]) # make sure the time units are consistent
                       * value(b.dispatchPeriod[dispatchPeriod].periodLength),
                       ) * m.thermalCapacity[generator]
                if dispatchPeriod != 1
                else Constraint.Skip
            )

    @b.Disjunct(m.thermalGenerators)
    def genShutdown(disj, generator):
        b = disj.parent_block()

        # operating limits
        ## NOTE: Reminder: thermalMin is a percentage of thermalCapacity
        @disj.Constraint(b.dispatchPeriods)
        def operating_limit_min(d, dispatchPeriod):
            return 0 <= b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]

        # Maximum operating limits
        @disj.Constraint(b.dispatchPeriods)
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalMin[generator]
            )

        # Ramp down constraints for generators shutting down
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        def ramp_down_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                <= max(
                    value(m.thermalMin[generator]),
                    value(m.rampDownRates[generator]) # make sure the time units are consistent
                    * b.dispatchPeriod[dispatchPeriod].periodLength,
                )
                * m.thermalCapacity[generator]
                if dispatchPeriod != 1
                else Constraint.Skip
            )

    @b.Disjunct(m.thermalGenerators)
    def genOff(disj, generator):
        b = disj.parent_block()

        # operating limits
        ## NOTE: Reminder: thermalMin is a percentage of thermalCapacity
        @disj.Constraint(b.dispatchPeriods)
        def operating_limit_max(disj, dispatchPeriod):
            return b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator] <= 0

        # Maximum quickstart reserve constraint
        ## NOTE: maxQuickstartReserve is a percentage of thermalCapacity
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        def max_quickstart_reserve(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].quickstartReserve[generator]
                <= m.maxQuickstartReserve[generator] * m.thermalCapacity[generator]
            )

    @b.Disjunction(m.thermalGenerators)
    def genStatus(disj, generator):
        return [
            disj.genOn[generator],
            disj.genStartup[generator],
            disj.genShutdown[generator],
            disj.genOff[generator],
        ]

    # Generators cannot be committed unless they are operational or just installed
    @b.LogicalConstraint(m.thermalGenerators)
    def commit_active_gens_only(b, generator):
        return lor(
            b.genOn[generator].indicator_var,
            b.genStartup[generator].indicator_var,
            b.genShutdown[generator].indicator_var,
        ).implies(
            lor(
                i_p.genOperational[generator].indicator_var,
                i_p.genInstalled[generator].indicator_var,
                i_p.genExtended[generator].indicator_var,
            )
        )


def add_commitment_constraints(
    b,
    comm_per,
):
    """Add commitment-associated disjunctions and constraints to representative period block."""
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    @b.Expression(doc="Total renewable surplus/deficit for commitment block")
    def renewableSurplusCommitment(b):
        return sum(
            # [ESR WIP: Q: Commenting the commitment period since I
            # don't think we need to include it.]
            # pyo.units.convert(m.dispatchPeriodLength, to_units=u.hr)
            b.dispatchPeriod[disp_per].renewableSurplusDispatch # in MW
            for disp_per in b.dispatchPeriods
        )

    # Define total operating costs for commitment block
    ## TODO: Replace this constraint with expressions using bounds transform
    ## NOTE: expressions are stored in gtep_cleanup branch
    ## costs considered need to be re-assessed and account for missing data
    
    # [ESR WIP: The fixed costs for thermal and renewable generators
    # are included in the dispatch stage.]
    @b.Expression()
    def operatingCostCommitment(b):
        return (
            sum(
                ## FIXME: update test objective value when this changes; ready to uncomment
                # (m.dispatchPeriodLength / 60) *
                # [ESR WIP: This term includes the op cost for each
                # 15-min dispatch period.]
                b.dispatchPeriod[disp_per].operatingCostDispatch # in $
                for disp_per in b.dispatchPeriods
            )
            + sum(
                m.startupCost[gen] # in $
                * b.genStartup[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
        )

    @b.Expression(doc="Total curtailment for commitment block in MW")
    def renewableCurtailmentCommitment(b):
        return sum(
            b.dispatchPeriod[disp_per].renewableCurtailmentDispatch
            for disp_per in b.dispatchPeriods
        )


def commitment_period_rule(b, commitment_period):
    """Create commitment period block.

    :param b: commitment period block
    :param commitment_period: corresponding commitment period label
    """
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    b.commitmentPeriod = commitment_period
    b.commitmentPeriodLength = Param(within=PositiveReals, default=1, units=u.hr)
    b.dispatchPeriods = RangeSet(m.numDispatchPeriods[r_p.currentPeriod])
    b.carbonTax = Param(default=0)
    b.dispatchPeriod = Block(b.dispatchPeriods)

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

    # [WIP: Corrected to be in the block "b", not in "m". Also,
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
            b.renewableCapacityExpected[renewableGen] = m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"][commitment_period - 1] * units_renewable_capacity
    
    ## TODO: Redesign load scaling and allow nature of it as argument

    # Demand at each bus
    temp_scale = 3
    temp_scale = 10

    scale_loads = True
    if scale_loads:
        for load_n in m.load_buses:
            m.loads[load_n] = (
                temp_scale
                * (
                    1
                    + (temp_scale + i_p.investmentStage) / (temp_scale + len(m.stages))
                )
            ) * m.md.data["elements"]["load"][load_n]["p_load"]["values"][commitment_period - 1]
    else:
        for load_n in m.load_buses:
            m.loads[load_n] = m.md.data["elements"]["load"][load_n]["p_load"]["values"][commitment_period - 1]

    ## TODO: This feels REALLY inelegant and bad.
    ## TODO: Something weird happens if I say periodLength has a unit
    for period in b.dispatchPeriods:
        b.dispatchPeriod[period].periodLength = Param(within=PositiveReals, default=1)
        add_dispatch_variables(b.dispatchPeriod[period], period)

    add_commitment_variables(b, commitment_period)
    add_commitment_constraints(b, commitment_period)

    for period in b.dispatchPeriods:
        add_dispatch_constraints(b.dispatchPeriod[period], period)


def add_representative_period_variables(b, rep_per):
    m = b.model()
    i_p = b.parent_block()

    b.renewableSurplusRepresentative = Var(within=Reals, initialize=0, units=u.USD)


def add_representative_period_constraints(b, rep_per):
    m = b.model()
    i_p = b.parent_block()

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_shutdown(b, commitmentPeriod, thermalGen):
        req_shutdown_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"])
        )
        return (
            atmost(
                req_shutdown_periods - 1,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genShutdown[thermalGen]
                    .indicator_var
                    for j in range(min(req_shutdown_periods, commitmentPeriod - 1))
                ],
            ).land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genShutdown[thermalGen]
                .indicator_var
            )
            # | b.commitmentPeriod[commitmentPeriod-1].genOn.indicator_var)
            .implies(
                b.commitmentPeriod[commitmentPeriod]
                .genShutdown[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_off_after_shutdown(b, commitmentPeriod, thermalGen):
        req_shutdown_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"])
        )
        return (
            atleast(
                req_shutdown_periods,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genShutdown[thermalGen]
                    .indicator_var
                    for j in range(min(req_shutdown_periods, commitmentPeriod - 1))
                ],
            )
            .land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genShutdown[thermalGen]
                .indicator_var
            )
            .implies(
                b.commitmentPeriod[commitmentPeriod].genOff[thermalGen].indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_startup(b, commitmentPeriod, thermalGen):
        req_startup_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"])
        )
        return (
            atmost(
                req_startup_periods - 1,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genStartup[thermalGen]
                    .indicator_var
                    for j in range(min(req_startup_periods, commitmentPeriod - 1))
                ],
            ).land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genStartup[thermalGen]
                .indicator_var
            )
            # | b.commitmentPeriod[commitmentPeriod-1].genOn.indicator_var)
            .implies(
                b.commitmentPeriod[commitmentPeriod]
                .genStartup[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_on_after_startup(b, commitmentPeriod, thermalGen):
        req_startup_periods = ceil(
            1 / float(m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"])
        )
        return (
            atleast(
                req_startup_periods,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genStartup[thermalGen]
                    .indicator_var
                    for j in range(min(req_startup_periods, commitmentPeriod - 1))
                ],
            )
            .land(
                b.commitmentPeriod[commitmentPeriod - 1]
                .genStartup[thermalGen]
                .indicator_var
            )
            .implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_uptime(b, commitmentPeriod, thermalGen):
        return (
            atmost(
                int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"]) - 1,
                [
                    b.commitmentPeriod[commitmentPeriod - j - 1]
                    .genOn[thermalGen]
                    .indicator_var
                    for j in range(
                        min(
                            int(
                                m.md.data["elements"]["generator"][thermalGen][
                                    "min_up_time"
                                ]
                            ),
                            commitmentPeriod - 1,
                        )
                    )
                ],
            )
            .land(
                b.commitmentPeriod[commitmentPeriod - 1].genOn[thermalGen].indicator_var
            )
            .implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
            )
            if commitmentPeriod
            != 1  # int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"])+1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_shutdown_after_uptime(b, commitmentPeriod, thermalGen):
        return (
            (
                atleast(
                    int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"]),
                    [
                        b.commitmentPeriod[commitmentPeriod - j - 1]
                        .genOn[thermalGen]
                        .indicator_var
                        for j in range(
                            min(
                                int(
                                    m.md.data["elements"]["generator"][thermalGen][
                                        "min_up_time"
                                    ]
                                ),
                                commitmentPeriod - 1,
                            )
                        )
                    ],
                ).land(
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .genOn[thermalGen]
                    .indicator_var
                )
            ).implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
                | b.commitmentPeriod[commitmentPeriod]
                .genShutdown[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_downtime(b, commitmentPeriod, thermalGen):
        return (
            (
                atmost(
                    int(m.md.data["elements"]["generator"][thermalGen]["min_down_time"])
                    - 1,
                    [
                        b.commitmentPeriod[commitmentPeriod - j - 1]
                        .genOff[thermalGen]
                        .indicator_var
                        for j in range(
                            min(
                                int(
                                    m.md.data["elements"]["generator"][thermalGen][
                                        "min_down_time"
                                    ]
                                ),
                                commitmentPeriod - 1,
                            )
                        )
                    ],
                ).land(
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .genOff[thermalGen]
                    .indicator_var
                )
            ).implies(
                b.commitmentPeriod[commitmentPeriod].genOff[thermalGen].indicator_var
            )
            if commitmentPeriod
            != 1  # >= int(m.md.data["elements"]["generator"][thermalGen]["min_down_time"])+1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_start_after_downtime(b, commitmentPeriod, thermalGen):
        return (
            (
                atleast(
                    int(
                        m.md.data["elements"]["generator"][thermalGen]["min_down_time"]
                    ),
                    [
                        b.commitmentPeriod[commitmentPeriod - j - 1]
                        .genOff[thermalGen]
                        .indicator_var
                        for j in range(
                            min(
                                int(
                                    m.md.data["elements"]["generator"][thermalGen][
                                        "min_down_time"
                                    ]
                                ),
                                commitmentPeriod - 1,
                            )
                        )
                    ],
                ).land(
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .genOff[thermalGen]
                    .indicator_var
                )
            ).implies(
                b.commitmentPeriod[commitmentPeriod].genOff[thermalGen].indicator_var
                | b.commitmentPeriod[commitmentPeriod]
                .genStartup[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )


def representative_period_rule(
    b,
    representative_period,
):
    """Create representative period block.

    :b: Representative period block
    :representative_period: corresponding representative period label
    """
    m = b.model()
    i_s = b.parent_block()
    
    b.currentPeriod = representative_period

    b.commitmentPeriods = RangeSet(m.numCommitmentPeriods[representative_period])
    b.commitmentPeriod = Block(b.commitmentPeriods, rule=commitment_period_rule)

    add_representative_period_variables(b, representative_period)
    add_representative_period_constraints(b, representative_period)


def investment_stage_rule(
    b,
    investment_stage,
):
    """Creates investment stage block.

    :b: Investment block
    :investment_stage: ID for current investment stage
    """
    m = b.parent_block()

    b.year = m.years[investment_stage - 1]
    print(f'b.year = {b.year}')

    ##########
    # [ESR WIP: Save lists with all relevant costs (fixed and variable
    # operating costs, fuel costs, and investment costs) for thermal
    # and renewable generators. Please refer to gtep_data_processing
    # for more details about the preprocessing of this data. NOTES:
    # The "capex" in the investment costs already include the interest
    # rate for each generator. Also, note that this data only covers
    # three years: 2025, 2030, and 2035. If more investment years are
    # needed, more data should be included in the data file for data
    # processing.

    # [ESR WIP: Assume we have two types of generators: thermal "CT"
    # (with gas fuel) and renewable "PV" (with "sun" as fuel).]

    gen_thermal_type = 'CT'
    gen_renewable_type = 'PV'

    m.genThermalInvCost = []
    m.genThermalFuelCost = []
    m.genThermalFixOpCost = []
    m.genThermalVarOpCost = []
    m.genRenewableInvCost = []
    m.genRenewableFuelCost = []
    m.genRenewableFixOpCost = []
    m.genRenewableVarOpCost = []
    for index, row in m.mc.gen_data_target.iterrows():
        if row['Unit Type'].startswith(gen_thermal_type):
            m.genThermalInvCost.append(row[f'capex_{b.year}']) # in $/kW
            m.genThermalFixOpCost.append(row[f'fixed_ops_{b.year}']) # in $/kW-yr
            m.genThermalVarOpCost.append(row[f'var_ops_{b.year}']) # $/MWh
            m.genThermalFuelCost.append(row[f'fuel_costs_{b.year}'])
                    
        elif row['Unit Type'].startswith(gen_renewable_type):
            m.genRenewableInvCost.append(row[f'capex_{b.year}']) # in $/kW
            m.genRenewableFixOpCost.append(row[f'fixed_ops_{b.year}']) # in $/kW-yr
            m.genRenewableVarOpCost.append(row[f'var_ops_{b.year}']) # $/MWh
            m.genRenewableFuelCost.append(row[f'fuel_costs_{b.year}'])
            
        else:
            continue

    # [ESR WIP: Update data for fixed and variable costs here since
    # they depend on the investment year. Also, convert the units to
    # be consistent.]
    units_fixed_cost = u.USD / (u.kW * u.year)
    units_var_cost = u.USD / (u.MW * u.hr)
    units_inv_cost = u.USD / u.kW
    units_fuel_cost = u.USD / (u.MW * u.hr)
    for gen in m.generators:
        if m.md.data["elements"]["generator"][gen]["generator_type"] == "thermal":
            m.fixedCost[gen] = pyo.units.convert(m.genThermalFixOpCost[0] * units_fixed_cost,
                                                 to_units= u.USD / (u.MW * u.hr)) 

            m.varCost[gen] = m.genThermalVarOpCost[0] * units_var_cost

            m.generatorInvestmentCost[gen] = pyo.units.convert(
                m.genThermalInvCost[0] * units_inv_cost,
                to_units=u.USD/u.MW
            )

            # [ESR WIP: Add fuel costs from preprocessed
            # data. Consider this cost is for Natural Gas generators,
            # not coal.]
            m.fuelCost[gen] = m.genThermalFuelCost[0] * units_fuel_cost
            
        else:
            # For renewable

            m.fixedCost[gen] = pyo.units.convert(m.genRenewableFixOpCost[0] * units_fixed_cost, 
                                                 to_units= u.USD / (u.MW * u.hr)) 
            m.varCost[gen] = m.genRenewableVarOpCost[0] * units_var_cost

            m.generatorInvestmentCost[gen] = pyo.units.convert(
                m.genRenewableInvCost[0] * units_inv_cost,
                to_units=u.USD/u.MW
            )

    # Final (converted) units are:
    # fixed cost = $/ MWh
    # var cost = $/MWh
    # inv cost = $/Mw
    # fuel cost = $/MWh
    
    # Cost per MW of curtailed renewable energy (Original) NOTE: what
    # should this be valued at?  This being both curtailment and load
    # shed.

    # [ESR WIP: Recalculate "curtailmentCost" and
    # "loadShedCostperCurtailment" (formerly "loadShedCost") since
    # they depend on "fuelCost". NOTE: These were originally defined
    # as parameters in the function model_data_reference after
    # "fuelCost" was defined.]
    m.curtailmentCost = 2 * max(value(m.fuelCost[gen]) for gen in m.thermalGenerators)
    m.loadShedCostperCurtailment = 1000 * m.curtailmentCost

    ##########
    
    b.representativePeriods = [
        p
        for p in m.representativePeriods
        # if m.representativePeriodStage[p] == investment_stage
    ]
    add_investment_variables(b, investment_stage)   

    b.representativePeriod = Block(
        b.representativePeriods, rule=representative_period_rule
    )
    b.maxThermalInvestment = Param(m.regions, default=1000, units=u.MW)
    b.maxRenewableInvestment = Param(m.regions, default=1000, units=u.MW)

    add_investment_constraints(b, investment_stage)


def create_objective_function(m):
    """Creates objective function.  Total cost is operating cost plus
    expansion cost plus penalty cost (penalties include generation deficits,
    renewable quota deficits, and curtailment)
    :param m: Pyomo GTEP model.
    """
    if len(m.stages) > 1:
        m.operatingCost = sum(
            m.investmentStage[stage].operatingCostInvestment for stage in m.stages
        )
        m.expansionCost = sum(
            m.investmentStage[stage].expansionCost for stage in m.stages
        )
        m.penaltyCost = sum(
            m.deficitPenalty[stage]
            * m.investmentFactor[stage]
            * m.investmentStage[stage].quotaDeficit
            + m.investmentStage[stage].renewableCurtailmentInvestment
            for stage in m.stages
        )

    ##### units problem
    @m.Objective()
    def total_cost_objective_rule(m):
        if len(m.stages) > 1:
            return (
                m.operatingCost
                + m.expansionCost
                + m.penaltyCost
            )
        else:
            return (
                m.investmentStage[1].operatingCostInvestment
                + m.investmentStage[1].expansionCost
                + m.deficitPenalty[1]
                * m.investmentFactor[1]
                * m.investmentStage[1].quotaDeficit
                + m.investmentStage[1].renewableCurtailmentInvestment
            )


def model_set_declaration(m, stages, rep_per=["a", "b"], com_per=2, dis_per=2):
    """
    Creates Pyomo Sets necessary (convenient) for solving the GTEP model.

    :m: Pyomo model object
    :stages: Number of stages in investment horizon
    """

    m.buses = Set(
        initialize=m.md.data["elements"]["bus"].keys(), doc="Individual buses"
    )

    m.regions = Set(
        initialize=(
            m.md.data["elements"]["bus"][bus]["area"]
            for bus in m.md.data["elements"]["bus"]
        ),
        doc="Regions / clusters of buses",
    )

    ## TODO: Right now, this means that branches can only be specified entirely as standard
    ## or as dc ... not mix-and-match
    if len(m.md.data["elements"]["branch"]) == 0:
        m.md.data["elements"]["branch"] = m.md.data["elements"]["dc_branch"]

    m.transmission = {
        branch: {
            "from_bus": m.md.data["elements"]["branch"][branch]["from_bus"],
            "to_bus": m.md.data["elements"]["branch"][branch]["to_bus"],
            "reactance": m.md.data["elements"]["branch"][branch]["reactance"],
        }
        for branch in m.md.data["elements"]["branch"]
    }

    m.generators = Set(
        initialize=m.md.data["elements"]["generator"].keys(),
        doc="All generators",
    )

    m.thermalGenerators = Set(
        within=m.generators,
        initialize=(
            gen
            for gen in m.generators
            if m.md.data["elements"]["generator"][gen]["generator_type"] == "thermal"
        ),
        doc="Thermal generators; subset of all generators",
    )

    m.renewableGenerators = Set(
        within=m.generators,
        initialize=(
            gen
            for gen in m.generators
            if m.md.data["elements"]["generator"][gen]["generator_type"] == "renewable"
        ),
        doc="Renewable generators; subset of all generators",
    )

    # [ESR WIP: Add set for transmission lines, relevant in
    # model_data_references.]
    m.lines =  Set(
        initialize=m.transmission.keys(), doc="Individual transmission lines"
    )

    ## NOTE: will want to cover baseline generator types in IDAES

    if m.md.data["elements"].get("storage"):
        m.storage = Set(
            initialize=(batt for batt in m.md.data["elements"]["storage"]),
            doc="Potential storage units",
        )

    ## TODO: make sure time units are both definable and consistent without being forced

    m.stages = RangeSet(stages, doc="Set of planning periods")

    m.representativePeriods = Set(
        initialize=rep_per,
        doc="Set of representative periods for each planning period",
    )
    
def model_data_references(m):
    """Creates and labels data for GTEP model; ties input data
    to model directly.
    :param m: Pyomo model object
    """

    m.thermalCapacity = Param(
        m.thermalGenerators,
        initialize={thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_max"]
                    for thermalGen in m.thermalGenerators},
        mutable=True,
        units=u.MW,
        doc="Maximum output of each thermal generator"
    )

    m.lifetimes = Param(
        m.generators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["lifetime"]
                    for gen in m.generators},
        mutable=True,
        units=u.year,
        doc="Lifetime of each generator"
    )

    m.thermalMin = Param(
        m.thermalGenerators,
        initialize={thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_min"]
                    for thermalGen in m.thermalGenerators},
        mutable=True,
        units=u.MW,
        doc="Minimum output of each thermal generator"
    )

    # [WIP: Rename since the name was repeated in the
    # commitment_period_rule function. Check if this is correct.]    
    m.renewableCapacityNameplate = Param(
        m.renewableGenerators,
        initialize={renewableGen: (0 if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float
                                   else max(m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"]))
                    for renewableGen in m.renewableGenerators},
        mutable=True,
        units=u.MW,
        doc="Maximum output of each renewable generator"
    )

    # TODO: WHAT HAVE I DONE HERE I HATE IT and JSC made it worse...

    # [ESR WIP: Take only the value for renewable capacity when using
    # max() to avoid errors.]
    m.renewableCapacityValue = Param(
        m.renewableGenerators,
        initialize={renewableGen: (0 if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float
                                   else min(m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"])
                                   / max(1, value(m.renewableCapacityNameplate[renewableGen])))
                    for renewableGen in m.renewableGenerators},
        mutable=True,
        units=u.dimensionless,
        doc="Fraction of generation capacity that can be reliably counted toward planning reserve"
    )

    # [ESR WIP: From case data, the value is divided by 100, which is
    # the per units conversion.]
    m.transmissionCapacity = Param(
        m.lines, 
        initialize={transmissionLine: m.md.data["elements"]["branch"][transmissionLine]["rating_long_term"]
                    for transmissionLine in m.lines},
        # mutable=True,
        units=u.MW,
        doc="Long term thermal rating of each transmission line"
    )
    
    m.spinningReserveFraction = Param(
        m.thermalGenerators,
        initialize={thermalGen: m.md.data["elements"]["generator"][thermalGen]["spinning_reserve_frac"]
                    for thermalGen in m.thermalGenerators},
        mutable=True,
        units=u.dimensionless,
        doc="Maximum fraction maximum thermal generation output that can be supplied as spinning reserve"
    )
    
    m.quickstartReserveFraction = Param(
        m.thermalGenerators,
        initialize={thermalGen: m.md.data["elements"]["generator"][thermalGen]["quickstart_reserve_frac"]
                    for thermalGen in m.thermalGenerators},
        # mutable=True,
        units=u.dimensionless,
        doc="Maximum fraction of maximum thermal generation output that can be supplied as quickstart reserve"
    )

    # [ESR WIP: When creating a Param for loads, an error occurs since
    # the load at each bus is a dictionary. To avoid this, I
    # initialized a m.loads parameter with a value of 0 and scaled it
    # with the right value in commitment_period_rule. I also created a
    # new set for the buses that have loads.]
    m.load_buses = pyo.Set(initialize=[i for i in m.md.data["elements"]["load"]])
    m.loads = pyo.Param(
        m.buses,
        initialize={load_n: 0 for load_n in m.buses},
        mutable=True,
        units=u.MW,
        doc="Demand at each bus"
    )
    # m.loads = {
    #     m.md.data["elements"]["load"][load_n]["bus"]: m.md.data["elements"]["load"][load_n]["p_load"]
    #     # for load_n in m.md.data["elements"]["load"]
    #     for load_n in m.load_buses
    # }

    ## NOTE: lazy fixing for dc_branch and branch... but should be an ok lazy fix

    # [ESR WIP: Commented for now since it is not use in this case but
    # might be used in the future when considering ACOPF]
    # m.lossRate = Param(
    #     m.transmission,
    #     initialize={branch: (m.md.data["elements"]["branch"][branch].get("loss_rate") or 0)
    #                 for branch in m.transmission},
    #     mutable=True,
    #     # units=,
    #     doc="Per-distance-unit multiplicative loss rate for each transmission line"
    # )


    ## NOTE: lazy fixing for dc_branch and branch... but should be an ok lazy fix
    m.distance = Param(
        m.transmission,
        initialize={branch: (m.md.data["elements"]["branch"][branch].get("distance") or 0)
                    for branch in m.transmission},
        mutable=True,
        units=u.m,
        doc="Distance between terminal buses for each transmission line"
    )

    # JSC TODO: Add cost of investment in each new branch to input data. Currently
    # selected 0 to ensure investments will be selected if needed
    m.branchInvestmentCost = Param(
        m.transmission,
        initialize={branch: (m.md.data["elements"]["branch"][branch].get("capital_cost") or 0)
                    for branch in m.transmission},
        mutable=True,
        units=u.USD,
        doc="Investment cost for each new branch"
    )

    # JSC TODO: Add branch capital multiplier to input data.
    m.branchCapitalMultiplier = Param(
        m.transmission,
        initialize={branch: (m.md.data["elements"]["branch"][branch].get("capital_multiplier") or 1)
                    for branch in m.transmission},
        mutable=True,
        units=u.dimensionless
    )
    
    m.branchExtensionMultiplier = Param(
        m.transmission,
        initialize={branch: (m.md.data["elements"]["branch"][branch].get("extension_multiplier") or 1)
                    for branch in m.transmission},
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life extension for each generator expressed as a fraction of initial investment cost"
    )

    ## TODO: These should go into each stage -- check where these
    ## values should come from
    m.peakLoad = Param(m.stages, default=0, units=u.MW)
    m.reserveMargin = Param(m.stages, default=0, units=u.MW)
    m.renewableQuota = Param(m.stages, default=0, units=u.MW)
    m.weights = Param(m.representativePeriods, default=1)
    m.investmentFactor = Param(m.stages, default=1, mutable=True, units=u.dimensionless)
    m.deficitPenalty = Param(m.stages, default=1, units=u.USD / u.MW)

    # (Original) NOTE: Lazy approx for NPV. [TODO: don't lazily approx
    # NPV, add it into unit handling and calculate from actual time
    # frames]

    # [ESR WIP: Commented since it is already included in the costs we
    # have from preprocessing stage.]
    # for stage in m.stages:
    #     m.investmentFactor[stage] *= 1 / ((1.04) ** (5 * stage))  

    # # [ESR WIP: Commented for now but depends on the type of data we
    # # are using for generators.]
    # m.startFuel = Param(
    #     m.generators,
    #     initialize={gen: m.md.data["elements"]["generator"][gen]["start_fuel"]
    #                 for gen in m.generators},
    #     mutable=True,
    #     # units=
    #     doc="Amount of fuel required to be consumed for startup process for each generator"
    # )

    # [ESR WIP: Original fuel cost. This is re-defined in the function
    # investment_stage_rule with values from preprocessed data.]
    m.fuelCost = Param(
        m.thermalGenerators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["fuel_cost"] if "RTS-GMLC" in m.md.data["system"]["name"]
                    else m.md.data["elements"]["generator"][gen]["p_cost"]["values"][1]
                    for gen in m.thermalGenerators},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Cost per unit of fuel at each generator"
    )
    
    m.emissionsFactor = Param(
        m.generators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["emissions_factor"]
                    for gen in m.generators},
        mutable=True,
        units=u.dimensionless,
        doc="Full lifecycle CO_2 emission factor for each generator"
    )

    # [ESR WIP: Include start-up cost only in thermal generators assuming a natural gas plant.]
    m.startupCost = Param(
        m.thermalGenerators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["non_fuel_startup_cost"]
                    for gen in m.thermalGenerators},
        mutable=True,
        units=u.USD,
        doc="Flat startup cost for thermal generators"
    )

    # (Arbitrary) multiplier corresponds to depreciation schedules for
    # individual technologies; higher values are indicative of slow
    # depreciation
    m.capitalMultiplier = Param(
        m.generators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["capital_multiplier"]
                    for gen in m.generators},
        mutable=True,
        units=u.dimensionless,
        doc="(Arbitrary) multiplier for new generator investments"
    )

    m.extensionMultiplier = Param(
        m.generators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["extension_multiplier"]
                    for gen in m.generators},
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life extension for each generator expressed as a fraction of initial investment cost"
    )

    # [ESR WIP: Replace original generator investment costs with costs
    # from preprocessed data. These are fixed to 0 here but re-defined
    # in the function investment_stage_rule.]
    m.generatorInvestmentCost = Param(
        m.generators,
        initialize={gen:1 for gen in m.generators},
        mutable=True,
        units=u.USD / u.MW,
        doc="Investment cost for all generators"
    )

    m.minOperatingReserve = Param(
        m.regions,
        initialize={region: m.md.data["system"]["min_operating_reserve"] for region in m.regions},
        mutable=True,
        units=u.dimensionless,
        doc="Minimum operating reserve as a fraction of load within a region"
    )

    m.minSpinningReserve = Param(
        m.regions,
        initialize={region: m.md.data["system"]["min_spinning_reserve"] for region in m.regions},
        mutable=True,
        units=u.dimensionless,
        doc="Minimum spinning reserve as a fraction of load within a region"
    )

    m.maxSpinningReserve = Param(
        m.thermalGenerators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["max_spinning_reserve"]
                    for gen in m.thermalGenerators},
        mutable=True,
        units=u.dimensionless,
        doc="Maximum spinning reserve available for each generator as a fraction maximum generator output"
    )

    m.maxQuickstartReserve = Param(
        m.thermalGenerators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["max_quickstart_reserve"]
                    for gen in m.thermalGenerators},
        # mutable=True,
        units=u.dimensionless,
        doc="Maximum quickstart reserve available for each generator as a fraction maximum generator output"
    )
    
    m.rampUpRates = Param(
        m.thermalGenerators,
        initialize={thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"]
                    for thermalGen in m.thermalGenerators},
        # mutable=True,
        units=u.MW / u.minutes,
        doc="Ramp up rates for each generator as a fraction of maximum generator output"
    )

    m.rampDownRates = Param(
        m.thermalGenerators,
        initialize={thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"]
                    for thermalGen in m.thermalGenerators},
        # mutable=True,
        units=u.MW / u.minutes,
        doc="Ramp down rates for each generator as a fraction of maximum generator output"
    )

    # for gen in m.thermalGenerators:
    #     print(f'm.rampUpRates[{gen}]={value(m.rampUpRates[gen])}')
    #     print(f'm.rampDownRates[{gen}]={value(m.rampDownRates[gen])}')
    
    # Matching for each generator to the region containing the bus at which the generator
    # is located
    m.gensAtRegion = {
        region: gen
        for region in m.regions
        for gen in m.generators
        if m.md.data["elements"]["bus"][m.md.data["elements"]["generator"][gen]["bus"]][
            "area"
        ]
        == region
    }

    #[ESR WIP: Declare fixed and operating costs here to avoid
    #multiple declarations of the same parameter. Set the value to 1
    #for now and then update it.]
    m.fixedCost = Param(
        m.generators,
        initialize={gen:1 for gen in m.generators},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Fixed operating costs"
    )
    m.varCost = Param(
        m.generators,
        initialize={gen:1 for gen in m.generators},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Variable costs"
    )

    # [ESR WIP: Declare and initialize curtailment and load shed costs
    # as parameters. These are re-calculated in
    # investment_stage_rule. Also, note that th original
    # "loadShedCost" was renamed "loadShedCostperCurtailment" to avoid
    # its repetition. ]
    m.curtailmentCost = Param(
        initialize=1,
        units=u.USD / (u.MW * u.hr),
        mutable=True,
        doc="Curtailment cost"
    )
    m.loadShedCostperCurtailment = Param(
        initialize=1000,
        units=u.USD / (u.MW * u.hr),
        mutable=True,
    )





def model_create_investment_stages(m, stages):
    """Creates investment blocks and linking constraints for GTEP model.
    Largely manages retirements and links operational units in a given investment stage
    to operational + installed - retired in the previous investment stage.

    :m: Pyomo model object
    :stages: Number of investment stages in planning horizon
    """

    # [ESR WIP: Add investment years]
    m.years = [2025, 2030, 2035]
    
    m.investmentStage = Block(m.stages, rule=investment_stage_rule)

    # Retirement/extension relationships over investment periods -- C&P'd
    # from the paper.  These are okay.
    # if len(m.stages) > 1:

    #     @m.Constraint(m.stages, m.thermalGenerators)
    #     def gen_retirement(m, stage, gen):
    #         return sum(
    #             m.investmentStage[t_2]
    #             .genInstalled[gen]
    #             .indicator_var.get_associated_binary()
    #             for t_2 in m.stages
    #             if t_2 <= stage - m.lifetimes[gen]
    #         ) <= sum(
    #             m.investmentStage[t_1]
    #             .genRetired[gen]
    #             .indicator_var.get_associated_binary()
    #             + m.investmentStage[t_1]
    #             .genExtended[gen]
    #             .indicator_var.get_associated_binary()
    #             for t_1 in m.stages
    #             if t_1 <= stage
    #         )

    # Linking generator investment status constraints
    @m.Constraint(m.stages, m.thermalGenerators)
    def gen_stats_link(m, stage, gen):
        return (
            m.investmentStage[stage]
            .genOperational[gen]
            .indicator_var.get_associated_binary()
            == m.investmentStage[stage - 1]
            .genOperational[gen]
            .indicator_var.get_associated_binary()
            + m.investmentStage[stage - 1]
            .genInstalled[gen]
            .indicator_var.get_associated_binary()
            - m.investmentStage[stage - 1]
            .genRetired[gen]
            .indicator_var.get_associated_binary()
            if stage != 1
            else Constraint.Skip
        )

    # Renewable generation (in MW) retirement relationships
    if len(m.stages) > 1:
        # [ESR WIP: Consider one stage = one year. This way lifetimes
        # units are consistent in the if statement below.]

        @m.Constraint(m.stages, m.renewableGenerators)
        def renewable_retirement(m, stage, gen):
            return sum(
                m.investmentStage[t_2].renewableInstalled[gen]
                for t_2 in m.stages
                # [ESR WIP: Only take the value of it so avoid errors
                # but always make sure the units makes sense.]
                if t_2 <= stage - value(m.lifetimes[gen])
            ) <= sum(
                m.investmentStage[t_1].renewableRetired[gen]
                + m.investmentStage[t_1].renewableExtended[gen]
                for t_1 in m.stages
                if t_1 <= stage
            )

    # Total renewable generation (in MW) operational at a given stage
    # is equal to what was operational and/or installed in the previous stage
    # less what was retired in the previous stage
    @m.Constraint(m.stages, m.renewableGenerators)
    def renewable_stats_link(m, stage, gen):
        return (
            m.investmentStage[stage].renewableOperational[gen]
            == m.investmentStage[stage - 1].renewableOperational[gen]
            + m.investmentStage[stage - 1].renewableInstalled[gen]
            - m.investmentStage[stage - 1].renewableRetired[gen]
            if stage != 1
            else Constraint.Skip
        )

    # If a gen is online at time t, it must have been online or installed at time t-1
    @m.LogicalConstraint(m.stages, m.thermalGenerators)
    def consistent_operation(m, stage, gen):
        return (
            m.investmentStage[stage]
            .genOperational[gen]
            .indicator_var.implies(
                m.investmentStage[stage - 1].genOperational[gen].indicator_var
                | m.investmentStage[stage - 1].genInstalled[gen].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a gen is online at time t, it must be online, extended, or retired at time t+1
    @m.LogicalConstraint(m.stages, m.thermalGenerators)
    def consistent_operation_future(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genOperational[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genOperational[gen].indicator_var
                | m.investmentStage[stage].genExtended[gen].indicator_var
                | m.investmentStage[stage].genRetired[gen].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # Retirement in period t-1 implies disabled in period t
    @m.LogicalConstraint(m.stages, m.thermalGenerators)
    def full_retirement(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genRetired[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genDisabled[gen].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a gen is disabled at time t-1, it must stay disabled  at time t
    @m.LogicalConstraint(m.stages, m.thermalGenerators)
    def consistent_disabled(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genDisabled[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genDisabled[gen].indicator_var
                | m.investmentStage[stage].genInstalled[gen].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a gen is extended at time t-1, it must stay extended or be retired at time t
    @m.LogicalConstraint(m.stages, m.thermalGenerators)
    def consistent_extended(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genExtended[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genExtended[gen].indicator_var
                | m.investmentStage[stage].genRetired[gen].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # Installation in period t-1 implies operational in period t
    @m.LogicalConstraint(m.stages, m.thermalGenerators)
    def full_investment(m, stage, gen):
        return (
            m.investmentStage[stage - 1]
            .genInstalled[gen]
            .indicator_var.implies(
                m.investmentStage[stage].genOperational[gen].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a branch is online at time t, it must have been online or installed at time t-1
    @m.LogicalConstraint(m.stages, m.transmission)
    def consistent_branch_operation(m, stage, branch):
        return (
            m.investmentStage[stage]
            .branchOperational[branch]
            .indicator_var.implies(
                m.investmentStage[stage - 1].branchOperational[branch].indicator_var
                | m.investmentStage[stage - 1].branchInstalled[branch].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a branch is online at time t, it must be online, extended, or retired at time t+1
    @m.LogicalConstraint(m.stages, m.transmission)
    def consistent_branch_operation_future(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchOperational[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchOperational[branch].indicator_var
                | m.investmentStage[stage].branchExtended[branch].indicator_var
                | m.investmentStage[stage].branchRetired[branch].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # Retirement in period t-1 implies disabled in period t
    @m.LogicalConstraint(m.stages, m.transmission)
    def full_branch_retirement(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchRetired[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchDisabled[branch].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a branch is disabled at time t-1, it must stay disabled or be installed at time t
    @m.LogicalConstraint(m.stages, m.transmission)
    def consistent_branch_disabled(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchDisabled[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchDisabled[branch].indicator_var
                | m.investmentStage[stage].branchInstalled[branch].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a branch is extended at time t-1, it must stay extended or be retired at time t
    @m.LogicalConstraint(m.stages, m.transmission)
    def consistent_branch_extended(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchExtended[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchExtended[branch].indicator_var
                | m.investmentStage[stage].branchRetired[branch].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # Installation in period t-1 implies operational in period t
    @m.LogicalConstraint(m.stages, m.transmission)
    def full_branch_investment(m, stage, branch):
        return (
            m.investmentStage[stage - 1]
            .branchInstalled[branch]
            .indicator_var.implies(
                m.investmentStage[stage].branchOperational[branch].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )
