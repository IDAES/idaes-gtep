# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from pyomo.environ import units as u

# from pyomo.gdp import *

from egret.data.model_data import ModelData
from egret.model_library.transmission.tx_utils import scale_ModelData_to_pu
from pyomo.common.timing import TicTocTimer
from pyomo.repn.linear import LinearRepnVisitor
import json
import numpy as np
import re

import math


from math import ceil
from gtep.config_options import (
    _get_model_config,
    _add_common_configs,
    _add_investment_configs,
)


# Define what a USD is for pyomo units purposes
# This will be set to a base year and we will do NPV calculations
# based on automatic pyomo unit transformations
u.load_definitions_from_strings(["USD = [currency]"])


####################################
########## New Work Here ###########
####################################

## TODO: Egret features


def data_update(investment_stage, storage_object, target_storage_object):
    pass


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
        num_reps=3,
        len_reps=24,
        num_commit=24,
        num_dispatch=4,
    ):
        """Initialize generation & expansion planning model object.

        :param stages: integer number of investment periods
        :param formulation: Egret stuff, to be filled
        :param data: full set of model data
        :param num_reps: integer number of representative periods per investment period
        :param len_reps: (for now integer) length of each representative period (in hours)
        :param num_commit: integer number of commitment periods per representative period
        :param num_dispatch: integer number of dispatch periods per commitment period
        :return: Pyomo model for full GTEP
        """

        self.stages = stages
        self.formulation = formulation
        self.data = data
        self.num_reps = num_reps
        self.len_reps = len_reps
        self.num_commit = num_commit
        self.num_dispatch = num_dispatch
        self.config = _get_model_config()
        self.timer = TicTocTimer()

        _add_common_configs(self.config)
        _add_investment_configs(self.config)

    def create_model(self):
        """Create concrete Pyomo model object associated with the ExpansionPlanningModel"""

        self.timer.tic("Creating GTEP Model")
        m = ConcreteModel()
        m.config = self.config
        m.rng = np.random.default_rng(seed=123186)

        ## TODO: checks for active/built/inactive/unbuilt/etc. gen
        ## NOTE: scale_ModelData_to_pu doesn't account for expansion data -- does it need to?
        if self.data is None:
            raise
        elif type(self.data.representative_data) is list:
            # If self.data is a list, it is a list of data for
            # representative periods
            m.data_list = self.data.representative_data
            ##TEXAS: testing this for proper scaling
            # for data in m.data_list:
            #     scale_ModelData_to_pu(data)
            m.md = m.data_list[0]
            m.data = self.data
        else:
            # If self.data is an Egret model data object, representative periods will just copy it unchanged
            m.data_list = None
            m.md = scale_ModelData_to_pu(self.data)
            m.formulation = self.formulation

        model_set_declaration(
            m, self.stages, rep_per=[i for i in range(1, self.num_reps + 1)]
        )
        m.representativePeriodLength = Param(
            m.representativePeriods, within=PositiveReals, default=24, units=u.hr
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
        m.commitmentPeriodLength = Param(within=PositiveReals, default=1, units=u.hr)
        # TODO: index by dispatch period? Certainly index by commitment period
        m.dispatchPeriodLength = Param(within=PositiveReals, default=0.25, units=u.hr)

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


def add_investment_variables(b, investment_stage):
    """Add variables to investment stage block.

    :param b: Investment block
    :param investment_stage: Investment stage index or name
    :return: None
    """

    m = b.model()
    b.investmentStage = investment_stage

    # Thermal generator disjuncts (operational, installed, retired, disabled, extended)
    @b.Disjunct(m.thermalGenerators)
    def genOperational(disj, gen):
        return

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
    
    """ Energy Storage Investment States"""
    # Energy Storage (Battery) disjuncts. For now mimicking thermal generators
    @b.Disjunct(m.batteryStorageSystems)
    def batOperational(disj, bat):
        return

    @b.Disjunct(m.batteryStorageSystems)
    def batInstalled(disj, bat):
        return

    @b.Disjunct(m.batteryStorageSystems)
    def batRetired(disj, bat):
        return

    @b.Disjunct(m.batteryStorageSystems)
    def batDisabled(disj, bat):
        return

    @b.Disjunct(m.batteryStorageSystems)
    def batExtended(disj, bat):
        return

    @b.Disjunction(m.batteryStorageSystems)
    def batInvestStatus(disj, bat):
        return [
            disj.batOperational[bat],
            disj.batInstalled[bat],
            disj.batRetired[bat],
            disj.batDisabled[bat],
            disj.batExtended[bat],
        ]
    

    if m.config["transmission"]:
        # Line disjuncts. For now mimicking thermal generator disjuncts, though different states may need to be defined
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
        m.renewableGenerators, within=NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableInstalled = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableRetired = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableExtended = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0, units=u.MW
    )
    b.renewableDisabled = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0, units=u.MW
    )

    # Track and accumulate costs and penalties
    b.quotaDeficit = Var(within=NonNegativeReals, initialize=0, units=u.MW * u.hr)
    # b.expansionCost = Var(within=Reals, initialize=0, units=u.USD)
    b.renewableCurtailmentInvestment = Var(
        within=NonNegativeReals, initialize=0, units=u.USD
    )


def add_investment_constraints(b, investment_stage):
    """Add standard inequalities (i.e., those not involving disjunctions) to investment stage block."""

    m = b.model()

    for gen in m.thermalGenerators:
        if (
            m.md.data["elements"]["generator"][gen]["in_service"] == False
            and investment_stage == 1
        ):
            b.genOperational[gen].indicator_var.fix(False)
            # b.genDisabled[gen].binary_indicator_var.fix(1)
        elif (
            m.md.data["elements"]["generator"][gen]["in_service"] == True
            and investment_stage == 1
        ):
            b.genInstalled[gen].indicator_var.fix(True)
    
    """ Energy Storage: Fixing In-Service batteries initial investment state based on input"""
    for bat in m.batteryStorageSystems:
        if (
            m.md.data["elements"]["storage"][bat]["in_service"] == False
            and investment_stage == 1
        ):
            b.batDisabled[bat].indicator_var.fix(True)
        elif (
            m.md.data["elements"]["storage"][bat]["in_service"] == True
            and investment_stage == 1
        ):
            b.batInstalled[bat].indicator_var.fix(True)
            # Also initialize storage level 

    for branch in m.transmission:
            b.genOperational[gen].indicator_var.fix(True)
            # b.genInstalled[gen].binary_indicator_var.fix(1)
    # for gen in m.thermalGenerators:
    #     if (
    #         m.md.data["elements"]["generator"][gen]["lifetime"] == 1
    #         and investment_stage == 2
    #     ):
    #         b.genRetired[gen].indicator_var.fix(True)
    for gen in m.renewableGenerators:
        if (
            m.md.data["elements"]["generator"][gen]["in_service"] == False
            and investment_stage == 1
        ):
            # print(gen)
            b.renewableOperational[gen].fix(0)
            #b.renewableDisabled[gen].fix(m.renewableCapacity[gen])
        elif (
            m.md.data["elements"]["generator"][gen]["in_service"] == True
            and investment_stage == 1
        ):
            b.renewableOperational[gen].fix(m.renewableCapacity[gen])

    if m.config["transmission"]:
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
                m.renewableCapacity[gen]
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
    #             m.renewableCapacity[gen]
    #             * b.genInstalled[gen].indicator_var.get_associated_binary()
    #             for gen in m.renewableGenerators & m.gensAtRegion[region]
    #         )
    #         <= b.maxRenewableInvestment[region]
    #         if m.renewableGenerators & m.gensAtRegion[region]
    #         else Constraint.Skip
    #     )

    ## NOTE: The following constraints can be split into rep_per and invest_stage components if desired

    # Operating costs for investment period
    @b.Expression()
    def operatingCostInvestment(b):
        operatingCostRepresentative = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                operatingCostRepresentative += (
                    m.weights[rep_per]
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .operatingCostCommitment
                )
                
        return m.investmentFactor[investment_stage] * operatingCostRepresentative
    
    """ Energy Storage Cost """            
    @b.Expression()
    def storageCostInvestment(b):
        storageCostRepresentative = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                storageCostRepresentative += (
                    m.weights[rep_per]
                    * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .storageCostCommitment
                    )
                
        return m.investmentFactor[investment_stage] * storageCostRepresentative

    # Investment costs for investment period
    ## FIXME: investment cost definition needs to be revisited AND possibly depends on
    ## data format.  It is _rare_ for these values to be defined at all, let alone consistently.
    @b.Expression()
    def investment_cost(b):
        baseline_cost = sum(
                m.generatorInvestmentCost[gen]
                * m.capitalMultiplier[gen]
                * b.genInstalled[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
                ) + sum(
                m.generatorInvestmentCost[gen]
                * m.capitalMultiplier[gen]
                * b.renewableInstalled[gen]
                for gen in m.renewableGenerators
            ) + sum(
                m.generatorInvestmentCost[gen]
                * m.extensionMultiplier[gen]
                * b.genExtended[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            ) + sum(
                m.generatorInvestmentCost[gen]
                * m.extensionMultiplier[gen]
                * b.renewableExtended[gen]
                for gen in m.renewableGenerators
            ) + sum(
                m.batteryInvestmentCost[bat]
                * m.batteryCapitalMultiplier[bat]
                * b.batInstalled[bat].indicator_var.get_associated_binary()
                for bat in m.batteryStorageSystems
            ) + sum(
                m.batteryInvestmentCost[bat]
                * m.batteryExtensionMultiplier[bat]
                * b.batExtended[bat].indicator_var.get_associated_binary()
                for bat in m.batteryStorageSystems
            ) + sum(
                m.generatorInvestmentCost[gen]
                * m.retirementMultiplier[gen]
                * b.renewableRetired[gen]
                for gen in m.renewableGenerators
            ) + sum(
                m.generatorInvestmentCost[gen]
                * m.retirementMultiplier[gen]
                * b.genRetired[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            ) + sum(
                m.branchInvestmentCost[branch]
                * m.branchCapitalMultiplier[branch]
                * b.branchInstalled[branch].indicator_var.get_associated_binary()
                for branch in m.transmission
            ) + sum(
                m.branchInvestmentCost[branch]
                * m.branchExtensionMultiplier[branch]
                * b.branchExtended[branch].indicator_var.get_associated_binary()
                for branch in m.transmission
            )
        return m.investmentFactor[investment_stage] * baseline_cost

    # Curtailment penalties for investment period
    @b.Constraint()
    def renewable_curtailment_cost(b):
        renewableCurtailmentRep = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                renewableCurtailmentRep += (
                    m.weights[rep_per]
                    * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .renewableCurtailmentCommitment
                )
        return (
            b.renewableCurtailmentInvestment
            == m.investmentFactor[investment_stage] * renewableCurtailmentRep
        )
    
    """ 
    # Initial, untested attempt for enforcing identical storage level at 
    # beginning and end of representative periods
    # Need to update to use init and end batteryChargeLevel?
    """
    # @b.Constraint(b.representativePeriods, m.batteryStorageSystems)
    # def consistent_battery_charge_level_commitment(b, rep_per, bat):

    #     return (
            
    #             b.representativePeriod[rep_per]
    #             .commitmentPeriod[
    #                 b.representativePeriod[rep_per]
    #                 .commitmentPeriods.first()
    #                 ]
    #                 .dispatchPeriod[
    #                     b.representativePeriod[rep_per]
    #                     .commitmentPeriod[
    #                         b.representativePeriod[rep_per]
    #                         .commitmentPeriods.first()
    #                         ]
    #                         .dispatchPeriods.first()
    #                     ]
    #                     .batteryChargeLevel[bat]
    #               == 
    #               b.representativePeriod[rep_per]
    #               .commitmentPeriod[
    #                   b.representativePeriod[rep_per]
    #                   .commitmentPeriods.last()
    #                   ]
    #                   .dispatchPeriod[
    #                       b.representativePeriod[rep_per]
    #                       .commitmentPeriod[
    #                           b.representativePeriod[rep_per]
    #                           .commitmentPeriods.last()
    #                           ]
    #                           .dispatchPeriods.last()
    #                       ]
    #                       .batteryChargeLevel[bat]
    #     )
    # @b.Constraint()
    # def renewable_curtailment_cost(b):
    #     renewableCurtailmentRep = 0
    #     for rep_per in b.representativePeriods:
    #         for com_per in b.representativePeriod[rep_per].commitmentPeriods:
    #             renewableCurtailmentRep += (
    #                 m.weights[rep_per]
    #                 * m.commitmentPeriodLength
    #                 * b.representativePeriod[rep_per]
    #                 .commitmentPeriod[com_per]
    #                 .renewableCurtailmentCommitment
    #             )
    #     return (
    #         b.renewableCurtailmentInvestment
    #         == m.investmentFactor[investment_stage] * renewableCurtailmentRep
    #     )

    ## NOTE: Constraint (13) in the reference paper
    # Minimum per-stage renewable generation requirement
    if m.config["include_investment"]:

        @b.Constraint()
        def renewable_generation_requirement(b):
            renewableSurplusRepresentative = 0
            ## TODO: preprocess loads for the appropriate sum here
            ed = 0
            for rep_per in b.representativePeriods:
                for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                    renewableSurplusRepresentative += (
                        m.weights[rep_per]
                        * b.representativePeriod[rep_per]
                        .commitmentPeriod[com_per]
                        .renewableSurplusCommitment
                    )
            return (
                renewableSurplusRepresentative + b.quotaDeficit
                >= m.renewableQuota[investment_stage] * ed
            )


def add_dispatch_variables(b, dispatch_period):
    """Add dispatch-associated variables to representative period block."""

    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    # Define bounds on thermal generator active generation
    def thermal_generation_limits(b, thermalGen):
        return (0, m.thermalCapacity[thermalGen])

    b.thermalGeneration = Var(
        m.thermalGenerators,
        domain=NonNegativeReals,
        bounds=thermal_generation_limits,
        initialize=0,
        units=u.MW,
    )
    
    """ Battery Parameters """
    def battery_capacity_limits(b, bat):
        return (m.minBatteryChargeLevel[bat], m.batteryCapacity[bat]) # The lower bound should be > 0 - data input

    # TODO: Note that this does not fix initial battery capacity at the first dispatch period - need to adjust constraint
    def init_battery_capacity(b, bat):
        return m.initBatteryChargeLevel[bat]

    b.batteryChargeLevel = Var(
        m.batteryStorageSystems,
        domain=NonNegativeReals,
        bounds=battery_capacity_limits,
        initialize=init_battery_capacity,
        units=u.MW,
    )
    
    # Define bounds on charging/discharging capability. Note that constraints 
    # enforce that there are min & max charge/discharge levels if the bat is in
    # the charging or discharging state
    def battery_charge_limits(b, bat):
        return(0, m.chargeMax[bat])
    
    def battery_discharge_limits(b, bat):
        return(0, m.dischargeMax[bat])

    b.batteryCharged = Var(
        m.batteryStorageSystems,
        domain=NonNegativeReals,
        bounds=battery_charge_limits,
        initialize=0,
        units=u.MW,
    )
    
    b.batteryDischarged = Var(
        m.batteryStorageSystems,
        domain=NonNegativeReals,
        bounds=battery_discharge_limits,
        initialize=0,
        units=u.MW,
    )

    # Define bounds on renewable generator active generation
    def renewable_generation_limits(b, renewableGen):
        return (0, m.renewableCapacity[renewableGen])

    b.renewableGeneration = Var(
        m.renewableGenerators,
        domain=NonNegativeReals,
        bounds=renewable_generation_limits,
        initialize=0,
        units=u.MW,
    )

    # Fix hydro when we don't have a time series
    for gen in m.renewableGenerators:
        if m.md.data["elements"]["generator"][gen]['fuel'] == 'H':
            b.renewableGeneration[gen].fix(m.renewableCapacity[gen])
            # print(m.renewableCapacity[gen])
            # import sys
            # sys.exit()

    # Define bounds on renewable generator curtailment
    def curtailment_limits(b, renewableGen):
        return (0, m.renewableCapacity[renewableGen])

    b.renewableCurtailment = Var(
        m.renewableGenerators,
        domain=NonNegativeReals,
        bounds=curtailment_limits,
        initialize=0,
        units=u.MW,
    )

    # Per generator surplus
    @b.Expression(m.renewableGenerators)
    def renewableGenerationSurplus(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] - b.renewableCurtailment[renewableGen]
        )

    # Per generator curtailment cost
    @b.Expression(m.renewableGenerators)
    def renewableCurtailmentCost(b, renewableGen):
        return b.renewableCurtailment[renewableGen] * m.curtailmentCost

    # Per generator cost
    ## TEXAS: added varCost below
    @b.Expression(m.thermalGenerators)
    def generatorCost(b, gen):
        return b.thermalGeneration[gen] * i_p.fuelCost[gen]
        return b.thermalGeneration[gen] * (i_p.fuelCost[gen] + i_p.varCost[gen])

    # * b.dispatchLength

    # Load shed per bus
    b.loadShed = Var(m.buses, domain=NonNegativeReals, initialize=0, units=u.MW * u.hr)

    # Per bus load shed cost
    @b.Expression(m.buses)
    def loadShedCost(b, bus):
        return b.loadShed[bus] * m.loadShedCost
    
    """ Per-Battery Operational cost variables"""
    @b.Expression(m.batteryStorageSystems)
    def batteryChargingCost(b, bat):
        return b.batteryCharged[bat] * m.chargingCost[bat] 

    # JSC addn Per Battery Discharging Cost
    @b.Expression(m.batteryStorageSystems)
    def batteryDischargingCost(b, bat):
        return b.batteryDischarged[bat] * m.dischargingCost[bat] 


    # Track total dispatch values and costs
    b.renewableSurplusDispatch = sum(b.renewableGenerationSurplus.values())

    b.generationCostDispatch = sum(b.generatorCost.values())

    b.loadShedCostDispatch = sum(b.loadShedCost.values())

    b.curtailmentCostDispatch = sum(b.renewableCurtailmentCost.values())

    b.operatingCostDispatch = (
        b.generationCostDispatch + b.loadShedCostDispatch + b.curtailmentCostDispatch
    )
    
    """ Per-Battery Operational costs """
    b.chargingCostDispatch = sum(b.batteryChargingCost.values())
    
    b.dischargingCostDispatch = sum(b.batteryDischargingCost.values())
    
    b.storageCostDispatch = (
        b.chargingCostDispatch + b.dischargingCostDispatch
        )
    b.operatingCostDispatch = b.generationCostDispatch + b.loadShedCostDispatch
    # + b.curtailmentCostDispatch
    # )

    b.renewableCurtailmentDispatch = sum(
        b.renewableCurtailment[gen] for gen in m.renewableGenerators
    )

    # Define bounds on transmission line capacity - restrictions on flow over
    # uninvested lines are enforced in a disjuction below
    ## TEXAS
    def power_flow_limits(b, branch):
        return (-m.transmissionCapacity[branch] * 8, m.transmissionCapacity[branch])

    # NOTE: this is an abuse of units and needs to be fixed for variable temporal resolution
    b.powerFlow = Var(
        m.transmission,
        domain=Reals,
        bounds=power_flow_limits,
        initialize=0,
        units=u.MW,
    )

    @b.Disjunct(m.transmission)
    def branchInUse(disj, branch):
        b = disj.parent_block()

        # Voltage angle
        def bus_angle_bounds(disj, bus):
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

        ##FIXME: we need this for all buses all the time
        disj.busAngle = Var(
            disj.branch_buses, domain=Reals, initialize=0, bounds=bus_angle_bounds
        )

        # Voltage angle
        def delta_bus_angle_bounds(disj, bus):
            return (-math.pi / 6, math.pi / 6)

        # Rule for maximum bus angle discrepancy
        def delta_bus_angle_rule(disj):
            fb = m.transmission[branch]["from_bus"]
            tb = m.transmission[branch]["to_bus"]
            return disj.busAngle[tb] - disj.busAngle[fb]

        ##FIXME: we can just add this as a constraint rather than a variable, this is weird
        # @KyleSkolfield - I think this var is unused and commented it out, can we delete?
        disj.deltaBusAngle = Var(
            domain=Reals, bounds=delta_bus_angle_bounds, rule=delta_bus_angle_rule
        )

        if m.config["flow_model"] == "DC":

            @disj.Constraint()
            def dc_power_flow(disj):
                fb = m.transmission[branch]["from_bus"]
                tb = m.transmission[branch]["to_bus"]
                reactance = m.md.data["elements"]["branch"][branch]["reactance"]
                if (
                    m.md.data["elements"]["branch"][branch]["branch_type"]
                    == "transformer"
                ):
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
                )

    @b.Disjunct(m.transmission)
    def branchNotInUse(disj, branch):

        # Fixing power flow to 0 and not creating bus angle variables for
        # branches that are not in use.
        @disj.Constraint()
        def dc_power_flow(disj):
            return b.powerFlow[branch] == 0

        return

    # Branches are either in-use or not. This disjunction may provide the
    # basis for transmission switching in the future.
    @b.Disjunction(m.transmission)
    def branchInUseStatus(disj, branch):
        return [disj.branchInUse[branch], disj.branchNotInUse[branch]]

    if m.config["transmission"]:
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

        ##FIXME: this logic isn't true.  remove when con fig fixes switching.
        ##FIXME: replace with disabled/retired \implies not in use
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

    # Define bounds on thermal generator spinning reserve supply
    def spinning_reserve_limits(b, thermalGen):
        return (
            0,
            m.spinningReserveFraction[thermalGen] * m.thermalCapacity[thermalGen],
        )

    b.spinningReserve = Var(
        m.thermalGenerators,
        domain=NonNegativeReals,
        bounds=spinning_reserve_limits,
        initialize=0,
        units=u.MW * u.hr,
    )

    # Define bounds on thermal generator quickstart reserve supply
    def quickstart_reserve_limits(b, thermalGen):
        return (
            0,
            m.quickstartReserveFraction[thermalGen] * m.thermalCapacity[thermalGen],
        )

    b.quickstartReserve = Var(
        m.thermalGenerators,
        domain=NonNegativeReals,
        bounds=quickstart_reserve_limits,
        initialize=0,
        units=u.MW * u.hr,
    )


def add_dispatch_constraints(b, disp_per):
    """Add dispatch-associated inequalities to representative period block."""
    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    # for key in m.loads.keys():
    #     m.loads[key] *= max(0, m.rng.normal(0.5, 0.2))

    # Energy balance constraint
    @b.Constraint(m.buses)
    def flow_balance(b, bus):
        balance = 0
        load = c_p.loads.get(bus) or 0
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
        # JSC addn
        batts = [
            bat for bat in m.batteryStorageSystems 
            if m.md.data["elements"]["storage"][bat]["bus"] == bus
        ]
        balance -= sum(b.powerFlow[i] for i in end_points)
        balance += sum(b.powerFlow[i] for i in start_points)
        balance += sum(b.thermalGeneration[g] for g in gens if g in m.thermalGenerators)
        """ Battery Storage added to flow balance constraint """
        balance += sum(b.batteryDischarged[bt] for bt in batts)
        balance -= sum(b.batteryCharged[bt] for bt in batts)
        balance += sum(
            b.renewableGeneration[g] for g in gens if g in m.renewableGenerators
        )
        balance -= load
        balance += b.loadShed[bus]
        return balance == 0

    # Capacity factor constraint
    # NOTE: In comparison to reference work, this is *per renewable generator*
    # JKS - charging costs from non-colocated plants?
    @b.Constraint(m.renewableGenerators)
    def capacity_factor(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] + b.renewableCurtailment[renewableGen]
            == m.renewableCapacity[renewableGen]
        )

    ## TODO: (@jkskolf) add renewableExtended to this and anywhere else
    @b.Constraint(m.renewableGenerators)
    def operational_renewables_only(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen]
            <= i_p.renewableInstalled[renewableGen]
            + i_p.renewableOperational[renewableGen]
            + i_p.renewableExtended[renewableGen]
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

        # Minimum operating Limits
        @disj.Constraint(b.dispatchPeriods)
        def operating_limit_min(d, dispatchPeriod):
            return (
                m.thermalMin[generator]
                <= b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
            )

        # Maximum operating limits
        ##FIXME: don't need this constraint
        @disj.Constraint(b.dispatchPeriods)
        def operating_limit_max(d, dispatchPeriod):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                + b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.thermalCapacity[generator]
            )

        # Ramp up limit constraints for fully on generators
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
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

        # Ramp down limit constraints for fully on generators
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        def ramp_down_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                <= m.rampDownRates[generator]
                * b.dispatchPeriod[dispatchPeriod].periodLength
                * m.thermalCapacity[generator]
                if dispatchPeriod != 1
                else Constraint.Skip
            )

        # Maximum spinning reserve constraint
        ##NOTE: maxSpiningReserve is a percentage of thermalCapacity
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        def max_spinning_reserve(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].spinningReserve[generator]
                <= m.maxSpinningReserve[generator] * m.thermalCapacity[generator]
            )

        ##FIXME: add quick start reserve = 0

    @b.Disjunct(m.thermalGenerators)
    def genStartup(disj, generator):
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

        # Ramp up constraints for generators starting up
        ## TODO: is this max necessary? I would like to remove
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        def ramp_up_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                <= max(
                    m.thermalMin[generator],
                    m.rampUpRates[generator]
                    * b.dispatchPeriod[dispatchPeriod].periodLength,
                )
                ##FIXME: I don't think this parenthesis is correct -- thermal capacity should go inside with the second term.  Or do I need to make this two constraints?
                * m.thermalCapacity[generator]
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

        ## RMA:
        ## We may need to turn off ramp down constraints for feasibility purposes
        ## We will need to think about this for future work, but commenting this out 
        ## is probably fine for the purposes of this paper 

        # Ramp down constraints for generators shutting down
        ## FIXME: uncomment out this stuff
        # @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        # def ramp_down_limits(disj, dispatchPeriod, generator):
        #     return (
        #         b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
        #         - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
        #         <= max(
        #             m.thermalMin[generator],
        #             m.rampDownRates[generator]
        #             * b.dispatchPeriod[dispatchPeriod].periodLength,
        #         )
        #         ##FIXME: I don't think this parenthesis is correct -- thermal capacity should go inside with the second term.  Or do I need to make this two constraints?
        #         * m.thermalCapacity[generator]
        #         if dispatchPeriod != 1
        #         else Constraint.Skip
        #     )

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
        ##FIXME: This isn't needed.  instead we need to set spinning reserve to 0.
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
            
    """
    Create constraints within disjunctions on battery storage commitment (charging/discharging/off)
    """
    
    """
    Battery Discharging Constraints
    """
    @b.Disjunct(m.batteryStorageSystems)
    def batDischarging(disj, bat):
        # operating limits
        ## NOTE: Reminder: thermalMin is a percentage of thermalCapacity
        b = disj.parent_block()

        # Minimum operating Limits if storage unit is on
        @disj.Constraint(b.dispatchPeriods)
        def discharge_limit_min(d, disp_per):
            return (
                m.dischargeMin[bat] # Assuming dischargeMin is an absolute value (MW)
                <= b.dispatchPeriod[disp_per].batteryDischarged[bat]
            )

        # Maximum operating limits
        @disj.Constraint(b.dispatchPeriods)
        def discharge_limit_max(d, disp_per):
            return (
                b.dispatchPeriod[disp_per].batteryDischarged[bat]
                <= m.dischargeMax[bat]
            )
        

        # Ramp up limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def discharge_ramp_up_limits(disj, disp_per):
            return (
                b.dispatchPeriod[disp_per].batteryDischarged[bat]
                - b.dispatchPeriod[disp_per - 1].batteryDischarged[bat]
                <= m.batteryDischargingRampUpRates[bat] # battery ramp rates are currently absolute values
                if disp_per != 1
                else Constraint.Skip
            )

        # Ramp down limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def discharge_ramp_down_limits(disj, disp_per):
            return (
                b.dispatchPeriod[disp_per - 1].batteryDischarged[bat]
                - b.dispatchPeriod[disp_per].batteryDischarged[bat]
                <= m.batteryDischargingRampDownRates[bat] # battery ramp rates are currently absolute values
                if disp_per != 1
                else Constraint.Skip 
            )
        
        # Force no charge when discharging
        @disj.Constraint(b.dispatchPeriods)
        def no_charge(disj, disp_per):
            return b.dispatchPeriod[disp_per].batteryCharged[bat] <= 0
        
        # Batteries that are charging both gain and lose energy
        @disj.Constraint(b.dispatchPeriods)
        def discharging_battery_storage_balance(disj, disp_per):
            return ( 
                    b.dispatchPeriod[disp_per].batteryChargeLevel[bat] == 
                    m.batteryRetentionRate[bat]*b.dispatchPeriod[disp_per-1].batteryChargeLevel[bat] -
                    b.dispatchPeriod[disp_per].batteryDischarged[bat]
                    if disp_per != 1
                    else Constraint.Skip
                ) 

    """
    Battery Charging Constraints
    """
    @b.Disjunct(m.batteryStorageSystems)
    def batCharging(disj, bat):
        b = disj.parent_block()

        @disj.Constraint(b.dispatchPeriods)
        def charge_limit_min(d, disp_per):
            return (
                m.chargeMin[bat] # Assuming chargeMin is an absolute value (MW)
                <= b.dispatchPeriod[disp_per].batteryCharged[bat]
            )

        # Maximum operating limits
        @disj.Constraint(b.dispatchPeriods)
        def charge_limit_max(d, disp_per):
            return (
                b.dispatchPeriod[disp_per].batteryCharged[bat]
                <= m.chargeMax[bat]
            )
        

        # Ramp up limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def charge_ramp_up_limits(disj, disp_per):
            return (
                b.dispatchPeriod[disp_per].batteryCharged[bat]
                - b.dispatchPeriod[disp_per - 1].batteryCharged[bat]
                <= m.batteryChargingRampUpRates[bat] # battery ramp rates are currently absolute values
                if disp_per != 1
                else Constraint.Skip
            )

        # Ramp down limit constraints for fully on bats
        @disj.Constraint(b.dispatchPeriods)
        def charge_ramp_down_limits(disj, disp_per):
            return (
                b.dispatchPeriod[disp_per - 1].batteryCharged[bat]
                - b.dispatchPeriod[disp_per].batteryCharged[bat]
                <= m.batteryChargingRampDownRates[bat] # battery ramp rates are currently absolute values
                if disp_per != 1
                else Constraint.Skip 
            )

        @disj.Constraint(b.dispatchPeriods)
        def no_discharge(disj, disp_per):
            return b.dispatchPeriod[disp_per].batteryDischarged[bat] <= 0
        
        # Batteries that are charging both gain and lose energy
        @disj.Constraint(b.dispatchPeriods)
        def charging_battery_storage_balance(disj, disp_per):
            return ( 
                    b.dispatchPeriod[disp_per].batteryChargeLevel[bat] == 
                    m.batteryRetentionRate[bat]*b.dispatchPeriod[disp_per-1].batteryChargeLevel[bat] +
                    m.batteryChargingEfficiency[bat]*b.dispatchPeriod[disp_per].batteryCharged[bat]
                    if disp_per != 1
                    else Constraint.Skip 
                ) # @JKS Evaluate if we need charging efficiency in this eqn and/or in flow balance


    """
    Battery Off Constraints
    """
    @b.Disjunct(m.batteryStorageSystems)
    def batOff(disj, bat):
        b = disj.parent_block()

        # If battery is off, it is not discharging in terms of sending energy 
        # to the grid
        @disj.Constraint(b.dispatchPeriods)
        def no_discharge(disj, disp_per):
            return b.dispatchPeriod[disp_per].batteryDischarged[bat] == 0

        # Batteries that are off cannot charge
        @disj.Constraint(b.dispatchPeriods)
        def no_charge(disj, disp_per):
            return b.dispatchPeriod[disp_per].batteryCharged[bat] == 0
        
        # Batteries that are off still lose energy, and none goes to the grid
        @disj.Constraint(b.dispatchPeriods)
        def off_batteries_lose_storage(disj, disp_per):
            return ( 
                    b.dispatchPeriod[disp_per].batteryChargeLevel[bat] == 
                    m.batteryRetentionRate[bat]*b.dispatchPeriod[disp_per-1].batteryChargeLevel[bat]
                    if disp_per != 1
                    else Constraint.Skip
                )
    
    # Batteries are exclusively either Charging, Discharging, or Off
    @b.Disjunction(m.batteryStorageSystems)
    def batStatus(disj, bat):
        return [
            disj.batCharging[bat],
            disj.batDischarging[bat],
            disj.batOff[bat],
        ]
    

    # bats cannot be committed unless they are operational or just installed
    @b.LogicalConstraint(m.batteryStorageSystems)
    def commit_active_batts_only(b, bat):
        return lor(
            b.batCharging[bat].indicator_var,
            b.batDischarging[bat].indicator_var
        ).implies(
            lor(
                i_p.batOperational[bat].indicator_var,
                i_p.batInstalled[bat].indicator_var,
                i_p.batExtended[bat].indicator_var,
            )
        )
                    


def add_commitment_constraints(
    b,
    comm_per,
):
    
    ## RMA:
    ## Here's where we'll fix commitment status to Off for thermal generators in outaged areas
    ## And we'll fix renewableGeneration to 0 for the dispatch periods here

    ## Note: I am guessing that the may 20 date is the .05 cutoff and the may 24 date is the .2 cutoff
    # we're going to lazily hard code this for now because oh lawd deadlines
    target_month = 5
    target_day = 20
    if r_p.month == target_month and r_p.day == target_day:
        # okay now for hourly
        current_hour_bus_outage_list = m.data.bus_hours[m.data.bus_hours["hour"] == b.commitmentPeriod - 1]
        bus_outages = current_hour_bus_outage_list.to_dict('list')["Bus Number"]
        for gen in m.thermalGenerators:
            if m.md.data["elements"]["generator"][gen]["bus"] in bus_outages:
                b.genOff[gen].indicator_var.fix(True)
        for gen in m.renewableGenerators:
            if m.md.data["elements"]["generator"][gen]["bus"] in bus_outages:
                b.dispatchPeriod[1].renewableGeneration[gen].fix(0)
        



def add_commitment_constraints(b, comm_per):
    """Add commitment-associated disjunctions and constraints to representative period block."""
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    # Define total renewable surplus/deficit for commitment block
    @b.Expression()
    def renewableSurplusCommitment(b):
        return sum(
            m.dispatchPeriodLength * b.dispatchPeriod[disp_per].renewableSurplusDispatch
            for disp_per in b.dispatchPeriods
        )

    # Define total operating costs for commitment block
    ## TODO: Replace this constraint with expressions using bounds transform
    ## NOTE: expressions are stored in gtep_cleanup branch
    ## costs considered need to be re-assessed and account for missing data

    # fixed cost units are WEIRD
    fixed_cost_coefs = 1000 / (5 * 8760)

    @b.Expression()
    def operatingCostCommitment(b):
        return (
            sum(
                b.dispatchPeriod[disp_per].operatingCostDispatch
                for disp_per in b.dispatchPeriods
            )
            + sum(
                i_p.fixedCost[gen]
                * b.commitmentPeriodLength
                * (
                    b.genOn[gen].indicator_var.get_associated_binary()
                    + b.genShutdown[gen].indicator_var.get_associated_binary()
                    + b.genStartup[gen].indicator_var.get_associated_binary()
                )
                for gen in m.thermalGenerators
            )
            ## FIXME: how do we do assign fixed operating costs to renewables; flat per location or per MW
            ## TEXAS: doing something wacky with those momentarily
            ## TEXAS FOR REAL WHAT LSAKDFJALWKERJ
            # + sum(
            #     i_p.fixedCost[gen]
            #     * b.commitmentPeriodLength
            #     * (
            #         i_p.renewableOperational[gen]
            #         + i_p.renewableInstalled[gen]
            #         + i_p.renewableExtended[gen]
            #     )
            #     # * m.renewableCapacity[gen]
            #     for gen in m.renewableGenerators
            # )
            + sum(
                m.startupCost[gen]
                * b.genStartup[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
        )
    
    # Define total storage costs for commitment block
    ## TODO: Replace this constraint with expressions using bounds transform
    ## NOTE: expressions are stored in gtep_cleanup branch
    ## costs considered need to be re-assessed and account for missing data
    """ Compute Battery Storage cost per dispatch period"""
    @b.Expression()
    def storageCostCommitment(b):
        return (
            sum(
                ## FIXME: update test objective value when this changes; ready to uncomment
                # (m.dispatchPeriodLength / 60) *
                b.dispatchPeriod[disp_per].storageCostDispatch
                for disp_per in b.dispatchPeriods
            )
        )

    # Define total curtailment for commitment block
    @b.Expression()
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
    m.renewableCapacity = {
        renewableGen: (
            0
            if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float
            else m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"][
                commitment_period - 1
            ]
        )
        for renewableGen in m.renewableGenerators
    }
    
    ## TEXAS: solar is too small what's up with that?
    for gen in m.renewableGenerators:
        if m.md.data['elements']['generator'][gen]['fuel'] == 'S':
            m.renewableCapacity[gen] *= 10

    ## TODO: Redesign load scaling and allow nature of it as argument
    # Demand at each bus
    b.load_scaling = r_p.load_scaling[r_p.load_scaling["hour"] == b.commitmentPeriod]
    # print(b.load_scaling.head())

    if m.config["scale_texas_loads"]:
        false_loads = []
        for load in m.md.data['elements']['load']:
            if type(m.md.data['elements']['load'][load]) == float:
                false_loads.append(load)
        for load in false_loads:
            del m.md.data['elements']['load'][load]
            # del m.loads[load]
        # print(m.loads)
        b.loads = {
            m.md.data["elements"]["load"][load_n]["bus"]: m.md.data["elements"]["load"][
                load_n
            ]["p_load"]["values"][commitment_period - 1]
            * b.load_scaling[m.md.data["elements"]["load"][load_n]["zone"]].iloc[0]
            for load_n in m.md.data["elements"]["load"]
        }
        # print(m.loads)
        for key, val in b.loads.items():
            # print(f"{key=}")
            # print(f"{val=}")
            b.loads[key] *= 1
            # for i, v in enumerate(val['values']):
            #     val['values'][i] *= 1/300
        # print(sum(m.loads.values()))
        # for key, val in m.thermalCapacity.items():
        #     if m.md.data["elements"]["generator"][key]["fuel"] == 'G':
        #         m.thermalCapacity[key] *= 1/10
                

    # if m.config["scale_loads"]:
    #     temp_scale = 3
    #     temp_scale = 10

    #     m.loads = {
    #         m.md.data["elements"]["load"][load_n]["bus"]: (
    #             temp_scale
    #             * (
    #                 1
    #                 + (temp_scale + i_p.investmentStage) / (temp_scale + len(m.stages))
    #             )
    #         )
    #         * m.md.data["elements"]["load"][load_n]["p_load"]["values"][
    #             commitment_period - 1
    #         ]
    #         for load_n in m.md.data["elements"]["load"]
    #     }

    # else:
    #     m.loads = {
    #         m.md.data["elements"]["load"][load_n]["bus"]: m.md.data["elements"]["load"][
    #             load_n
    #         ]["p_load"]["values"][commitment_period - 1]
    #         for load_n in m.md.data["elements"]["load"]
    #     }

    ## TODO: This feels REALLY inelegant and bad.
    ## TODO: Something weird happens if I say periodLength has a unit
    for period in b.dispatchPeriods:
        b.dispatchPeriod[period].periodLength = Param(within=PositiveReals, default=1)
        add_dispatch_variables(b.dispatchPeriod[period], period)

    ## TODO: if commitment is neglected but dispatch is still desired, pull something different here? or simply don't enforce linked commitment constraints?
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
    if m.config["include_commitment"]:
        ##FIXME this needs to be updated for variable length commitment periods
        ## do this by (pre) processing the set of commitment periods for req_shutdown_periods
        @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
        def consistent_commitment_shutdown(b, commitmentPeriod, thermalGen):
            req_shutdown_periods = ceil(
                1
                / float(
                    m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"]
                )
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
                1
                / float(
                    m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"]
                )
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
                    b.commitmentPeriod[commitmentPeriod]
                    .genOff[thermalGen]
                    .indicator_var
                )
                if commitmentPeriod != 1
                else LogicalConstraint.Skip
            )

        @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
        def consistent_commitment_startup(b, commitmentPeriod, thermalGen):
            req_startup_periods = ceil(
                1
                / float(m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"])
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
                1
                / float(m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"])
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
                    int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"])
                    - 1,
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
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .genOn[thermalGen]
                    .indicator_var
                )
                .implies(
                    b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
                )
                if commitmentPeriod
                != 1  # int(m.md.data["elements"]["generator"][thermalGen]["min_up_time"])+1
                else LogicalConstraint.Skip
            )

        ##FIXME: Is this constraint necessary?
        # @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
        # def consistent_commitment_shutdown_after_uptime(
        #     b, commitmentPeriod, thermalGen
        # ):
        #     return (
        #         (
        #             atleast(
        #                 int(
        #                     m.md.data["elements"]["generator"][thermalGen][
        #                         "min_up_time"
        #                     ]
        #                 ),
        #                 [
        #                     b.commitmentPeriod[commitmentPeriod - j - 1]
        #                     .genOn[thermalGen]
        #                     .indicator_var
        #                     for j in range(
        #                         min(
        #                             int(
        #                                 m.md.data["elements"]["generator"][thermalGen][
        #                                     "min_up_time"
        #                                 ]
        #                             ),
        #                             commitmentPeriod - 1,
        #                         )
        #                     )
        #                 ],
        #             ).land(
        #                 b.commitmentPeriod[commitmentPeriod - 1]
        #                 .genOn[thermalGen]
        #                 .indicator_var
        #             )
        #         ).implies(
        #             b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
        #             | b.commitmentPeriod[commitmentPeriod]
        #             .genShutdown[thermalGen]
        #             .indicator_var
        #         )
        #         if commitmentPeriod != 1
        #         else LogicalConstraint.Skip
        #     )

        @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
        def consistent_commitment_downtime(b, commitmentPeriod, thermalGen):
            return (
                (
                    atmost(
                        int(
                            m.md.data["elements"]["generator"][thermalGen][
                                "min_down_time"
                            ]
                        )
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
                    b.commitmentPeriod[commitmentPeriod]
                    .genOff[thermalGen]
                    .indicator_var
                )
                if commitmentPeriod
                != 1  # >= int(m.md.data["elements"]["generator"][thermalGen]["min_down_time"])+1
                else LogicalConstraint.Skip
            )

    # TODO: The inter-commitment linking charge constraints are very ugly and 
    # repetitive. Can we make a variable for the particular periods we need
    # for cleaner code?
    """ Link battery charge level in consecutive commitment periods """
    @b.Constraint(b.commitmentPeriods, m.batteryStorageSystems)
    def consistent_battery_charge_level_commitment(b, commitmentPeriod, bat):
        if commitmentPeriod != 1:
            return (
                    # Charge Level in last dispatch period of previous commitment period
                    # less losses from inefficient retention
                    m.batteryRetentionRate[bat] * 
                    (
                        b.commitmentPeriod[commitmentPeriod - 1]
                        .dispatchPeriod
                        [
                            b.commitmentPeriod[commitmentPeriod - 1]
                            .dispatchPeriods.last()
                        ]
                        .batteryChargeLevel[bat] 
                    ) +
                    # Amount charged in first dispatch period of new commitment period
                    (
                        m.batteryChargingEfficiency[bat] * 
                        (
                            b.commitmentPeriod[commitmentPeriod]
                            .dispatchPeriod
                            [
                                b.commitmentPeriod[commitmentPeriod]
                                .dispatchPeriods.first()
                            ]
                            .batteryCharged[bat] 
                        )
                    ) - 
                    # Amount discharged in first dispatch period of new commitment period
                    (
                        b.commitmentPeriod[commitmentPeriod]
                        .dispatchPeriod
                        [
                            b.commitmentPeriod[commitmentPeriod]
                            .dispatchPeriods.first()
                        ]
                        .batteryDischarged[bat] 
                    )
                    == 
                    # Charge level in first dispatch period of new commitment period
                    b.commitmentPeriod[commitmentPeriod].dispatchPeriod
                    [
                        b.commitmentPeriod[commitmentPeriod]
                        .dispatchPeriods.first()
                        ]
                    .batteryChargeLevel[bat]
                )
        
        else: 
            # Initial value for each representative period.
            # Will constraints linking the representative period force 
            # a small amount of charging to offset the retention drop?
            return (
                # Initial charge level (data input)
                m.initBatteryChargeLevel[bat]
                +
                # Amount charged in first dispatch period of new commitment period
                (
                    m.batteryChargingEfficiency[bat] * 
                    (
                        b.commitmentPeriod[commitmentPeriod]
                        .dispatchPeriod
                        [
                            b.commitmentPeriod[commitmentPeriod]
                            .dispatchPeriods.first()
                        ]
                        .batteryCharged[bat] 
                    )
                ) - 
                # Amount discharged in first dispatch period of new commitment period
                (
                    b.commitmentPeriod[commitmentPeriod]
                    .dispatchPeriod
                    [
                        b.commitmentPeriod[commitmentPeriod]
                        .dispatchPeriods.first()
                    ]
                    .batteryDischarged[bat] 
                )
                
                
                 == 
                 (
                     b.commitmentPeriod[commitmentPeriod].dispatchPeriod
                     [
                         b.commitmentPeriod[commitmentPeriod]
                         .dispatchPeriods.first()
                     ]
                     .batteryChargeLevel[bat]
                 )
            )
        
        
    """ Link battery charge level in consecutive charging commitment periods """
    """ Ramp Up """
    @b.Constraint(b.commitmentPeriods, m.batteryStorageSystems)
    def consistent_battery_charge_charge_ramp_up_commitment(b, commitmentPeriod, bat):

        return (
            (
                b.commitmentPeriod[commitmentPeriod]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod]
                    .dispatchPeriods.first()
                ]
                .batteryCharged[bat] 
            )
            -
            (
                b.commitmentPeriod[commitmentPeriod - 1]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .dispatchPeriods.last()
                ]
                .batteryCharged[bat] 
            )
            <= 
            m.batteryChargingRampUpRates[bat]
            
            
        if commitmentPeriod != 1 and 
        b.commitmentPeriod[commitmentPeriod].batCharging[bat] and 
        b.commitmentPeriod[commitmentPeriod-1].batCharging[bat]
        
                
        else Constraint.Skip
        
              ) 
 
    
    """ Link battery charge level in consecutive charging commitment periods """
    """ Ramp Down """
    @b.Constraint(b.commitmentPeriods, m.batteryStorageSystems)
    def consistent_battery_charge_charge_ramp_down_commitment(b, commitmentPeriod, bat):

        return (
            (
                b.commitmentPeriod[commitmentPeriod - 1]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .dispatchPeriods.last()
                ]
                .batteryCharged[bat]
            )
            -
            (
                b.commitmentPeriod[commitmentPeriod]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod]
                    .dispatchPeriods.first()
                ]
                .batteryCharged[bat] 
            )
            <= 
            m.batteryChargingRampDownRates[bat]
            
            
        if commitmentPeriod != 1 and 
        b.commitmentPeriod[commitmentPeriod].batCharging[bat] and 
        b.commitmentPeriod[commitmentPeriod-1].batCharging[bat]
        
                
        else Constraint.Skip
        
              ) 

    
    """ Link battery discharge level in consecutive discharging commitment periods """
    """ Ramp Up """
    @b.Constraint(b.commitmentPeriods, m.batteryStorageSystems)
    def consistent_battery_discharge_discharge_ramp_up_commitment(b, commitmentPeriod, bat):

        return (
            (
                b.commitmentPeriod[commitmentPeriod]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod]
                    .dispatchPeriods.first()
                ]
                .batteryDischarged[bat] 
            )
            -
            (
                b.commitmentPeriod[commitmentPeriod - 1]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .dispatchPeriods.last()
                ]
                .batteryDischarged[bat] 
            )
            <= 
            m.batteryDischargingRampUpRates[bat]
            
            
        if commitmentPeriod != 1 and 
        b.commitmentPeriod[commitmentPeriod].batDischarging[bat] and 
        b.commitmentPeriod[commitmentPeriod-1].batDischarging[bat]
        
                
        else Constraint.Skip
        
             ) 
   
    
    """ Link battery discharge level in consecutive discharging commitment periods """
    """ Ramp Down """
    @b.Constraint(b.commitmentPeriods, m.batteryStorageSystems)
    def consistent_battery_discharge_discharge_ramp_down_commitment(b, commitmentPeriod, bat):

        return (
            (
                b.commitmentPeriod[commitmentPeriod - 1]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod - 1]
                    .dispatchPeriods.last()
                ]
                .batteryDischarged[bat] 
            )
            -
            (
                b.commitmentPeriod[commitmentPeriod]
                .dispatchPeriod
                [
                    b.commitmentPeriod[commitmentPeriod]
                    .dispatchPeriods.first()
                ]
                .batteryDischarged[bat] 
            )
            <= 
            m.batteryDischargingRampDownRates[bat]
            
            
        if commitmentPeriod != 1 and 
        b.commitmentPeriod[commitmentPeriod].batDischarging[bat] and 
        b.commitmentPeriod[commitmentPeriod-1].batDischarging[bat]
        
                
        else Constraint.Skip
        
             ) 
        
        
# def representative_period_rule(
#     b,
#     representative_period,
# ):
        ##FIXME: is this constraint necessary?
        # @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
        # def consistent_commitment_start_after_downtime(b, commitmentPeriod, thermalGen):
        #     return (
        #         (
        #             atleast(
        #                 int(
        #                     m.md.data["elements"]["generator"][thermalGen][
        #                         "min_down_time"
        #                     ]
        #                 ),
        #                 [
        #                     b.commitmentPeriod[commitmentPeriod - j - 1]
        #                     .genOff[thermalGen]
        #                     .indicator_var
        #                     for j in range(
        #                         min(
        #                             int(
        #                                 m.md.data["elements"]["generator"][thermalGen][
        #                                     "min_down_time"
        #                                 ]
        #                             ),
        #                             commitmentPeriod - 1,
        #                         )
        #                     )
        #                 ],
        #             ).land(
        #                 b.commitmentPeriod[commitmentPeriod - 1]
        #                 .genOff[thermalGen]
        #                 .indicator_var
        #             )
        #         ).implies(
        #             b.commitmentPeriod[commitmentPeriod]
        #             .genOff[thermalGen]
        #             .indicator_var
        #             | b.commitmentPeriod[commitmentPeriod]
        #             .genStartup[thermalGen]
        #             .indicator_var
        #         )
        #         if commitmentPeriod != 1
        #         else LogicalConstraint.Skip
        #     )


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
    b.load_scaling = i_s.load_scaling[
        (i_s.load_scaling["month"] == b.month) & (i_s.load_scaling["day"] == b.day)
    ]

    b.currentPeriod = representative_period
    if m.config["include_commitment"] or m.config["include_redispatch"]:
        b.commitmentPeriods = RangeSet(m.numCommitmentPeriods[representative_period])
        b.commitmentPeriod = Block(b.commitmentPeriods, rule=commitment_period_rule)

        add_representative_period_variables(b, representative_period)
        add_representative_period_constraints(b, representative_period)


def investment_stage_rule(b, investment_stage):
    """Creates investment stage block.

    :b: Investment block
    :investment_stage: ID for current investment stage
    """
    m = b.parent_block()

    b.year = m.years[investment_stage - 1]
    if m.config["scale_texas_loads"]:
        b.load_scaling = m.data.load_scaling[m.data.load_scaling["year"] == b.year]

        kw_to_mw_option = 1000
        other_option = 1
        ##TEXAS: lmao this is garbage; generalize this
        if investment_stage == 1:
            b.fixedCost = Param(m.generators, initialize=m.fixedCost1)
            b.varCost = Param(m.generators, initialize=m.varCost1)
            b.fuelCost = Param(m.generators, initialize=m.fuelCost1)
            thermalInvestmentCost = {
                gen: other_option
                * m.thermalCapacity[gen]
                * m.md.data["elements"]["generator"][gen]["capex1"]
                for gen in m.thermalGenerators
            }
            renewableInvestmentCost = {
                gen: other_option
                * m.renewableCapacity[gen]
                * m.md.data["elements"]["generator"][gen]["capex1"]
                for gen in m.renewableGenerators
            }
            m.generatorInvestmentCost = thermalInvestmentCost | renewableInvestmentCost
            print("gen investment cost")
            print(sum(m.generatorInvestmentCost.values()))
        elif investment_stage == 2:
            b.fixedCost = Param(m.generators, initialize=m.fixedCost2)
            b.varCost = Param(m.generators, initialize=m.varCost2)
            b.fuelCost = Param(m.generators, initialize=m.fuelCost2)
            thermalInvestmentCost = {
                gen: other_option
                * m.thermalCapacity[gen]
                * m.md.data["elements"]["generator"][gen]["capex2"]
                for gen in m.thermalGenerators
            }
            renewableInvestmentCost = {
                gen: other_option
                * m.renewableCapacity[gen]
                * m.md.data["elements"]["generator"][gen]["capex2"]
                for gen in m.renewableGenerators
            }
            m.generatorInvestmentCost = thermalInvestmentCost | renewableInvestmentCost
        else:
            b.fixedCost = Param(m.generators, initialize=m.fixedCost3)
            b.varCost = Param(m.generators, initialize=m.varCost3)
            b.fuelCost = Param(m.generators, initialize=m.fuelCost3)
            thermalInvestmentCost = {
                gen: other_option
                * m.thermalCapacity[gen]
                * m.md.data["elements"]["generator"][gen]["capex3"]
                for gen in m.thermalGenerators
            }
            renewableInvestmentCost = {
                gen: other_option
                * m.renewableCapacity[gen]
                * m.md.data["elements"]["generator"][gen]["capex3"]
                for gen in m.renewableGenerators
            }
            m.generatorInvestmentCost = thermalInvestmentCost | renewableInvestmentCost

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
    """ Added Battery Storage Cost to Objective Function """
    
    if len(m.stages) > 1:
        m.operatingCost = sum(
            m.investmentStage[stage].operatingCostInvestment for stage in m.stages
        )
        m.storageCost = sum(
            m.investmentStage[stage].storageCostInvestment for stage in m.stages
        ) 
        m.expansionCost = sum(
            m.investmentStage[stage].investment_cost for stage in m.stages
        )
        m.penaltyCost = sum(
            m.deficitPenalty[stage]
            * m.investmentFactor[stage]
            * m.investmentStage[stage].quotaDeficit
            + m.investmentStage[stage].renewableCurtailmentInvestment
            for stage in m.stages
        )

    @m.Objective()
    def total_cost_objective_rule(m):
        if len(m.stages) > 1:
            return m.operatingCost + m.storageCost + m.expansionCost + m.penaltyCost
        else:
            return (
                m.investmentStage[1].operatingCostInvestment
                + m.investmentStage[1].storageCostInvestment # JSC Addn
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
        initialize=m.md.data["elements"]["generator"].keys(), doc="All generators"
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

    """ Hard-coded a test battery. Data inputs from csv file had issues. """
    ## NOTE: will want to cover baseline generator types in IDAES
    # This should be updated for battery. @JKS is this using the 
    # built-in structure from EGRET or just a placeholder?
    if m.md.data["elements"].get("storage"):
        m.storage = Set(
            initialize=(ess for ess in m.md.data["elements"]["storage"]),
            doc="Potential storage units",
        )
        
    else:
        # TODO: assign and modify data below with better parameters. 
        # Currently, not all are used.
        m.md.data["elements"]["storage"] = {
            "test_battery": {
                "name": "ideas_spelled_wrong",
                "bus": 3,
                "generator": None,
                "storage_type": "battery",
                "energy_capacity": 100,
                "initial_state_of_charge": 5, 
                "end_state_of_charge": 5,
                "minimum_state_of_charge": 5,
                "charge_efficiency": 1, 
                "discharge_efficiency": 1, 
                "max_discharge_rate": 20, 
                "min_discharge_rate": 2, 
                "max_charge_rate": 100, 
                "min_charge_rate": 1, 
                "initial_charge_rate": 0,
                "initial_discharge_rate": 0,
                "charge_cost": 0, 
                "discharge_cost": 0, 
                "retention_rate_60min": 1, # This has been verified to work at levels below 1; currently set to 1 for testing other storage components
                "ramp_up_input_60min": 1, 
                "ramp_down_input_60min": 1, 
                "ramp_up_output_60min": 2,
                "ramp_down_output_60min": 2, 
                "in_service": True,
                "capital_multiplier": 1, 
                "extension_multiplier": 1}} # Thermal generator fuel costs are on [0.5,1.5]; renewables have no fuel cost. What should go here?
        
        m.storage = Set(
            initialize=(ess for ess in m.md.data["elements"]["storage"]),
            doc="Potential storage units",
        )
        
    m.batteryStorageSystems = Set(
        within=m.storage,
        initialize=(
            batt for batt in m.storage
            if m.md.data["elements"]["storage"][batt]["storage_type"] == "battery"),
        doc="Batteries; subset of all energy storage systems")

    ## TODO: make sure time units are both definable and consistent without being forced

    m.stages = RangeSet(stages, doc="Set of planning periods")

    m.representativePeriods = Set(
        initialize=rep_per, doc="Set of representative periods for each planning period"
    )


def model_data_references(m):
    """Creates and labels data for GTEP model; ties input data
    to model directly.
    :param m: Pyomo model object
    """

    # Maximum output of each thermal generator
    m.thermalCapacity = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_max"]
        for thermalGen in m.thermalGenerators
    }

    print(sum(m.thermalCapacity.values()))

    # Lifetime of each generator; needs units
    m.lifetimes = {
        gen: m.md.data["elements"]["generator"][gen]["lifetime"] for gen in m.generators
    }

    # Minimum output of each thermal generator
    m.thermalMin = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_min"]
        for thermalGen in m.thermalGenerators
    }
    
    """ Battery Storage properties read-in from data """
    m.batteryCapacity = {
        bat: m.md.data["elements"]["storage"][bat]["energy_capacity"] 
        for bat in m.batteryStorageSystems
        } # maximum storage capacity
    
    m.initBatteryChargeLevel = {
        bat: m.md.data["elements"]["storage"][bat]["initial_state_of_charge"] 
        for bat in m.batteryStorageSystems
        } # initial storage capacity
    
    m.minBatteryChargeLevel = {
        bat: m.md.data["elements"]["storage"][bat]["minimum_state_of_charge"] 
        for bat in m.batteryStorageSystems
        } # minimum storage capacity
    
    m.chargingCost = {
        bat: m.md.data["elements"]["storage"][bat]["charge_cost"] 
        for bat in m.batteryStorageSystems
        } # cost to charge per unit electricity
    
    m.dischargingCost = {
        bat: m.md.data["elements"]["storage"][bat]["discharge_cost"] 
        for bat in m.batteryStorageSystems
        } # cost to discharge per unit electricity
    
    m.dischargeMin = {
        bat: m.md.data["elements"]["storage"][bat]["min_discharge_rate"] 
        for bat in m.batteryStorageSystems
        } # minimum amount to discharge per dispatch period when discharging
    
    m.dischargeMax = {
        bat: m.md.data["elements"]["storage"][bat]["max_discharge_rate"] 
        for bat in m.batteryStorageSystems
        } # maximum amount to discharge per dispatch period when discharging
    
    m.chargeMin = {
        bat: m.md.data["elements"]["storage"][bat]["min_charge_rate"] 
        for bat in m.batteryStorageSystems
        } # minimum amount to charge per dispatch period when charging
    
    m.chargeMax = {
        bat: m.md.data["elements"]["storage"][bat]["max_charge_rate"] 
        for bat in m.batteryStorageSystems
        } # maximum amount to charge per dispatch period when charging
    
    m.batteryDischargingRampUpRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_up_output_60min"] 
        for bat in m.batteryStorageSystems
        } # maximum amount of ramp up between dispatch periods when discharging. 
    # Notice that default EGRET naming convention assumes dispatch periods are 60 minutes
    
    m.batteryDischargingRampDownRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_down_output_60min"] 
        for bat in m.batteryStorageSystems
        } # maximum amount of ramp down between dispatch periods when discharging.
    
    m.batteryChargingRampUpRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_up_input_60min"] 
        for bat in m.batteryStorageSystems
        } # maximum amount of ramp up between dispatch periods when charging.
    
    m.batteryChargingRampDownRates = {
        bat: m.md.data["elements"]["storage"][bat]["ramp_down_input_60min"] 
        for bat in m.batteryStorageSystems
        } # maximum amount of ramp down between dispatch periods when charging.
    
    m.batteryDischargingEfficiency  = {
        bat: m.md.data["elements"]["storage"][bat]["discharge_efficiency"] 
        for bat in m.batteryStorageSystems
        } # proportion of energy discharged that is not lost to technological 
    # inefficiencies with in dispatch periods and which is usable in the flow balance
    
    m.batteryChargingEfficiency  = {
        bat: m.md.data["elements"]["storage"][bat]["charge_efficiency"] 
        for bat in m.batteryStorageSystems
        } # proportion of energy charged that is not lost to technological 
    # inefficiencies within dispatch periods and which is usable in the flow balance
    
    m.batteryRetentionRate  = {
        bat: m.md.data["elements"]["storage"][bat]["retention_rate_60min"] 
        for bat in m.batteryStorageSystems
        } # proportion of energy discharged that is not lost to technological 
    # inefficiencies between dispatch periods and which is usable in the flow balance
    
    # (Arbitrary) multiplier for new battery investments corresponds to depreciation schedules
    # for individual technologies; higher values are indicative of slow depreciation
    m.batteryCapitalMultiplier = {
        bat: m.md.data["elements"]["storage"][bat]["capital_multiplier"]
        for bat in m.batteryStorageSystems
    }

    # Cost of life extension for each battery, expressed as a fraction of initial investment cost
    m.batteryExtensionMultiplier = {
        bat: m.md.data["elements"]["storage"][bat]["extension_multiplier"]
        for bat in m.batteryStorageSystems
    }
    
    m.batteryInvestmentCost = {
        bat: 0
        for bat in m.batteryStorageSystems
        } # Future not real cost: idealized DoE 10-yr targets or something
    

    # Maximum output of each renewable generator
    m.renewableCapacity = {
        renewableGen: (
            0
            if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float
            else max(
                m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"]
            )
        )
        for renewableGen in m.renewableGenerators
    }

    #print(sum(m.renewableCapacity.values()))

    # A fraction of renewableCapacity representing fraction of capacity
    # that can be reliably counted toward planning reserve requirement
    # TODO: WHAT HAVE I DONE HERE I HATE IT 
    m.renewableCapacityValue = {
        renewableGen: (
            0
            if type(m.md.data["elements"]["generator"][renewableGen]["p_max"]) == float
            else min(
                m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"]
            )
            / max(1, m.renewableCapacity[renewableGen])
        )
        for renewableGen in m.renewableGenerators
    }

    # Long term thermal rating of each transmission line
    m.transmissionCapacity = {
        transmissionLine: m.md.data["elements"]["branch"][transmissionLine][
            "rating_long_term"
        ]
        for transmissionLine in m.transmission.keys()
    }

    # Maximum fraction of a thermal generator's maximum output that can be
    # supplied as spinning reserve
    m.spinningReserveFraction = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen][
            "spinning_reserve_frac"
        ]
        for thermalGen in m.thermalGenerators
    }

    # Maximum fraction of a thermal generator's maximum output that can be
    # supplied as quickstart reserve
    m.quickstartReserveFraction = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen][
            "quickstart_reserve_frac"
        ]
        for thermalGen in m.thermalGenerators
    }

    # Demand at each bus
    m.loads = {
        m.md.data["elements"]["load"][load_n]["bus"]: m.md.data["elements"]["load"][
            load_n
        ]["p_load"]
        for load_n in m.md.data["elements"]["load"]
    }
    # for key, val in m.loads.items():
    #     for i, v in enumerate(val['values']):
    #         val['values'][i] *= 1/10

    ## NOTE: lazy fixing for dc_branch and branch... but should be an ok lazy fix
    # Per-distance-unit multiplicative loss rate for each transmission line
    m.lossRate = {
        branch: (m.md.data["elements"]["branch"][branch].get("loss_rate") or 0)
        for branch in m.transmission
    }

    ## NOTE: lazy fixing for dc_branch and branch... but should be an ok lazy fix
    # Distance between terminal buses for each transmission line
    m.distance = {
        branch: (m.md.data["elements"]["branch"][branch].get("distance") or 0)
        for branch in m.transmission
    }

    # TODO: Add cost of investment in each new branch to input data. Currently
    # selected 0 to ensure investments will be selected if needed
    m.branchInvestmentCost = {
        branch: (m.md.data["elements"]["branch"][branch].get("capital_cost") or 0)
        for branch in m.transmission
    }

    # TODO: Add branch capital multiplier to input data.
    m.branchCapitalMultiplier = {
        branch: (m.md.data["elements"]["branch"][branch].get("capital_multiplier") or 1)
        for branch in m.transmission
    }

    # Cost of life extension for each generator, expressed as a fraction of initial investment cost
    m.branchExtensionMultiplier = {
        branch: (
            m.md.data["elements"]["branch"][branch].get("extension_multiplier") or 1
        )
        for branch in m.transmission
    }

    ## TODO: These should go into each stage -- check where these values should come from
    m.peakLoad = Param(m.stages, default=0, units=u.MW)
    m.reserveMargin = Param(m.stages, default=0, units=u.MW)
    m.renewableQuota = Param(m.stages, default=0, units=u.MW)
    m.weights = Param(m.representativePeriods, initialize=m.data.representative_weights, default=5 * 365 / 4)
    m.investmentFactor = Param(m.stages, default=1, mutable=True)
    ## NOTE: Lazy approx for NPV
    ## TODO: don't lazily approx NPV, add it into unit handling and calculate from actual time frames
    # for stage in m.stages:
    #     m.investmentFactor[stage] *= 1 / ((1.04) ** (5 * stage))
    m.fixedOperatingCost = Param(m.generators, default=1, units=u.USD / u.hr)
    m.deficitPenalty = Param(m.stages, default=1, units=u.USD / (u.MW * u.hr))

    # Amount of fuel required to be consumed for startup process for each generator
    m.startFuel = {
        gen: m.md.data["elements"]["generator"][gen]["start_fuel"]
        for gen in m.generators
    }

    # TEXAS: make this a list per investment stage or whatever
    fuelCost1 = {}
    fuelCost2 = {}
    fuelCost3 = {}
    # Cost per unit of fuel at each generator
    if m.config["scale_texas_loads"]:
        for gen in m.thermalGenerators:
            fuelCost1[gen] = m.md.data["elements"]["generator"][gen]["fuel_cost1"]
            fuelCost2[gen] = m.md.data["elements"]["generator"][gen]["fuel_cost2"]
            fuelCost3[gen] = m.md.data["elements"]["generator"][gen]["fuel_cost3"]
    elif "RTS-GMLC" in m.md.data["system"]["name"]:
        for gen in m.thermalGenerators:
            fuelCost[gen] = m.md.data["elements"]["generator"][gen]["fuel_cost"]
    else:
        for gen in m.thermalGenerators:
            fuelCost[gen] = m.md.data["elements"]["generator"][gen]["p_cost"]["values"][
                1
            ]

    m.fuelCost1 = Param(
        m.thermalGenerators, initialize=fuelCost1, units=u.USD / (u.MW * u.hr)
    )
    m.fuelCost2 = Param(
        m.thermalGenerators, initialize=fuelCost2, units=u.USD / (u.MW * u.hr)
    )
    m.fuelCost3 = Param(
        m.thermalGenerators, initialize=fuelCost3, units=u.USD / (u.MW * u.hr)
    )

    fixedCost1 = {}
    fixedCost2 = {}
    fixedCost3 = {}
    varCost1 = {}
    varCost2 = {}
    varCost3 = {}
    if m.config["scale_texas_loads"]:
        for gen in m.generators:
            fixedCost1[gen] = m.md.data["elements"]["generator"][gen]["fixed_ops1"]
            fixedCost2[gen] = m.md.data["elements"]["generator"][gen]["fixed_ops2"]
            fixedCost3[gen] = m.md.data["elements"]["generator"][gen]["fixed_ops3"]
            varCost1[gen] = m.md.data["elements"]["generator"][gen]["var_ops1"]
            varCost2[gen] = m.md.data["elements"]["generator"][gen]["var_ops2"]
            varCost3[gen] = m.md.data["elements"]["generator"][gen]["var_ops3"]

    m.fixedCost1 = Param(
        m.generators, initialize=fixedCost1, units=u.USD / (u.MW * u.hr)
    )
    m.fixedCost2 = Param(
        m.generators, initialize=fixedCost2, units=u.USD / (u.MW * u.hr)
    )
    m.fixedCost3 = Param(
        m.generators, initialize=fixedCost3, units=u.USD / (u.MW * u.hr)
    )
    m.varCost1 = Param(m.generators, initialize=varCost1, units=u.USD / (u.MW * u.hr))
    m.varCost2 = Param(m.generators, initialize=varCost2, units=u.USD / (u.MW * u.hr))
    m.varCost3 = Param(m.generators, initialize=varCost3, units=u.USD / (u.MW * u.hr))

    # Cost per MW of curtailed renewable energy
    # NOTE: what should this be valued at?  This being both curtailment and load shed.
    # TODO: update valuations
    m.curtailmentCost = Param(
        initialize=100,
        units=u.USD / (u.MW * u.hr),
    )
    m.loadShedCost = Param(
        initialize=5000,
        units=u.USD / (u.MW * u.hr),
    )

    # Full lifecycle CO_2 emission factor for each generator
    m.emissionsFactor = {
        gen: m.md.data["elements"]["generator"][gen]["emissions_factor"]
        for gen in m.generators
    }

    # Flat startup cost for each generator
    if "RTS-GMLC" in m.md.data["system"]["name"]:
        startupCost = {
            gen: m.md.data["elements"]["generator"][gen]["non_fuel_startup_cost"]
            for gen in m.thermalGenerators
        }
        m.startupCost = Param(m.thermalGenerators, initialize=startupCost, units=u.USD)
    else:
        startupCost = {
            gen: m.md.data["elements"]["generator"][gen]["startup_cost"]
            for gen in m.generators
        }
        m.startupCost = Param(m.generators, initialize=startupCost, units=u.USD)

    # (Arbitrary) multiplier for new generator investments corresponds to depreciation schedules
    # for individual technologies; higher values are indicative of slow depreciation
    m.capitalMultiplier = {
        gen: m.md.data["elements"]["generator"][gen]["capital_multiplier"]
        for gen in m.generators
    }

    # Cost of life extension for each generator, expressed as a fraction of initial investment cost
    m.extensionMultiplier = {
        gen: m.md.data["elements"]["generator"][gen]["extension_multiplier"]
        for gen in m.generators
    }
    if m.config["scale_texas_loads"]:
        m.extensionMultiplier = {gen: 0.06 for gen in m.generators}
    thermal_retirement = {gen: 0.1 for gen in m.thermalGenerators}
    renewable_retirement = {gen: 1 for gen in m.renewableGenerators}
    m.retirementMultiplier = thermal_retirement | renewable_retirement

    # Cost of investment in each new generator
    m.generatorInvestmentCost = {
        # gen: m.md.data["elements"]["generator"][gen]["investment_cost"]
        # for gen in m.generators]
        gen: 0
        for gen in m.generators
    }

    # Minimum operating reserve, expressed as a fraction of load within a region
    m.minOperatingReserve = {
        region: m.md.data["system"]["min_operating_reserve"] for region in m.regions
    }

    # Minimum spinning reserve, expressed as a fraction of load within a region
    m.minSpinningReserve = {
        region: m.md.data["system"]["min_spinning_reserve"] for region in m.regions
    }

    # Maximum spinning reserve available for each generator; expressed as a fraction
    # maximum generator output
    m.maxSpinningReserve = {
        gen: m.md.data["elements"]["generator"][gen]["max_spinning_reserve"]
        for gen in m.thermalGenerators
    }

    # Maximum quickstart reserve available for each generator; expressed as a fraction
    # maximum generator output
    m.maxQuickstartReserve = {
        gen: m.md.data["elements"]["generator"][gen]["max_quickstart_reserve"]
        for gen in m.thermalGenerators
    }

    # Ramp up rates for each generator; expressed as a fraction of maximum generator output
    m.rampUpRates = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"]
        for thermalGen in m.thermalGenerators
    }

    # Ramp down rates for each generator; expressed as a fraction of maximum generator output
    m.rampDownRates = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"]
        for thermalGen in m.thermalGenerators
    }

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


def model_create_investment_stages(m, stages):
    """Creates investment blocks and linking constraints for GTEP model.
    Largely manages retirements and links operational units in a given investment stage
    to operational + installed - retired in the previous investment stage.

    :m: Pyomo model object
    :stages: Number of investment stages in planning horizon
    """

    ## NOTE: temporary years handling for texas case study
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
    if m.config["include_investment"]:

        # # Linking generator investment status constraints
        # @m.Constraint(m.stages, m.thermalGenerators)
        # def gen_stats_link(m, stage, gen):
        #     return (
        #         m.investmentStage[stage]
        #         .genOperational[gen]
        #         .indicator_var.get_associated_binary()
        #         == m.investmentStage[stage - 1]
        #         .genOperational[gen]
        #         .indicator_var.get_associated_binary()
        #         + m.investmentStage[stage - 1]
        #         .genInstalled[gen]
        #         .indicator_var.get_associated_binary()
        #         - m.investmentStage[stage - 1]
        #         .genRetired[gen]
        #         .indicator_var.get_associated_binary()
        #         if stage != 1
        #         else Constraint.Skip
        #     )

    # TODO: Do we need these for branches and storage? Would guess yes, but branches seemed to work without it?
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

    """ Battery investment stage state change logic """
    @m.Constraint(m.stages, m.batteryStorageSystems)
    def bat_stats_link(m, stage, bat):
        return (
            m.investmentStage[stage]
            .batOperational[bat]
            .indicator_var.get_associated_binary()
            == m.investmentStage[stage - 1]
            .batOperational[bat]
            .indicator_var.get_associated_binary()
            + m.investmentStage[stage - 1]
            .batInstalled[bat]
            .indicator_var.get_associated_binary()
            - m.investmentStage[stage - 1]
            .batRetired[bat]
            .indicator_var.get_associated_binary()
            if stage != 1
            else Constraint.Skip
        )
    
    @m.Constraint(m.stages, m.transmission)
    def branch_stats_link(m, stage, branch):
        return (
            m.investmentStage[stage]
            .branchOperational[branch]
            .indicator_var.get_associated_binary()
            == m.investmentStage[stage - 1]
            .branchOperational[branch]
            .indicator_var.get_associated_binary()
            + m.investmentStage[stage - 1]
            .branchInstalled[branch]
            .indicator_var.get_associated_binary()
            - m.investmentStage[stage - 1]
            .branchRetired[branch]
            .indicator_var.get_associated_binary()
            if stage != 1
            else Constraint.Skip
        )

    # Renewable generation (in MW) retirement relationships
    if len(m.stages) > 1:
        if len(m.stages) > 1:
            ##FIXME Rewrite as logic
            @m.Constraint(m.stages, m.thermalGenerators)
            def gen_retirement(m, stage, gen):
                return sum(
                    m.investmentStage[t_2]
                    .genOperational[gen]
                    .indicator_var.get_associated_binary()
                    + m.investmentStage[t_2]
                    .genInstalled[gen]
                    .indicator_var.get_associated_binary()
                    for t_2 in m.stages
                    if t_2 <= stage - m.lifetimes[gen]
                ) <= sum(
                    m.investmentStage[t_1]
                    .genRetired[gen]
                    .indicator_var.get_associated_binary()
                    + m.investmentStage[t_1]
                    .genExtended[gen]
                    .indicator_var.get_associated_binary()
                    for t_1 in m.stages
                    if t_1 <= stage
                )

        # Renewable generation (in MW) retirement relationships
        # if len(m.stages) > 1:

        #     @m.Constraint(m.stages, m.renewableGenerators)
        #     def renewable_retirement(m, stage, gen):
        #         return sum(
        #             m.investmentStage[t_2].renewableInstalled[gen]
        #             for t_2 in m.stages
        #             if t_2 <= stage - m.lifetimes[gen]
        #         ) <= sum(
        #             m.investmentStage[t_1].renewableRetired[gen]
        #             + m.investmentStage[t_1].renewableExtended[gen]
        #             for t_1 in m.stages
        #             if t_1 <= stage
        #         )

        # Total renewable generation (in MW) operational at a given stage
        # is equal to what was operational and/or installed in the previous stage
        # less what was retired in the previous stage
        @m.Constraint(m.stages, m.renewableGenerators)
        def renewable_stats_link(m, stage, gen):
            return (
                m.investmentStage[stage].renewableOperational[gen]
                == m.investmentStage[stage - 1].renewableOperational[gen]
                + m.investmentStage[stage - 1].renewableInstalled[gen]
                + m.investmentStage[stage - 1].renewableExtended[gen]
                - m.investmentStage[stage - 1].renewableRetired[gen]
                if stage != 1
                else Constraint.Skip
            )

        # @m.Constraint(m.stages, m.renewableGenerators)
        # def renewable_more_stats_link(m, stage, gen):
        #     return (
        #         m.investmentStage[stage].renewableDisabled[gen]
        #         == m.investmentStage[stage - 1].renewableDisabled[gen]
        #         + m.investmentStage[stage - 1].renewableRetired[gen]
        #         - m.investmentStage[stage - 1].renewableInstalled[gen]
        #         if stage != 1
        #         else Constraint.Skip
        #     )

        # @m.Constraint(m.stages, m.renewableGenerators)
        # def renewable_capacity_enforcement(m, stage, gen):
        #     return m.investmentStage[stage].renewableOperational[gen] + m.investmentStage[stage].renewableInstalled[gen] <= m.renewableCapacity[gen]

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
        ##FIXME Disabling is permanent.  Re investment is a "new" unit.  Remove the "or"
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
    
    """ Battery investment stage logic Pt 2"""
    # If a bat is online at time t, it must have been online or installed at time t-1
    @m.LogicalConstraint(m.stages, m.batteryStorageSystems)
    def consistent_battery_operation(m, stage, bat):
        return (
            m.investmentStage[stage]
            .batOperational[bat]
            .indicator_var.implies(
                m.investmentStage[stage - 1].batOperational[bat].indicator_var
                | m.investmentStage[stage - 1].batInstalled[bat].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a bat is online at time t, it must be online, extended, or retired at time t+1
    @m.LogicalConstraint(m.stages, m.batteryStorageSystems)
    def consistent_operation_battery_future(m, stage, bat):
        return (
            m.investmentStage[stage - 1]
            .batOperational[bat]
            .indicator_var.implies(
                m.investmentStage[stage].batOperational[bat].indicator_var
                | m.investmentStage[stage].batExtended[bat].indicator_var
                | m.investmentStage[stage].batRetired[bat].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # Retirement in period t-1 implies disabled in period t
    @m.LogicalConstraint(m.stages, m.batteryStorageSystems)
    def full_battery_retirement(m, stage, bat):
        return (
            m.investmentStage[stage - 1]
            .batRetired[bat]
            .indicator_var.implies(
                m.investmentStage[stage].batDisabled[bat].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a bat is disabled at time t-1, it must stay disabled  at time t
    @m.LogicalConstraint(m.stages, m.batteryStorageSystems)
    def consistent_battery_disabled(m, stage, bat):
        return (
            m.investmentStage[stage - 1]
            .batDisabled[bat]
            .indicator_var.implies(
                m.investmentStage[stage].batDisabled[bat].indicator_var
                | m.investmentStage[stage].batInstalled[bat].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # If a bat is extended at time t-1, it must stay extended or be retired at time t
    @m.LogicalConstraint(m.stages, m.batteryStorageSystems)
    def consistent_battery_extended(m, stage, bat):
        return (
            m.investmentStage[stage - 1]
            .batExtended[bat]
            .indicator_var.implies(
                m.investmentStage[stage].batExtended[bat].indicator_var
                | m.investmentStage[stage].batRetired[bat].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

    # Installation in period t-1 implies operational in period t
    @m.LogicalConstraint(m.stages, m.batteryStorageSystems)
    def full_battery_investment(m, stage, bat):
        return (
            m.investmentStage[stage - 1]
            .batInstalled[bat]
            .indicator_var.implies(
                m.investmentStage[stage].batOperational[bat].indicator_var
            )
            if stage != 1
            else LogicalConstraint.Skip
        )

        if m.config["transmission"]:
            # If a branch is online at time t, it must have been online or installed at time t-1
            @m.LogicalConstraint(m.stages, m.transmission)
            def consistent_branch_operation(m, stage, branch):
                return (
                    m.investmentStage[stage]
                    .branchOperational[branch]
                    .indicator_var.implies(
                        m.investmentStage[stage - 1]
                        .branchOperational[branch]
                        .indicator_var
                        | m.investmentStage[stage - 1]
                        .branchInstalled[branch]
                        .indicator_var
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
