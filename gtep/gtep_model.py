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


# Define what a USD is for pyomo units purposes
# This will be set to a base year and we will do NPV calculations
# based on automatic pyomo unit transformations
u.load_definitions_from_strings(["USD = [currency]"])

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
    """A generalized generation and transmission expansion planning model.
    """
    def __init__(
        self,
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
        self.timer = TicTocTimer()

    def create_model(self):
        """Create concrete Pyomo model object associated with the ExpansionPlanningModel
        """
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
        m.dispatchPeriodLength = Param(within=PositiveReals, default=15, units=u.min)

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

    def report_large_coefficients(self, outfile, magnitude_cutoff):
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
    """Add continuous variables to investment stage block.
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

    # Renewable generator MW values (operational, installed, retired, extended)
    b.renewableOperational = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0
    )
    b.renewableInstalled = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0
    )
    b.renewableRetired = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0
    )
    b.renewableExtended = Var(
        m.renewableGenerators, within=NonNegativeReals, initialize=0
    )

    # Track and accumulate costs and penalties
    b.quotaDeficit = Var(within=NonNegativeReals, initialize=0, units=u.MW)
    b.operatingCostInvestment = Var(within=Reals, initialize=0, units=u.USD)
    b.expansionCost = Var(within=Reals, initialize=0, units=u.USD)
    b.renewableCurtailmentInvestment = Var(
        within=NonNegativeReals, initialize=0, units=u.USD
    )


def add_investment_constraints(
    b,
    investment_stage,
):
    """Add standard inequalities (i.e., those not involving disjunctions) to investment stage block.
    """

    m = b.model()

    ## TODO: Fix var value rather than add constraint
    @b.LogicalConstraint(m.thermalGenerators)
    def thermal_uninvested(b, gen):
        if m.md.data["elements"]["generator"][gen]["in_service"] == False:
            return exactly(0, b.genOperational[gen].indicator_var)
        else:
            return LogicalConstraint.Skip

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

    ## NOTE: Constraint (13) in the reference paper
    # Minimum per-stage renewable generation requirement
    @b.Constraint()
    def renewable_generation_requirement(b):
        renewableSurplusRepresentative = 0
        ## TODO: preprocess loads for the appropriate sum here
        ed = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                renewableSurplusRepresentative += (
                    m.weights[rep_per]
                    * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .renewableSurplusCommitment
                )
        return (
            renewableSurplusRepresentative + b.quotaDeficit
            >= m.renewableQuota[investment_stage] * ed
        )

    # Operating costs for investment period
    @b.Constraint()
    def operating_cost_investment(b):
        operatingCostRepresentative = 0
        for rep_per in b.representativePeriods:
            for com_per in b.representativePeriod[rep_per].commitmentPeriods:
                operatingCostRepresentative += (
                    m.weights[rep_per]
                    * m.commitmentPeriodLength
                    * b.representativePeriod[rep_per]
                    .commitmentPeriod[com_per]
                    .operating_cost_commitment
                )
        return (
            b.operatingCostInvestment
            == m.investmentFactor[investment_stage] * operatingCostRepresentative
        )

    # Investment costs for investment period
    ## NOTE: investment cost definition needs to be revisited AND possibly depends on
    ## data format.  It is _rare_ for these values to be defined at all, let alone consistently.
    @b.Constraint()
    def investment_cost(b):
        return b.expansionCost == m.investmentFactor[investment_stage] * (
            sum(
                m.generatorInvestmentCost[gen] * m.capitalMultiplier[gen]
                # * m.thermalCapacity[gen]
                * b.genInstalled[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.capitalMultiplier[gen]
                * m.renewableCapacity[gen]
                * b.renewableInstalled[gen]
                for gen in m.renewableGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen] * m.extensionMultiplier[gen]
                # * m.thermalCapacity[gen]
                * b.genExtended[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
            + sum(
                m.generatorInvestmentCost[gen]
                * m.extensionMultiplier[gen]
                * m.renewableCapacity[gen]
                * b.renewableExtended[gen]
                for gen in m.renewableGenerators
            )
        )

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


def add_dispatch_variables(
    b,
    dispatch_period,
):
    """Add dispatch-associated variables to representative period block.
    """

    m = b.model()
    c_p = b.parent_block()

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

    # Define bounds on transmission line capacity
    def power_flow_limits(b, transmissionLine):
        return (
            -m.transmissionCapacity[transmissionLine],
            m.transmissionCapacity[transmissionLine],
        )

    # NOTE: this is an abuse of units and needs to be fixed for variable temporal resolution
    b.powerFlow = Var(
        m.transmission,
        domain=Reals,
        bounds=power_flow_limits,
        initialize=0,
        units=u.MW,
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
        units=u.MW,
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
        units=u.MW,
    )

    # Track total renewable surplus/deficit for other expressions
    b.renewableSurplusDispatch = Var(domain=Reals, initialize=0, units=u.MW)
    b.operatingCostDispatch = Var(domain=Reals, initialize=0, units=u.USD)
    b.renewableCurtailmentDispatch = Var(
        domain=NonNegativeReals, initialize=0, units=u.MW
    )

    # Load shed per bus
    b.loadShed = Var(m.buses, domain=NonNegativeReals, initialize=0)

    # Voltage angle
    def bus_angle_bounds(b, bus):
        return (-30, 30)

    b.busAngle = Var(m.buses, domain=Reals, initialize=0, bounds=bus_angle_bounds)

    # Voltage angle difference
    def delta_bus_angle_bounds(b, line):
        return (-60, 60)

    b.deltaBusAngle = Var(
        m.transmission, domain=Reals, initialize=0, bounds=delta_bus_angle_bounds
    )


def add_dispatch_constraints(
    b,
    disp_per,
):
    """Add dispatch-associated inequalities to representative period block.
    """
    m = b.model()
    c_p = b.parent_block()
    r_p = c_p.parent_block()
    i_p = r_p.parent_block()

    ## TODO: what do we do when reactance isn't supplied in the dataset?
    @b.Constraint(m.transmission)
    def dc_power_flow(b, line):
        fb = m.transmission[line]["from_bus"]
        tb = m.transmission[line]["to_bus"]
        reactance = m.md.data["elements"]["branch"][line]["reactance"]
        if m.md.data["elements"]["branch"][line]["branch_type"] == "transformer":
            reactance *= m.md.data["elements"]["branch"][line]["transformer_tap_ratio"]
            shift = m.md.data["elements"]["branch"][line]["transformer_phase_shift"]
        else:
            shift = 0
        return b.powerFlow[line] == (-1 / reactance) * (
            b.busAngle[tb] - b.busAngle[fb] + shift
        )

    # Energy balance constraint
    ## FIXME: curtailment handling is problematic; see below
    @b.Constraint(m.buses)
    def flow_balance(b, bus):
        balance = 0
        load = m.loads.get(bus) or 0
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
        # balance -= sum(
        #    b.renewableCurtailment[g] for g in gens if g in m.renewableGenerators
        # )
        balance -= load
        balance += b.loadShed[bus]
        return balance == 0

    # Capacity factor constraint
    # NOTE: In comparison to reference work, this is *per renewable generator*
    @b.Constraint(m.renewableGenerators)
    def capacity_factor(b, renewableGen):
        return (
            b.renewableGeneration[renewableGen] + b.renewableCurtailment[renewableGen]
            == m.renewableCapacity[renewableGen]
        )

    ## FIXME
    # NOTE: I believe this makes curtailment crazy expensive
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

    # Define total renewable surplus/deficit for dispatch block
    @b.Constraint()
    def renewable_surplus_dispatch(b):
        return b.renewableSurplusDispatch == sum(
            b.renewableGeneration[gen] - b.renewableCurtailment[gen]
            for gen in m.renewableGenerators
        )

    # Define total renewable curtailment for dispatch block
    @b.Constraint()
    def renewable_curtailment_dispatch(b):
        return b.renewableCurtailmentDispatch == sum(
            b.renewableCurtailment[gen] for gen in m.renewableGenerators
        )

    # Define total operating cost for dispatch block
    @b.Constraint()
    ## FIXME ... because of above mentioned curtailment issues?
    def operating_cost_dispatch(b):
        b.generationCost = Expression(
            expr=sum(
                b.thermalGeneration[gen] * m.fuelCost[gen]
                for gen in m.thermalGenerators
            )
        )
        b.loadshedwhatever = Expression(
            expr=sum(b.loadShed[bus] * m.loadShedCost for bus in m.buses)
        )
        return b.operatingCostDispatch == sum(
            b.thermalGeneration[gen] * m.fuelCost[gen] for gen in m.thermalGenerators
        ) + sum(b.loadShed[bus] * m.loadShedCost for bus in m.buses)
        # ) + sum(
        #     b.renewableCurtailment[gen] * m.curtailmentCost
        #     for gen in m.renewableGenerators
        # )


def add_commitment_variables(b, commitment_period):
    """Add variables and disjuncts to commitment period block.
    """
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

        # Ramp down constraints for generators shutting down
        @disj.Constraint(b.dispatchPeriods, m.thermalGenerators)
        def ramp_down_limits(disj, dispatchPeriod, generator):
            return (
                b.dispatchPeriod[dispatchPeriod - 1].thermalGeneration[generator]
                - b.dispatchPeriod[dispatchPeriod].thermalGeneration[generator]
                <= max(
                    m.thermalMin[generator],
                    m.rampDownRates[generator]
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
        def operating_limit_max(d, dispatchPeriod):
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

    # Track total renewable surplus/deficit for future expressions
    b.renewableSurplusCommitment = Var(within=Reals, initialize=0, units=u.MW)

    # Track total operating cost for commitment stage
    b.operatingCostCommitmentzzz = Var(within=Reals, initialize=0, units=u.USD)

    # Track total curtailment for objective function
    b.renewableCurtailmentCommitment = Var(
        within=NonNegativeReals, initialize=0, units=u.MW
    )


def add_commitment_constraints(
    b,
    comm_per,
):
    """Add commitment-associated disjunctions and constraints to representative period block.
    """
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    # Define total renewable surplus/deficit for commitment block
    @b.Constraint()
    def renewable_surplus_commitment(b):
        return b.renewableSurplusCommitment == sum(
            m.dispatchPeriodLength * b.dispatchPeriod[disp_per].renewableSurplusDispatch
            for disp_per in b.dispatchPeriods
        )

    # Define total operating costs for commitment block
    ## TODO: Replace this constraint with expressions using bounds transform
    ## NOTE: expressions are stored in gtep_cleanup branch
    ## FIXME: costs blowing up by dispatch length?
    ## costs considered need to be re-assessed and account for missing data
    @b.Expression()
    def operating_cost_commitment(b):
        b.dispatch_operating_cost = Expression(
            expr=sum(
                # (60 / m.dispatchPeriodLength) *
                b.dispatchPeriod[disp_per].operatingCostDispatch
                for disp_per in b.dispatchPeriods
            )
        )
        b.fixed_operating_cost = Expression(
            expr=sum(
                m.fixedOperatingCost[gen]
                # * m.thermalCapacity[gen]
                * (
                    b.genOn[gen].indicator_var.get_associated_binary()
                    + b.genShutdown[gen].indicator_var.get_associated_binary()
                    + b.genStartup[gen].indicator_var.get_associated_binary()
                )
                for gen in m.thermalGenerators
            )
            + sum(
                m.fixedOperatingCost[gen]
                # * m.renewableCapacity[gen]
                for gen in m.renewableGenerators
            )
        )
        b.startup_operating_cost = Expression(
            expr=sum(
                # m.commitmentPeriodLength
                # * m.thermalCapacity[gen]
                # * (
                # m.startFuel[gen] * m.fuelCost[gen]
                # + m.startFuel[gen] * b.carbonTax * m.emissionsFactor[gen]
                # +m.startupCost[gen]
                # )
                m.startupCost[gen]
                * b.genStartup[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
        )
        return (
            sum(
                # (60 / m.dispatchPeriodLength) *
                b.dispatchPeriod[disp_per].operatingCostDispatch
                for disp_per in b.dispatchPeriods
            )
            + sum(
                m.fixedOperatingCost[gen]
                # * m.thermalCapacity[gen]
                * (
                    b.genOn[gen].indicator_var.get_associated_binary()
                    + b.genShutdown[gen].indicator_var.get_associated_binary()
                    + b.genStartup[gen].indicator_var.get_associated_binary()
                )
                for gen in m.thermalGenerators
            )
            + sum(
                m.fixedOperatingCost[gen]
                # * m.renewableCapacity[gen]
                for gen in m.renewableGenerators
            )
            + sum(
                # m.commitmentPeriodLength
                # * m.thermalCapacity[gen]
                # * (
                # m.startFuel[gen] * m.fuelCost[gen]
                # + m.startFuel[gen] * b.carbonTax * m.emissionsFactor[gen]
                # +m.startupCost[gen]
                # )
                m.startupCost[gen]
                * b.genStartup[gen].indicator_var.get_associated_binary()
                for gen in m.thermalGenerators
            )
        )

    # Define total curtailment for commitment block
    @b.Constraint()
    def renewable_curtailment_commitment(b):
        return b.renewableCurtailmentCommitment == sum(
            b.dispatchPeriod[disp_per].renewableCurtailmentDispatch
            for disp_per in b.dispatchPeriods
        )


def commitment_period_rule(b, commitment_period):
    """Create commitment period block.

    :b: commitment period block
    :commitment_period: corresponding commitment period label
    """
    m = b.model()
    r_p = b.parent_block()
    i_p = r_p.parent_block()

    b.commitmentPeriod = commitment_period
    b.dispatchPeriods = RangeSet(m.numDispatchPeriods[r_p.currentPeriod])
    b.carbonTax = Param(default=0)
    b.dispatchPeriod = Block(b.dispatchPeriods)

    # update properties for this time period!!!!
    if m.data_list:
        m.md = m.data_list[i_p.representativePeriods.index(r_p.currentPeriod)]

    # Maximum output of each renewable generator
    m.renewableCapacity = {
        renewableGen: m.md.data["elements"]["generator"][renewableGen]["p_max"][
            "values"
        ][commitment_period - 1]
        for renewableGen in m.renewableGenerators
    }

    ## TODO: Redesign load scaling and allow nature of it as argument
    # Demand at each bus
    scale_loads = True
    if scale_loads:
        m.loads = {
            m.md.data["elements"]["load"][load_n]["bus"]: (
                3 + (3 + i_p.investmentStage) / (3 + len(m.stages))
            )
            * m.md.data["elements"]["load"][load_n]["p_load"]["values"][
                commitment_period - 1
            ]
            for load_n in m.md.data["elements"]["load"]
        }
    else:
        m.loads = {
            m.md.data["elements"]["load"][load_n]["bus"]: m.md.data["elements"]["load"][
                load_n
            ]["p_load"]["values"][commitment_period - 1]
            for load_n in m.md.data["elements"]["load"]
        }

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
    def shutdown_commitment(b, commitmentPeriod, thermalGen):
        return (
            b.commitmentPeriod[commitmentPeriod - 1]
            .genShutdown[thermalGen]
            .indicator_var.implies(
                b.commitmentPeriod[commitmentPeriod].genOff[thermalGen].indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def startup_commitment(b, commitmentPeriod, thermalGen):
        return (
            b.commitmentPeriod[commitmentPeriod - 1]
            .genStartup[thermalGen]
            .indicator_var.implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_activity(b, commitmentPeriod, thermalGen):
        return (
            b.commitmentPeriod[commitmentPeriod - 1]
            .genOn[thermalGen]
            .indicator_var.implies(
                b.commitmentPeriod[commitmentPeriod].genOn[thermalGen].indicator_var
                | b.commitmentPeriod[commitmentPeriod]
                .genShutdown[thermalGen]
                .indicator_var
            )
            if commitmentPeriod != 1
            else LogicalConstraint.Skip
        )

    @b.LogicalConstraint(b.commitmentPeriods, m.thermalGenerators)
    def consistent_commitment_inactivity(b, commitmentPeriod, thermalGen):
        return (
            b.commitmentPeriod[commitmentPeriod - 1]
            .genOff[thermalGen]
            .indicator_var.implies(
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
    :m: Pyomo GTEP model.
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

    @m.Objective()
    def total_cost_objective_rule(m):
        if len(m.stages) > 1:
            return m.operatingCost + m.expansionCost + m.penaltyCost
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
    :m: Pyomo model object
    """

    # Maximum output of each thermal generator
    m.thermalCapacity = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_max"]
        for thermalGen in m.thermalGenerators
    }

    # Lifetime of each generator; needs units
    m.lifetimes = {
        gen: m.md.data["elements"]["generator"][gen]["lifetime"] for gen in m.generators
    }

    # Minimum output of each thermal generator
    m.thermalMin = {
        thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_min"]
        for thermalGen in m.thermalGenerators
    }

    # Maximum output of each renewable generator
    m.renewableCapacity = {
        renewableGen: max(
            m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"]
        )
        for renewableGen in m.renewableGenerators
    }

    # A fraction of renewableCapacity representing fraction of capacity
    # that can be reliably counted toward planning reserve requirement
    # TODO: WHAT HAVE I DONE HERE I HATE IT
    m.renewableCapacityValue = {
        renewableGen: min(
            m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"]
        )
        / max(1, m.renewableCapacity[renewableGen])
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

    ## TODO: These should go into each stage -- check where these values should come from
    m.peakLoad = Param(m.stages, default=0, units=u.MW)
    m.reserveMargin = Param(m.stages, default=0, units=u.MW)
    m.renewableQuota = Param(m.stages, default=0, units=u.MW)
    m.weights = Param(m.representativePeriods, default=1)
    m.investmentFactor = Param(m.stages, default=1, mutable=True)
    ## NOTE: Lazy approx for NPV
    ## TODO: don't lazily approx NPV, add it into unit handling and calculate from actual time frames
    for stage in m.stages:
        m.investmentFactor[stage] *= 1 / ((1.04) ** (5 * stage))
    m.fixedOperatingCost = Param(m.generators, default=1, units=u.USD)
    m.deficitPenalty = Param(m.stages, default=1, units=u.USD / u.MW)

    # Amount of fuel required to be consumed for startup process for each generator
    m.startFuel = {
        gen: m.md.data["elements"]["generator"][gen]["start_fuel"]
        for gen in m.generators
    }

    # Cost per unit of fuel at each generator
    if "RTS-GMLC" in m.md.data["system"]["name"]:
        m.fuelCost = {
            gen: m.md.data["elements"]["generator"][gen]["fuel_cost"]
            for gen in m.thermalGenerators
        }
    else:
        m.fuelCost = {
            gen: m.md.data["elements"]["generator"][gen]["p_cost"]["values"][1]
            for gen in m.thermalGenerators
        }

    # Cost per MW of curtailed renewable energy
    # NOTE: what should this be valued at?  This being both curtailment and load shed.
    # TODO: update valuations
    m.curtailmentCost = 2 * max(m.fuelCost.values())
    m.loadShedCost = 1000 * m.curtailmentCost

    # Full lifecycle CO_2 emission factor for each generator
    m.emissionsFactor = {
        gen: m.md.data["elements"]["generator"][gen]["emissions_factor"]
        for gen in m.generators
    }

    # Flat startup cost for each generator
    if "RTS-GMLC" in m.md.data["system"]["name"]:
        m.startupCost = {
            gen: m.md.data["elements"]["generator"][gen]["non_fuel_startup_cost"]
            for gen in m.thermalGenerators
        }
    else:
        m.startupCost = {
            gen: m.md.data["elements"]["generator"][gen]["startup_cost"]
            for gen in m.generators
        }

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

    # Cost of investment in each new generator
    m.generatorInvestmentCost = {
        gen: m.md.data["elements"]["generator"][gen]["investment_cost"]
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

    m.investmentStage = Block(m.stages, rule=investment_stage_rule)

    # Retirement/extension relationships over investment periods -- C&P'd
    # from the paper.  These are okay.
    if len(m.stages) > 1:

        @m.Constraint(m.stages, m.thermalGenerators)
        def gen_retirement(m, stage, gen):
            return sum(
                m.investmentStage[t_2]
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

        @m.Constraint(m.stages, m.renewableGenerators)
        def renewable_retirement(m, stage, gen):
            return sum(
                m.investmentStage[t_2].renewableInstalled[gen]
                for t_2 in m.stages
                if t_2 <= stage - m.lifetimes[gen]
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

    # If a gen is disabled at time t-1, it must stay disabled or be installed at time t
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
