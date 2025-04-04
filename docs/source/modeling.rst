.. _Modeling:

Modeling 
=========

The `ExpansionPlanningModel()` class implements a modular and flexible
Generalized Disjunctive Programming (GDP) formulation for power
infrastructure planning problems. This formulation establishes a
generalized model for transmission and generation expansion planning
and it is represented as a dispatch problem integrated within a unit
commitment and investment problem for a specified representative
period. in this formulation, each of these problems corresponds to a
distinct stage with the commitment stage serving as the primary link
between them.

Figure 1 illustrates a graphical representation of these stages, their
components (variables and constraints), and the connections within the
formulation.

.. figure:: ../images/model_stages.png
   :align: center

   Figure 1. Modular representation of the three stages included in
   the GTEP model formulation.

In the following sections, we will provide a detailed explanation of
each stage and its components.

Steps for Model Construction
----------------------------

To develop the Generation Transmission and Expansion Planning (GTEP)
model depicted in Figure 1, we divide the formulation into four key
stages:

1. **Data allocation stage**: Defines all components and parameters of
   the model and assigns their respective values. (Note that this
   stage is not represented in Figure 1.)
  
2. **Investment stage**: Estimates all investment variables and
   constraints, including representative costs for both thermal and
   renewable generators as well as transmission equipment (e.g.,
   transmission lines). Additionally, it encompasses a set of discrete
   decisions (represented as disjunctions through the GDP formulation)
   to determine the operational state of generation and transmission
   equipment (retired, disabled, etc.) that can potentially influence
   the investment costs.

3. **Commitment stage**: Determines the operational status of the
   generators via discrete decisions using disjunctions. It takes into
   account relevant operational costs and consistent constraints to
   ensure reliable power generation service and establishes operating
   limits considering the availability of the equipment from the
   dispatch stage.

4. **Dispatch stage**: Defines the status of transmission lines using
   discrete decisions, considering operational costs, operability
   limits (such as power flow, ramping, reserves, etc.), and
   consistent constraints to ensure the effective utilization of
   active power lines.

The GTEP model is constructed by following the steps outlined
below. The first four steps correspond directly to the stages
described above, and in some steps, the description is accompanied by
example code.


(1) **Data allocation stage**

(a) The initial step to create this model is to import the necessary
    libraries.
   
    .. code-block::

       import json
       import numpy as np
       import math
       from math import ceil
   
       # Import Pyomo components
       from pyomo.environ import * 
       from pyomo.environ import units as u
       from pyomo.common.timing import TicTocTimer
       from pyomo.repn.linear import LinearRepnVisitor
   
       from egret.data.model_data import ModelData
       from egret.model_library.transmission.tx_utils import scale_ModelData_to_pu
       from config_options import _get_model_config

       # Import IDAES-GTEP functions 
       from gtep.gtep_model import model_data_references
       from gtep.gtep_model import model_create_investment_stages
       from gtep.gtep_model import create_objective_function

(b) Continue the creation of the GTEP model by including a new Python
    `class`, which includes a set of parameters as arguments of the
    function.

    .. code-block::

       class ExpansionPlanningModel:
            """A generalized generation and
               transmission expansion planning model.
	    """

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
            
                self.stages = stages
                self.formulation = formulation
                self.data = data
                self.num_reps = num_reps
                self.len_reps = len_reps
                self.num_commit = num_commit
                self.num_dispatch = num_dispatch
                self.config = _get_model_config()
                self.timer = TicTocTimer()

(c) Inside the `ExpansionPlanningModel()` class, we begin by creating
    a concrete Pyomo model object using the `ConcreteModel`
    component. Within this model, we define several time-dependent
    parameters for each stage, including:

    (1) **Number of Commitment and Dispatch Periods**: These
        parameters specify the total number of periods allocated for
        the commitment and dispatch stages, respectively. These are
        key for determining the operational status of generators over
        time and how power flows are managed across the transmission
        network.
    (2) **Commitment and Dispatch Period Lengths**: These parameters
        define the duration of each commitment and dispatch period,
        respectively. These allow us to model the time intervals
        during which generators can change their status (on, off,
        etc.) and evaluate operational decisions made during power
        dispatch.

    .. code-block::

  	    def create_model(self):
            """Create concrete Pyomo model object associated with the ExpansionPlanningModel"""

            self.timer.tic("Creating GTEP Model")
            m = ConcreteModel()

            if self.data is None:
                raise
	    elif type(self.data) is list:
	        m.data_list = self.data
                m.md = scale_ModelData_to_pu(self.data[0])
        
            model_set_declaration(
                m, self.stages, 
                rep_per=[i for i in range(1, self.num_reps + 1)]
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
            m.dispatchPeriodLength = Param(
                within=PositiveReals, 
                default=15, 
                units=u.min
            )

(d) Create and label data in the model using a pre-defined
    function. This function ties the input data directly to the model
    by assigning values to all parameters. Refer to :ref:`Data` in
    this documentation to know more details about this procedure.

    .. code-block::

       model_data_references(m)

 
(2) **Investment stage**

(a) The investment stage is implemented as a Pyomo block using the
    `Block` component. This block includes all relevant cost variables
    and constraints associated with generation and transmission
    equipment including operating, expansion, maintenance, and penalty
    costs. Additionally, it incorporates a set of discrete decisions
    (represented as disjunctions) to establish the operational status
    of generation and transmission equipment, while considering a
    series of linking constraints to monitor the operational,
    installed, and retired states of the equipment. In the GTEP model,
    we include a pre-defined function that incorporates all these
    components.

    .. code-block::

       model_create_investment_stages(m, self.stages)


    It is during the calculation of the operating costs that we
    establish a connection with the commitment and dispatch stages. To
    facilitate this, we define a second block specifically dedicated
    to the nested commitment stage. Figure 2 illustrates these
    interconnections and how these relationships propagate throughout
    the formulation, while Table 1 presents some of these key
    components along with their definition.
    
    .. figure:: ../images/model_investment_stage.png
       :align: center

       Figure 2. Modular representation of the investment stage.

    .. table:: Table 1: Variables and constraints in the investment stage
       :widths: 20 15 35

       ==================================== =========================== ===========================================================================
       Component                            Type                        Definition
       ==================================== =========================== ===========================================================================
       `representativePeriods`              Set                         Representative periods for the investment stage
       `commitmentPeriods`                  Set                         Commitment periods in the investment stage
       `dispatchPeriods`                    Set                         Dispatch periods in the investment stage
       `renewableSurplusRepresentative`     Variable                    Minimum pr-stage renewable generation requirement
       `renewableOperational`               Variable
       `renewableInstalled`                 Variable
       `renewableRetired`                   Variable
       `renewableExtended`                  Variable
       `quotaDeficit`                       Variable                   
       `operatingCostInvestment`            Variable
       `expansionCost`                      Variable
       `renewableCurtailmentInvestment`     Variable
       `renewableCapacity`                  Parameter                   Maximum output of each renewable
       `maxThermalInvestment`               Parameter                   Maximum investment for thermal generators
       `maxRenewableInvestment`             Parameter                   Maximum investment for renewable generators
       `renewable_generation_requirement`   Constraint                  Minimum per-stage renewable generation requirement
       `operatingCostInvestment`            Constraint                  Operating costs for investment period
       `investment_cost`                    Constraint                  Investment cost for each investment period
       `renewable_curtailment_cost`         Constraint                  Curtailment penalties for investment period
       `investStatus`                       Disjunction                 Discrete decision for the generation investment status. The alternatives
                                                                        (as disjuncts) are: `genOPerational`,
                                                                        `genInstalled`, `genRetired`, `genDisabled`, and `genExtended`
       `branchinvestStatus`                 Disjunction                 Discrete decision for the transmission equipment investment status. The
                                                                        alternatives (as disjuncts) are: `branchOPerational`, `branchInstalled`,
									`branchRetired`, `branchDisabled`, and `branchExtended`
       `gen_stats_link`                     Linking constraint          Generation investment status that depends on the state of the equipment:
                                                                        operational, installed, and retired
       `renewable_stats_link`               Linking constraint          Total renewable generation for each stage that depends on the state of the
                                                                        equipment (operational, installed, and retired)
       `renewable retirement`               Linking constraint          Renewable retirement that depends on the installed, retired, or extended
                                                                        and lifetime of renewable equipment
	`consistent_extended`               Logical constraint          If a generator is extended at time `t`, it must stay extended or retired at
	                                                                time `t`
       `consistent_branch_operation`        Logical constraint          If a branch is online at time `t`, it must have been on or installed at
                                                                        time t-1
       `consistent_branch_operation_future` Logical constraint          If a branch is online at time `t`, it must be on, extended, or retired at
                                                                        time `t+1`
       `consistent_branch_extended`         Logical constraint          If a branch is extended at time `t-1`, it must stay extended or be retired
                                                                        at time `t`
       `consistent_operation`               Logical constraint          If a generator is online at time `t`, it should be online or installed at
                                                                        time `t-1`
       `consistent_operation_future`        Logical constraint          If a generator is online at time `t`, it can be online, extended, or
                                                                        retired at time `t+1`
       `consistent_branch_disabled`         Logical constraint          If a branch is disabled at time `t-1`, it must stay disabled or be
                                                                        installed at time `t`
       `full_retirement`                    Logical constraint          If a generator is retired at time `t-1`, it should be disable at time `t`
       `full_investment`                    Logical constraint          Installation in period `t-1` implies operational in period `t`
       `full_branch_retirement`             Logical constraint          Retirement in period `t-1` implies disabled in period `t`
       `full_branch_investment`             Logical constraint          Installation in period `t-1` implies operational in period `t`
       ==================================== =========================== ===========================================================================
       
(3) **Commitment stage**

(a) Within the commitment block defined in step 2(a), we incorporate
    the commitment periods and introduce relevant components to ensure
    the operation of transmission and generation equipment aligns with
    their designated statuses. During this stage, we also introduce
    cost and operational constraints that are contingent upon the
    commitment decisions and the outcomes of the dispatch stage.

    Figure 3 illustrates the various components utilized in this
    stage, highlighting the terms that connect to the dispatch
    stage. Additionally, Table 3 provides a detailed overview of all
    the components included in the stage model formulation.
    
    .. figure:: ../images/model_commitment_stage.png
       :align: center

       Figure 3. Modular representation of the commitment stage.

    .. table:: Table 2: Components in the commitment stage
       :widths: 20 15 35

       ============================================= =================== =============================================================================
       Component                                     Type                Description
       ============================================= =================== =============================================================================
       `dispatchPeriods`                             Set                 Number of dispatch periods
       `carbonTax`                                   Parameter           Carbon tax
       `renewableSurplusCommitment`                  Expression          Total renewable surplus/deficit for commitment
       `operatingCostCommitment`                     Expression          Define total operating costs for commitment
       `renewableCurtailmentCommitment`              Expression          Define total curtailment for commitment
       `genStatus`                                   Disjunction         Discrete decision to set generation equipment status. The alternatives
                                                                         are: `genOn`, `genStartup`, `genShutDown`, and `genOff`.
       `consistent_commitment_shutdown`              Logical constraint
       `consistent_commitment_off_after_shutdown`    Logical constraint
       `consistent_commitment_startup`               Logical constraint
       `consistent_commitment_on_after_startup`      Logical Constraint
       `consistent_commitment_uptime`                Logical Constraint
       `consistent_commitment_shutdown_after_uptime` Logical Constraint
       `consistent_commitment_downtime`              Logical Constraint 
       `consistent_commitment_start_after_downtime`  Logical Constraint
       `commit_active_gens_only`                     Logical constraint  Generators cannot be committed unless they are operational or just installed
       ============================================= =================== =============================================================================

 
(4) **Dispatch Stage**

(a) This stage incorporates all key variables and constraints
    associated to the dispatch problem, including operational limits
    for generation and curtailment, as well as definitions for load
    shedding and surplus dispatch. Additionally, this stage outlines
    the costs associated with these variables.

    Figure 4 illustrates some of these components, while Table 3
    provides a more detailed description of their definition.
    
    .. figure:: ../images/model_dispatch_stage.png
       :align: center

       Figure 4. Modular representation of the dispatch stage.

  
    .. table:: Table 3: Components in dispatch stage
       :widths: 20 15 35

       ================================ =================== =========================================================================================
       Component                        Type                Description
       ================================ =================== =========================================================================================
       `thermalGeneration`              Variable            Thermal generation power during the dispatch stage
       `renewableGeneration`            Variable            Renewable generation power during the dispatch stage
       `renewableCurtailment`           Variable            Amount of curtailed renewable power
       `spinningReserve`                Variable            
       `quickstartReserve`              Variable
       `powerFlow`                      Variable             Power flow
       `renewableGenerationSurplus`     Expression           Power surplus per renewable generator
       `renewableCurtailmentCost`       Expression           Cost of curtailment per renewable generator
       `generatorCost`                  Expression           Thermal generator cost based on the fuel cost
       `loadShedCost`                   Expression           Load shed cost per bus
       `flow_balance`                   Constraint           Energy power balance
       `branchInUseStatus`              Disjunction          Discrete decision to determine the status of the transmission equipment. The alternatives
                                                             are: `branchInUse` (which includes the DC OPF equations) and `branchNotInUse` (which
							     fixes power flow to zero).
       `must_use_active_branches`       Logical constraint   If a branch is in use, it must be active
       `cannot_use_inactive_branches`   Logical constraint   If a branch is not in use, it must be inactive
       ================================ =================== =========================================================================================


(5) Once the data and blocks for the different stages are
    incorporated, we define the total cost as the objective function
    using the Pyomo `Objective` component. The total cost includes not
    only investment costs, but also operating, expansion, and penalty
    costs that connect all the stages in the formulation. The
    objective function is integrated within the following pre-defined
    function.

    .. code-block::

       create_objective_function(m)


.. currentmodule:: gtep.gtep_model

.. automodule:: gtep
    :members:

.. automodule:: gtep.gtep_model
    :members: 

