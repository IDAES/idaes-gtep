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

"""Sets and Parameters in the Generation and Transmission Expansion
Planning (GTEP) Model

"""

import logging

import pyomo.environ as pyo
from pyomo.environ import units as u
from collections import defaultdict

import gtep.model_library.storage as stor

logger = logging.getLogger("gtep.model_library.components")


def add_model_sets(m, stages, rep_per=["a", "b"], com_per=2, dis_per=2):
    """
    Creates Pyomo Sets necessary (convenient) for solving the GTEP model.

    :m: Pyomo model object
    :stages: Number of stages in investment horizon
    """

    # [TODO: make sure time units are both definable and consistent
    # without being forced.]
    m.stages = pyo.RangeSet(stages, doc="Set of planning periods")

    m.representativePeriods = pyo.Set(
        initialize=rep_per,
        doc="Set of representative periods for each planning period",
    )

    m.buses = pyo.Set(
        initialize=m.md.data["elements"]["bus"].keys(), doc="Individual buses"
    )

    m.regions = pyo.Set(
        initialize=(
            m.md.data["elements"]["bus"][bus]["area"]
            for bus in m.md.data["elements"]["bus"]
        ),
        doc="Regions / clusters of buses",
    )

    # NOTE: Right now, this means that branches can only be specified
    # entirely as standard or as dc ... not mix-and-match.
    if len(m.md.data["elements"]["branch"]) == 0:
        m.md.data["elements"]["branch"] = m.md.data["elements"]["dc_branch"]

    m.lines = pyo.Set(
        initialize=m.md.data["elements"]["branch"].keys(),
        doc="Individual transmission lines",
    )

    _toBranches = defaultdict(list)
    _fromBranches = defaultdict(list)
    for branch, data in m.md.data["elements"]["branch"].items():
        _toBranches[data["to_bus"]].append(branch)
        _fromBranches[data["from_bus"]].append(branch)

    m.branchByToBus = pyo.Set(m.buses, initialize=_toBranches)
    m.branchByFromBus = pyo.Set(m.buses, initialize=_fromBranches)

    m.generators = pyo.Set(
        initialize=m.md.data["elements"]["generator"].keys(),
        doc="All generators",
    )

    m.generatorsByBus = pyo.Set(
        m.buses,
        initialize={
            bus: [
                gen
                for gen in m.generators
                if m.md.data["elements"]["generator"][gen]["bus"] == bus
            ]
            for bus in m.buses
        },
    )

    m.thermalGenerators = pyo.Set(
        within=m.generators,
        initialize=(
            gen
            for gen in m.generators
            if m.md.data["elements"]["generator"][gen]["generator_type"] == "thermal"
        ),
        doc="Thermal generators; subset of all generators",
    )

    m.thermalGeneratorsByBus = pyo.Set(
        m.buses, initialize=lambda m, b: m.generatorsByBus[b] & m.thermalGenerators
    )

    if m.config["advanced_hydro"]:

        m.hydroGenerators = pyo.Set(
            within=m.generators,
            initialize=(
                gen
                for gen in m.generators
                if m.md.data["elements"]["generator"][gen]["unit_type"] == "HYDRO"
            ),
            doc="Hydropower generators; subset of all generators",
        )

        m.hydroGeneratorsByBus = pyo.Set(
            m.buses, initialize=lambda m, b: m.generatorsByBus[b] & m.hydroGenerators
        )

        m.renewableGenerators = pyo.Set(
            within=m.generators,
            initialize=(
                gen
                for gen in m.generators
                if (
                    m.md.data["elements"]["generator"][gen]["generator_type"]
                    == "renewable"
                    and m.md.data["elements"]["generator"][gen]["unit_type"] != "HYDRO"
                )
            ),
            doc="Renewable generators; subset of all generators",
        )

    else:

        m.renewableGenerators = pyo.Set(
            within=m.generators,
            initialize=(
                gen
                for gen in m.generators
                if m.md.data["elements"]["generator"][gen]["generator_type"]
                == "renewable"
            ),
            doc="Renewable generators; subset of all generators",
        )

    m.renewableGeneratorsByBus = pyo.Set(
        m.buses, initialize=lambda m, b: m.generatorsByBus[b] & m.renewableGenerators
    )

    m.load_buses = pyo.Set(initialize=[i for i in m.md.data["elements"]["load"]])

    # NOTE: We will want to cover baseline generator types in IDAES
    # This should be updated for battery. @JKS is this using the
    # built-in structure from EGRET or just a placeholder?
    if m.config["storage"]:
        m.storage = pyo.Set(
            initialize=(stor for stor in m.md.data["elements"]["storage"]),
            doc="Potential storage units",
        )

    m.storageByBus = pyo.Set(
        m.buses,
        initialize={
            bus: (
                [
                    batt
                    for batt in m.storage
                    if m.md.data["elements"]["storage"][batt]["bus"] == bus
                ]
                if m.config["storage"]
                else []
            )
            for bus in m.buses
        },
    )


def add_model_parameters(m, years):
    """Creates and labels all the parameters in the GTEP model. This
    method ties input data directly to the model.

    :param m: Pyomo model object

    """

    # Add investment years. This data is coming from the model data
    # processing object.
    m.years = years

    # Add power-related parameters
    m.thermalCapacity = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_max"]
            for thermalGen in m.thermalGenerators
        },
        mutable=True,
        units=u.MW,
        doc="Maximum output of each thermal generator",
    )
    if m.config["flow_model"] == "ACR" or m.config["flow_model"] == "ACP":
        m.thermalReactiveMax = pyo.Param(
            m.thermalGenerators,
            initialize={
                thermalGen: m.md.data["elements"]["generator"][thermalGen]["q_max"]
                for thermalGen in m.thermalGenerators
            },
            mutable=True,
            units=u.MVAR,
            doc="Maximum reactive output of each thermal generator",
        )

        m.thermalReactiveMin = pyo.Param(
            m.thermalGenerators,
            initialize={
                thermalGen: m.md.data["elements"]["generator"][thermalGen]["q_min"]
                for thermalGen in m.thermalGenerators
            },
            mutable=True,
            units=u.MVAR,
            doc="Minimum reactive output of each thermal generator",
        )

    if m.config["advanced_hydro"]:
        m.hydroCapacity = pyo.Param(
            m.hydroGenerators,
            initialize={
                gen: max(m.md.data["elements"]["generator"][gen]["p_max"]["values"])
                for gen in m.hydroGenerators
            },
            domain=pyo.NonNegativeReals,
            mutable=True,
            units=u.MW,
            doc="Maximum output of each hydropower generator",
        )

    m.generatorLifetime = pyo.Param(
        m.generators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["lifetime"]
            for gen in m.generators
        },
        mutable=True,
        units=u.year,
        doc="Lifetime of each generator",
    )

    m.branchLifetime = pyo.Param(
        m.lines,
        initialize={line: 3 for line in m.lines for line in m.lines},
        mutable=True,
        units=u.year,
        doc="Lifetime of each branch",
    )

    m.thermalMin = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen]["p_min"]
            for thermalGen in m.thermalGenerators
        },
        mutable=True,
        units=u.MW,
        doc="Minimum output of each thermal generator",
    )

    m.renewableCapacityNameplate = pyo.Param(
        m.renewableGenerators,
        initialize={
            renewableGen: (
                m.md.data["elements"]["generator"][renewableGen]["p_max"]
                if type(m.md.data["elements"]["generator"][renewableGen]["p_max"])
                == float
                else max(
                    [
                        max(
                            m.data_list[i].data["elements"]["generator"][renewableGen][
                                "p_max"
                            ]["values"]
                        )
                        for i in range(len(m.data_list))
                    ]
                )
            )
            for renewableGen in m.renewableGenerators
        },
        # mutable is set to False because it is used to set the
        # capacity factor and should not change
        mutable=False,
        units=u.MW,
        doc="Maximum output of each renewable generator",
    )

    # NOTE: From case data, the value is divided by 100, which is the
    # per units conversion.
    m.transmissionCapacity = pyo.Param(
        m.lines,
        initialize={
            transmissionLine: m.md.data["elements"]["branch"][transmissionLine][
                "rating_long_term"
            ]
            for transmissionLine in m.lines
        },
        units=u.MW,
        doc="Long term thermal rating of each transmission line",
    )

    m.from_bus = pyo.Param(
        m.lines,
        initialize={
            branch: m.md.data["elements"]["branch"][branch]["from_bus"]
            for branch in m.lines
        },
        within=m.buses,
        doc="Start bus for each transmission line",
    )

    m.to_bus = pyo.Param(
        m.lines,
        initialize={
            branch: m.md.data["elements"]["branch"][branch]["to_bus"]
            for branch in m.lines
        },
        within=m.buses,
        doc="End bus for each transmission line",
    )

    m.reactance = pyo.Param(
        m.lines,
        initialize={
            branch: m.md.data["elements"]["branch"][branch]["reactance"]
            for branch in m.lines
        },
        domain=pyo.Reals,
        doc="Reactance for each transmission line (ohms)",
    )

    m.spinningReserveFraction = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen][
                "spinning_reserve_frac"
            ]
            for thermalGen in m.thermalGenerators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Maximum fraction of maximum thermal generation output that can be supplied as spinning reserve",
    )

    m.quickstartReserveFraction = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen][
                "quickstart_reserve_frac"
            ]
            for thermalGen in m.thermalGenerators
        },
        units=u.dimensionless,
        doc="Maximum fraction of maximum thermal generation output that can be supplied as quickstart reserve",
    )

    # Initialize the m.loads parameter with a value of 0 and scaled it
    # with the right value in the commitment stage.
    m.loads = pyo.Param(
        m.buses,
        initialize={load_n: 0 for load_n in m.buses},
        mutable=True,
        units=u.MW,
        doc="Demand at each bus",
    )

    # [TODO: Fixing for dc_branch and branch, but we should revisit
    # this.]
    m.distance = pyo.Param(
        m.lines,
        initialize={
            branch: (m.md.data["elements"]["branch"][branch].get("distance") or 0)
            for branch in m.lines
        },
        mutable=True,
        units=u.m,
        doc="Distance between terminal buses for each transmission line",
    )

    # [TODO: These should go into each stage. Check where these values
    # should come from.]
    m.peakLoad = pyo.Param(m.stages, default=0, units=u.MW)
    m.reserveMargin = pyo.Param(m.stages, default=0, units=u.MW)
    m.renewableQuota = pyo.Param(m.stages, default=0, units=u.MW)

    # Map each representative period to its corresponding date, then
    # use that date to retrieve the appropriate representative
    # weight. This keeps the model indexed by representative period
    # number.
    representative_dates_dict = {
        i: date for i, date in zip(m.representativePeriods, m.data.representative_dates)
    }
    m.representativeDate = pyo.Param(
        m.representativePeriods,
        initialize=representative_dates_dict,
        within=pyo.Any,
        mutable=False,
        doc="Representative date associated with each representative period",
    )
    weights_dict = {
        i: m.data.representative_weights_dict[representative_dates_dict[i]]
        for i in m.representativePeriods
    }
    m.weights = pyo.Param(
        m.representativePeriods,
        initialize=weights_dict,
        mutable=False,
        doc="Representative period weights indexed by representative period",
    )

    m.investmentFactor = pyo.Param(
        m.stages, default=1, mutable=True, units=u.dimensionless
    )
    m.deficitPenalty = pyo.Param(m.stages, default=1, units=u.USD / u.MW)

    m.emissionsFactor = pyo.Param(
        m.generators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["emissions_factor"]
            for gen in m.generators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Full lifecycle CO_2 emission factor for each generator",
    )

    # Include start-up cost only in thermal generators assuming a
    # natural gas plant.
    m.startupCost = pyo.Param(
        m.thermalGenerators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["non_fuel_startup_cost"]
            for gen in m.thermalGenerators
        },
        mutable=True,
        units=u.USD,
        doc="Flat startup cost for thermal generators",
    )

    # (Arbitrary) multiplier corresponds to depreciation schedules for
    # individual technologies; higher values are indicative of slow
    # depreciation
    m.capitalMultiplier = pyo.Param(
        m.generators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["capital_multiplier"]
            for gen in m.generators
        },
        mutable=True,
        units=u.dimensionless,
        doc="(Arbitrary) multiplier for new generator investments",
    )

    m.extensionMultiplier = pyo.Param(
        m.generators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["extension_multiplier"]
            for gen in m.generators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life extension for each generator expressed as a fraction of initial investment cost",
    )

    # [BLN TODO: Check what value should be used here]
    m.retirementMultiplier = pyo.Param(
        m.generators,
        initialize={
            gen: (0.1 if gen in m.thermalGenerators else 1.0) for gen in m.generators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life retirement for each generator expressed as a fraction of initial investment cost",
    )

    m.minOperatingReserve = pyo.Param(
        m.regions,
        initialize={
            region: m.md.data["system"]["min_operating_reserve"] for region in m.regions
        },
        mutable=True,
        units=u.dimensionless,
        doc="Minimum operating reserve as a fraction of load within a region",
    )

    m.minSpinningReserve = pyo.Param(
        m.regions,
        initialize={
            region: m.md.data["system"]["min_spinning_reserve"] for region in m.regions
        },
        mutable=True,
        units=u.dimensionless,
        doc="Minimum spinning reserve as a fraction of load within a region",
    )

    m.maxSpinningReserve = pyo.Param(
        m.thermalGenerators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["max_spinning_reserve"]
            for gen in m.thermalGenerators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Maximum spinning reserve available for each generator as a fraction maximum generator output",
    )

    m.maxQuickstartReserve = pyo.Param(
        m.thermalGenerators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["max_quickstart_reserve"]
            for gen in m.thermalGenerators
        },
        units=u.dimensionless,
        doc="Maximum quickstart reserve available for each generator as a fraction maximum generator output",
    )

    m.rampUpRates = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_up_rate"]
            for thermalGen in m.thermalGenerators
        },
        units=u.MW / u.minutes,
        doc="Ramp up rates for each generator",
    )

    m.rampDownRates = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"]
            for thermalGen in m.thermalGenerators
        },
        units=u.MW / u.minutes,
        doc="Ramp down rates for each generator",
    )

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

    # Initialize curtailment and load shed costs as parameters and
    # re-calculate them during the investment stage.
    m.curtailmentCost = pyo.Param(
        initialize=0,
        units=u.USD / (u.MW * u.hr),
        mutable=True,
        doc="Curtailment cost",
    )
    m.loadShedCostperCurtailment = pyo.Param(
        initialize=30000,
        units=u.USD / (u.MW * u.hr),
        mutable=True,
    )

    # Add cost parameters
    add_gen_branch_cost_parameters(m)
    if m.config["storage"] == True:
        stor.add_storage_params(m)

    # Add legacy parameters.These parameters are commented in the
    # original model. Keep here to check if we should include them in
    # future version of the model.
    """
    # NOTE: Commented for now since it is not use in this case but
    # might be used in the future when considering ACOPF.
    # m.lossRate = pyo.Param(
    #     m.lines,
    #     initialize={branch: (m.md.data["elements"]["branch"][branch].get("loss_rate") or 0)
    #                 for branch in m.lines},
    #     mutable=True,
    #     # units=,
    #     doc="Per-distance-unit multiplicative loss rate for each transmission line"
    # )

    # Simple approximation for NPV. NOTE: This is commented for now
    # since it is already included in the costs we have from
    # preprocessing stage.
    for stage in m.stages:
        m.investmentFactor[stage] *= 1 / ((1.04) ** (5 * stage))

    # NOTE: Commented for now but depends on the type of data we
    # are using for generators.
    m.startFuel = pyo.Param(
        m.generators,
        initialize={gen: m.md.data["elements"]["generator"][gen]["start_fuel"]
                    for gen in m.generators},
        mutable=True,
        doc="Amount of fuel required to be consumed for startup process for each generator"
    )
    """


def add_gen_branch_cost_parameters(m):
    """Adds cost parameters for transmission lines, generators, and
    thermal generators fuel costs for stage years. The costs
    parameters support multi-year planning with year-indexed values.

    The cost data source depends on configuration values:

        - If include_cost_data=True and
          use_cost_from_data_files=False, it adds costs from
          pre-converted spreadsheet data from m.mc.gen_cost_data.

        - If include_cost_data=True and use_cost_from_data_files=True,
          it adds {cost_type}_{year} data from input data files
          'gen.csv' and 'branch.csv' files.

    The required units for relevant costs are:
        fixed_ops_cost_<year>     : USD/MWh
        var_ops_cost_<year>       : USD/MWh
        capex_<year>              : USD/MW
        fuel_cost_<year>          : USD/MWh
        fuel_reactive_cost_<year> : USD/MWh

    """

    m.years_set = pyo.Set(initialize=m.years)

    def get_branch_cost_init(m, years, cost_type, default_value=0):
        """Creates branch cost initialization dictionary.

        :m:             Pyomo model object
        :years:         List of years to initialize costs for
        :cost_type:     String indicating which cost type to get
                        ('capex', 'fixed_ops_cost', or 'var_ops_cost')
        :default_value: Default value to use when cost data isn't
                        available

        returns: Dictionary with (branch, year) keys and cost values

        """

        init_dict = {}

        if m.mc is not None and hasattr(m.mc, "branch_cost_from_csv_data"):

            # Use cost data from 'branch.csv' file with the style
            # {cost_type}_{year}.w
            for index, row in m.mc.branch_cost_from_csv_data.iterrows():
                branch = row["UID"]
                if branch in m.lines:
                    for year in years:
                        # Build the column name based on cost_type and
                        # year.
                        column_name = f"{cost_type}_{year}"
                        init_dict[(branch, year)] = float(row[column_name])
        else:

            # Use capital/investment cost data from 'branch.csv' file
            # from column 'capital_cost'.
            if cost_type == "capex":
                for branch in m.lines:
                    capital_cost = m.md.data["elements"]["branch"][branch].get(
                        "capital_cost", default_value
                    )
                    for year in years:
                        init_dict[(branch, year)] = capital_cost
            else:

                # Use default value if no other value is available.
                for branch in m.lines:
                    for year in years:
                        init_dict[(branch, year)] = default_value

        return init_dict

    m.branchInvestmentCost = pyo.Param(
        m.lines,
        m.years_set,
        initialize=get_branch_cost_init(m, m.years, "capex", 0),
        mutable=True,
        units=u.USD / u.MW,
        doc="Investment cost for each branch by year",
    )

    m.branchFixedCost = pyo.Param(
        m.lines,
        m.years_set,
        initialize=get_branch_cost_init(m, m.years, "fixed_ops_cost", 0),
        mutable=True,
        units=u.USD / u.MW,
        doc="Fixed cost for each branch by year",
    )

    m.branchVariableCost = pyo.Param(
        m.lines,
        m.years,
        initialize=get_branch_cost_init(m, m.years, "var_ops_cost", 0),
        mutable=True,
        units=u.USD / u.MW,
        doc="Variable cost for each branch by year",
    )

    # Branch multipliers remain constant across years
    m.branchCapitalMultiplier = pyo.Param(
        m.lines,
        initialize={
            branch: (
                m.md.data["elements"]["branch"][branch].get("capital_multiplier") or 1
            )
            for branch in m.lines
        },
        mutable=True,
        units=u.dimensionless,
    )

    m.branchExtensionMultiplier = pyo.Param(
        m.lines,
        initialize={
            branch: (
                m.md.data["elements"]["branch"][branch].get("extension_multiplier") or 1
            )
            for branch in m.lines
        },
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life extension for each branch expressed as a fraction of initial investment cost",
    )

    def get_generator_cost_init(m, years, cost_type, default_value=0):
        """Create generator cost initialization dictionary.

        :m:             Pyomo model object
        :years:         List of years to initialize costs for
        :cost_type:     String indicating which cost type to get
                        ('capex', 'fixed_ops_cost', or 'var_ops_cost')
        :default_value: Default value to use when cost data isn't available

        returns:        Dictionary with (generator, year) keys and cost
        values

        """
        init_dict = {}

        if cost_type == "fuel_cost":
            set_gens = m.thermalGenerators
        else:
            set_gens = m.generators

        if m.mc is not None and hasattr(m.mc, "gen_cost_from_csv_data"):
            # Use cost data from 'gen.csv' with the style
            # '{cost_type}_{year}`.
            for index, row in m.mc.gen_cost_from_csv_data.iterrows():
                gen = row["GEN UID"]
                if gen in set_gens:
                    for year in years:
                        column_name = f"{cost_type}_{year}"
                        init_dict[(gen, year)] = float(row[column_name])

        elif m.mc is not None and hasattr(m.mc, "gen_cost_data"):
            # Use standardized cost data based on the generator type.
            cost_df = m.mc.gen_cost_data.copy()
            cost_df["Unit Type"] = cost_df["Unit Type"].astype(str).str.upper()

            def get_cost_row_by_unit_type(gen):
                """Gets cost row based on generator type."""

                gen_data = m.md.data["elements"]["generator"][gen]
                unit_type = str(gen_data.get("unit_type", "")).strip().upper()
                generator_type = str(gen_data.get("generator_type", "")).strip().lower()

                if generator_type == "renewable":
                    cost_unit_type = "PV"
                elif generator_type == "thermal":
                    cost_unit_type = "CT"
                else:
                    cost_unit_type = "CT"

                cost_rows = cost_df[cost_df["Unit Type"] == cost_unit_type]
                return cost_rows.iloc[0]

            for gen in set_gens:
                cost_row = get_cost_row_by_unit_type(gen)
                for year in years:
                    column_name = f"{cost_type}_{year}"
                    init_dict[(gen, year)] = float(cost_row[column_name])
        else:

            # Use default value when no cost data is available.
            for gen in set_gens:
                for year in years:
                    if cost_type == "fuel_cost":
                        init_dict[(gen, year)] = (
                            m.fuelCostperMMBTU[gen] * m.heatRate[gen]
                        )
                    else:
                        init_dict[(gen, year)] = default_value

        return init_dict

    m.generatorInvestmentCost = pyo.Param(
        m.generators,
        m.years_set,
        initialize=get_generator_cost_init(m, m.years, "capex", 0),
        mutable=True,
        units=u.USD / u.MW,
        doc="Investment cost for all generators by year",
    )

    m.generatorFixedCost = pyo.Param(
        m.generators,
        m.years_set,
        initialize=get_generator_cost_init(m, m.years, "fixed_ops_cost", 0),
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Fixed operating costs for all generators by year",
    )

    m.generatorVariableCost = pyo.Param(
        m.generators,
        m.years_set,
        initialize=get_generator_cost_init(m, m.years, "var_ops_cost", 0),
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Variable costs for all generators by year",
    )

    # Fuel costs are calculated from heat rate and fuel price, but can
    # be overridden by cost data.
    m.fuelCostperMMBTU = pyo.Param(
        m.thermalGenerators,
        initialize={
            gen: (
                m.md.data["elements"]["generator"][gen]["fuel_cost"]
                if "RTS-GMLC" in m.md.data["system"]["name"]
                else m.md.data["elements"]["generator"][gen]["p_cost"]["values"][1]
            )
            for gen in m.thermalGenerators
        },
        mutable=True,
        units=u.USD / u.MMBTU,
        doc="Fuel cost per MMBTU at each generator",
    )

    m.heatRate = pyo.Param(
        m.thermalGenerators,
        initialize={
            gen: (
                m.md.data["elements"]["generator"][gen]["heat_rate"]
                if "RTS-GMLC" in m.md.data["system"]["name"]
                else m.md.data["elements"]["generator"][gen]["heat_rate"]
            )
            for gen in m.thermalGenerators
        },
        mutable=True,
        units=u.MMBTU / (u.MW * u.hr),
        doc="Heat rate for each thermal generator",
    )

    m.fuelCost = pyo.Param(
        m.thermalGenerators,
        m.years_set,
        initialize=get_generator_cost_init(m, m.years, "fuel_cost", 0),
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Cost per unit of fuel at each generator by year",
    )

    m.fuelCostReactive = pyo.Param(
        m.thermalGenerators,
        m.years_set,
        initialize=get_generator_cost_init(m, m.years, "fuel_cost", 0),
        mutable=True,
        units=u.USD / (u.MVAR * u.hr),
        doc="Cost per unit of fuel at each generator (reactive power) by year",
    )
