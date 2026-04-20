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

import pyomo.environ as pyo
from pyomo.environ import units as u

import gtep.model_library.storage as stor


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

    m.transmission = {
        branch: {
            "from_bus": m.md.data["elements"]["branch"][branch]["from_bus"],
            "to_bus": m.md.data["elements"]["branch"][branch]["to_bus"],
            "reactance": m.md.data["elements"]["branch"][branch]["reactance"],
        }
        for branch in m.md.data["elements"]["branch"]
    }

    m.generators = pyo.Set(
        initialize=m.md.data["elements"]["generator"].keys(),
        doc="All generators",
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

    m.renewableGenerators = pyo.Set(
        within=m.generators,
        initialize=(
            gen
            for gen in m.generators
            if m.md.data["elements"]["generator"][gen]["generator_type"] == "renewable"
        ),
        doc="Renewable generators; subset of all generators",
    )

    m.lines = pyo.Set(
        initialize=m.transmission.keys(), doc="Individual transmission lines"
    )

    m.load_buses = pyo.Set(initialize=[i for i in m.md.data["elements"]["load"]])

    # NOTE: We will want to cover baseline generator types in IDAES
    # This should be updated for battery. @JKS is this using the
    # built-in structure from EGRET or just a placeholder?
    if m.md.data["elements"].get("storage"):
        m.storage = pyo.Set(
            initialize=(ess for ess in m.md.data["elements"]["storage"]),
            doc="Potential storage units",
        )


def add_model_parameters(m, num_commit, num_dispatch, duration_dispatch):
    """Creates and labels all the parameters in the GTEP model. This
    method ties input data directly to the model.

    :param m: Pyomo model object

    """

    # Add investment years. [TODO: Make sure this value comes from a
    # configuration arg and not hardcoded values.]
    m.years = [2025, 2030, 2035]

    # Add parameters related to the representative periods for the
    # different stages
    m.representativePeriodLength = pyo.Param(
        m.representativePeriods, within=pyo.PositiveReals, default=24, units=u.hr
    )
    m.numCommitmentPeriods = pyo.Param(
        m.representativePeriods,
        within=pyo.PositiveIntegers,
        default=2,
        initialize=num_commit,
    )
    m.numDispatchPeriods = pyo.Param(
        m.representativePeriods,
        within=pyo.PositiveIntegers,
        default=2,
        initialize=num_dispatch,
    )
    m.commitmentPeriodLength = pyo.Param(
        within=pyo.PositiveReals, default=1, units=u.hr
    )

    # [TODO: Index by dispatch period? Certainly index by
    # commitment period.]
    m.dispatchPeriodLength = pyo.Param(
        within=pyo.PositiveReals, initialize=duration_dispatch, units=u.minutes
    )

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

    m.lifetimes = pyo.Param(
        m.generators,
        initialize={
            gen: m.md.data["elements"]["generator"][gen]["lifetime"]
            for gen in m.generators
        },
        mutable=True,
        units=u.year,
        doc="Lifetime of each generator",
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
        mutable=True,
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
        m.transmission,
        initialize={
            branch: (m.md.data["elements"]["branch"][branch].get("distance") or 0)
            for branch in m.transmission
        },
        mutable=True,
        units=u.m,
        doc="Distance between terminal buses for each transmission line",
    )

    # Initialize investment costs in each new transmission
    # line. Currently selected the value of 0 to ensure investments
    # will be selected, if needed.
    m.branchInvestmentCost = pyo.Param(
        m.transmission,
        initialize={
            branch: (m.md.data["elements"]["branch"][branch].get("capital_cost") or 0)
            for branch in m.transmission
        },
        mutable=True,
        units=u.USD,
        doc="Investment cost for each new branch",
    )

    # [JSC TODO: Add branch capital multiplier to input data.]
    m.branchCapitalMultiplier = pyo.Param(
        m.transmission,
        initialize={
            branch: (
                m.md.data["elements"]["branch"][branch].get("capital_multiplier") or 1
            )
            for branch in m.transmission
        },
        mutable=True,
        units=u.dimensionless,
    )

    m.branchExtensionMultiplier = pyo.Param(
        m.transmission,
        initialize={
            branch: (
                m.md.data["elements"]["branch"][branch].get("extension_multiplier") or 1
            )
            for branch in m.transmission
        },
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life extension for each generator expressed as a fraction of initial investment cost",
    )

    # [TODO: These should go into each stage. Check where these values
    # should come from.]
    m.peakLoad = pyo.Param(m.stages, default=0, units=u.MW)
    m.reserveMargin = pyo.Param(m.stages, default=0, units=u.MW)
    m.renewableQuota = pyo.Param(m.stages, default=0, units=u.MW)
    m.weights = pyo.Param(m.representativePeriods, default=1)
    m.investmentFactor = pyo.Param(
        m.stages, default=1, mutable=True, units=u.dimensionless
    )
    m.deficitPenalty = pyo.Param(m.stages, default=1, units=u.USD / u.MW)

    # Initialize fuel cost. This is re-calculated during the
    # investment stage with values from preprocessing data.
    m.fuelCost = pyo.Param(
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
        units=u.USD / (u.MW * u.hr),
        doc="Cost per unit of fuel at each generator",
    )

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

    # Initialize generator investment costs to 0 and re-defined them
    # during investment stage using data from preprocessing.
    m.generatorInvestmentCost = pyo.Param(
        m.generators,
        initialize={gen: 1 for gen in m.generators},
        mutable=True,
        units=u.USD / u.MW,
        doc="Investment cost for all generators",
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
        # units=u.MW / u.minutes,
        units=u.dimensionless,
        doc="Ramp up rates for each generator as a fraction of maximum generator output",
    )

    m.rampDownRates = pyo.Param(
        m.thermalGenerators,
        initialize={
            thermalGen: m.md.data["elements"]["generator"][thermalGen]["ramp_down_rate"]
            for thermalGen in m.thermalGenerators
        },
        # units=u.MW / u.minutes,
        units=u.dimensionless,
        doc="Ramp down rates for each generator as a fraction of maximum generator output",
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

    # Initialize fixed and variable costs and update values during
    # investment stage.
    m.fixedCost = pyo.Param(
        m.generators,
        initialize={gen: 1 for gen in m.generators},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Fixed operating costs",
    )
    m.varCost = pyo.Param(
        m.generators,
        initialize={gen: 1 for gen in m.generators},
        mutable=True,
        units=u.USD / (u.MW * u.hr),
        doc="Variable costs",
    )

    # Initialize curtailment and load shed costs as parameters and
    # re-calculate them during the investment stage.
    m.curtailmentCost = pyo.Param(
        initialize=1,
        units=u.USD / (u.MW * u.hr),
        mutable=True,
        doc="Curtailment cost",
    )
    m.loadShedCostperCurtailment = pyo.Param(
        initialize=1000,
        units=u.USD / (u.MW * u.hr),
        mutable=True,
    )

    if m.config["storage"] == True:
        stor.add_storage_params(m)

    # Add legacy parameters.These parameters are commented in the
    # original model. Keep here to check if we should include them in
    # future version of the model.
    """
    # NOTE: Commented for now since it is not use in this case but
    # might be used in the future when considering ACOPF.
    # m.lossRate = pyo.Param(
    #     m.transmission,
    #     initialize={branch: (m.md.data["elements"]["branch"][branch].get("loss_rate") or 0)
    #                 for branch in m.transmission},
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


def add_model_cost_parameters(m, year):
    """This method saves lists with all relevant costs (fixed and
    variable operating costs, fuel costs, and investment costs) for
    thermal and renewable generators. Refer to gtep_data_processing
    script for more details about the preprocessing of this data.

    Note that the "capex" in the investment costs already include the
    interest rate for each generator. Additionally, this data only
    covers three years: 2025, 2030, and 2035. If more investment years
    are needed, more data should be included in the data file for data
    processing.

    """

    # [WIP: Assume we have two types of generators: thermal "CT" (with
    # gas fuel) and renewable "PV" (with "sun" as fuel).]
    gen_thermal_type = "CT"
    gen_renewable_type = "PV"

    m.genThermalInvCost = []
    m.genThermalFuelCost = []
    m.genThermalFixOpCost = []
    m.genThermalVarOpCost = []
    m.genRenewableInvCost = []
    m.genRenewableFuelCost = []
    m.genRenewableFixOpCost = []
    m.genRenewableVarOpCost = []

    if m.mc is not None:
        for index, row in m.mc.gen_data_target.iterrows():
            if row["Unit Type"].startswith(gen_thermal_type):
                m.genThermalInvCost.append(row[f"capex_{year}"])  # in $/kW
                m.genThermalFixOpCost.append(row[f"fixed_ops_{year}"])  # in $/kW-yr
                m.genThermalVarOpCost.append(row[f"var_ops_{year}"])  # $/MWh
                m.genThermalFuelCost.append(row[f"fuel_costs_{year}"])

            elif row["Unit Type"].startswith(gen_renewable_type):
                m.genRenewableInvCost.append(row[f"capex_{year}"])  # in $/kW
                m.genRenewableFixOpCost.append(row[f"fixed_ops_{year}"])  # in $/kW-yr
                m.genRenewableVarOpCost.append(row[f"var_ops_{year}"])  # $/MWh
                m.genRenewableFuelCost.append(row[f"fuel_costs_{year}"])

            else:
                continue
    else:
        # TODO: Check what the default costs should be
        print(
            "Cost data was not provided in m.mc instance (check DataProcessing for more details). Setting costs parameters to random values for now."
        )
        m.genThermalInvCost.append(1)  # in $/kW
        m.genThermalFixOpCost.append(1)  # in $/kW-yr
        m.genThermalVarOpCost.append(1)  # $/MWh
        m.genThermalFuelCost.append(1)
        m.genRenewableInvCost.append(1)  # in $/kW
        m.genRenewableFixOpCost.append(1)  # in $/kW-yr
        m.genRenewableVarOpCost.append(1)  # $/MWh
        m.genRenewableFuelCost.append(1)

    # Update data for fixed and variable costs (previously defined
    # with random default values in add_model_parameters) since they
    # depend on the investment year. Also, convert the units to be
    # consistent.
    units_fixed_cost = u.USD / (u.kW * u.year)
    units_var_cost = u.USD / (u.MW * u.hr)
    units_inv_cost = u.USD / u.kW
    units_fuel_cost = u.USD / (u.MW * u.hr)
    for gen in m.generators:
        if m.md.data["elements"]["generator"][gen]["generator_type"] == "thermal":
            m.fixedCost[gen] = pyo.units.convert(
                m.genThermalFixOpCost[0] * units_fixed_cost,
                to_units=u.USD / (u.MW * u.hr),
            )
            m.varCost[gen] = m.genThermalVarOpCost[0] * units_var_cost

            m.generatorInvestmentCost[gen] = pyo.units.convert(
                m.genThermalInvCost[0] * units_inv_cost, to_units=u.USD / u.MW
            )

            # Add fuel costs from preprocessed data. Consider this
            # cost is for Natural Gas generators, not coal.
            m.fuelCost[gen] = m.genThermalFuelCost[0] * units_fuel_cost

        else:

            # For renewable
            m.fixedCost[gen] = pyo.units.convert(
                m.genRenewableFixOpCost[0] * units_fixed_cost,
                to_units=u.USD / (u.MW * u.hr),
            )
            m.varCost[gen] = m.genRenewableVarOpCost[0] * units_var_cost

            m.generatorInvestmentCost[gen] = pyo.units.convert(
                m.genRenewableInvCost[0] * units_inv_cost, to_units=u.USD / u.MW
            )

    # Final (converted) units are:
    # fixed cost = $/MWh
    # var cost = $/MWh
    # inv cost = $/Mw
    # fuel cost = $/MWh

    # Cost per MW of curtailed renewable energy. This equation
    # re-calculates curtailment and load shed costa since they depend
    # on the re-defined parameter "fuelCost".
    m.curtailmentCost = 2 * max(
        pyo.value(m.fuelCost[gen]) for gen in m.thermalGenerators
    )
    m.loadShedCostperCurtailment = 1000 * m.curtailmentCost
