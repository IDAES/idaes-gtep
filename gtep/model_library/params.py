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

"""Parameters in the Generation and Transmission Expansion Planning
(GTEP) Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u

import gtep.model_library.storage as stor


def model_data_references(m):
    """Creates and labels all the parameters in the GTEP model. This
    method ties input data directly to the model.

    :param m: Pyomo model object

    """

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

    # [ESR WIP: Rename since the name was repeated in the
    # commitment_period_rule function. Check if this is correct.]
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

    # TODO: WHAT HAVE I DONE HERE I HATE IT and JSC made it worse...

    # [ESR WIP: Take only the value for renewable capacity when using
    # max() to avoid errors.]
    # BLN: Pretty sure this should be removed but double check commented constraint using this
    """ m.renewableCapacityValue = pyo.Param(
        m.renewableGenerators,
        initialize={
            renewableGen: (
                0
                if type(m.md.data["elements"]["generator"][renewableGen]["p_max"])
                == float
                else min(
                    m.md.data["elements"]["generator"][renewableGen]["p_max"]["values"]
                )
                / max(1, pyo.value(m.renewableCapacityNameplate[renewableGen]))
            )
            for renewableGen in m.renewableGenerators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Fraction of generation capacity that can be reliably counted toward planning reserve",
    )
 """
    # [ESR WIP: From case data, the value is divided by 100, which is
    # the per units conversion.]
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
        # mutable=True,
        units=u.dimensionless,
        doc="Maximum fraction of maximum thermal generation output that can be supplied as quickstart reserve",
    )

    # [ESR WIP: When creating a Param for loads, an error occurs since
    # the load at each bus is a dictionary. To avoid this, I
    # initialized a m.loads parameter with a value of 0 and scaled it
    # with the right value in commitment_period_rule. I also created a
    # new set for the buses that have loads.]
    m.loads = pyo.Param(
        m.buses,
        initialize={load_n: 0 for load_n in m.buses},
        mutable=True,
        units=u.MW,
        doc="Demand at each bus",
    )

    # [ESR WIP: Commented for now since it is not use in this case but
    # might be used in the future when considering ACOPF]
    # m.lossRate = pyo.Param(
    #     m.transmission,
    #     initialize={branch: (m.md.data["elements"]["branch"][branch].get("loss_rate") or 0)
    #                 for branch in m.transmission},
    #     mutable=True,
    #     # units=,
    #     doc="Per-distance-unit multiplicative loss rate for each transmission line"
    # )

    ## NOTE: lazy fixing for dc_branch and branch... but should be an ok lazy fix
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

    # TODO: Add cost of investment in each new branch to input data. Currently
    # selected 0 to ensure investments will be selected if needed
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

    # JSC TODO: Add branch capital multiplier to input data.
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

    ## TODO: These should go into each stage -- check where these
    ## values should come from
    m.peakLoad = pyo.Param(m.stages, default=0, units=u.MW)
    m.reserveMargin = pyo.Param(m.stages, default=0, units=u.MW)
    m.renewableQuota = pyo.Param(m.stages, default=0, units=u.MW)
    m.weights = pyo.Param(m.representativePeriods, default=1)
    m.investmentFactor = pyo.Param(
        m.stages, default=1, mutable=True, units=u.dimensionless
    )
    m.deficitPenalty = pyo.Param(m.stages, default=1, units=u.USD / u.MW)

    # (Original) NOTE: Lazy approx for NPV. [TODO: don't lazily approx
    # NPV, add it into unit handling and calculate from actual time
    # frames]

    # [ESR WIP: Commented since it is already included in the costs we
    # have from preprocessing stage.]
    # for stage in m.stages:
    #     m.investmentFactor[stage] *= 1 / ((1.04) ** (5 * stage))

    # # [ESR WIP: Commented for now but depends on the type of data we
    # # are using for generators.]
    # m.startFuel = pyo.Param(
    #     m.generators,
    #     initialize={gen: m.md.data["elements"]["generator"][gen]["start_fuel"]
    #                 for gen in m.generators},
    #     mutable=True,
    #     # units=
    #     doc="Amount of fuel required to be consumed for startup process for each generator"
    # )

    # [ESR WIP: Original fuel cost. This is re-defined in the function
    # investment_stage_rule with values from preprocessed data.]
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

    # [ESR WIP: Include start-up cost only in thermal generators assuming a natural gas plant.]
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

    # BLN: TODO: Check what value should be used here
    m.retirementMultiplier = pyo.Param(
        m.generators,
        initialize={
            gen: (0.1 if gen in m.thermalGenerators else 1.0) for gen in m.generators
        },
        mutable=True,
        units=u.dimensionless,
        doc="Cost of life retirement for each generator expressed as a fraction of initial investment cost",
    )

    # [ESR WIP: Replace original generator investment costs with costs
    # from preprocessed data. These are fixed to 0 here but re-defined
    # in the function investment_stage_rule.]
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

    
    # [ESR WIP: Declare fixed and operating costs here to avoid
    # multiple declarations of the same parameter. Set the value to 1
    # for now and updated in function investment_stage_rule.]
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

    # [ESR WIP: Declare and initialize curtailment and load shed costs
    # as parameters. These are re-calculated in
    # investment_stage_rule. Also, note that the original
    # "loadShedCost" was renamed "loadShedCostperCurtailment" to avoid
    # repetition. ]
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
