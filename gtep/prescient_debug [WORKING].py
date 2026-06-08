from egret.data.model_data import ModelData
from egret.parsers import rts_gmlc_parser as parser
import pandas as pd
import os


def _read_branches(base_dir: str, elements: dict, bus_id_to_name: dict) -> None:

    # add the branches
    elements["branch"] = {}
    branch_df = pd.read_csv(os.path.join(base_dir, "branch.csv"))
    for idx, row in branch_df.iterrows():

        branch_dict = {
            "from_bus": bus_id_to_name[str(row["From Bus"])],
            "to_bus": bus_id_to_name[str(row["To Bus"])],
            "in_service": True,
            "resistance": float(row["R"]),
            "reactance": float(row["X"]),
            "charging_susceptance": float(row["B"]),
            "rating_long_term": float(row["Cont Rating"]) or None,
            "rating_short_term": float(row["LTE Rating"]) or None,
            "rating_emergency": float(row["STE Rating"]) or None,
            "angle_diff_min": -90,
            "angle_diff_max": 90,
            "pf": None,
            "qf": None,
            "pt": None,
            "qt": None,
        }

        TAP = float(row["Tr Ratio"])
        if TAP != 0.0:
            branch_dict["branch_type"] = "transformer"
            branch_dict["transformer_tap_ratio"] = TAP
            branch_dict["transformer_phase_shift"] = 0.0
        else:
            branch_dict["branch_type"] = "line"

        name = str(row["UID"])
        elements["branch"][name] = branch_dict
    branch_df = None


def _read_generators(base_dir: str, elements: dict, bus_id_to_name: dict) -> None:
    from math import isnan

    # add the generators
    elements["generator"] = {}
    RENEWABLE_TYPES = {"WIND", "HYDRO", "RTPV", "PV", "ROR"}
    gen_df = pd.read_csv(os.path.join(base_dir, "gen.csv"))
    for idx, row in gen_df.iterrows():
        # if this is storage we need to handle it differently
        if row["Fuel"] == "Storage":
            continue

        # NOTE: for now, Egret doesn't handle CSP -- not clear how to model
        if row["Unit Type"] == "CSP":
            continue

        name = str(row["GEN UID"])
        bus_name = bus_id_to_name[str(row["Bus ID"])]
        gen_dict = {
            "bus": bus_name,
            "in_service": True,
            "mbase": 100.0,
            "pg": float(row["MW Inj"]),
            "qg": float(row["MVAR Inj"]),
            "p_min": float(row["PMin MW"]),
            "p_max": float(row["PMax MW"]),
            "q_min": float(row["QMin MVAR"]),
            "q_max": float(row["QMax MVAR"]),
            "ramp_q": float(row["Ramp Rate MW/Min"]),
            "fuel": str(row["Fuel"]),
            "unit_type": str(row["Unit Type"]),
            "area": elements["bus"][bus_name]["area"],
            "zone": elements["bus"][bus_name]["zone"],
        }

        # Remove optional values if not present
        for key in ("p_min", "p_max", "q_min", "q_max", "ramp_q"):
            if isnan(gen_dict[key]):
                del gen_dict[key]

        UNIT_TYPE = str(row["Unit Type"])
        if UNIT_TYPE in RENEWABLE_TYPES:
            gen_dict["generator_type"] = "renewable"
            # ROR is treated as HYDRO by Egret
            if UNIT_TYPE == "ROR":
                gen_dict["unit_type"] = "HYDRO"
        elif UNIT_TYPE == "SYNC_COND":
            ## TODO: should we have a flag for these?
            gen_dict["generator_type"] = "thermal"
        else:
            gen_dict["generator_type"] = "thermal"

        elements["generator"][name] = gen_dict

        # after this is only really needed for thermal units
        if UNIT_TYPE in RENEWABLE_TYPES:
            continue

        # Gen cost
        ## round as in RTS-GMLC Prescient/topysp.py
        pmax = float(row["PMax MW"])

        # There can be any number of 'Output_pct_<i>' columns.
        # Stop at the first one that doesn't exist or doesn't hold a number
        def valid_output_pcts():
            for i in range(50):
                try:
                    val = float(row[f"Output_pct_{i}"])
                    if isnan(val):
                        return
                    yield (i, val)
                except:
                    return

        x = {i: round(val * pmax, 1) for i, val in valid_output_pcts()}
        fuel_field_count = len(x)

        if fuel_field_count > 0:
            ## /1000. from the RTS-GMLC MATPOWER writer --
            ## heat rates are in BTU/kWh, 10^6 BTU == 1 MMBTU, 10^3 kWh == 1 MWh, so MMBTU/MWh == 10^3/10^6 * BTU/kWh
            f = {}
            f[0] = (float(row["HR_avg_0"]) * 1000.0 / 1000000.0) * x[0]
            for i in range(1, fuel_field_count):
                f[i] = (
                    (
                        (x[i] - x[i - 1])
                        * (float(row[f"HR_incr_{i}"]) * 1000.0 / 1000000.0)
                    )
                ) + f[i - 1]

            F_COEFF = [
                (x[i], round(f[i], 2))
                for i in range(fuel_field_count)
                if (
                    ((i == 0) or (x[i - 1], f[i - 1]) != (x[i], f[i]))
                    and (x[i], f[i]) != (0.0, 0.0)
                )
            ]
            if F_COEFF == []:
                F_COEFF = [(pmax, 0.0)]
            gen_dict["p_fuel"] = {"data_type": "fuel_curve", "values": F_COEFF}

        # UC Data
        MIN_DN_TIME = float(row["Min Down Time Hr"])

        # Startup types and costs, from hot to cold
        startup_heat = (
            float(row["Start Heat Hot MBTU"]),
            float(row["Start Heat Warm MBTU"]),
            float(row["Start Heat Cold MBTU"]),
        )
        startup_time = (
            float(row["Start Time Hot Hr"]),
            float(row["Start Time Warm Hr"]),
            float(row["Start Time Cold Hr"]),
        )

        # Arrange fuel requirements from hottest to coldest, ignoring missing values.
        startup_fuel = []
        for i in range(3):
            # Skip blank values
            if isnan(startup_time[i]) or isnan(startup_heat[i]):
                continue

            t = max(startup_time[i], MIN_DN_TIME)
            f = startup_heat[i]

            # For entries with matching times, use to the colder data
            if len(startup_fuel) > 0 and startup_fuel[-1][0] == t:
                startup_fuel[-1] = (t, f)
            else:
                startup_fuel.append((t, f))

        # If the warmest fuel requirement has a time longer than the minimum
        # down time, extend that warmest requirement down to minimum down time.
        if len(startup_fuel) > 0 and startup_fuel[0][0] > MIN_DN_TIME:
            startup_fuel[0] = (MIN_DN_TIME, startup_fuel[0][1])

        gen_dict["startup_fuel"] = startup_fuel
        fixed_startup_cost = float(row["Non Fuel Start Cost $"])
        if not isnan(fixed_startup_cost):
            gen_dict["non_fuel_startup_cost"] = fixed_startup_cost
        gen_dict["shutdown_cost"] = 0.0

        gen_dict["agc_capable"] = True
        gen_dict["p_min_agc"] = gen_dict["p_min"]
        gen_dict["p_max_agc"] = gen_dict["p_max"]

        ramp_q = gen_dict["ramp_q"]
        gen_dict["ramp_agc"] = ramp_q
        gen_dict["ramp_up_60min"] = 60.0 * ramp_q
        gen_dict["ramp_down_60min"] = 60.0 * ramp_q

        gen_dict["fuel_cost"] = float(row["Fuel Price $/MMBTU"])

        # these assumptions are the same as prescient-rtsgmlc
        gen_dict["startup_capacity"] = gen_dict["p_min"]
        gen_dict["shutdown_capacity"] = gen_dict["p_min"]
        gen_dict["min_up_time"] = float(row["Min Up Time Hr"])
        gen_dict["min_down_time"] = MIN_DN_TIME

        elements["generator"][name] = gen_dict
    gen_df = None


data_path = "./gtep/data/123_Bus_Resil_Week"
# data_path = "./gtep/data/123_Bus_Coal"

### from egret _create_rtsgmlc_skeleton ###
model_data = ModelData.empty_model_data_dict()

elements = model_data["elements"]
system = model_data["system"]

system["name"] = "RTS-GMLC"
system["baseMVA"] = 100.0

bus_id_to_name = parser._read_buses_and_areas(
    data_path, elements, system
)  # no error either 123
parser._read_branches(
    data_path, elements, bus_id_to_name
)  # no error 123_resil #error 123_coal
_read_generators(
    data_path, elements, bus_id_to_name
)  # error 123_resil (ramp_q key error)

pass
