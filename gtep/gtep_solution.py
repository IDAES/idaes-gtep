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

# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield, Thom R. Edwards
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

import os
import logging
import json
import pandas as pd
import numpy as np
from pathlib import Path
from collections import namedtuple, defaultdict

import pyomo.environ as pyo
import pyomo.gdp as gdp
from pyomo.environ import units as u
from pyomo.core.base.param import IndexedParam
from pyomo.core.base.expression import ScalarExpression, IndexedExpression

import squarify

import plotly.graph_objects as go

logger = logging.getLogger(__name__)


# [TODO] inject units into plots
class ExpansionPlanningSolution:
    """A class that stores the solution to the ExpansionPlanningModel
    class for writing and visualization."""

    def __init__(self, data_path):
        self.gen_df = pd.read_csv(f"{data_path}/gen.csv")
        self.gen_types = {
            gen_type: self.gen_df[self.gen_df["Unit Type"] == gen_type]["PMax MW"].sum()
            for gen_type in set(self.gen_df["Unit Type"])
        }

    def load_from_file(self):
        pass

    def save_results_in_json_files(self, gtep_model, dir_name):

        folder_name = dir_name

        os.makedirs(folder_name, exist_ok=True)
        print(
            f"\n Creating the directory '{folder_name}' to save the results. Working on it ..."
        )
        m = gtep_model.model

        valid_names = ["Inst", "Oper", "Disa", "Ext", "Ret"]
        renewable_investments = {}
        dispatchable_investments = {}
        load_shed = {}
        power_flow = {}
        generation = {}
        curtailment = {}
        reserves = {}
        charging = {}
        discharging = {}
        for var in m.component_objects(pyo.Var, descend_into=True):
            for index in var:
                if "Shed" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        load_shed[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Reserve" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        reserves[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Flow" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        power_flow[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Generation" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        generation[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Curtailment" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        curtailment[var.name + "." + str(index)] = pyo.value(var[index])
                elif "storageCharged" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        charging[var.name + "." + str(index)] = pyo.value(var[index])
                elif "storageDischarge" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        discharging[var.name + "." + str(index)] = pyo.value(var[index])
                for name in valid_names:
                    if name in var.name:
                        if pyo.value(var[index]) >= 0.001:
                            renewable_investments[var.name + "." + str(index)] = (
                                pyo.value(var[index])
                            )
        for var in m.component_objects(gdp.Disjunct, descend_into=True):
            for index in var:
                for name in valid_names:
                    if name in var.name:
                        if pyo.value(var[index].indicator_var) == True:
                            dispatchable_investments[var.name + "." + str(index)] = (
                                pyo.value(var[index].indicator_var)
                            )

        costs = {}
        for exp in m.component_objects(pyo.Expression, descend_into=True):
            if "Cost" in exp.name or "cost" in exp.name:
                if type(exp) is ScalarExpression:
                    costs[exp.name] = pyo.value(exp)
            if type(exp) is IndexedExpression:
                for e in exp:
                    costs[exp[e].name] = pyo.value(exp[e])

        loads = {}
        for param in m.component_objects(pyo.Param, descend_into=True):
            if "commitment" in param.name and "loads" in param.name:
                if type(param) is IndexedParam:
                    for p in param:
                        loads[param[p].name] = pyo.value(param[p])

        # Output file names
        output_files = {
            "renewable_investments": renewable_investments,
            "dispatchable_investments": dispatchable_investments,
            "load_shed": load_shed,
            "costs": costs,
            "flows": power_flow,
            "generation": generation,
            "curtailment": curtailment,
            "loads": loads,
            "reserves": reserves,
            "charging": charging,
            "discharging": discharging,
        }

        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        for name, data in output_files.items():
            filename = f"{folder_name}/{name}.json"
            with open(filename, "w") as fil:
                json.dump(data, fil)

        print(" -> The following files have been created in '{}':".format(folder_name))
        for name in output_files:
            print(f" - {folder_name}/{name}.json")

    def read_json(self, filepath):
        # Read a json file
        json_filepath = Path(filepath)
        with open(json_filepath, "r") as fobj:
            json_read = json.loads(fobj.read())

        return json_read

    def to_dict(self, dict_in):
        """Converts a flat dictionary with dot-separated keys into a
        nested dictionary: {time_key: {state: {gen_key: value}}}

        Ignores entries where the second part of the key is 'branch'.

        """

        ignore_this = "branch"
        out_dict = {}

        for key, val in dict_in.items():
            # split the name to figure out depth
            split_name = key.split(".")

            if ignore_this not in split_name[1]:
                # set toplevel defaults
                out_dict.setdefault(split_name[0], {})

                # split things by a predefined prefix
                out_dict[split_name[0]].setdefault(split_name[1], {})

                # specific_key = split_name[1].split(subsplit_key, 1)[1]
                specific_key = "".join(split_name[2:])
                out_dict[split_name[0]][split_name[1]][specific_key] = val

        return out_dict

    def create_plots(self, case_json, results_path, data_path, plot_type="all"):
        """This function reads a solution .json file, uses gen.csv and
        candidate_generators_initial_list.csv to map generator UIDs to
        unit types and PMax MW, and generates a stacked bar plot of
        generation mix by investment year.

        """

        plots_dir = os.path.join(results_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        print(f" Created the subdirectory '{plots_dir}' to save the plots.")

        GenerationType = namedtuple("GenerationType", ["label", "color"])
        GENERATION_TYPES = {
            "CC": GenerationType("Gas CC", "#20b2aa"),
            "CT": GenerationType("Gas CT", "#6e8b3d"),
            "PV": GenerationType("Solar", "#ffb90f"),
            "NUC": GenerationType("Nuclear", "#39FF14"),
            "STEAM": GenerationType("Steam", "#b0b0b0"),
            "THERMAL": GenerationType("Thermal", "#e25822"),
            "COAL": GenerationType("Coal", "#333333"),
            "WIND": GenerationType("Wind", "#4f94cd"),
            "DR": GenerationType("Demand Response", "#a020f0"),
            "HYDRO": GenerationType("Hydro", "#00bfff"),
            "BATTERY": GenerationType("Battery", "#25ccff"),
        }

        def get_gen_arrays(gen_case_json, results_path, data_path, GENERATION_TYPES):

            # Read gen and candidate_gen .csv files for GEN UID to map for
            # Unit Type and PMax MW
            gen_df = pd.read_csv(f"{data_path}/gen.csv")
            gen_uid_to_type = {
                row["GEN UID"]: row["Unit Type"].upper() for _, row in gen_df.iterrows()
            }
            gen_uid_to_pmax = {
                row["GEN UID"]: float(row["PMax MW"]) for _, row in gen_df.iterrows()
            }
            gen_cand_df = pd.read_csv(
                f"{data_path}/candidate_generators_initial_list.csv"
            )
            gen_cand_uid_to_type = {
                row["GEN UID"]: row["Unit Type"].upper()
                for _, row in gen_cand_df.iterrows()
            }
            gen_cand_uid_to_pmax = {
                row["GEN UID"]: float(row["PMax MW"])
                for _, row in gen_cand_df.iterrows()
            }

            # Read storage.csv for storage units
            storage_df = pd.read_csv(f"{data_path}/storage.csv")
            storage_uid_to_type = {
                row["name"]: row["storage_type"].upper()
                for _, row in storage_df.iterrows()
            }
            storage_uid_to_pmax = {
                row["name"]: float(row.get("energy_capacity", 0))
                for _, row in storage_df.iterrows()
            }

            # Read and process .json. These names are based on the saved
            # .json files from function save_results_in_json_files
            if gen_case_json == "renewables":
                json_file = f"{results_path}/renewable_investments.json"
            elif gen_case_json == "dispatchables":
                json_file = f"{results_path}/dispatchable_investments.json"
            else:
                print("WARNING: Case not debugged")

            dict_in = self.to_dict(self.read_json(json_file))
            time_keys = list(dict_in.keys())

            # Collect all generator keys
            states_set = set()
            keys_set = set()
            for this_time_key in time_keys:
                states_set.update(dict_in[this_time_key].keys())
                for this_state in dict_in[this_time_key].keys():
                    keys_set.update(dict_in[this_time_key][this_state].keys())

            # Map generator keys to canonical unit types and PMax MW
            gens_keys_to_type = {}
            gens_keys_to_pmax = {}
            for this_key in list(keys_set):
                if this_key in gen_cand_uid_to_type:
                    unit_type = gen_cand_uid_to_type[this_key]
                    pmax = gen_cand_uid_to_pmax[this_key]
                elif this_key in gen_uid_to_type:
                    unit_type = gen_uid_to_type.get(this_key)
                    pmax = gen_uid_to_pmax.get(this_key)
                elif this_key in storage_uid_to_type:
                    unit_type = storage_uid_to_type[this_key]
                    pmax = storage_uid_to_pmax[this_key]
                else:
                    unit_type = None
                    pmax = 0

                # Make unit_type uppercase to ensure case-insensitive
                # matching
                unit_type_upper = unit_type.upper() if unit_type else None

                # Only use if in GENERATION_TYPES
                if unit_type_upper and unit_type_upper in GENERATION_TYPES:
                    gens_keys_to_type[this_key] = unit_type_upper
                    gens_keys_to_pmax[this_key] = pmax
                else:
                    raise ValueError(
                        f"[ERROR] Generator or storage '{this_key}' has unknown or unsupported unit type '{unit_type}'."
                    )

            # After building gens_keys_to_type
            unique_types = set(gens_keys_to_type.values())
            gen_types_sorted = sorted(unique_types)  # Alphabetical order

            # Read the DAY_AHEAD .csv file with the year column and get
            # unique years in order of appearance
            time_periods_df = pd.read_csv(f"{data_path}/DAY_AHEAD_renewables.csv")
            time_periods = (
                time_periods_df["Year"].drop_duplicates().astype(str).tolist()
            )

            # Build gen_mix using PMax MW and solution values
            gen_mix = {tp: {gt: 0.0 for gt in gen_types_sorted} for tp in time_periods}
            for tp in time_periods:
                for k, val in dict_in.items():
                    for state, gen_dict in val.items():
                        for gen, value in gen_dict.items():
                            unit_type = gens_keys_to_type.get(gen)
                            if gen_case_json == "dispatchables":
                                pmax = gens_keys_to_pmax.get(gen, 0.0)
                                if unit_type in gen_mix[tp]:
                                    gen_mix[tp][unit_type] += pmax * value
                            elif gen_case_json == "renewables":
                                if unit_type in gen_mix[tp]:
                                    gen_mix[tp][unit_type] += value

            gen_mix_arrays = {
                k: np.array([gen_mix[stage].get(k, 0.0) for stage in gen_mix.keys()])
                for k in gen_types_sorted
            }
            # print('gen_mix_arrays:', gen_mix_arrays)

            return gen_mix, gen_mix_arrays, time_periods

        # Define multiple functions to create interactive Plotly plots
        # for the generation mix: a stack plot, a pie chart, and a
        # treemap. The user can select which one to use by setting up
        # the plot_type option. By default, this function will plot
        # all if no value is given.
        def plotly_stackplot_gen_mix(
            time_periods, gen_mix_arrays, GENERATION_TYPES, results_path, case_json
        ):
            """This function creates an interactive Plotly stacked bar
            chart of generation mix by investment year and saves it as
            an HTML file.

            """

            # Prepare the bottom (cumulative sum) for stacking
            fig = go.Figure()
            bottom = np.zeros(len(time_periods))

            for gen_class, mix_array in gen_mix_arrays.items():
                if gen_class in GENERATION_TYPES:
                    component_label = GENERATION_TYPES[gen_class].label
                    component_color = GENERATION_TYPES[gen_class].color
                    fig.add_bar(
                        x=time_periods,
                        y=mix_array,
                        name=component_label,
                        marker_color=component_color,
                    )

            fig.update_layout(
                barmode="relative",
                title="Generation Mix",
                xaxis_title="Investment Year",
                yaxis_title="Nameplate Capacity [MW]",
                legend=dict(
                    yanchor="middle",
                    y=0.5,
                    xanchor="left",
                    x=1.02,
                    font=dict(size=14),
                ),
                width=1200,
                height=600,
                plot_bgcolor="white",
            )
            fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="lightgray")
            fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgray")

            # Save as an interactive HTML
            plot_path = (
                f"{results_path}/plots/{case_json}_gen_mix_summary_interactive.html"
            )
            fig.write_html(f"{plot_path}")
            print(f" -> Saved interactive stack plot for generation mix to {plot_path}")

        def plotly_treemap_gen_mix(
            gen_mix, GENERATION_TYPES, results_path, case_json, small_pct_threshold=5
        ):
            """This function creates and saves an interactive Plotly
            treemap of generation mix for each time period as an HTML
            file.

            """
            for tp, mix in gen_mix.items():
                filtered_mix = {
                    k: v for k, v in mix.items() if v > 0 and k in GENERATION_TYPES
                }
                if not filtered_mix:
                    continue

                total = sum(filtered_mix.values())
                sorted_items = sorted(filtered_mix.items(), key=lambda x: -x[1])
                sizes = [v for k, v in sorted_items]
                pcts = [v / total * 100 for k, v in sorted_items]
                labels = []
                colors = []
                customdata = []

                for (k, v), pct in zip(sorted_items, pcts):
                    label = f"{GENERATION_TYPES[k].label}<br>{int(v)} MW<br>{pct:.1f}%"
                    labels.append(label if pct >= small_pct_threshold else "")
                    colors.append(GENERATION_TYPES[k].color)
                    customdata.append(
                        f"{GENERATION_TYPES[k].label} ({int(v)} MW, {pct:.1f}%)"
                    )

                # Plotly treemap
                fig = go.Figure(
                    go.Treemap(
                        labels=[GENERATION_TYPES[k].label for k, v in sorted_items],
                        parents=[""] * len(sorted_items),
                        values=sizes,
                        marker=dict(colors=colors),
                        textinfo="label+value+percent entry",
                        hovertext=customdata,
                        hoverinfo="text",
                    )
                )

                fig.update_layout(
                    title=f"Generation Mix Treemap - {tp}",
                    width=900,
                    height=600,
                    margin=dict(t=50, l=25, r=25, b=25),
                )

                # Save as an interactive HTML
                plot_path = (
                    f"{results_path}/plots/{case_json}_treemap_{tp}_interactive.html"
                )
                fig.write_html(f"{plot_path}")
                print(f" -> Saved interactive treemap for {tp} to {plot_path}")

        def plotly_pie_gen_mix(gen_mix, GENERATION_TYPES, results_path, case_json):
            """This method creates and saves an interactive Plotly pie
            chart of generation mix for each time period as an HTML
            file.

            """
            for tp, mix in gen_mix.items():
                filtered_mix = {
                    k: v for k, v in mix.items() if v > 0 and k in GENERATION_TYPES
                }
                if not filtered_mix:
                    continue

                sizes = [filtered_mix[k] for k in filtered_mix]
                total = sum(sizes)
                labels = [
                    f"{GENERATION_TYPES[k].label} ({int(filtered_mix[k])} MW, {filtered_mix[k]/total*100:.1f}%)"
                    for k in filtered_mix
                ]
                colors = [GENERATION_TYPES[k].color for k in filtered_mix]

                # Plotly pie chart
                fig = go.Figure(
                    go.Pie(
                        labels=[GENERATION_TYPES[k].label for k in filtered_mix],
                        values=sizes,
                        marker=dict(colors=colors, line=dict(color="white", width=1)),
                        textinfo="label+percent",
                        hoverinfo="label+value+percent",
                        pull=[0.05]
                        * len(sizes),  # Slightly "explode" all slices for separation
                        hole=0,  # 0 for pie, >0 for donut
                    )
                )

                fig.update_layout(
                    title=f"Generation Mix Pie Chart - {tp}",
                    width=700,
                    height=700,
                    margin=dict(t=50, l=25, r=25, b=25),
                    showlegend=True,
                )

                # Save as an interactive HTML
                plot_path = (
                    f"{results_path}/plots/{case_json}_pie_leader_{tp}_interactive.html"
                )
                fig.write_html(f"{plot_path}")
                print(f" -> Saved interactive pie chart for {tp} to {plot_path}")

        # Create gen_mix dictionary and arrays needed for the plots
        if case_json == "combined":
            gen_mix_ren, gen_mix_arrays_ren, time_periods_ren = get_gen_arrays(
                "renewables", results_path, data_path, GENERATION_TYPES
            )
            gen_mix_disp, gen_mix_arrays_disp, time_periods_disp = get_gen_arrays(
                "dispatchables", results_path, data_path, GENERATION_TYPES
            )

            # Check that time_periods are the same
            if time_periods_ren != time_periods_disp:
                raise ValueError(
                    "Time periods for renewables and dispatchables do not match!"
                )
            time_periods = time_periods_ren  # or use time_periods_disp too

            # Get the union of all time periods and all types
            all_time_periods = sorted(
                set(gen_mix_ren.keys()) | set(gen_mix_disp.keys())
            )
            all_types = sorted(
                set(
                    t
                    for mix in [gen_mix_ren, gen_mix_disp]
                    for v in mix.values()
                    for t in v
                )
            )

            # Merge gen_mix
            gen_mix = {}
            for tp in all_time_periods:
                gen_mix[tp] = {}
                for gt in all_types:
                    val_ren = gen_mix_ren.get(tp, {}).get(gt, 0.0)
                    val_disp = gen_mix_disp.get(tp, {}).get(gt, 0.0)
                    gen_mix[tp][gt] = val_ren + val_disp

            # Merge gen_mix_arrays
            gen_mix_arrays = {
                gt: np.array([gen_mix[tp].get(gt, 0.0) for tp in all_time_periods])
                for gt in all_types
            }

        else:
            gen_mix, gen_mix_arrays, time_periods = get_gen_arrays(
                case_json, results_path, data_path, GENERATION_TYPES
            )

        if plot_type == "stackplot":
            plotly_stackplot_gen_mix(
                time_periods, gen_mix_arrays, GENERATION_TYPES, results_path, case_json
            )
        elif plot_type == "treemap":
            plotly_treemap_gen_mix(gen_mix, GENERATION_TYPES, results_path, case_json)
        elif plot_type == "piechart":
            plotly_pie_gen_mix(gen_mix, GENERATION_TYPES, results_path, case_json)
        elif plot_type == "all":
            plotly_stackplot_gen_mix(
                time_periods, gen_mix_arrays, GENERATION_TYPES, results_path, case_json
            )
            plotly_treemap_gen_mix(gen_mix, GENERATION_TYPES, results_path, case_json)
            plotly_pie_gen_mix(gen_mix, GENERATION_TYPES, results_path, case_json)
        else:
            raise ValueError(
                f"Plot type '{plot_type}' is not supported. Please choose between 'stackplot', 'treemap', or 'piechart'."
            )

    def create_stackgraph_and_metrics(self, results_path, rep_days):

        try:
            import ujson as json
        except ImportError:
            import json

        with open(f"{results_path}/generation.json", "r") as F:
            gen_data = json.load(F)

        with open(f"{results_path}/loads.json", "r") as f:
            loads_data = json.load(f)

        with open(f"{results_path}/load_shed.json", "r") as f:
            load_shed_data = json.load(f)

        with open(f"{results_path}/reserves.json", "r") as f:
            reserves_data = json.load(f)

        with open(f"{results_path}/charging.json", "r") as f:
            charging_data = json.load(f)

        with open(f"{results_path}/discharging.json", "r") as f:
            discharging_data = json.load(f)

        # Note that these are the same colors used in the stack plots
        # and pie charts above
        def darken_color(hex_color, percent=0.2):
            """Darken a hex color by a given percent (0.2 = 20%)"""
            hex_color = hex_color.lstrip("#")
            rgb = [int(hex_color[i : i + 2], 16) for i in (0, 2, 4)]
            darker_rgb = [max(0, int(c * (1 - percent))) for c in rgb]
            return "#" + "".join(f"{c:02x}" for c in darker_rgb)

        GEN_TYPES = {
            "nuclear": "#39FF14",
            "coal": "#333333",
            "hydro": "#00bfff",
            "cc_gas": "#20b2aa",
            "ct_gas": "#6e8b3d",
            "battery_discharge": "#25ccff",
            "wind": "#4f94cd",
            "solar": "#ffb90f",
            "thermal_other": "#e25822",
            "steam": "#b0b0b0",
            "dr": "#a020f0",
            "ES4": "#a0522d",
            "battery_charge": "#25ccff",
            # Candidates: 20% darker than original, same pattern for all
            "hydro-c": darken_color("#00bfff"),
            "gas_cc-c": darken_color("#20b2aa"),
            "gas_ct-c": darken_color("#6e8b3d"),
            "battery-c": darken_color("#7b9095"),
            "wind-c": darken_color("#4f94cd"),
            "pv-c": darken_color("#ffb90f"),
            "steam-c": darken_color("#b0b0b0"),
            "ES4-c": darken_color("#a0522d"),
        }
        GEN_TYPE_HATCHES = {
            # No hatch pattern for "original" types
            "coal": "",
            "cc_gas": "",
            "ct_gas": "",
            "solar": "",
            "wind": "",
            "thermal_other": "",
            "dr": "",
            "hydro": "",
            "nuclear": "",
            "ES4": "",
            "battery_discharge": "",
            "battery_charge": "",
            # Candidates get a hatch pattern
            "hydro-c": "////",
            "gas_cc-c": "////",
            "gas_ct-c": "////",
            "battery-c": "////",
            "wind-c": "////",
            "pv-c": "////",
            "steam-c": "////",
            "ES4-c": "////",
        }
        GEN_TYPE_ALIASES = {
            "coal": "Coal",
            "cc_gas": "CC",
            "ct_gas": "CT",
            "dr": "DR",
            "solar": "Solar",
            "thermal_other": "Thermal",
            "wind": "Wind",
            "nuclear": "Nuclear",
            "hydro": "Hydro",
            "ES4": "ES4",
            "steam": "Steam",
            "gas_cc-c": "CC (Candidate)",
            "gas_ct-c": "CT (Candidate)",
            "pv-c": "Solar (Candidate)",
            "wind-c": "Wind (Candidate)",
            "hydro-c": "Hydro (Candidate)",
            "ES4-c": "ES4 (Candidate)",
            "steam-c": "Steam (Candidate)",
            "battery_charge": "Battery Charging",
            "battery_discharge": "Battery Discharging",
        }

        generation = {}
        for g, val in gen_data.items():
            c = list(pyo.ComponentUID(g)._cids)
            _, (stage,) = c.pop(0)
            if stage not in generation:
                generation[stage] = {}
            stage_dict = generation[stage]

            _, (period,) = c.pop(0)
            if period not in stage_dict:
                stage_dict[period] = {}
            period_dict = stage_dict[period]

            _, (commitment,) = c.pop(0)
            if commitment not in period_dict:
                period_dict[commitment] = {}
            commitment_dict = period_dict[commitment]

            _, (dispatch,) = c.pop(0)
            if dispatch not in commitment_dict:
                commitment_dict[dispatch] = dict.fromkeys(GEN_TYPES, 0)
            dispatch_dict = commitment_dict[dispatch]

            gen_name = c[-1][0]
            _type = None
            for gt in GEN_TYPES:
                if gen_name.endswith(gt):
                    _type = gt
                    break
            if _type is None:
                raise RuntimeError(f"Cannot map generator name '{gen_name}' to type")
            dispatch_dict[_type] += val

        # Add battery charging data to generation structure
        for g, val in charging_data.items():
            c = list(pyo.ComponentUID(g)._cids)
            _, (stage,) = c.pop(0)
            if stage not in generation:
                generation[stage] = {}
            stage_dict = generation[stage]

            _, (period,) = c.pop(0)
            if period not in stage_dict:
                stage_dict[period] = {}
            period_dict = stage_dict[period]

            _, (commitment,) = c.pop(0)
            if commitment not in period_dict:
                period_dict[commitment] = {}
            commitment_dict = period_dict[commitment]

            _, (dispatch,) = c.pop(0)
            if dispatch not in commitment_dict:
                commitment_dict[dispatch] = dict.fromkeys(GEN_TYPES, 0)
            dispatch_dict = commitment_dict[dispatch]

            dispatch_dict["battery_charge"] -= val

        # Add battery discharging data to generation structure
        # Per request, plot discharge as negative (below x-axis)
        for g, val in discharging_data.items():
            c = list(pyo.ComponentUID(g)._cids)
            _, (stage,) = c.pop(0)
            if stage not in generation:
                generation[stage] = {}
            stage_dict = generation[stage]

            _, (period,) = c.pop(0)
            if period not in stage_dict:
                stage_dict[period] = {}
            period_dict = stage_dict[period]

            _, (commitment,) = c.pop(0)
            if commitment not in period_dict:
                period_dict[commitment] = {}
            commitment_dict = period_dict[commitment]

            _, (dispatch,) = c.pop(0)
            if dispatch not in commitment_dict:
                commitment_dict[dispatch] = dict.fromkeys(GEN_TYPES, 0)
            dispatch_dict = commitment_dict[dispatch]

            dispatch_dict["battery_discharge"] += val

        print("\n[DEBUG] Storage summary from JSON inputs")

        total_charging = sum(charging_data.values())
        total_discharging = sum(discharging_data.values())

        print(f"[DEBUG] Total charging (raw): {total_charging:,.3f}")
        print(f"[DEBUG] Total discharging (raw): {total_discharging:,.3f}")

        charging_by_suffix = defaultdict(float)
        for g, val in charging_data.items():
            name = g.split(".")[-1]
            if name.endswith("_battery"):
                charging_by_suffix["battery"] += val
            elif name.endswith("_ps"):
                charging_by_suffix["ps"] += val
            else:
                charging_by_suffix["other"] += val

        discharging_by_suffix = defaultdict(float)
        for g, val in discharging_data.items():
            name = g.split(".")[-1]
            if name.endswith("_battery"):
                discharging_by_suffix["battery"] += val
            elif name.endswith("_ps"):
                discharging_by_suffix["ps"] += val
            else:
                discharging_by_suffix["other"] += val

        # print("[DEBUG] Charging by storage type suffix:")
        # for k, v in charging_by_suffix.items():
        #     print(f"    {k}: {v:,.3f}")

        # print("[DEBUG] Discharging by storage type suffix:")
        # for k, v in discharging_by_suffix.items():
        #     print(f"    {k}: {v:,.3f}")

        time_periods = [
            (s, p, c, d)
            for s in generation
            for p in generation[s]
            for c in generation[s][p]
            for d in generation[s][p][c]
        ]
        times = list(range(len(time_periods)))

        loads = {}
        for g, val in loads_data.items():
            c = list(pyo.ComponentUID(g)._cids)
            _, (stage,) = c.pop(0)
            if stage not in loads:
                loads[stage] = {}
            stage_dict = loads[stage]
            _, (period,) = c.pop(0)
            if period not in stage_dict:
                stage_dict[period] = {}
            period_dict = stage_dict[period]
            _, (commitment,) = c.pop(0)
            if commitment not in period_dict:
                period_dict[commitment] = 0
            period_dict[commitment] += val  # Sum all buses for this time period

        loads_trace = []
        for s, p, c, d in time_periods:
            try:
                total_load = loads[s][p][c]
            except KeyError:
                total_load = 0
            loads_trace.append(total_load)

        # Build load_shed dict: sum all buses for each (stage, period, commitment)
        load_shed = {}
        for g, val in load_shed_data.items():
            c = list(pyo.ComponentUID(g)._cids)
            _, (stage,) = c.pop(0)
            if stage not in load_shed:
                load_shed[stage] = {}
            stage_dict = load_shed[stage]
            _, (period,) = c.pop(0)
            if period not in stage_dict:
                stage_dict[period] = {}
            period_dict = stage_dict[period]
            _, (commitment,) = c.pop(0)
            if commitment not in period_dict:
                period_dict[commitment] = 0
            period_dict[commitment] += val  # Sum all buses for this time period

        # Build load_shed_trace to match time_periods (repeat for each dispatch)
        load_shed_trace = []
        for s, p, c, d in time_periods:
            try:
                total_shed = load_shed[s][p][c]
            except KeyError:
                total_shed = 0
            load_shed_trace.append(total_shed)
        # print(load_shed_trace)

        HATCH_TO_PATTERN = {
            "": "",  # solid fill
            "....": ".",  # dots
            "////": "/",  # slashes
            "xxxx": "x",  # crosshatch
        }
           
        def plotly_stackgraph(
            times,
            time_periods,
            generation,
            GEN_TYPES,
            GEN_TYPE_ALIASES,
            GEN_TYPE_HATCHES,
            HATCH_TO_PATTERN,
            results_path,
        ):
            """This function creates an interactive Plotly stackgraph
            for given representative days.  Each bar represents one
            hour in one representative day. The x-axis is labeled with
            the representative day and hour (shown at hour 0 and 12).

            """

            n_hours_per_day = 24
            n_rep_days = len(rep_days)
            n_points = n_hours_per_day * n_rep_days

            # Convert the rep_days strings to pandas Timestamps for
            # formatting
            rep_days_dt = [pd.to_datetime(d) for d in rep_days]

            # Build x-axis labels and tick positions: For each hour
            # in each representative day, create a label.  Only show
            # the label for hour 0 and hour 12 of each day, leave
            # others blank for clarity.
            x_labels = []
            tickvals = []
            ticktext = []
            for i, day in enumerate(rep_days_dt):
                for h in range(n_hours_per_day):
                    idx = i * n_hours_per_day + h  # Position in the x-axis
                    if h == 0:
                        label = day.strftime("%b-%d 00:00")
                        x_labels.append(label)
                        tickvals.append(idx)
                        ticktext.append(label)
                    elif h == 12:
                        label = day.strftime("%b-%d 12:00")
                        x_labels.append(label)
                        tickvals.append(idx)
                        ticktext.append(label)
                    else:
                        x_labels.append("")

            # The x-axis for the bars is just integer positions (0 to n_points-1)
            times = list(range(n_points))

            # Prepare traces for each generator type
            traces = []
            for name, color in GEN_TYPES.items():
                label = GEN_TYPE_ALIASES.get(name, name)
                # One value per hour, for all representative days
                values = np.array(
                    [generation[s][p][c][d][name] for s, p, c, d in time_periods]
                )
                hatch = GEN_TYPE_HATCHES.get(name, "")
                pattern_shape = HATCH_TO_PATTERN.get(hatch, "")
                # Use lower opacity for candidate types (those with a
                # hatch)
                opacity = 0.7 if hatch else 1.0

                traces.append(
                    go.Bar(
                        x=times,  # integer positions for each hour
                        y=values,
                        name=label,
                        marker_color=color,
                        marker_pattern_shape=pattern_shape,
                        opacity=opacity,
                        marker_line_width=0,  # remove white line
                    )
                )
            # Add load shed as a stacked bar
            traces.append(
                go.Bar(
                    x=times,
                    y=load_shed_trace,
                    name="Load Shed",
                    marker_color="magenta",
                    opacity=0.7,
                    marker_line_width=0,
                )
            )
            fig = go.Figure(data=traces)
            fig.add_trace(
                go.Scatter(
                    x=times,
                    y=loads_trace,
                    mode="lines+markers",
                    name="Total Load",
                    line=dict(color="black", width=3, dash="dash"),
                    marker=dict(size=4, color="black"),
                    showlegend=True,
                )
            )
            # fig.add_trace(
            #     go.Scatter(
            #         x=times,
            #         y=load_shed_trace,
            #         mode="lines+markers",
            #         name="Load Shed",
            #         line=dict(color="red", width=3, dash="dot"),
            #         marker=dict(size=4, color="magenta"),
            #         showlegend=True,
            #     )
            # )
            fig.update_layout(
                barmode="relative",
                bargap=0,  # remove white spacing between bars
                title="Generation Mix (Representative Days)",
                xaxis=dict(
                    # title="Hours",
                    title="Representative Days (labeled every 12 hours)",
                    tickvals=tickvals,  # show ticks at hour 0 and 12 of each rep day
                    ticktext=ticktext,  # show corresponding label
                    showgrid=True,
                    gridcolor="gray",
                    gridwidth=0.7,
                    linecolor="black",
                    mirror=True,
                ),
                yaxis=dict(
                    title="Nameplate Capacity [MW]",
                    showgrid=True,
                    gridcolor="gray",
                    gridwidth=0.7,
                    linecolor="black",
                    mirror=True,
                ),
                legend=dict(
                    yanchor="middle",
                    y=0.5,
                    xanchor="left",
                    x=1.02,
                    font=dict(size=14),
                    title="Generation Type",
                ),
                width=1200,
                height=600,
                plot_bgcolor="white",
                paper_bgcolor="white",
            )

            # Add vertical lines to visually separate each
            # representative day
            for i in range(1, n_rep_days):
                fig.add_vline(
                    x=i * n_hours_per_day,
                    line=dict(color="gray", width=1, dash="dot"),
                    opacity=0.5,
                )

            # Add a little space above the tallest bar
            all_series = {
                name: np.array(
                    [generation[s][p][c][d][name] for s, p, c, d in time_periods]
                )
                for name in GEN_TYPES
            }

            positive_stack = np.sum(
                [np.clip(vals, 0, None) for vals in all_series.values()],
                axis=0,
            )
            negative_stack = np.sum(
                [np.clip(vals, None, 0) for vals in all_series.values()],
                axis=0,
            )

            ymin = negative_stack.min() if len(negative_stack) else 0
            ymax = positive_stack.max() if len(positive_stack) else 0

            if loads_trace:
                ymax = max(ymax, max(loads_trace))

            lower = ymin * 1.25 if ymin < 0 else -1
            upper = ymax * 1.25 if ymax > 0 else 1

            fig.update_yaxes(
                range=[lower, upper],
                zeroline=True,
                zerolinewidth=2,
                zerolinecolor="black",
            )

            # Save as interactive HTML
            plot_path = f"{results_path}/plots/stackgraph_generators_interactive.html"
            fig.write_html(f"{plot_path}")
            print(f" -> Saved interactive stackgraph to {plot_path}")

            # print("\n[DEBUG] First 10 plotted values by type:")
            # for name in ["battery_charge", "battery_discharge"]:
            #     values = np.array(
            #         [generation[s][p][c][d][name] for s, p, c, d in time_periods]
            #     )
            #     print(f"[DEBUG] {name}: {values[:10]}")

        def calculate_metrics(
                folder_name, rep_days, time_periods, loads_trace, generation, target_day, target_hour
        ):
            """
            Reads the generation.json file and calculates the total generation by generator type.
            
            :param folder_name: Directory containing the generation.json file.
            :return: Dictionary with total generation by type.
            """
            gen_types = [
                "cc_gas", "ct_gas", "coal", "nuclear", "thermal_other", "hydro", "solar",
                "wind", "battery_discharge", "steam", "dr", "ES4", "battery_charge",
                "hydro-c", "gas_cc-c", "gas_ct-c", "battery-c", "wind-c", "pv-c", "steam-c", "ES4-c"
            ]
            
            total_gen_by_type = defaultdict(float)
            file_path = os.path.join(folder_name, "generation.json")
            if not os.path.exists(file_path):
                print(f"[WARNING] File not found: {file_path}")
                return {}

            with open(file_path, "r") as f:
                data = json.load(f)
                
                for key, value in data.items():
                    # The generator type is always at the end after the last '.'
                    gen_name = key.split('.')[-1]
                    for gen_type in gen_types:
                        if gen_name.endswith(gen_type):
                            total_gen_by_type[gen_type] += value
                            break
                        
            def mw_to_gwh(total_mw, hours_per_period=1):
                """
                Converts total MW (sum over all periods) to GWh.
                :param total_mw: Total MW (sum of MW for each period)
                :param hours_per_period: Duration of each period in hours (default 1)
                :return: Total GWh
                """
                total_mwh = total_mw * hours_per_period
                total_gwh = total_mwh / 1000
                return total_gwh

            total_gen_all_types = sum(total_gen_by_type.values())
            total_gen_gwh = mw_to_gwh(total_gen_all_types, hours_per_period=1)
            print(f"Total generation (GWh): {total_gen_gwh}")

            print("Total generation by generator type:")
            for gen_type, total in total_gen_by_type.items():
                total_per_type_gwh = mw_to_gwh(total, hours_per_period=1)
                print(f"  {gen_type}_(GWh): {total_per_type_gwh}")
                
            # Read charging.json and discharging.json, sums all _battery keys
            for file_type in ["charging", "discharging"]:
                file_path = os.path.join(folder_name, f"{file_type}.json")
                if not os.path.exists(file_path):
                    print(f"[WARNING] File not found: {file_path}")
                    continue

                with open(file_path, "r") as f:
                    data = json.load(f)
                    
                total_battery = sum(
                    value for key, value in data.items() if key.endswith("_battery") or key.endswith("_ps")
                )
                total_battery_gwh = mw_to_gwh(total_battery, hours_per_period=1)
                print(f"Total battery {file_type} (GWh): {total_battery_gwh}")

            # Find the index for desired day and time
            try:
                day_idx = rep_days.index(target_day)
            except ValueError:
                raise ValueError(f"Target day {target_day} not found in rep_days!")

            target_idx = day_idx * 24 + target_hour
            target_time_period = time_periods[target_idx]
            
            print(f"Results for representative day: {target_day}, hour: {target_hour} (index {target_idx})")

            # Total load
            total_load_gw = loads_trace[target_idx]/1000
            print(f"Total load (GW): {total_load_gw:.2f}")
            
            # Generation per type
            print("Generation by type:")
            for gen_type in [
                    "cc_gas", "ct_gas", "coal", "nuclear", "thermal_other", "hydro", "solar",
                    "wind", "battery_discharge", "steam", "dr", "ES4", "battery_charge",
                    "hydro-c", "gas_cc-c", "gas_ct-c", "battery-c", "wind-c", "pv-c", "steam-c", "ES4-c"
            ]:
                value = generation[target_time_period[0]][target_time_period[1]][target_time_period[2]][target_time_period[3]].get(gen_type, 0)
                value_gw = value/1000
                print(f"  {gen_type} (GW): {value_gw:.2f}")


        plotly_stackgraph(
            times,
            time_periods,
            generation,
            GEN_TYPES,
            GEN_TYPE_ALIASES,
            GEN_TYPE_HATCHES,
            HATCH_TO_PATTERN,
            results_path,
        )

        calculate_metrics(
            results_path,
            rep_days=rep_days,
            time_periods=time_periods,
            loads_trace=loads_trace,
            generation=generation,
            target_day="2034-07-12 00:00",
            target_hour=19,
        )
    
    def create_html_report(self, results_path, plot_type):

        def html_results_tab(total_cost):
            return f"""
            <div id="Results" class="tabcontent">
            <h2>Results</h2>
            <table>
            <tr><th>Total Cost</th></tr>
            <tr><td>{total_cost}</td></tr>
            </table>
            </div>
            """

        def html_plots_tab(plot_files, plot_type):
            html = """
            <div id="Plots" class="tabcontent">
            <h2>Plots</h2>
            """

            plot_map = {
                "stackplot": ("gen_mix_summary_interactive.html", "Stack Plot"),
                "treemap": ("treemap_2034_interactive.html", "Treemap"),
                "piechart": ("pie_leader_2034_interactive.html", "Pie Chart"),
            }
            categories = ["Dispatchables", "Renewables"]

            for cat in categories:
                cat_lower = cat.lower()
                html += f"<h3>{cat}</h3><ul>\n"
                if plot_type == "all":
                    for _, (fname, label) in plot_map.items():
                        html += f'<li><a href="plots/{cat_lower}_{fname}" target="_blank">{label} ({cat})</a></li>\n'
                elif plot_type in plot_map:
                    fname, label = plot_map[plot_type]
                    html += f'<li><a href="plots/{cat_lower}_{fname}" target="_blank">{label} ({cat})</a></li>\n'
                else:
                    html += "<li>Plot type not supported</li>\n"
                html += "</ul>\n"

            # Add Stackgraph section and plot
            html += "<h3>Stackgraph</h3><ul>\n"
            html += '<li><a href="plots/stackgraph_generators_interactive.html" target="_blank">Stackgraph</a></li>\n'
            html += "</ul>\n"

            html += "</div>"
            return html

        # Read total cost from costs.json
        costs_file = f"{results_path}/costs.json"
        try:
            with open(costs_file, "r") as f:
                costs = json.load(f)
            total_cost = None
            for k in costs:
                if "total" in k.lower():
                    total_cost = costs[k]
                    break
            if total_cost is None and costs:
                total_cost = list(costs.values())[0]
        except Exception as e:
            print(f"[WARNING] Could not read total cost from {costs_file}: {e}")
            total_cost = "N/A"

        # Manually specify plot files
        plot_files = [
            ("Stack Plot (Dispatchables)", "plots/dispatchables_gen_mix_summary.png"),
            ("Stack Plot (Renewables)", "plots/renewables_gen_mix_summary.png"),
            ("Treemap (Dispatchables)", "plots/dispatchables_treemap_2034.png"),
            ("Treemap (Renewables)", "plots/renewables_treemap_2034.png"),
            ("Pie Chart (Dispatchables)", "plots/dispatchables_pie_leader_2034.png"),
            ("Pie Chart (Renewables)", "plots/renewables_pie_leader_2034.png"),
            ("Stackgraph", "plots/stackgraph_generators.png"),
        ]

        html = f"""
        <html>
        <head>
        <title>Expansion Planning Report</title>
        <style>
        body {{ font-family: Arial, sans-serif; }}
        .tab {{
        overflow: hidden;
        border-bottom: 1px solid #ccc;
        background-color: #f1f1f1;
        }}
        .tab button {{
        background-color: inherit;
        float: left;
        border: none;
        outline: none;
        cursor: pointer;
        padding: 14px 16px;
        transition: 0.3s;
        font-size: 17px;
        }}
        .tab button:hover {{ background-color: #ddd; }}
        .tab button.active {{ background-color: #ccc; }}
        .tabcontent {{
        display: none;
        padding: 20px;
        border: 1px solid #ccc;
        border-top: none;
        }}
        table, th, td {{
        border: 1px solid #ccc;
        border-collapse: collapse;
        padding: 8px;
        }}
        </style>
        <script>
        function openTab(evt, tabName) {{
        var i, tabcontent, tablinks;
        tabcontent = document.getElementsByClassName("tabcontent");
        for (i = 0; i < tabcontent.length; i++) {{
        tabcontent[i].style.display = "none";
        }}
        tablinks = document.getElementsByClassName("tablinks");
        for (i = 0; i < tablinks.length; i++) {{
        tablinks[i].className = tablinks[i].className.replace(" active", "");
        }}
        document.getElementById(tabName).style.display = "block";
        evt.currentTarget.className += " active";
        }}
        </script>
        </head>
        <body>
        <h1>Expansion Planning Report</h1>
        <div class="tab">
        <button class="tablinks" onclick="openTab(event, 'Results')" id="defaultOpen">Results</button>
        <button class="tablinks" onclick="openTab(event, 'Plots')">Plots</button>
        </div>
        {html_results_tab(total_cost)}
        {html_plots_tab(plot_files, plot_type)}
        <script>
        document.getElementById("defaultOpen").click();
        </script>
        </body>
        </html>
        """

        html_file = os.path.join(results_path, "gtep_report.html")
        with open(html_file, "w") as f:
            f.write(html)
        print(f" HTML report written to {html_file}")

