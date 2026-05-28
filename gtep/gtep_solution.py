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
from pyomo.core.base.expression import ScalarExpression, IndexedExpression

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import squarify

from gtep.gtep_model import ExpansionPlanningModel
from gtep.tutorial_helper_fns import to_dict, plot_binaries
import calendar

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
        for var in m.component_objects(pyo.Var, descend_into=True):
            for index in var:
                if "Shed" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        load_shed[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Flow" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        power_flow[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Generation" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        generation[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Curtailment" in var.name:
                    if pyo.value(var[index]) >= 0.001:
                        curtailment[var.name + "." + str(index)] = pyo.value(var[index])
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

        renewable_investment_name = folder_name + "/renewable_investments.json"
        dispatchable_investment_name = folder_name + "/dispatchable_investments.json"
        load_shed_name = folder_name + "/load_shed.json"
        costs_name = folder_name + "/costs.json"
        flow_name = folder_name + "/flows.json"
        generation_name = folder_name + "/generation.json"
        curtailment_name = folder_name + "/curtailment.json"

        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        with open(renewable_investment_name, "w") as fil:
            json.dump(renewable_investments, fil)
        with open(dispatchable_investment_name, "w") as fil:
            json.dump(dispatchable_investments, fil)
        with open(load_shed_name, "w") as fil:
            json.dump(load_shed, fil)
        with open(costs_name, "w") as fil:
            json.dump(costs, fil)
        with open(flow_name, "w") as fil:
            json.dump(power_flow, fil)
        with open(generation_name, "w") as fil:
            json.dump(generation, fil)
        with open(curtailment_name, "w") as fil:
            json.dump(curtailment, fil)

        print(f""" -> The following files have been created in '{folder_name}':
         - {renewable_investment_name}
         - {dispatchable_investment_name}
         - {load_shed_name}
         - {costs_name}
         - {flow_name}
         - {generation_name}
         - {curtailment_name}
        """)

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

    def create_plots(self, case_json, results_path, data_path):
        """This function reads a solution .json file, uses gen.csv and
        candidate_generators_initial_list.csv to map generator UIDs to
        unit types and PMax MW, and generates a stacked bar plot of
        generation mix by investment year.

        """

        plots_dir = os.path.join(results_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        print(
            f" Created the subdirectory '{plots_dir}' to save the plots."
        )

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
            "BATTERY": GenerationType("Battery", "#7b9095"),
        }

        # Read gen and candidate_gen .csv files for GEN UID to map for
        # Unit Type and PMax MW
        gen_df = pd.read_csv(f"{data_path}/gen.csv")
        gen_uid_to_type = {
            row["GEN UID"]: row["Unit Type"].upper() for _, row in gen_df.iterrows()
        }
        gen_uid_to_pmax = {
            row["GEN UID"]: float(row["PMax MW"]) for _, row in gen_df.iterrows()
        }

        gen_cand_df = pd.read_csv(f"{data_path}/candidate_generators_initial_list.csv")
        gen_cand_uid_to_type = {
            row["GEN UID"]: row["Unit Type"].upper()
            for _, row in gen_cand_df.iterrows()
        }
        gen_cand_uid_to_pmax = {
            row["GEN UID"]: float(row["PMax MW"]) for _, row in gen_cand_df.iterrows()
        }

        # Read and process .json. These names are based on the saved
        # .json files from function save_results_in_json_files
        if case_json == "renewables":
            json_file = f"{results_path}/renewable_investments.json"
        elif case_json == "dispatchables":
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
            else:
                unit_type = gen_uid_to_type.get(this_key)
                pmax = gen_uid_to_pmax.get(this_key)

            # Make unit_type uppercase to ensure case-insensitive
            # matching
            if unit_type:
                unit_type_upper = unit_type.upper()
            else:
                unit_type_upper = None

            # Only use if in GENERATION_TYPES
            if unit_type and unit_type in GENERATION_TYPES:
                gens_keys_to_type[this_key] = unit_type_upper
                gens_keys_to_pmax[this_key] = pmax
            else:
                raise ValueError(
                    f"[ERROR] Generator '{this_key}' has unknown or unsupported unit type '{unit_type}'."
                )

        # After building gens_keys_to_type
        unique_types = set(gens_keys_to_type.values())

        gen_types_sorted = sorted(unique_types)  # Alphabetical order

        # Read the DAY_AHEAD .csv file with the year column and get
        # unique years in order of appearance
        time_periods_df = pd.read_csv(f"{data_path}/DAY_AHEAD_renewables.csv")
        time_periods = time_periods_df["Year"].drop_duplicates().astype(str).tolist()

        # Build gen_mix using PMax MW and solution values
        gen_mix = {tp: {gt: 0.0 for gt in gen_types_sorted} for tp in time_periods}
        for tp in time_periods:
            for k, val in dict_in.items():
                for state, gen_dict in val.items():
                    for gen, value in gen_dict.items():
                        unit_type = gens_keys_to_type.get(gen)
                        if case_json == "dispatchables":
                            pmax = gens_keys_to_pmax.get(gen, 0.0)
                            if unit_type in gen_mix[tp]:
                                gen_mix[tp][unit_type] += pmax * value
                        elif case_json == "renewables":
                            if unit_type in gen_mix[tp]:
                                gen_mix[tp][unit_type] += value

        gen_mix_arrays = {
            k: np.array([gen_mix[stage].get(k, 0.0) for stage in gen_mix.keys()])
            for k in gen_types_sorted
        }

        # print('gen_mix_arrays:', gen_mix_arrays)

        # Plot 1: Create stack plot with gen mix
        width = 1
        fig, ax = plt.subplots(figsize=(16, 8))
        bottom = np.zeros(len(time_periods))
        for gen_class, mix_array in gen_mix_arrays.items():
            if gen_class in GENERATION_TYPES:
                component_label = GENERATION_TYPES[gen_class].label
                component_color = GENERATION_TYPES[gen_class].color
                p = ax.bar(
                    time_periods,
                    mix_array,
                    width,
                    label=component_label,
                    color=component_color,
                    bottom=bottom,
                    edgecolor="#FFFFFF",
                    linewidth=0.5,
                )
                bottom += mix_array

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        ax.set_title("Generation Mix")
        ax.set_ylabel("Nameplate Capacity [MW]")
        ax.set_xlabel("Investment Year")
        plt.savefig(f"{results_path}/plots/{case_json}_gen_mix_summary.png")
        plt.close()
        print(
            f" -> Saved stack plot for generation mix to {results_path}/{case_json}_gen_mix_summary.png"
        )

        # Plot 2: Create treemaps
        small_pct_threshold = 5
        for tp, mix in gen_mix.items():
            filtered_mix = {
                k: v for k, v in mix.items() if v > 0 and k in GENERATION_TYPES
            }
            if not filtered_mix:
                continue

            total = sum(filtered_mix.values())
            # Prepare lists for squarify
            sorted_items = sorted(filtered_mix.items(), key=lambda x: -x[1])
            sizes = [v for k, v in sorted_items]
            pcts = [v / total * 100 for k, v in sorted_items]
            labels = []
            legend_labels = []
            legend_colors = []

            for (k, v), pct in zip(sorted_items, pcts):
                if pct >= small_pct_threshold:
                    labels.append(
                        f"{GENERATION_TYPES[k].label}\n{int(v)} MW\n{pct:.1f}%"
                    )
                else:
                    # No label for small rectangles, but adding them
                    # in a legend box
                    labels.append("")
                    legend_labels.append(
                        f"{GENERATION_TYPES[k].label} ({int(v)} MW, {pct:.1f}%)"
                    )
                    legend_colors.append(GENERATION_TYPES[k].color)

            colors = [GENERATION_TYPES[k].color for k, v in sorted_items]

            plt.figure(figsize=(10, 7))
            squarify.plot(
                sizes=sizes,
                label=labels,
                color=colors,
                alpha=0.85,
                text_kwargs={"fontsize": 12, "weight": "bold"},
            )
            plt.title(f"{tp}", fontsize=16)
            plt.axis("off")

            # Add legend for small rectangles
            if legend_labels:
                handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in legend_colors]
                plt.legend(
                    handles,
                    legend_labels,
                    title=f"Small Rectangles (<{small_pct_threshold}%)",
                    loc="center left",
                    bbox_to_anchor=(1.05, 0.5),
                    fontsize=12,
                    title_fontsize=13,
                    frameon=True,
                )

            plt.tight_layout()
            plt.savefig(
                f"{results_path}/plots/{case_json}_treemap_{tp}.png",
                dpi=150,
                bbox_inches="tight",
            )
            plt.close()
            print(
                f" -> Saved treemap for {tp} to {results_path}/{case_json}_treemap_{tp}.png"
            )

        # Plot 3: Create pie chart
        for tp, mix in gen_mix.items():
            filtered_mix = {
                k: v for k, v in mix.items() if v > 0 and k in GENERATION_TYPES
            }
            if not filtered_mix:
                continue

            sizes = [filtered_mix[k] for k in filtered_mix]
            labels = [
                f"{GENERATION_TYPES[k].label} ({int(filtered_mix[k])} MW, {filtered_mix[k]/sum(sizes)*100:.1f}%)"
                for k in filtered_mix
            ]
            colors = [GENERATION_TYPES[k].color for k in filtered_mix]

            # Add a small explode for each slice to create space
            # between slices
            explode = [0.03] * len(sizes)

            fig, ax = plt.subplots(figsize=(10, 10))
            wedges, _ = ax.pie(
                sizes,
                labels=None,
                colors=colors,
                startangle=90,
                explode=explode,
                wedgeprops={"edgecolor": "none"},  # No lines between slices
                shadow=False,
            )

            total = sum(sizes)
            for i, (wedge, label) in enumerate(zip(wedges, labels)):
                size = sizes[i]
                angle = (wedge.theta2 + wedge.theta1) / 2
                x = np.cos(np.deg2rad(angle))
                y = np.sin(np.deg2rad(angle))
                # Position label outside the pie
                # print(size/total)
                if size / total < 0.0001:
                    label_dist = 1.8
                elif size / total < 0.0005:
                    label_dist = 1.6
                elif size / total < 0.001:
                    label_dist = 1.5
                elif size / total < 0.005:
                    label_dist = 1.4
                elif size / total < 0.01:
                    label_dist = 1.2
                elif i % 2 == 0:
                    label_dist = 1.6
                else:
                    label_dist = 1.75

                label_x = label_dist * x
                label_y = label_dist * y

                ax.annotate(
                    label,
                    xy=(x, y),
                    xytext=(label_x, label_y),
                    ha="center",
                    va="center",
                    fontsize=12,
                    weight="bold",
                    arrowprops=dict(arrowstyle="-", color="gray", lw=1.5),
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        fc="white",
                        ec="gray",
                        lw=0.8,
                        alpha=0.7,
                    ),
                )

            ax.axis("equal")
            plt.tight_layout()
            plt.savefig(
                f"{results_path}/plots/{case_json}_pie_leader_{tp}.png",
                dpi=150,
                bbox_inches="tight",
            )
            plt.close()
            print(
                f" -> Saved pie chart with leader lines for {tp} to {results_path}/{case_json}_pie_leader_{tp}.png"
            )

    def create_stackgraph(self, results_path, data_path):

        # Read the DAY_AHEAD .csv file with the year column and get
        # unique years in order of appearance
        years_df = pd.read_csv(f"{data_path}/DAY_AHEAD_renewables.csv")
        years_val = years_df["Year"].drop_duplicates().astype(str).tolist()

        try:
            import ujson as json
        except ImportError:
            import json

        with open(f"{results_path}/generation.json", "r") as F:
            gen_data = json.load(F)

        # CC	CC	CC	tab20	1
        # CT	CT	CT	tab20	3
        # COAL	COAL	CO	tab20	5
        # NUCLEAR	NUCLEAR	NU	tab20	2
        # PV	PV	PV	tab20	9
        # WIND	WIND	WI	tab20	11
        # THERMAL	THERMAL	TH	tab20	13
        # HYDRO	HYDRO	HY	tab20	20
        # BATT	BATT	BA	tab20	15
        # ES4	ES4	ES	tab20	17
        # PS	PS	PS	tab20	19
        # Load Shed	Load Shed	SL	tab20	7

        # [ESR: Comment out original colors for now]
        # tab20 = plt.get_cmap("tab20")
        # GEN_TYPES = {
        #     "coal": tab20(5),
        #     "cc_gas": tab20(1),
        #     "ct_gas": tab20(3),
        #     "dr": tab20(19),
        #     "solar": tab20(9),
        #     "thermal_other": tab20(13),
        #     "wind": tab20(11),
        #     "gas_cc-c": tab20(0),
        #     "gas_ct-c": tab20(2),
        #     "pv-c": tab20(8),
        #     "wind-c": tab20(10),
        # }

        # Note that these are the same colors used in the stack plots
        # and pie charts above
        GEN_TYPES = {
            "coal": "#333333",
            "cc_gas": "#20b2aa",
            "ct_gas": "#6e8b3d",
            "dr": "#a020f0",
            "solar": "#ffb90f",
            "thermal_other": "#e25822",
            "wind": "#4f94cd",
            # "hydro": "#00bfff",
            # "battery": "#7b9095",
            # "nuclear": "#39FF14",
            # "steam": "#b0b0b0",
            # Candidates use the same color as their base type
            "gas_cc-c": "#20b2aa",
            "gas_ct-c": "#6e8b3d",
            "pv-c": "#ffb90f",
            "wind-c": "#4f94cd",
            # "hydro-c": "#00bfff",
            # "battery-c": "#7b9095",
            # "steam-c": "#b0b0b0",
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
            # Candidates get a hatch pattern
            "gas_cc-c": "....",
            "gas_ct-c": "....",
            "pv-c": "////",
            "wind-c": "xxxx",
            # "hydro-c": "....",
        }
        GEN_TYPE_ALIASES = {
            "coal": "Coal",
            "cc_gas": "CC",
            "ct_gas": "CT",
            "dr": "DR",
            "solar": "Solar",
            "thermal_other": "Thermal",
            "wind": "Wind",
            # "hydro": "Hydro",
            "gas_cc-c": "CC (Candidate)",
            "gas_ct-c": "CT (Candidate)",
            "pv-c": "Solar (Candidate)",
            "wind-c": "Wind (Candidate)",
            # "hydro-c": "Hydro (Candidate)",
        }

        generation = {}
        for g, val in gen_data.items():
            c = list(pyo.ComponentUID(g)._cids)
            _, (stage,) = c.pop(0)
            if stage not in generation:
                generation[stage] = {}
            stage = generation[stage]
            _, (period,) = c.pop(0)
            if period not in stage:
                stage[period] = {}
            period = stage[period]
            _, (committment,) = c.pop(0)
            if committment not in period:
                period[committment] = {}
            committment = period[committment]
            _, (dispatch,) = c.pop(0)
            if dispatch not in committment:
                committment[dispatch] = dict.fromkeys(GEN_TYPES, 0)
            dispatch = committment[dispatch]
            gen_name = c[-1][0]
            _type = None
            for gt in GEN_TYPES:
                if gen_name.endswith(gt):
                    _type = gt
            if _type is None:
                raise RuntimeError(f"Cannot map generator name '{gen_name}' to type")
            dispatch[_type] += val

        time_periods = [
            (s, p, c, d)
            for s in generation
            for p in generation[s]
            for c in generation[s][p]
            for d in generation[s][p][c]
        ]
        times = list(range(len(time_periods)))

        fig, ax = plt.subplots(figsize=(18, 8))
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        bottom = np.zeros(len(time_periods))
        width = 1

        for name, color in GEN_TYPES.items():
            hatch = GEN_TYPE_HATCHES.get(name, "")
            label = GEN_TYPE_ALIASES.get(name, name)
            values = np.array(
                [generation[s][p][c][d][name] for s, p, c, d in time_periods]
            )
            p = ax.bar(
                times,
                values,
                width,
                # label=name,
                label=label,
                color=color,
                bottom=bottom,
                # edgecolor="#FFFFFF",
                edgecolor=None,  # no white lines between bars
                # linewidth=0.5,
                linewidth=0.05,
                hatch=hatch,
            )
            bottom += values

        # [TODO: We are only plotting for one year, but we should
        # extend this for the case when we have multiple years.]
        year = int(years_val[0])
        month_starts = []
        day = 0
        for month in range(1, 13):
            month_starts.append(day)
            days_in_month = calendar.monthrange(year, month)[1]
            day += days_in_month
        month_labels = [calendar.month_abbr[m] for m in range(1, 13)]

        # Add vertical lines at month boundaries
        for ms in month_starts:
            ax.axvline(ms, color='gray', linestyle=':', linewidth=1, alpha=0.5)

        # Set x-ticks and labels
        ax.set_xticks(month_starts)
        ax.set_xticklabels(month_labels, fontsize=12)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14, title="Generation Type")
        ax.set_title(f"Generation Mix for Investment Year {year}", fontsize=20)
        ax.set_ylabel('Nameplate Capacity [MW]', fontsize=16)
        ax.set_xlabel('Months', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.grid(axis='both', color='gray', linestyle=':', linewidth=0.7, alpha=0.5)

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.10)
        
        plt.tight_layout()
        plt.savefig(
            f"{results_path}/plots/stackgraph_generators.png", dpi=150, bbox_inches="tight"
        )
        plt.close()
        print(f" -> Saved stackgraph to {results_path}/stackgraph_generators.png")


    def create_html_report(self, results_path):

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

        def html_plots_tab(plot_files):
            html = """
            <div id="Plots" class="tabcontent">
            <h2>Plots</h2>
            """

            # Group plots by category
            renewables = []
            dispatchables = []
            stackgraph = []
            for label, fname in plot_files:
                label_lower = label.lower()
                if "renewable" in label_lower:
                    renewables.append((label, fname))
                elif "dispatchable" in label_lower:
                    dispatchables.append((label, fname))
                elif "stackgraph" in label_lower:
                    stackgraph.append((label, fname))

            if dispatchables:
                html += "<h3>Dispatchables</h3><ul>\n"
                for label, fname in dispatchables:
                    html += f'<li><a href="{fname}" target="_blank">{label}</a></li>\n'
                html += "</ul>\n"

            if renewables:
                html += "<h3>Renewables</h3><ul>\n"
                for label, fname in renewables:
                    html += f'<li><a href="{fname}" target="_blank">{label}</a></li>\n'
                html += "</ul>\n"

            if stackgraph:
                html += "<h3>Stackgraph</h3><ul>\n"
                for label, fname in stackgraph:
                    html += f'<li><a href="{fname}" target="_blank">{label}</a></li>\n'
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
        {html_plots_tab(plot_files)}
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


