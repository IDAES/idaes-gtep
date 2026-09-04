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
import csv
import json
import logging
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import plotly.graph_objects as go

from collections import namedtuple, defaultdict
from pathlib import Path

import pyomo.environ as pyo
import pyomo.gdp as gdp
from pyomo.environ import units as u
from pyomo.core.base.param import IndexedParam
from pyomo.core.base.expression import ScalarExpression, IndexedExpression

from gtep.gtep_model import ExpansionPlanningModel

logger = logging.getLogger(__name__)


class ExpansionPlanningSolution:
    """This class stores the solution to the ExpansionPlanningModel
    class for writing and visualization.

    """

    def __init__(self, data_path):
        self.gen_df = pd.read_csv(f"{data_path}/gen.csv")
        self.storage_csv_path = os.path.join(data_path, "storage.csv")
        if os.path.exists(self.storage_csv_path):
            self.storage_df = pd.read_csv(self.storage_csv_path)
        self._get_generation_types()

    def _get_generation_types(self):
        """This method returns generation type labels and colors used
        in plots and stackgraphs.

        """

        GenerationType = namedtuple("GenerationType", ["label", "color"])
        tab20 = plt.get_cmap("tab20")

        def darken_color(hex_color, percent=0.2):
            hex_color = hex_color.lstrip("#")
            rgb = [int(hex_color[i : i + 2], 16) for i in (0, 2, 4)]
            darker_rgb = [max(0, int(c * (1 - percent))) for c in rgb]
            return "#" + "".join(f"{c:02x}" for c in darker_rgb)

        gas_cc = mcolors.to_hex(tab20(1))
        gas_ct = mcolors.to_hex(tab20(3))
        coal = mcolors.to_hex(tab20(5))
        nuclear = mcolors.to_hex(tab20(2))
        solar = mcolors.to_hex(tab20(9))
        rt_solar = mcolors.to_hex(tab20(8))
        wind = mcolors.to_hex(tab20(11))
        thermal = mcolors.to_hex(tab20(13))
        steam = mcolors.to_hex(tab20(14))
        hydro = mcolors.to_hex(tab20(19))
        storage = mcolors.to_hex(tab20(15))
        other = mcolors.to_hex(tab20(0))

        self.gen_types = {
            "CC": GenerationType("Gas CC", gas_cc),
            "CT": GenerationType("Gas CT", gas_ct),
            "COAL": GenerationType("Coal", coal),
            "NUC": GenerationType("Nuclear", nuclear),
            "PV": GenerationType("Solar", solar),
            "RTPV": GenerationType("RT Solar", rt_solar),
            "WIND": GenerationType("Wind", wind),
            "THERMAL": GenerationType("Thermal", thermal),
            "STEAM": GenerationType("Steam", steam),
            "HYDRO": GenerationType("Hydro", hydro),
            "BATTERY": GenerationType("Storage", storage),
            "PS": GenerationType("Pumped Storage", storage),
            "OTHER": GenerationType("Other", other),
        }
        gen_types_candidate = {
            unit
            + "-c": GenerationType(
                gentype.label + " Candidate", darken_color(gentype.color)
            )
            for unit, gentype in self.gen_types.items()
        }
        for k, v in gen_types_candidate.items():
            self.gen_types[k] = v

    def load_from_model(self, gtep_model):
        """This method loads the results from the solved model
        and the metadata into the solution object.

        This method stores solver results, model dimensions, input
        data, and selected commitment/investment expression values
        from a solved ExpansionPlanningModel.

        """
        # Check that the input is a GTEP model object.
        if type(gtep_model) is not ExpansionPlanningModel:
            logger.warning(
                f"Solutions must be loaded from ExpansionPlanningModel objects, not %s"
                % type(gtep_model)
            )
            raise ValueError

        # Check that the model has solver results.
        if gtep_model.results is None:
            raise ValueError(
                "ExpansionPlanningSolution objects loaded from model must have a results component."
            )

        # Store solver results.
        self.results = gtep_model.results

        # Store model dimensions and formulation metadata.
        self.stages = gtep_model.stages
        self.formulation = gtep_model.formulation
        self.data = gtep_model.data
        self.num_reps = gtep_model.num_reps
        self.num_commit = gtep_model.num_commit
        self.num_dispatch = gtep_model.num_dispatch

        # Store selected expression values for validation/reporting.
        self.expressions = {
            expr.name: pyo.value(expr)
            for expr in gtep_model.model.component_data_objects(pyo.Expression)
            if ("Commitment" in expr.name) or ("Investment" in expr.name)
        }

    def _to_dict(self) -> dict:
        """This method converts solution results into a
        JSON-friendly dictionary.

        This method stores solver summary information, primal variable
        values, selected expression values, and nested result trees
        for downstream validation, reporting, or serialization.

        """

        # Store top-level solver results and expression values.
        results_dict = {
            "solution_loader": self.results.solution_loader,
            "termination_condition": self.results.termination_condition,
            "best_feasible_objective": self.results.best_feasible_objective,
            "best_objective_bound": self.results.best_objective_bound,
            "wallclock_time": self.results.wallclock_time,
            "expressions": self.expressions,
        }

        # Convert termination condition to a JSON-friendly dictionary.
        results_dict["termination_condition"] = {
            "value": self.results.termination_condition.value,
            "name": self.results.termination_condition.name,
        }

        # Store flat primal variable values.
        results_dict["solution_loader"] = {"primals": {}}
        for key, val in self.results.solution_loader.get_primals()._dict.items():
            tmp_key = key

            results_dict["solution_loader"]["primals"][tmp_key] = {
                "name": val[0].name,
                "value": val[0].value,
                "bounds": val[0].bounds,
            }

            # Add binary flag when applicable.
            if val[0].is_binary():
                results_dict["solution_loader"]["primals"][tmp_key]["is_binary"] = val[
                    0
                ].is_binary()

            # Add units when available.
            if val[0].get_units() is not None:
                results_dict["solution_loader"]["primals"][tmp_key]["units"] = (
                    val[0].get_units().name
                )
            else:
                results_dict["solution_loader"]["primals"][tmp_key]["units"] = val[
                    0
                ].get_units()

        # Initialize nested result trees.
        results_dict["primals_tree"] = {}
        results_dict["expressions_tree"] = {}
        for key, val in self.results.solution_loader.get_primals()._dict.items():
            # Split variable name to define nesting depth.
            split_name = val[0].name.split(".")

            tmp_dict = {
                "name": val[0].name,
                "value": val[0].value,
                "bounds": val[0].bounds,
            }

            # Add binary flag when applicable.
            if val[0].is_binary():
                tmp_dict["is_binary"] = val[0].is_binary()

            # Add units when available.
            if val[0].get_units() is not None:
                tmp_dict["units"] = val[0].get_units().name
            else:
                tmp_dict["units"] = val[0].get_units()

            # Add primal variable to nested dictionary.
            def nested_set(this_dict, key, val):
                if len(key) > 1:
                    if key[1] == "binary_indicator_var":
                        this_dict[key[0]] = val
                    else:
                        this_dict.setdefault(key[0], {})
                        nested_set(this_dict[key[0]], key[1:], val)
                else:
                    this_dict[key[0]] = val

            nested_set(results_dict["primals_tree"], split_name, tmp_dict)

        for key, val in self.expressions.items():
            # Split expression name to define nesting depth.
            split_name = key.split(".")

            tmp_dict = {"value": val}

            # Add expression to nested dictionary.
            def nested_set(this_dict, key, val):
                if len(key) > 1:
                    this_dict.setdefault(key[0], {})
                    nested_set(this_dict[key[0]], key[1:], val)
                else:
                    this_dict[key[0]] = val

            nested_set(results_dict["expressions_tree"], split_name, tmp_dict)

        # Store nested expression and primal trees on the solution
        # object.
        self.expressions_tree = results_dict["expressions_tree"]
        self.primals_tree = results_dict["primals_tree"]

        # Build final output dictionary.
        out_dict = {
            "data": self.data.representative_data[0].data,
            "results": results_dict,
        }

        return out_dict

    def read_json(self, filepath):
        # Read a json file
        json_filepath = Path(filepath)
        with open(json_filepath, "r") as fobj:
            json_read = json.loads(fobj.read())

        return json_read

    def to_nested_dict(self, dict_in):
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

    def create_results_directory(self, dir_name):
        """This method creates a directory to save model results.

        :param dir_name: Name or path of the directory where results
                         will be saved. Defaults to "results".
        :return: Directory name/path.

        """
        os.makedirs(dir_name, exist_ok=True)
        logger.info(f"\nCreating the directory '{dir_name}' to save the results. ")

        return dir_name

    def save_model_config_to_csv(self, gtep_model, dir_name):
        """Save model configuration settings to a CSV file.

        :param gtep_model: Expansion planning model object.
        :param dir_name: Directory where the CSV file is written.
        """
        config_csv_path = f"{dir_name}/model_config.csv"

        with open(config_csv_path, mode="w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["config_key", "config_value", "value_type"])

            for key, value in sorted(gtep_model.config.items()):
                writer.writerow([key, repr(value), type(value).__name__])

        logger.info(f">>>Saved model configuration to: {config_csv_path}")

    def save_results_in_json_files(self, gtep_model, dir_name, value_threshold=1e-3):
        """This method saves the model results to JSON files.

        Outputs include investments, load shed, costs, flows,
        generation, curtailment, loads, reserves, and storage
        charge/discharge. Creates the results directory if needed.

        Only variable values greater than or equal to the argument
        `value_threshold` are saved to avoid writing near-zero
        numerical values. The default threshold is 1e-3.

        :param gtep_model: Solved expansion planning model object.
        :param dir_name: Directory where JSON files are written.
        :param value_threshold: Minimum variable value to save in the
                                JSON outputs. Defaults to ``1e-3``.

        """

        folder_name = dir_name
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
                # Save only values above value_threshold to avoid
                # writing near-zero numerical values. The threshold is
                # configurable with a default value of 1e-3.
                if "Shed" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        load_shed[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Reserve" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        reserves[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Flow" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        power_flow[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Generation" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        generation[var.name + "." + str(index)] = pyo.value(var[index])
                elif "Curtailment" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        curtailment[var.name + "." + str(index)] = pyo.value(var[index])
                elif "storageCharged" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        charging[var.name + "." + str(index)] = pyo.value(var[index])
                elif "storageDischarge" in var.name:
                    if pyo.value(var[index]) >= value_threshold:
                        discharging[var.name + "." + str(index)] = pyo.value(var[index])
                for name in valid_names:
                    if name in var.name:
                        if pyo.value(var[index]) >= value_threshold:
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

        # Loads are currently read through Prescient, which maps them
        # to buses and stores them as indexed parameters. If a future
        # workflow loads non-Prescient scalar load parameters, this
        # logic may need to be updated.
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

        logger.info(
            f"The following files have been created in the directory '{folder_name}':"
        )
        for name in output_files:
            logger.info(f" - {folder_name}/{name}.json")

    def create_plots(
        self, case_json, results_path, data_path, plot_type="all", savefig=True
    ):
        """This method creates generation-mix plots from saved
        solution JSON files. These plots visualize the total amount
        of generation capacity available based on investment decisions.

        It reads investment results for renewable, dispatchable, or
        combined generation and maps generator/storage IDs to unit
        types using `gen.csv` and, when available, `storage.csv`. It
        creates interactive Plotly treemap and/or pie chart HTML.
        The figure(s) are returned and also saved to HTML
        file(s) in the results `plots` directory if `savefig=True`.

        :param case_json: Results group to plot. Options are
                          "renewables", "dispatchables", or
                          "combined".
        :param results_path: Directory containing saved JSON result
                             files.
        :param data_path: Directory containing input data files.
        :param plot_type: Plot type to generate. Options are
                          "treemap", "piechart", or "all".
                          Defaults to "all".
        :param savefig:   Whether to write the figure(s) to file

        """
        if plot_type not in ["treemap", "piechart", "all"]:
            raise ValueError(
                f"Plot type '{plot_type}' is not supported. "
                "Please choose between 'treemap', 'piechart', or 'all."
            )

        if savefig:
            plots_dir = os.path.join(results_path, "plots")
            if not os.path.exists(plots_dir):
                os.makedirs(plots_dir)
                logger.info(
                    f"\nCreated the subdirectory '{plots_dir}' to save the plots."
                )

        def get_gen_arrays(gen_case_json, results_path, data_path, gen_types):
            """This function builds generation-mix dictionaries used
            by plotting functions.

            It reads generator and optional storage data, loads the
            selected investment JSON file, maps assets to generation
            types, and returns generation capacity by time period and
            unit type.

            :param gen_case_json: Results group to process, either
                                  "renewables" or "dispatchables".
            :param results_path: Directory containing saved JSON result
                                 files.
            :param data_path: Directory containing input data files.
            :param gen_types: Mapping of supported unit types to
                                     plot labels and colors.
            :return: Tuple containing gen_mix, gen_mix_arrays,
                     and time_periods.

            """

            # Map generators IDs to Unit Type and PMax
            gen_uid_to_type = {
                row["GEN UID"]: row["Unit Type"].upper()
                for _, row in self.gen_df.iterrows()
            }
            gen_uid_to_pmax = {
                row["GEN UID"]: float(row["PMax MW"])
                for _, row in self.gen_df.iterrows()
            }

            # Read and process saved JSON files for renewables and
            # dispatchables units.
            json_files = {
                "renewables": os.path.join(results_path, "renewable_investments.json"),
                "dispatchables": os.path.join(
                    results_path, "dispatchable_investments.json"
                ),
            }

            if gen_case_json not in json_files:
                raise ValueError(
                    f"Unsupported gen_case_json '{gen_case_json}'. "
                    "Choose 'renewables' or 'dispatchables'."
                )

            dict_in = self.to_nested_dict(self.read_json(json_files[gen_case_json]))

            # Collect all generator/storage keys that appear in the
            # investment results. For this, we loop over investment
            # stages/time keys, over investment-state dictionaries,
            # and at the end collect each generator/storage ID. Use
            # set() to remove duplicates.
            time_keys = list(dict_in.keys())
            keys_set = set()
            for time_key in time_keys:
                for state_dict in dict_in[time_key].values():
                    for asset in state_dict:
                        keys_set.add(asset)

            # Map generator keys to unit types and PMax MW
            gens_keys_to_type = {}
            gens_keys_to_pmax = {}
            for this_key in list(keys_set):
                if this_key in gen_uid_to_type:
                    unit_type = gen_uid_to_type.get(this_key)
                    pmax = gen_uid_to_pmax.get(this_key)
                else:
                    unit_type = None
                    pmax = 0

                # Make unit_type uppercase to ensure case-insensitive
                # matching
                unit_type_upper = unit_type.upper() if unit_type else None

                # Only use if in gen_types dictionary.
                if unit_type_upper and unit_type_upper in gen_types:
                    gens_keys_to_type[this_key] = unit_type_upper
                    gens_keys_to_pmax[this_key] = pmax
                else:
                    raise ValueError(
                        f"[ERROR] Generator or storage '{this_key}' has unknown or unsupported unit type '{unit_type}'."
                    )

            # After building gens_keys_to_type
            unique_types = set(gens_keys_to_type.values())
            gen_types_sorted = sorted(unique_types)  # Alphabetical order

            # Get the modeled year(s) (in order of appearance) from
            # the DAY_AHEAD_renewables.csv time-series.
            time_periods_df = pd.read_csv(f"{data_path}/DAY_AHEAD_renewables.csv")
            time_periods = (
                time_periods_df["Year"].drop_duplicates().astype(str).tolist()
            )

            # Build generation mix by time period and generation type
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
                gen_type: np.array(
                    [gen_mix[tp].get(gen_type, 0.0) for tp in time_periods]
                )
                for gen_type in gen_types_sorted
            }

            return gen_mix, gen_mix_arrays, time_periods

        # Define functions to create a treemap and pie chart
        # interactive Plotly plots for the generation mix. The user
        # can select which one to use by setting up the plot_type
        # option.
        def plotly_treemap_gen_mix(
            gen_mix, gen_types, results_path, case_json, small_pct_threshold=5
        ):
            """This function creates interactive Plotly treemap plots
            of generation mix.

            One HTML treemap is created and shows the share of
            generation capacity by the unit type defined in case_json.

            :param gen_mix: Dictionary of generation mix by time
                            period and unit type.
            :param gen_types: Mapping of unit types to plot labels and
                              colors.
            :param results_path: Directory where plot files are saved.
            :param case_json: Results group name used in output file
                              names. Options are "renewables",
                              "dispatchables", or "combined".
            :param small_pct_threshold: Minimum percentage used for
                                        displaying labels. Defaults to
                                        5.

            """

            for tp, mix in gen_mix.items():
                filtered_mix = {
                    k: v for k, v in mix.items() if v > 0 and k in gen_types
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
                    label = f"{gen_types[k].label}<br>{int(v)} MW<br>{pct:.1f}%"
                    labels.append(label if pct >= small_pct_threshold else "")
                    colors.append(gen_types[k].color)
                    customdata.append(f"{gen_types[k].label} ({int(v)} MW, {pct:.1f}%)")

                fig = go.Figure(
                    go.Treemap(
                        labels=[gen_types[k].label for k, v in sorted_items],
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

                fname = f"{results_path}/plots/treemap_{case_json}_{tp}.html"
                return fig, fname

        def plotly_pie_gen_mix(gen_mix, gen_types, results_path, case_json):
            """This function creates interactive Plotly pie charts of
            generation mix.

            One HTML pie chart is created for each unit type and shows
            the share of generation capacity by unit type.

            :param gen_mix: Dictionary of generation mix by time period
                            and unit type.
            :param gen_types: Mapping of unit types to plot labels and
                              colors.
            :param results_path: Directory where plot files are saved.
            :param case_json: Results group name used in output file
                              names.

            """

            for tp, mix in gen_mix.items():
                filtered_mix = {
                    k: v for k, v in mix.items() if v > 0 and k in gen_types
                }
                if not filtered_mix:
                    continue

                sizes = [filtered_mix[k] for k in filtered_mix]
                total = sum(sizes)
                labels = [
                    f"{gen_types[k].label} ({int(filtered_mix[k])} MW, {filtered_mix[k]/total*100:.1f}%)"
                    for k in filtered_mix
                ]
                colors = [gen_types[k].color for k in filtered_mix]

                # Plotly pie chart
                fig = go.Figure(
                    go.Pie(
                        labels=[gen_types[k].label for k in filtered_mix],
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

                fname = f"{results_path}/plots/piechart_{case_json}_{tp}.html"
                return fig, fname

        # ------------------------------------------------------------
        # The creation of plots under create_plots starts here, after
        # the plotting helper functions.

        if case_json != "combined":
            # Create gen_mix dictionary and arrays needed to plot
            # renewables and dispatchables types in separate plots.
            gen_mix, gen_mix_arrays, time_periods = get_gen_arrays(
                case_json, results_path, data_path, self.gen_types
            )
        else:
            # Combine gen_mix dictionary and arrays needed to plot
            # renewables and dispatchables types in the same plot.
            gen_mix_ren, gen_mix_arrays_ren, time_periods_ren = get_gen_arrays(
                "renewables", results_path, data_path, self.gen_types
            )
            gen_mix_disp, gen_mix_arrays_disp, time_periods_disp = get_gen_arrays(
                "dispatchables", results_path, data_path, self.gen_types
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

            # Merge gen_mix and gen_mix arrays
            gen_mix = {}
            for tp in all_time_periods:
                gen_mix[tp] = {}
                for gt in all_types:
                    val_ren = gen_mix_ren.get(tp, {}).get(gt, 0.0)
                    val_disp = gen_mix_disp.get(tp, {}).get(gt, 0.0)
                    gen_mix[tp][gt] = val_ren + val_disp

            gen_mix_arrays = {
                gt: np.array([gen_mix[tp].get(gt, 0.0) for tp in all_time_periods])
                for gt in all_types
            }

        figs, fnames = [], []
        if plot_type in ("treemap", "all"):
            fig, fname = plotly_treemap_gen_mix(
                gen_mix, self.gen_types, results_path, case_json
            )
            figs.append(fig), fnames.append(fname)
        if plot_type in ("piechart", "all"):
            fig, fname = plotly_pie_gen_mix(
                gen_mix, self.gen_types, results_path, case_json
            )
            figs.append(fig), fnames.append(fname)

        if savefig:
            for fig, fname in zip(figs, fnames):
                fig.write_html(fname)
                logger.info(f" -> Saved to {fname}")

        return figs[0] if len(figs) == 1 else figs

    def create_stackgraph(self, results_path, rep_days, savefig=True):
        """This method creates and saves an interactive stackgraph of
        dispatch results.

        This method reads saved JSON result files, organizes
        generation, storage charge/discharge if enabled, load, and
        load-shed values by stage, representative period, commitment
        period, and dispatch period, and creates a Plotly stacked bar
        chart with total load.

        :param results_path: Directory containing saved JSON result
                             files and where the plot will be saved.
        :param rep_days: List of representative day labels used for
                         formatting the x-axis.
        :param savefig: Whether to save the figure to file. Defaults
                            to `True`.

        """

        def load_json(name):
            with open(f"{results_path}/{name}.json", "r") as f:
                return json.load(f)

        gen_data = load_json("generation")
        loads_data = load_json("loads")
        load_shed_data = load_json("load_shed")
        reserves_data = load_json("reserves")
        charging_data = load_json("charging")
        discharging_data = load_json("discharging")

        # Use the generation-type mapping for stackgraph colors,
        # labels, and candidate hatch patterns.
        GEN_TYPES = {key: val.color for key, val in self.gen_types.items()}
        GEN_TYPE_ALIASES = {key: val.label for key, val in self.gen_types.items()}
        GEN_TYPE_HATCHES = {
            key: "/" if str(key).endswith("-c") else "" for key in GEN_TYPES
        }

        gen_uid_to_type = {
            row["GEN UID"]: (
                row["Unit Type"].upper() + "-c"
                if row["GEN UID"].endswith("-c")
                else row["Unit Type"].upper()
            )
            for _, row in self.gen_df.iterrows()
        }
        missing_unit_types = [
            unit_type
            for unit_type in gen_uid_to_type.values()
            if unit_type not in GEN_TYPES
        ]
        if missing_unit_types:
            raise ValueError(f"Unit type(s) {missing_unit_types} not recognized")

        # Create dictionary with generation by time index and
        # generator type for plotting. For this, parse each generation
        # variable name to recover its stage, representative period,
        # commitment period, dispatch period, and generator name.
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

            gen_name = str(c[-1][0])
            try:
                _type = gen_uid_to_type[gen_name]
            except KeyError:
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

        total_charging = sum(charging_data.values())
        total_discharging = sum(discharging_data.values())

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

        # Build load_shed dict: sum all buses for each (stage, period,
        # commitment)
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

        # Build load_shed_trace to match time_periods (repeat for each
        # dispatch)
        load_shed_trace = []
        for s, p, c, d in time_periods:
            try:
                total_shed = load_shed[s][p][c]
            except KeyError:
                total_shed = 0
            load_shed_trace.append(total_shed)

        # --------------------------------------------------------------
        # Build and save an interactive Plotly stackgraph of
        # generation results. The plot stacks generation by technology
        # over all representative-day hours, marks candidate resources
        # with a pattern and includes load shed and total load as a
        # dashed line. The x-axis is labeled at hour 0 and hour 12 for
        # each representative day.
        n_hours_per_day = 24
        n_rep_days = len(rep_days)
        n_points = n_hours_per_day * n_rep_days
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

        # The x-axis for the bars is just integer positions (0 to
        # n_points-1)
        times = list(range(n_points))

        # Prepare traces for each generator type
        traces = []
        for name, color in GEN_TYPES.items():
            label = GEN_TYPE_ALIASES.get(name, name)
            # One value per hour, for all representative days
            values = np.array(
                [generation[s][p][c][d][name] for s, p, c, d in time_periods]
            )
            pattern_shape = GEN_TYPE_HATCHES.get(name, "")
            # Use lower opacity for candidate types (those with a
            # hatch)
            opacity = 0.7 if pattern_shape else 1.0

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
        tab20 = plt.get_cmap("tab20")
        traces.append(
            go.Bar(
                x=times,
                y=load_shed_trace,
                name="Load Shed",
                marker_color=mcolors.to_hex(tab20(7)),
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
        fig.update_layout(
            barmode="relative",
            bargap=0,  # remove white spacing between bars
            title="Generation Mix (Representative Days)",
            xaxis=dict(
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

        if savefig:
            plot_path = f"{results_path}/plots/stackgraph_generators.html"
            fig.write_html(f"{plot_path}")
            logger.info(f" -> Saved interactive stackgraph to {plot_path}")
        return fig
