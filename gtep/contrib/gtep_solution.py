# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield, Thom R. Edwards
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
import logging

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import networkx as nx
import pandas as pd
import numpy as np
import re


from matplotlib.patches import Rectangle, RegularPolygon, PathPatch
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
from matplotlib.transforms import Affine2D
from matplotlib.colors import Normalize
import matplotlib.path as mpath

logger = logging.getLogger(__name__)


# [TODO] inject units into plots
class ExpansionPlanningSolution:
    def __init__(self):
        # PopPop the Power Optimization Possum says ౿ᓕ ̤Ꜥ·> --- "eat trash, heck __init__, hold me"
        pass

    def load_from_file(self):
        pass

    def load_from_model(self, gtep_model):
        if type(gtep_model) is not ExpansionPlanningModel:
            logger.warning(
                f"Solutions must be loaded from ExpansionPlanningModel objects, not %s"
                % type(gtep_model)
            )
            raise ValueError
        if gtep_model.results is None:
            raise ValueError(
                "ExpansionPlanningSolution objects loaded from model must have a results component."
            )
        self.results = gtep_model.results  # Highs results object
        self.stages = gtep_model.stages  # int
        self.formulation = gtep_model.formulation  # None (???)
        self.data = gtep_model.data  # ModelData object
        self.num_reps = gtep_model.num_reps  # int
        self.len_reps = gtep_model.len_reps  # int
        self.num_commit = gtep_model.num_commit  # int
        self.num_dispatch = gtep_model.num_dispatch  # int

        self.expressions = {
            expr.name: value(expr)
            for expr in gtep_model.model.component_data_objects(Expression)
            if ("Commitment" in expr.name) or ("Investment" in expr.name)
        }

        self.variables = {
            expr.name: value(expr)
            for expr in gtep_model.model.component_data_objects(Expression)
            if ("Commitment" in expr.name) or ("Investment" in expr.name)
        }

    def import_data_object(self, data_obj):
        self.data = data_obj.md

    def read_json(self, filepath):
        # read a json file and recover a solution primals
        json_filepath = Path(filepath)
        with open(json_filepath, "r") as fobj:
            json_read = json.loads(fobj.read())
        self.primals_tree = json_read["results"]["primals_tree"]

    def dump_json(self, filename="./gtep_solution_jscTest.json"):
        # def dump_json(self, filename="./gtep_solution.json"):

        dump_filepath = Path(filename)
        with open(dump_filepath, "w") as fobj:
            json.dump(self._to_dict(), fobj)

    def _to_dict(self) -> dict:

        results_dict = {
            "solution_loader": self.results.solution_loader,  # object
            "termination_condition": self.results.termination_condition,  # object
            "best_feasible_objective": self.results.best_feasible_objective,
            "best_objective_bound": self.results.best_objective_bound,
            "wallclock_time": self.results.wallclock_time,
            "expressions": self.expressions,
        }

        # "best_feasible_objective", "best_objective_bound", and "wallclock_time" are all numbers, dont need subhandlers

        # subhandle "termination_condition"
        results_dict["termination_condition"] = {
            "value": self.results.termination_condition.value,
            "name": self.results.termination_condition.name,
        }

        # subhandle "solution_loader"
        results_dict["solution_loader"] = {"primals": {}}

        for key, val in self.results.solution_loader.get_primals()._dict.items():
            tmp_key = key

            # handle binary vars by delving one layer in
            results_dict["solution_loader"]["primals"][tmp_key] = {
                "name": val[0].name,
                "value": val[0].value,
                "bounds": val[0].bounds,
            }

            # handle binary
            if val[0].is_binary():
                results_dict["solution_loader"]["primals"][tmp_key]["is_binary"] = val[
                    0
                ].is_binary()
            # handle units, sometimes they dont have anything
            if val[0].get_units() is not None:
                results_dict["solution_loader"]["primals"][tmp_key]["units"] = (
                    val[0].get_units().name
                )
            else:
                results_dict["solution_loader"]["primals"][tmp_key]["units"] = val[
                    0
                ].get_units()

        # renest "termination_condition" as a json-friendly dictionary
        # things are either vars (which have some sort of signifier in [] brackets) or are an attribute, which dont
        # the name variable will give it away
        results_dict["primals_tree"] = {}
        results_dict["expressions_tree"] = {}

        for key, val in self.results.solution_loader.get_primals()._dict.items():
            # split the name to figure out depth
            split_name = val[0].name.split(".")

            # start at the bottom and nest accordingly
            tmp_dict = {
                "name": val[0].name,
                "value": val[0].value,
                "bounds": val[0].bounds,
            }

            # handle binary
            if val[0].is_binary():
                tmp_dict["is_binary"] = val[0].is_binary()

            # handle units, sometimes they dont have anything
            if val[0].get_units() is not None:
                tmp_dict["units"] = val[0].get_units().name
            else:
                tmp_dict["units"] = val[0].get_units()

            # allocate the nested dictionary
            def nested_set(this_dict, key, val):
                if len(key) > 1:
                    # check if it's a binary var and pull up one layer
                    if key[1] == "binary_indicator_var":
                        this_dict[key[0]] = val
                    else:
                        this_dict.setdefault(key[0], {})
                        nested_set(this_dict[key[0]], key[1:], val)
                else:
                    this_dict[key[0]] = val

            nested_set(results_dict["primals_tree"], split_name, tmp_dict)

        pass

        for key, val in self.expressions.items():
            # split the name to figure out depth
            split_name = key.split(".")

            # start at the bottom and nest accordingly
            tmp_dict = {
                "value": val,
            }

            # allocate the nested dictionary
            def nested_set(this_dict, key, val):
                if len(key) > 1:
                    # check if it's a binary var and pull up one layer
                    this_dict.setdefault(key[0], {})
                    nested_set(this_dict[key[0]], key[1:], val)
                else:
                    this_dict[key[0]] = val

            nested_set(results_dict["expressions_tree"], split_name, tmp_dict)
        # split out expressions
        self.expressions_tree = results_dict["expressions_tree"]

        # mint the final dictionary to save
        out_dict = {"data": self.data.data, "results": results_dict}

        self.primals_tree = results_dict["primals_tree"]

        return out_dict

    def discover_level_relationships(self, dispatch_level_dict):
        list_of_keys = list(dispatch_level_dict.keys())

        relationships_dict = {}
        # go through each key and split them into categories and names
        # each name should have a handful of categories, which should be the same across a group of names
        for this_key in list_of_keys:
            # check if it has a bracketed relationship, and if it does go ahead, otherwise skip
            try:
                primal_category = this_key.split("[")[0]
                primal_name = this_key.split("[")[1].split("]")[0]
                relationships_dict.setdefault(primal_name, set())
                relationships_dict[primal_name].add(primal_category)

            except IndexError as iEx:
                print(
                    f'[WARNING] discover_level_relationships has encountered an error: Attempted to split out {this_key}, failed with error: "{iEx}". Assigning as axuilary.'
                )

        # convert sets to frozensets to be hashable
        for this_key in relationships_dict:
            relationships_dict[this_key] = frozenset(relationships_dict[this_key])

        # now go through each primal name and check for the groups who match
        matching_groups_dict = {}
        for this_primal_name, this_primal_set in relationships_dict.items():
            matching_groups_dict.setdefault(this_primal_set, set())
            matching_groups_dict[this_primal_set].add(this_primal_name)

        return matching_groups_dict

    def _level_relationship_dict_to_df_workhorse(
        self, level_key, timeseries_dict, keys_of_interest, vars_of_interest
    ):
        df_data_dict = {}
        units_dict = {}
        # set our defaults
        df_data_dict.setdefault(level_key, [])
        for this_koi in keys_of_interest:
            for this_voi in vars_of_interest:
                df_data_dict.setdefault(f"{this_koi}_{this_voi}_value", [])
                df_data_dict.setdefault(f"{this_koi}_{this_voi}_lower_bound", [])
                df_data_dict.setdefault(f"{this_koi}_{this_voi}_upper_bound", [])

        # dump data into dict to read into df
        for period_dict in timeseries_dict:
            df_data_dict[level_key].append(period_dict["period_number"])
            for this_koi in keys_of_interest:
                for this_voi in vars_of_interest:
                    # check if this is a variable by checking if it has a "value"
                    if "value" in period_dict["primals_by_name"][this_koi][this_voi]:
                        # if its an integer, cast it as a boolean for now
                        if (
                            "is_binary"
                            in period_dict["primals_by_name"][this_koi][this_voi]
                        ):
                            if period_dict["primals_by_name"][this_koi][this_voi][
                                "is_binary"
                            ]:
                                df_data_dict[f"{this_koi}_{this_voi}_value"].append(
                                    bool(
                                        round(
                                            period_dict["primals_by_name"][this_koi][
                                                this_voi
                                            ]["value"]
                                        )
                                    )  # have to cast to int because there are floating point errors
                                )
                                units_dict.setdefault(
                                    f"{this_koi}_{this_voi}_value",
                                    period_dict["primals_by_name"][this_koi][this_voi][
                                        "units"
                                    ],
                                )
                            else:
                                df_data_dict[f"{this_koi}_{this_voi}_value"].append(
                                    period_dict["primals_by_name"][this_koi][this_voi][
                                        "value"
                                    ]
                                )
                                units_dict.setdefault(
                                    f"{this_koi}_{this_voi}_value",
                                    period_dict["primals_by_name"][this_koi][this_voi][
                                        "units"
                                    ],
                                )
                        else:
                            df_data_dict[f"{this_koi}_{this_voi}_value"].append(
                                period_dict["primals_by_name"][this_koi][this_voi][
                                    "value"
                                ]
                            )
                            units_dict.setdefault(
                                f"{this_koi}_{this_voi}_value",
                                period_dict["primals_by_name"][this_koi][this_voi][
                                    "units"
                                ],
                            )
                        df_data_dict[f"{this_koi}_{this_voi}_lower_bound"].append(
                            period_dict["primals_by_name"][this_koi][this_voi][
                                "bounds"
                            ][0]
                        )
                        units_dict.setdefault(
                            f"{this_koi}_{this_voi}_value",
                            period_dict["primals_by_name"][this_koi][this_voi]["units"],
                        )
                        df_data_dict[f"{this_koi}_{this_voi}_upper_bound"].append(
                            period_dict["primals_by_name"][this_koi][this_voi][
                                "bounds"
                            ][1]
                        )
                        units_dict.setdefault(
                            f"{this_koi}_{this_voi}_value",
                            period_dict["primals_by_name"][this_koi][this_voi]["units"],
                        )

        # try to make a DF, and if not just pass back an empty
        try:
            data_df = pd.DataFrame(df_data_dict)
            # fix any Nones and make them NaNs
            data_df = data_df.fillna(value=np.nan)
            return data_df, units_dict
        except ValueError as vEx:
            print(
                f"[WARNING] _level_relationship_dict_to_df_workhorse attempted to create dataframe and failed: {vEx}"
            )
            return pd.DataFrame(), {}

    def _plot_workhorse_relational(
        self,
        level_key,
        df,
        keys,
        vars,
        parent_key_string,
        pretty_title="Selected Data",
        plot_bounds=False,
        save_dir=".",
        aspect_ratio=1,
    ):

        # figure out how big the plot needs to be
        gridspec_height = 2 * max(len(keys), len(vars))
        gridspec_width = 2
        fig_width_padding = 0
        fig_height_padding = 0
        max_figheight = 48
        total_periods = len(df[level_key])
        key_gridspec_div = floor(
            gridspec_height / len(keys)
        )  # number of gridspec heights a key plot can be
        var_gridspec_div = floor(
            gridspec_height / len(vars)
        )  # number of gridspec heights a var plot can be

        # to make things look nice, we dont want height or width to be more than twice the other
        fig_width = (total_periods * gridspec_width * 4) + fig_width_padding
        fig_width = min(max_figheight, fig_width)
        fig_height = (2 * gridspec_height) + fig_height_padding
        if fig_width / fig_height > aspect_ratio:
            fig_height = floor(fig_width / aspect_ratio)
        elif fig_height / fig_width > aspect_ratio:
            fig_width = floor(fig_height / aspect_ratio)

        # set up plot
        fig = plt.figure(
            figsize=(fig_width, fig_height), tight_layout=False
        )  # (32, 16) works will for 4 plots tall and about 6 periods wide per plot
        gs = fig.add_gridspec(gridspec_height, gridspec_width)
        # plot out the keys of interest
        ax_koi_list = []
        for ix_koi, this_koi in enumerate(keys):
            ax_koi = fig.add_subplot(
                gs[(ix_koi * key_gridspec_div) : ((ix_koi + 1) * key_gridspec_div), 0]
            )
            ax_koi_list.append(ax_koi)

            for iy, this_voi in enumerate(vars):
                ax_koi.plot(
                    df[level_key],
                    df[f"{this_koi}_{this_voi}_value"],
                    label=f"{this_koi}_{this_voi}",
                    marker="o",
                )
                if plot_bounds:
                    ax_koi.fill_between(
                        df[level_key],
                        df[f"{this_koi}_{this_voi}_lower_bound"],
                        df[f"{this_koi}_{this_voi}_upper_bound"],
                        alpha=0.1,
                    )

            ax_koi.set_ylabel("Value $[n]$")
            ax_koi.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax_koi.legend()

        # label axes
        ax_koi_list[-1].set_xlabel(f"{level_key} $[n]$")
        ax_koi_list[0].set_title(f"{pretty_title} by Type")

        # plot variables of interest
        ax_voi_list = []
        # plot generations and curtailmentsagainst each outher
        for ix_voi, this_voi in enumerate(vars):
            ax_voi = fig.add_subplot(
                gs[(ix_voi * var_gridspec_div) : ((ix_voi + 1) * var_gridspec_div), 1]
            )
            ax_voi_list.append(ax_voi)
            for this_koi in keys:
                ax_voi.plot(
                    df[level_key],
                    df[f"{this_koi}_{this_voi}_value"],
                    label=f"{this_koi}_{this_voi}",
                    marker="o",
                )
                if plot_bounds:
                    ax_voi.fill_between(
                        df[level_key],
                        df[f"{this_koi}_{this_voi}_lower_bound"],
                        df[f"{this_koi}_{this_voi}_upper_bound"],
                        alpha=0.1,
                    )

            ax_voi.set_ylabel("Value $[n]$")
            ax_voi.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax_voi.legend()

        # label axes
        ax_voi_list[-1].set_xlabel(f"{level_key} $[n]$")
        ax_voi_list[0].set_title(f"{pretty_title} by Category")

        fig.align_labels()
        fig.suptitle(f"{parent_key_string}")
        fig.savefig(
            f"{save_dir}{parent_key_string}_{pretty_title.replace(' ', '_')}.png"
        )
        plt.close()

    def _plot_workhose_binaries(
        self,
        level_key,
        df,
        keys,
        vars,
        parent_key_string,
        pretty_title="Selected Data",
        save_dir=".",
    ):

        fig = plt.figure(figsize=(32, 16), tight_layout=False)
        gs = fig.add_gridspec(1, 1)  # only need 1 plot for now
        # if all the variables are binaries, we can assume that the vars are all binaries and the keys are all categories
        total_height = len(vars)
        interstate_height = 1.0 / (len(keys) + 2)
        width = 1
        width_padding = 0.05
        ax_bins = fig.add_subplot(gs[:, :])
        ax_bins.set_ylim([-0.5, total_height - 0.5])  # set ylims to support bools
        ax_bins.set_xlim([0.5, len(df[level_key]) + 0.5])  # set xlims to support bools
        ax_bins.set_yticklabels([None] + list(vars))
        ax_bins.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax_bins.xaxis.set_major_locator(MaxNLocator(integer=True))
        for axline_ix in range(total_height):
            ax_bins.axhline(
                axline_ix + 0.5,
                color="grey",
                linewidth=3,
            )  # draw a seperator line between each level
        for axline_ix in range(len(df[level_key])):
            ax_bins.axvline(
                axline_ix + 0.5,
                color="grey",
                linewidth=3,
                linestyle="dotted",
                alpha=0.5,
            )  # draw a seperator line between each level

        for ix_key, this_koi in enumerate(keys):
            # make a dummy line to steal the color cycler and make a single item for the legend
            (line,) = ax_bins.plot(
                [None],
                [None],
                label=f"{this_koi}",
                linewidth=5,
            )
            for ix_var, this_voi in enumerate(vars):
                for tx, is_it_on in zip(
                    df[level_key], df[f"{this_koi}_{this_voi}_value"]
                ):
                    if is_it_on:
                        tmp_rect = plt.Rectangle(
                            [
                                tx - 0.5 + width_padding,
                                ((ix_var) + (interstate_height * (ix_key + 1))) - 0.5,
                            ],
                            width - (width_padding * 2),
                            interstate_height,
                            alpha=0.9,
                            edgecolor="black",
                            color=line.get_color(),
                        )
                        ax_bins.add_patch(tmp_rect)

        ax_bins.set_xlabel(f"{level_key} $[n]$")
        ax_bins.set_title("State Variable Time History")
        ax_bins.set_ylabel("Binary State")
        ax_bins.legend()

        fig.align_labels()
        fig.suptitle(f"{parent_key_string}")
        fig.savefig(
            f"{save_dir}{parent_key_string}_{pretty_title.replace(' ', '_')}.png"
        )
        plt.close()

    def _level_relationship_df_to_plot(
        self,
        level_key,
        df,
        keys,
        vars,
        parent_key_string,
        pretty_title="Selected Data",
        plot_bounds=False,
        save_dir=".",
        config={},
    ):

        # [HACK] hard coding the generator state order, to be fixed later
        config["order_gen_state"] = ["genOff", "genShutdown", "genStartup", "genOn"]
        config["order_gen_invest_state"] = [
            "genDisabled",
            "genRetired",
            "genExtended",
            "genInstalled",
            "genOperational",
        ]
        config["order_branch_invest_state"] = [
            "branchDisabled",
            "branchRetired",
            "branchExtended",
            "branchInstalled",
            "branchOperational",
        ]

        # check if ALL the possible things to look at are binaries
        all_binaries = True
        for ix, this_voi in enumerate(vars):
            for iy, this_koi in enumerate(keys):
                if not (df[f"{this_koi}_{this_voi}_value"].dtype == "bool"):
                    all_binaries = False
                    break
        if all_binaries:

            # check the config to see if we have any overrides
            if "order_gen_state" in config:
                # check that everything can be mapped over
                matched_config_override = True
                for item in vars:
                    if not item in config["order_gen_state"]:
                        matched_config_override = False
                        break
                if matched_config_override:
                    vars = config["order_gen_state"]
            if "order_gen_invest_state" in config:
                # check that everything can be mapped over
                matched_config_override = True
                for item in vars:
                    if not item in config["order_gen_invest_state"]:
                        matched_config_override = False
                        break
                if matched_config_override:
                    vars = config["order_gen_invest_state"]
            if "order_branch_invest_state" in config:
                # check that everything can be mapped over
                matched_config_override = True
                for item in vars:
                    if not item in config["order_branch_invest_state"]:
                        matched_config_override = False
                        break
                if matched_config_override:
                    vars = config["order_branch_invest_state"]

            self._plot_workhose_binaries(
                level_key,
                df,
                keys,
                vars,
                parent_key_string,
                pretty_title,
                save_dir,
            )

        else:
            self._plot_workhorse_relational(
                level_key,
                df,
                keys,
                vars,
                parent_key_string,
                pretty_title,
                plot_bounds,
                save_dir,
            )

    def _expressions_plot_workhorse(
        self,
        level_key,
        upper_level_dict,
        parent_key_string,
        save_dir="./",
        plot_bounds=False,
    ):
        # go through a commitment period and parse out the dispatch periods
        # slice out all keys pertaining to dispatchPeriod
        level_period_keys = [
            this_key for this_key in upper_level_dict.keys() if (level_key in this_key)
        ]

        # scan level period for keys that have values associated
        keys_of_vals_of_interest = []
        for this_key in upper_level_dict[level_period_keys[0]]:
            if "value" in upper_level_dict[level_period_keys[0]][this_key]:
                keys_of_vals_of_interest.append(this_key)

        print(keys_of_vals_of_interest)

        # if we found things that have values, make a dictionary and then plot
        # why do we need to do it by level first and then pivot the dict? Because we don't actually know what the "periods" numbers are, they might be arbitrary and in an arbitrary order
        if len(keys_of_vals_of_interest) > 0:
            # check all the periods and sort them out
            vals_dict = {}

            for this_key in upper_level_dict:
                level_period_number = int(
                    re.split("\[|\]", this_key.split(level_key)[1])[1]
                )
                vals_dict.setdefault(level_period_number, {})
                for this_val_key in keys_of_vals_of_interest:
                    vals_dict[level_period_number][this_val_key] = upper_level_dict[
                        this_key
                    ][this_val_key]["value"]

            print(vals_dict)

            # now pivot the dictionary to make a dataframe
            # make a dictionary where the keys are the top layer
            df_dict = {key: [] for key in keys_of_vals_of_interest}
            sorted_vals_period = sorted(vals_dict)
            df_dict["period_number"] = sorted_vals_period
            for this_val_period in sorted_vals_period:
                for this_key in keys_of_vals_of_interest:
                    df_dict[this_key].append(vals_dict[this_val_period][this_key])

            print(df_dict)
            expression_level_df = pd.DataFrame(df_dict)
            print(expression_level_df)

            # plot the DF
            # figure out how big the plot needs to be
            gridspec_height = 2 * len(keys_of_vals_of_interest)
            gridspec_width = 2
            fig_width_padding = 0
            fig_height_padding = 0
            max_figheight = 48
            total_periods = len(expression_level_df)
            key_gridspec_div = floor(
                gridspec_height / len(keys_of_vals_of_interest)
            )  # number of gridspec heights a key plot can be

            # to make things look nice, we dont want height or width to be more than twice the other
            fig_width = (total_periods * gridspec_width * 4) + fig_width_padding
            fig_width = min(max_figheight, fig_width)
            fig_height = (2 * gridspec_height) + fig_height_padding
            if fig_width / fig_height > 2:
                fig_height = floor(fig_width / 2)
            elif fig_height / fig_width > 2:
                fig_width = floor(fig_height / 2)

            # set up plot
            fig = plt.figure(
                figsize=(fig_width, fig_height), tight_layout=False
            )  # (32, 16) works will for 4 plots tall and about 6 periods wide per plot
            gs = fig.add_gridspec(gridspec_height, gridspec_width)
            # plot out the keys of interest

            pretty_title = "Expression"

            ax_koi_list = []
            for ix_koi, this_koi in enumerate(keys_of_vals_of_interest):
                ax_koi = fig.add_subplot(
                    gs[
                        (ix_koi * key_gridspec_div) : ((ix_koi + 1) * key_gridspec_div),
                        :,
                    ]
                )
                ax_koi_list.append(ax_koi)

                ax_koi.plot(
                    expression_level_df["period_number"],
                    expression_level_df[this_koi],
                    label=f"{this_koi}",
                    marker="o",
                )
                ax_koi.set_ylabel("Value $[n]$")
                ax_koi.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax_koi.legend()

            # label axes
            ax_koi_list[-1].set_xlabel(f"{level_key} $[n]$")
            ax_koi_list[0].set_title(f"{pretty_title} by Type")

            # JSC update - " ", "_" to ' ', '_' for compilation. Not sure if this is due to a version diff or what
            fig.align_labels()
            fig.suptitle(f"{parent_key_string}")
            fig.savefig(
                f"{save_dir}{parent_key_string}_{pretty_title.replace(' ', '_')}.png"
            )
            plt.close()

    def _level_plot_workhorse(
        self,
        level_key,
        upper_level_dict,
        parent_key_string,
        save_dir="./",
        plot_bounds=False,
    ):
        # go through a commitment period and parse out the dispatch periods
        level_timeseries = []
        # slice out all keys pertaining to dispatchPeriod
        level_period_keys = [
            this_key for this_key in upper_level_dict.keys() if (level_key in this_key)
        ]
        # aux_var_dict = {}
        for this_key in level_period_keys:
            level_period_dict = {}
            # cut out which dispatch period this is
            level_period_number = int(
                re.split("\[|\]", this_key.split(level_key)[1])[1]
            )
            # print(level_period_number)

            level_period_dict["period_number"] = level_period_number

            # pivot the dictionary to get the primals into categories
            primals_by_category = {}
            primals_by_name = {}

            # split key on brackets to get title
            for this_primal in upper_level_dict[this_key]:
                # check if it has a bracketed relationship, and if it does go ahead, otherwise skip
                tmp_save_primal = upper_level_dict[this_key][this_primal]
                try:

                    # check if this_primal is splitable on brackets

                    # based on the name of the primal, split out the "category" and the "name"
                    # Example: commitmentPeriod[1]
                    # - category: "commitmentPeriod"
                    # -     name: "1"
                    primal_category = this_primal.split("[")[0]
                    primal_name = this_primal.split("[")[1].split("]")[0]

                    # create one view that shares the categories, and one the shares the names
                    primals_by_category.setdefault(primal_category, {})
                    primals_by_name.setdefault(primal_name, {})

                    primals_by_category[primal_category][primal_name] = tmp_save_primal
                    primals_by_name[primal_name][primal_category] = tmp_save_primal

                except IndexError as iEx:
                    print(
                        f"[WARNING] _level_plot_workhorse has encountered an error: Attempted to split out {this_primal} from {this_key}, failed with error {iEx}. Skipping."
                    )
                    # aux_var_dict.setdefault(this_primal, [])
                    # aux_var_dict[this_primal].append(this_key)
            level_period_dict["primals_by_category"] = primals_by_category
            level_period_dict["primals_by_name"] = primals_by_name
            level_timeseries.append(level_period_dict)

        # # clean up the aux vars based on what we found
        # for aux_var_primal_name, aux_var_category_name in aux_var_dict.items():
        #     print(aux_var_primal_name, aux_var_category_name)

        # sort by the dispatch period number
        level_timeseries = sorted(level_timeseries, key=lambda x: x["period_number"])

        # discover the relationships at the dispatch level
        # ASSUMES that all the dispatch levels have the exact same underlying variables and relationships
        level_relationships = self.discover_level_relationships(
            upper_level_dict[level_period_keys[0]]
        )
        # print("LEVEL RELATIONSHIPS DEBUG")
        # print(level_relationships)
        # print("END LEVEL RELATIONSHIPS DEBUG")

        # plot relationships
        for vars_of_interest, keys_of_interest in level_relationships.items():

            # sort the vars and keys for consistency
            tmp_voi = sorted(vars_of_interest)
            tmp_koi = sorted(keys_of_interest)

            # make a df for debug and also easy tabularness for plots
            this_df_of_interest, this_df_units = (
                self._level_relationship_dict_to_df_workhorse(
                    level_key, level_timeseries, tmp_koi, tmp_voi
                )
            )

            # check if we got anything in the df
            if not this_df_of_interest.empty:

                # just ram the variable names together, that'll be fine, right?
                this_pretty_title = ", ".join(tmp_voi)

                # [HACK]
                # if we find powerflow, plot it as a network
                if "powerFlow" in tmp_voi:
                    self._plot_graph_workhorse(
                        this_df_of_interest,
                        "powerFlow",
                        parent_key_string,
                        units=this_df_units,
                        pretty_title=this_pretty_title,
                        save_dir=save_dir,
                    )

                    # [HACK] put this back one indent level when done
                    ## tab this back and forth to do the things

                    # plot it
                self._level_relationship_df_to_plot(
                    level_key,
                    this_df_of_interest,
                    tmp_koi,
                    tmp_voi,
                    parent_key_string,
                    pretty_title=this_pretty_title,
                    save_dir=save_dir,
                    plot_bounds=plot_bounds,
                )

    # JSC update - 'dc_branch' to 'branch'
    def _plot_graph_workhorse(
        self,
        df,
        value_key,
        parent_key_string,
        what_is_a_bus_called="branch",  #'dc_branch',
        units=None,
        pretty_title="Selected Data",
        save_dir=".",
    ):
        # testing networkx plots

        # preslice out data of interest
        cols_of_interest = [col for col in df.columns if f"{value_key}_value" in col]
        df_of_interest = df[cols_of_interest]
        df_max = df_of_interest.to_numpy().max()
        df_min = df_of_interest.to_numpy().min()
        # assume all the units are the same, and pull the first one
        units_str = ""
        if units[cols_of_interest[0]] is not None:
            units_str = f" [{units[cols_of_interest[0]]}]"

        # construct graph object
        G = nx.Graph()
        labels = {}
        # add nodes
        for item in self.data.data["elements"]["bus"]:
            G.add_node(item)
            labels[item] = item

        # do edges manually later

        # set up plot
        fig = plt.figure(
            figsize=(16, 8), tight_layout=False
        )  # (32, 16) works will for 4 plots tall and about 6 periods wide per plot
        ax_graph = fig.add_subplot()
        ax_graph.grid(False)

        # # add edges
        # for item in self.data.data['elements']['branch']:
        #     G.add_edge(self.data.data['elements']['branch'][item]['from_bus'], self.data.data['elements']['branch'][item]['to_bus'])

        # G = nx.path_graph(5)
        graph_node_position_dict = nx.kamada_kawai_layout(G)
        # graph_node_position_dict = nx.planar_layout(G)
        # graph_node_position_dict = nx.spectral_layout(G)
        nx.drawing.draw_networkx_nodes(
            G, graph_node_position_dict, node_size=1000, ax=ax_graph
        )
        nx.draw_networkx_labels(
            G,
            graph_node_position_dict,
            labels,
            font_size=18,
            font_color="whitesmoke",
            ax=ax_graph,
        )

        def draw_single_edge_flow(
            item,
            glyph_values_slice,
            ax_graph,
            cmap=cm.rainbow,
            norm=Normalize(vmin=None, vmax=None),
            glyph_type="custom",
        ):

            def generate_flow_glyphs(
                num_glyphs,
                spacing=0.05,
                glyph_type="triangle",
                glyph_rotation=0.0,
                verts=3,
            ):

                flow_glyphs = []
                for this_block_ix in range(num_glyphs):
                    # normalizing this patch to 1

                    ####
                    # rectangle version
                    ####
                    if glyph_type == "rectangle":
                        # anchor for rectangles are set to bottom left
                        glyph_anchor_coord = [this_block_ix / float(num_glyphs), -0.5]
                        # height is y, width is x
                        consistent_width = 1.0 / float(num_glyphs)
                        # apply scaling
                        x_nudge = consistent_width * (spacing)
                        # nudge the start forward a bit (by the nudge factor)
                        glyph_anchor_coord[0] += x_nudge
                        patch_width = consistent_width - x_nudge
                        patch_height = 1
                        flow_glyphs.append(
                            Rectangle(glyph_anchor_coord, patch_width, patch_height)
                        )

                    ####
                    # triangle version
                    ####
                    if glyph_type == "triangle":
                        # triangles need to be in the center and then given a size
                        glyph_anchor_coord = [
                            (this_block_ix + 0.5) / float(num_glyphs),
                            0,
                        ]
                        glyph_verts = 3
                        glyph_radius = (1.0 / float(num_glyphs)) / 2.0
                        # apply nudges
                        glyph_radius *= 1 - (spacing / 2.0)
                        flow_glyphs.append(
                            RegularPolygon(
                                glyph_anchor_coord,
                                glyph_verts,
                                radius=glyph_radius,
                                orientation=glyph_rotation,
                            )
                        )

                        yscale_transform = Affine2D().scale(sx=1, sy=0.5 / glyph_radius)
                        # rescale y to make it fit in a 1x1 box
                        flow_glyphs[-1].set_transform(yscale_transform)

                    ####
                    # n-gon version
                    ####
                    if glyph_type == "n-gon":
                        # triangles need to be in the center and then given a size
                        glyph_anchor_coord = [
                            (this_block_ix + 0.5) / float(num_glyphs),
                            0,
                        ]
                        glyph_verts = verts
                        glyph_radius = (1.0 / float(num_glyphs)) / 2.0
                        # apply nudges
                        glyph_radius *= 1 - (spacing)
                        flow_glyphs.append(
                            RegularPolygon(
                                glyph_anchor_coord,
                                glyph_verts,
                                radius=glyph_radius,
                                orientation=glyph_rotation,
                            )
                        )

                        yscale_transform = Affine2D().scale(sx=1, sy=0.5 / glyph_radius)
                        # rescale y to make it fit in a 1x1 box
                        flow_glyphs[-1].set_transform(yscale_transform)

                    ####
                    # custom_flow
                    ####
                    if glyph_type == "custom":
                        # anchor for rectangles are set to bottom left
                        glyph_anchor_coord = [this_block_ix / float(num_glyphs), -0.5]
                        # height is y, width is x
                        consistent_width = 1.0 / float(num_glyphs)
                        # apply scaling
                        x_nudge = consistent_width * (spacing)
                        patch_width = consistent_width - x_nudge
                        patch_height = 1
                        codes, verts = zip(
                            *[
                                (mpath.Path.MOVETO, glyph_anchor_coord),
                                (
                                    mpath.Path.LINETO,
                                    [
                                        glyph_anchor_coord[0],
                                        glyph_anchor_coord[1] + patch_height,
                                    ],
                                ),
                                (
                                    mpath.Path.LINETO,
                                    [
                                        glyph_anchor_coord[0] + patch_width * 0.7,
                                        glyph_anchor_coord[1] + patch_height,
                                    ],
                                ),  # go 70% of the width along the top
                                (
                                    mpath.Path.LINETO,
                                    [
                                        glyph_anchor_coord[0] + patch_width,
                                        glyph_anchor_coord[1] + patch_height * 0.5,
                                    ],
                                ),  # go the rest of the width and meet in the center
                                (
                                    mpath.Path.LINETO,
                                    [
                                        glyph_anchor_coord[0] + patch_width * 0.7,
                                        glyph_anchor_coord[1],
                                    ],
                                ),  # go back a bit and to the bottom to finish the wedge
                                (mpath.Path.LINETO, glyph_anchor_coord),
                            ]
                        )  # go to home

                        flow_glyphs.append(
                            PathPatch(mpath.Path(verts, codes), ec="none"),
                        )

                        rotation_transofrm = Affine2D().rotate_around(
                            glyph_anchor_coord[0] + patch_width * 0.5,
                            glyph_anchor_coord[1] + patch_height * 0.5,
                            glyph_rotation,
                        )
                        # rescale y to make it fit in a 1x1 box
                        flow_glyphs[-1].set_transform(rotation_transofrm)

                return flow_glyphs

            # make some blocks
            # weights_top = (np.random.randn(4)+1)/2.
            # weights_bot = (np.random.randn(4)+1)/2.
            # weights_top = np.array(range(num_blocks))/(num_blocks*2)
            # weights_bot = (np.array(range(num_blocks))+num_blocks)/(num_blocks*2)
            weights_top = glyph_values_slice
            weights_bot = glyph_values_slice

            top_flow_glyphs = generate_flow_glyphs(
                len(weights_top), glyph_type=glyph_type
            )
            top_facecolors = cmap(norm(weights_top))
            top_flow_collection = PatchCollection(
                top_flow_glyphs, facecolors=top_facecolors, edgecolors="grey", alpha=0.5
            )
            # bot_flow_glyphs = generate_flow_glyphs(len(weights_bot), glyph_type=glyph_type, glyph_rotation=(np.pi/2.)) # for squares
            bot_flow_glyphs = generate_flow_glyphs(
                len(weights_bot), glyph_type=glyph_type, glyph_rotation=(np.pi)
            )  # for custom
            bot_flow_glyphs = reversed(bot_flow_glyphs)
            bot_facecolors = cmap(norm(weights_bot))
            # bot_flow_collection = PatchCollection(bot_flow_glyphs, facecolors=bot_facecolors, edgecolors='grey', alpha=0.5) # [HACK]

            # scale and move top and bottom collections
            # top_base_transform = Affine2D().scale(sx=1, sy=0.9) + Affine2D().translate(0, 0.5) #+ ax_graph.transData  # [HACK]
            top_base_transform = Affine2D().scale(sx=1, sy=1.0) + Affine2D().translate(
                0, 0.0
            )  # + ax_graph.transData
            top_flow_collection.set_transform(top_base_transform)
            bot_base_transform = Affine2D().scale(sx=1, sy=0.9) + Affine2D().translate(
                0, -0.5
            )  # + ax_graph.transData
            # bot_base_transform = Affine2D().scale(sx=1, sy=0.9) + Affine2D().translate(0, -0.5) + ax_graph.transData
            # bot_flow_collection.set_transform(bot_base_transform) # [HACK]

            # combine collections and move to edge between nodes

            # attempt to rotate
            start_key = self.data.data["elements"][what_is_a_bus_called][item][
                "from_bus"
            ]
            end_key = self.data.data["elements"][what_is_a_bus_called][item]["to_bus"]
            start_pos = graph_node_position_dict[start_key]
            end_pos = graph_node_position_dict[end_key]
            node_distance = np.linalg.norm(end_pos - start_pos)
            rot_angle_rad = np.arctan2(
                (end_pos[1] - start_pos[1]), (end_pos[0] - start_pos[0])
            )

            along_edge_scale = 0.4
            away_from_edge_scale = 0.1
            # set up transformations
            # stretch to the distance between target nodes
            length_transform = Affine2D().scale(
                sx=node_distance * along_edge_scale, sy=1
            )
            # squish
            scale_transform = Affine2D().scale(sx=1, sy=away_from_edge_scale)
            # rotate
            rot_transform = Affine2D().rotate_deg(np.rad2deg(rot_angle_rad))
            # translate to the node start, then push it along the edge until it's apprximately centered and scaled nicely
            translate_transform = Affine2D().translate(
                start_pos[0]
                + (
                    np.cos(rot_angle_rad) * node_distance * 0.5 * (1 - along_edge_scale)
                ),
                start_pos[1]
                + (
                    np.sin(rot_angle_rad) * node_distance * 0.5 * (1 - along_edge_scale)
                ),
            )
            t2 = (
                length_transform
                + scale_transform
                + rot_transform
                + translate_transform
                + ax_graph.transData
            )

            top_flow_collection.set_transform(top_flow_collection.get_transform() + t2)
            # bot_flow_collection.set_transform(bot_flow_collection.get_transform() + t2) # [HACK]

            # add collection
            ax_graph.add_collection(top_flow_collection)
            # ax_graph.add_collection(bot_flow_collection)  # [HACK]

        # add edges
        # define edge colorbar
        cmap = cm.rainbow
        normalize = Normalize(vmin=df_min, vmax=df_max)
        cmappable = cm.ScalarMappable(norm=normalize, cmap=cmap)

        for item in self.data.data["elements"][what_is_a_bus_called]:

            # grab the keys we care about
            start_key = self.data.data["elements"][what_is_a_bus_called][item][
                "from_bus"
            ]
            end_key = self.data.data["elements"][what_is_a_bus_called][item]["to_bus"]
            start_pos = graph_node_position_dict[start_key]
            end_pos = graph_node_position_dict[end_key]
            # edge_key = f"branch_{start_key}_{end_key}_{value_key}_value"
            # alt_edge_key = f"branch_{end_key}_{start_key}_{value_key}_value"
            edge_key = f"{item}_{value_key}_value"
            alt_edge_key = f"{item}_{value_key}_value"

            # @KyleSkolfield is there a reason to not do this?
            branch_name_edge_key = item + "_powerFlow_value"

            # kind = 'triangle'
            # kind = 'rectangle'
            kind = "custom"
            glyph_values_slice = df[branch_name_edge_key].values

            try:
                glyph_values_slice = df[edge_key].values
            except KeyError as kex:
                print(f"Attempted to slice DF in network using {edge_key}, failed.")

            draw_single_edge_flow(
                item,
                glyph_values_slice,
                ax_graph,
                cmap=cmap,
                norm=normalize,
                glyph_type=kind,
            )

            # forward arrow
            ax_graph.arrow(
                start_pos[0],
                start_pos[1],
                (end_pos[0] - start_pos[0]),
                (end_pos[1] - start_pos[1]),
                color="black",
            )
            # backward arrow
            ax_graph.arrow(
                end_pos[0],
                end_pos[1],
                (start_pos[0] - end_pos[0]),
                (start_pos[1] - end_pos[1]),
                color="black",
            )

        # insert colorbar
        fig.colorbar(cmappable, ax=ax_graph, label=f"{value_key}{units_str}")
        # make some titles
        fig.suptitle(f"{parent_key_string}_{value_key}")

        # JSC update - " ", "_" to ' ', '_' for compilation. Not sure if this is due to a version diff or what
        # save
        fig.savefig(
            f"{save_dir}{parent_key_string}_{pretty_title.replace(' ', '_')}_graph.png"
        )
        pass

    def plot_levels(self, save_dir="."):

        # self._plot_graph_workhorse()

        # plot or represent primals trees
        for this_root_level_key in self.primals_tree:
            if "investmentStage" in this_root_level_key:
                # run the toplevel keys
                parent_key_string = f"{this_root_level_key}"
                self._level_plot_workhorse(
                    "investmentStage", self.primals_tree, this_root_level_key, save_dir
                )

                # run the representative period subkeys
                investment_level_cut = self.primals_tree[this_root_level_key]
                parent_key_string = f"{this_root_level_key}"
                self._level_plot_workhorse(
                    "representativePeriod",
                    investment_level_cut,
                    parent_key_string,
                    save_dir,
                )

                for this_inv_level_key in self.primals_tree[this_root_level_key].keys():
                    if "representativePeriod" in this_inv_level_key:
                        representative_level_cut = self.primals_tree[
                            this_root_level_key
                        ][this_inv_level_key]
                        parent_key_string = (
                            f"{this_root_level_key}_{this_inv_level_key}"
                        )
                        self._level_plot_workhorse(
                            "commitmentPeriod",
                            representative_level_cut,
                            parent_key_string,
                            save_dir,
                        )

                        for this_rep_level_key in self.primals_tree[
                            this_root_level_key
                        ][this_inv_level_key].keys():
                            if "commitmentPeriod" in this_rep_level_key:
                                commitment_level_cut = self.primals_tree[
                                    this_root_level_key
                                ][this_inv_level_key][this_rep_level_key]

                                parent_key_string = f"{this_root_level_key}_{this_inv_level_key}_{this_rep_level_key}"

                                self._level_plot_workhorse(
                                    "dispatchPeriod",
                                    commitment_level_cut,
                                    parent_key_string,
                                    save_dir,
                                )

        # # plot or represent expressions
        # self._expressions_plot_workhorse(
        #     "investmentStage", self.expressions_tree, 'investmentStage', save_dir
        # )
        # for this_root_level_key in self.expressions_tree:
        #     if "investmentStage" in this_root_level_key:
        #         investment_level_cut = self.expressions_tree[this_root_level_key]
        #         parent_key_string = f"{this_root_level_key}"
        #         self._expressions_plot_workhorse(
        #             "representativePeriod", investment_level_cut, parent_key_string, save_dir
        #         )

        #         for this_inv_level_key in self.expressions_tree[this_root_level_key].keys():
        #             if "representativePeriod" in this_inv_level_key:
        #                 representative_level_cut = self.expressions_tree[this_root_level_key][this_inv_level_key]
        #                 parent_key_string = f"{this_root_level_key}_{this_inv_level_key}"
        #                 self._expressions_plot_workhorse(
        #                     "commitmentPeriod", representative_level_cut, parent_key_string, save_dir
        #                 )

        #                 for this_rep_level_key in self.expressions_tree[
        #                     this_root_level_key
        #                 ][this_inv_level_key].keys():
        #                     if "commitmentPeriod" in this_rep_level_key:
        #                         commitment_level_cut = self.expressions_tree[
        #                             this_root_level_key
        #                         ][this_inv_level_key][this_rep_level_key]

        #                         parent_key_string = f"{this_root_level_key}_{this_inv_level_key}_{this_rep_level_key}"

        #                         self._expressions_plot_workhorse(
        #                             "dispatchPeriod",commitment_level_cut, parent_key_string, save_dir
        #                         )
        # pass
