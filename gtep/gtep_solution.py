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


import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import re

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
        # add in

    def read_json(self, filepath):
        # read a json file and recover a solution primals
        json_filepath = Path(filepath)
        with open(json_filepath, "r") as fobj:
            json_read = json.loads(fobj.read())
        self.primals_tree = json_read["results"]["primals_tree"]

    def dump_json(self, filename="./gtep_solution.json"):

        dump_filepath = Path(filename)
        with open(dump_filepath, "w") as fobj:
            json.dump(self._to_dict(), fobj)

    def _to_dict(self):

        results_dict = {
            "solution_loader": self.results.solution_loader,  # object
            "termination_condition": self.results.termination_condition,  # object
            "best_feasible_objective": self.results.best_feasible_objective,
            "best_objective_bound": self.results.best_objective_bound,
            "wallclock_time": self.results.wallclock_time,
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
        out_dict = {"data": self.data.data, "results": results_dict}

        self.primals_tree = results_dict["primals_tree"]

        return out_dict

    def plot_investment_level(self):

        # can grab bus names from
        pass

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
                    f'[WARNING] discover_level_relationships has encountered an error: Attempted to split out {this_key}, failed with error: "{iEx}". Skipping.'
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

    def _level_dict_to_df_workhorse(
        self, level_key, timeseries_dict, keys_of_interest, vars_of_interest
    ):
        df_data_dict = {}
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
                            else:
                                df_data_dict[f"{this_koi}_{this_voi}_value"].append(
                                    period_dict["primals_by_name"][this_koi][this_voi][
                                        "value"
                                    ]
                                )
                        else:
                            df_data_dict[f"{this_koi}_{this_voi}_value"].append(
                                period_dict["primals_by_name"][this_koi][this_voi][
                                    "value"
                                ]
                            )
                        df_data_dict[f"{this_koi}_{this_voi}_lower_bound"].append(
                            period_dict["primals_by_name"][this_koi][this_voi][
                                "bounds"
                            ][0]
                        )
                        df_data_dict[f"{this_koi}_{this_voi}_upper_bound"].append(
                            period_dict["primals_by_name"][this_koi][this_voi][
                                "bounds"
                            ][1]
                        )

        # try to make a DF, and if not just pass back an empty
        try:
            data_df = pd.DataFrame(df_data_dict)
            # fix any Nones and make them NaNs
            data_df = data_df.fillna(value=np.nan)
            return data_df
        except ValueError as vEx:
            print(
                f"[WARNING] _level_dict_to_df_workhorse attempted to create dataframe and failed: {vEx}"
            )
            return pd.DataFrame()

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
        if fig_width / fig_height > 2:
            fig_height = floor(fig_width / 2)
        elif fig_height / fig_height > 2:
            fig_height = floor(fig_height / 2)

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

    def _level_df_to_plot(
        self,
        level_key,
        df,
        keys,
        vars,
        parent_key_string,
        pretty_title="Selected Data",
        plot_bounds=False,
        save_dir=".",
    ):

        # check if ALL the possible things to look at are binaries
        all_binaries = True
        for ix, this_voi in enumerate(vars):
            for iy, this_koi in enumerate(keys):
                if not (df[f"{this_koi}_{this_voi}_value"].dtype == "bool"):
                    all_binaries = False
                    break
        if all_binaries:

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
                try:
                    # print(this_primal)
                    primal_category = this_primal.split("[")[0]
                    primal_name = this_primal.split("[")[1].split("]")[0]

                    primals_by_category.setdefault(primal_category, {})
                    primals_by_name.setdefault(primal_name, {})

                    primals_by_category[primal_category][primal_name] = (
                        upper_level_dict[this_key][this_primal]
                    )
                    primals_by_name[primal_name][primal_category] = upper_level_dict[
                        this_key
                    ][this_primal]
                except IndexError as iEx:
                    print(
                        f"[WARNING] _level_plot_workhorse has encountered an error: Attempted to split out {this_key}, failed with error {iEx}. Skipping."
                    )
            level_period_dict["primals_by_category"] = primals_by_category
            level_period_dict["primals_by_name"] = primals_by_name
            level_timeseries.append(level_period_dict)

        # sort by the dispatch period number
        level_timeseries = sorted(level_timeseries, key=lambda x: x["period_number"])

        # discover the relationships at the dispatch level
        # ASSUMES that all the dispatch levels have the exact same underlying variables and relationships
        level_relationships = self.discover_level_relationships(
            upper_level_dict[level_period_keys[0]]
        )

        for vars_of_interest, keys_of_interest in level_relationships.items():

            # sort the vars and keys for consistency
            tmp_voi = sorted(vars_of_interest)
            tmp_koi = sorted(keys_of_interest)

            # make a df for debug and also easy tabularness for plots
            this_df_of_interest = self._level_dict_to_df_workhorse(
                level_key, level_timeseries, tmp_koi, tmp_voi
            )

            # check if we got anything in the df
            if not this_df_of_interest.empty:

                # just ram the variable names together, that'll be fine, right?
                this_pretty_title = ", ".join(tmp_voi)

                # plot it
                self._level_df_to_plot(
                    level_key,
                    this_df_of_interest,
                    tmp_koi,
                    tmp_voi,
                    parent_key_string,
                    pretty_title=this_pretty_title,
                    save_dir=save_dir,
                    plot_bounds=plot_bounds,
                )

    def plot_levels(self, save_dir="."):

        # plot or represent some dispatch periods or something I don't know
        for this_root_level_key in self.primals_tree:
            if "investmentStage" in this_root_level_key:
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

        pass
        # things to cut on
        # renewableGeneration
        # renewableCurtailment
        # powerFlow
        # busAngle
        # thermalGeneration
        # spinningReserve
        # quickstartReserve

        # [
        #     "renewableGeneration[10_PV]",
        #     "renewableCurtailment[10_PV]",
        #     "renewableGeneration[2_RTPV]",
        #     "renewableCurtailment[2_RTPV]",
        #     "renewableGeneration[1_HYDRO]",
        #     "renewableCurtailment[1_HYDRO]",
        #     "renewableGeneration[4_WIND]",
        #     "renewableCurtailment[4_WIND]",
        #     "powerFlow[branch_2_3]",
        #     "busAngle[bus3]",
        #     "busAngle[bus2]",
        #     "powerFlow[branch_1_10]",
        #     "busAngle[bus10]",
        #     "busAngle[bus1]",
        #     "powerFlow[branch_3_4_1]",
        #     "busAngle[bus4]",
        #     "powerFlow[branch_4_10]",
        #     "powerFlow[branch_1_4]",
        #     "powerFlow[branch_1_2]",
        #     "powerFlow[branch_3_4_0]",
        #     "loadShed[bus1]",
        #     "thermalGeneration[4_CC]",
        #     "thermalGeneration[4_STEAM]",
        #     "loadShed[bus4]",
        #     "thermalGeneration[10_STEAM]",
        #     "loadShed[bus10]",
        #     "loadShed[bus2]",
        #     "thermalGeneration[3_CT]",
        #     "loadShed[bus3]",
        #     "quickstartReserve[3_CT]",
        #     "quickstartReserve[10_STEAM]",
        #     "quickstartReserve[4_CC]",
        #     "quickstartReserve[4_STEAM]",
        #     "spinningReserve[3_CT]",
        #     "spinningReserve[10_STEAM]",
        #     "spinningReserve[4_CC]",
        #     "spinningReserve[4_STEAM]",
        # ]
