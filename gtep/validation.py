from pyomo.environ import *
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_solution import ExpansionPlanningSolution
import re
import os
import shutil
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def populate_generators(
    data_input_path: str, sol_object: ExpansionPlanningSolution, data_output_path: str
):
    # load existing and candidate generators from initial prescient data
    # note that -c in name indicates candidate
    input_df = pd.read_csv(data_input_path + "/gen.csv")

    # pull final stage solution variables for thermal and renewable investments
    # for thermal:
    # generator should exist in future grid if the status in the final stage
    # is installed, operational, or extended

    def gen_name_filter(gen_name):
        return "gen" in gen_name and (
            "Ext" in gen_name or "Ope" in gen_name or "Ins" in gen_name
        )

    solution_dict = sol_object._to_dict()["results"]["primals_tree"]
    end_investment_stage = list(solution_dict.keys())[0]
    end_investment_solution_dict = {
        k: v["value"]
        for k, v in solution_dict[end_investment_stage].items()
        if gen_name_filter(k) and v["value"] > 0.5
    }
    end_investment_gens = [
        re.search(r"\[.*\]", k).group(0)[1:-1]
        for k in end_investment_solution_dict.keys()
    ]

    # for renewable:
    # total capacity should be installed + operational + extended values
    def renewable_name_filter(gen_name):
        return "renew" in gen_name and (
            "Ext" in gen_name or "Ope" in gen_name or "Ins" in gen_name
        )

    end_investment_renewable_dict = {
        k: v["value"]
        for k, v in solution_dict[end_investment_stage].items()
        if renewable_name_filter(k)
    }
    end_investment_renewable_gens = {
        re.search(r"\[.*\]", k).group(0)[1:-1]: 0
        for k in end_investment_renewable_dict.keys()
    }
    for k, v in end_investment_renewable_dict.items():
        end_investment_renewable_gens[re.search(r"\[.*\]", k).group(0)[1:-1]] += v
    for k, v in end_investment_renewable_gens.items():
        ## NOTE: (@jkskolf) this will break in pandas 3.0
        input_df["PMax MW"].mask(input_df["GEN UID"] == k, v, inplace=True)

    end_investment_gens += [k for k in end_investment_renewable_gens.keys()]
    # populate output dataframe
    output_df = input_df[input_df["GEN UID"].isin(end_investment_gens)]

    # TODO: (@jkskolf) should we update prices here? I think no, but ...
    if not os.path.exists(data_output_path):
        os.makedirs(data_output_path)
    output_df.to_csv(data_output_path + "/gen.csv", index=False)


def populate_transmission(data_input_path, sol_object, data_output_path):
    # load existing and candidate generators from initial prescient data
    # note that -c in name indicates candidate
    input_df = pd.read_csv(data_input_path + "/branch.csv")

    # pull final stage solution variables for transmission
    def branch_name_filter(gen_name):
        return "bran" in gen_name and (
            "Ext" in gen_name or "Ope" in gen_name or "Ins" in gen_name
        )

    solution_dict = sol_object._to_dict()["results"]["primals_tree"]
    end_investment_stage = list(solution_dict.keys())[0]
    end_investment_solution_dict = {
        k: v["value"]
        for k, v in solution_dict[end_investment_stage].items()
        if branch_name_filter(k) and v["value"] > 0.5
    }
    end_investment_branches = [
        re.search(r"\[.*\]", k).group(0)[1:-1]
        for k in end_investment_solution_dict.keys()
    ]
    output_df = input_df[input_df["UID"].isin(end_investment_branches)]

    if not os.path.exists(data_output_path):
        os.makedirs(data_output_path)
    output_df.to_csv(data_output_path + "/branch.csv", index=False)


def filter_pointers(data_input_path, data_output_path):
    # load initial timeseries pointers
    input_pointers_df = pd.read_csv(data_input_path + "/timeseries_pointers.csv")

    # load final generators
    # NOTE: must be run _after_ populate_generators and with the same data_output_path
    # to pull resulting generator objects
    output_generators_df = pd.read_csv(data_output_path + "/gen.csv")

    # keep generators that exist at the final investment stage and remove the rest
    # keep all non-generator timeseries pointers
    matching_gen_list = [gen for gen in output_generators_df["GEN UID"]]
    output_df = input_pointers_df[
        input_pointers_df["Object"].isin(matching_gen_list)
        | input_pointers_df["Category"]
        != "Generator"
    ]

    if not os.path.exists(data_output_path):
        os.makedirs(data_output_path)
    output_df.to_csv(data_output_path + "/timeseries_pointers.csv")


def clone_timeseries(data_input_path, data_output_path):
    file_list = os.listdir(data_input_path)
    file_list.remove("timeseries_pointers.csv")
    file_list.remove("gen.csv")
    file_list.remove("branch.csv")

    # @jkskolf, I don't think I like this ...
    for fname in file_list:
        shutil.copy(data_input_path + "/" + fname, data_output_path + "/" + fname)
