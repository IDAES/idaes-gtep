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

from pyomo.environ import *
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_solution import ExpansionPlanningSolution
import re
import os
import shutil
import logging
from collections import defaultdict

import pandas as pd

logger = logging.getLogger(__name__)


def safe_extract_variable_index(variable_name: str) -> str:
    """
    Takes a variable name and extracts the index (enclosed in square brackets).
    If no matches are found, returns the whole variable name.

    :param variable_name: Variable name.
    :type variable_name: str
    :returns: Index.
    """
    search_result = re.search(r"\[.*\]", variable_name)
    return search_result.group(0)[1:-1] if search_result else variable_name


def extract_end_variable_values(
    sol_object: ExpansionPlanningSolution,
    variable_type: str,
    element_statuses: list[str] = ["Extended", "Operational", "Installed"],
) -> dict:
    """
    Accesses the primal variables for the last investment stage in `sol_object`
    and collects the name and value for all variables of the given `variable_type`
    which describe statuses in `element_statuses`.

    :param sol_object: Solution object.
    :type sol_object: gtep.gtep_solution.ExpansionPlanningSolution
    :param variable_type: Type of variable to extract (e.g., `"gen"`).
    :type variable_type: str
    :param element_status: Variable statuses to check for; should be substrings of
        variable names. Defaults to `["Extended", "Operational", "Installed"]`
    :type element_statuses: list[str], optional
    :returns: Dictionary of the form {var_name: var_value}
    """
    ALLOWED_VARIABLE_TYPES = ["gen", "renewable", "branch"]
    if variable_type not in ALLOWED_VARIABLE_TYPES:
        raise ValueError(
            f"variable_type argument must be one of {ALLOWED_VARIABLE_TYPES}."
        )

    sol_dict = sol_object._to_dict()["results"]["primals_tree"]
    end_investment_stage = list(sol_dict.keys())[0]  # more robust way to do this?

    end_investment_values_dict = {
        name: var["value"]
        for name, var in sol_dict[end_investment_stage].items()
        if variable_type in name
        and any([status in name for status in element_statuses])
    }

    return end_investment_values_dict


def sum_variable_values_by_index(value_dict: dict[str, float]) -> dict[str, float]:
    """
    Takes an input dict of variables and values, groups by index, and computes the sum of values.

    :param value_dict: Dictionary of the form {var[idx]: value}.
    :type value_dict: dict[str, float]
    :returns: Dictionary of the form {idx: summed_value}.
    """

    result = defaultdict(float)
    for name, value in value_dict.items():
        index = safe_extract_variable_index(name)
        result[index] += value

    return dict(result)


def safe_write_dataframe_to_csv(dataframe: pd.DataFrame, directory: str, filename: str):
    """
    Writes a DataFrame to CSV, creating the directory if necessary.

    :param dataframe: DataFrame to write.
    :type dataframe: pandas.DataFrame
    :param directory: Directory to write in.
    :type directory: str
    :param filename: Name of file to write to.
    :type filename: str
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    dataframe.to_csv(os.path.join(directory, filename), index=False)


def populate_generators(
    data_input_path: str, sol_object: ExpansionPlanningSolution, data_output_path: str
):
    """
    Takes a set of input generator data and updates it to reflect the solution object.

    :param data_input_path:     Path to folder with prescient generator data. Must contain a file named `"gen.csv"`.
    :param sol_object:          Solution object run with the prescient data at `data_input_path`
    :param data_output_path:    Path to write the generator data to
    :type data_input_path:      str
    :type sol_object:           gtep.gtep_solution.ExpansionPlanningSolution
    :type data_output_path:     str
    """

    # load existing and candidate generators from initial prescient data
    # note that -c in name indicates candidate
    input_df = pd.read_csv(os.path.join(data_input_path, "gen.csv"))

    # get the sum by index of extended, operational, and installed variables for thermal gens during last investment period
    end_gen_idxs = sum_variable_values_by_index(
        extract_end_variable_values(sol_object, "gen")
    )
    end_gen_idx_list = [
        idx for idx, val in end_gen_idxs.items() if val > 0.5
    ]  # keep only gens which are "active"

    # get the sum by index of extended, operational, and installed variables for renewables during last investment period
    end_renew_idxs = sum_variable_values_by_index(
        extract_end_variable_values(sol_object, "renewable")
    )
    end_renew_idx_list = list(end_renew_idxs.keys())  # why no similar check here?

    # update renewables with values from last investment period
    input_df["PMax MW"] = (
        input_df["GEN UID"].map(end_renew_idxs).fillna(input_df["PMax MW"])
    )

    # define output dataframe
    output_df = input_df[
        input_df["GEN UID"].isin(end_gen_idx_list + end_renew_idx_list)
    ]

    # TODO: (@jkskolf) should we update prices here? I think no, but ...
    safe_write_dataframe_to_csv(output_df, data_output_path, "gen.csv")


def populate_transmission(
    data_input_path: str, sol_object: ExpansionPlanningSolution, data_output_path: str
):
    """
    Takes a set of input transmission data and updates it to reflect the solution object.

    :param data_input_path:     Path to folder with prescient transmission data. Must contain a file named `"branch.csv"`.
    :param sol_object:          Solution object run with the prescient data at `data_input_path`
    :param data_output_path:    Path to write the transmission data to
    :type data_input_path:      str
    :type sol_object:           gtep.gtep_solution.ExpansionPlanningSolution
    :type data_output_path:     str
    """

    # load existing and candidate generators from initial prescient data
    # note that -c in name indicates candidate
    input_df = pd.read_csv(os.path.join(data_input_path, "branch.csv"))

    # get the sum by index of extended, operational, and installed variables for branches during last investment period
    end_branch_idxs = sum_variable_values_by_index(
        extract_end_variable_values(sol_object, "branch")
    )
    end_branch_idx_list = [
        idx for idx, val in end_branch_idxs.items() if val > 0.5
    ]  # keep only branches which are active

    # define output dataframe
    output_df = input_df[input_df["UID"].isin(end_branch_idx_list)]

    safe_write_dataframe_to_csv(output_df, data_output_path, "branch.csv")


def filter_pointers(data_input_path: str, data_output_path: str):
    """
    Takes a set of input timeseries pointers and updates it to reflect the solution object.
    Must be run _after_ `populate_generators` and with the same `data_output_path`.

    :param data_input_path:     Path to folder with prescient timeseries pointers. Must contain a file named `"timeseries_pointers.csv"`.
    :param sol_object:          Solution object run with the prescient data at `data_input_path`
    :param data_output_path:    Path to write the timeseries pointers to. Must contain a file named `"gen.csv"`.
    :type data_input_path:      str
    :type sol_object:           gtep.gtep_solution.ExpansionPlanningSolution
    :type data_output_path:     str
    """

    # load initial timeseries pointers
    input_pointers_df = pd.read_csv(
        os.path.join(data_input_path, "timeseries_pointers.csv")
    )

    # load final generators
    output_generators_df = pd.read_csv(os.path.join(data_output_path, "gen.csv"))

    # keep generators that exist at the final investment stage and remove the rest
    # keep all non-generator timeseries pointers
    matching_gen_list = [gen for gen in output_generators_df["GEN UID"]]
    output_df = input_pointers_df[
        input_pointers_df["Object"].isin(matching_gen_list)
        | input_pointers_df["Category"]
        != "Generator"
    ]

    safe_write_dataframe_to_csv(output_df, data_output_path, "timeseries_pointers.csv")


def copy_prescient_inputs(data_input_path: str, data_output_path: str):
    """
    Copies all files at `data_input_path` to `data_output_path`, except for:
        - timeseries_pointers.csv
        - gen.csv
        - branch.csv

    These files are instead handled by other functions in this module (namely,
    `filter_pointers`, `populate_generators`, and `populate_transmission`).

    :param data_input_path:     Path to folder with files to copy.
    :param data_output_path:    Path to write the files to.
    :type data_input_path:      str
    :type data_output_path:     str
    """

    if not os.path.exists(data_output_path):
        os.makedirs(data_output_path)

    file_list = os.listdir(data_input_path)
    file_list.remove("timeseries_pointers.csv")
    file_list.remove("gen.csv")
    file_list.remove("branch.csv")

    # @jkskolf, I don't think I like this ...
    for fname in file_list:
        from_file = os.path.join(data_input_path, fname)
        to_file = os.path.join(data_output_path, fname)
        shutil.copy(from_file, to_file)
