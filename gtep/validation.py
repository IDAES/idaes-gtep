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

import re
from pathlib import Path
import shutil
from collections import defaultdict
from numbers import Number

from gtep.gtep_solution import ExpansionPlanningSolution

import pandas as pd


def safe_extract_variable_index(variable_name: str) -> str:
    """
    Takes a variable name and extracts the index (enclosed in square brackets).
    If no matches are found, returns the whole variable name.

    :param variable_name:   Variable name.
    :type variable_name:    str
    :returns:               Index, not including the square brackets.
    """
    # TODO: possibly modify once solution object has sets
    search_result = re.search(r"\[.*\]", variable_name)
    return search_result.group(0)[1:-1] if search_result else variable_name


def extract_primals_last_investment_stage(
    sol_object: ExpansionPlanningSolution,
) -> dict:
    """
    Accesses the primal variables for the last investment stage in `sol_object`.

    :param sol_object:  Solution object.
    :type sol_object:   gtep.gtep_solution.ExpansionPlanningSolution
    :returns:           Dictionary of the form {var_name: var_data}, where var_data is a dict
    """
    # TODO: modify once solution object has sets
    sol_dict = sol_object._to_dict()["results"]["primals_tree"]
    end_investment_stage = list(sol_dict.keys())[0]  # more robust way to do this?
    return sol_dict[end_investment_stage]


def extract_variable_values(
    primals_dict: dict,
    variable_type: str,
    element_statuses: list[str] = ["Extended", "Operational", "Installed"],
) -> dict:
    """
    Collects the name and value for all variables of the given `variable_type`
    which describe statuses in `element_statuses`.

    :param primals_dict:        Primal variables, in the form {var_name: var_value}.
    :type primals_dict:         dict
    :param variable_type:       Type of variable to extract (e.g., `"gen"`).
    :type variable_type:        str
    :param element_statuses:    Variable statuses to check for; should be substrings of
                                    variable names. Defaults to
                                    `["Extended", "Operational", "Installed"]`
    :type element_statuses:     list[str], optional
    :returns:                   Dictionary of the form {var_name: var_value}
    """
    # TODO: modify once solution object has sets
    end_investment_values_dict = {
        name: var["value"]
        for name, var in primals_dict.items()
        if variable_type in name
        and any([status in name for status in element_statuses])
    }

    return end_investment_values_dict


def sum_variable_values_by_index(value_dict: dict[str, Number]) -> dict[str, Number]:
    """
    Takes an input dict of variables and values, groups by index, and computes the sum of values.
    This function expects the variable values to be non-bool Numbers (e.g., int, float) and
    raises an error if a different type is passed in via `value_dict`.

    :param value_dict:  Dictionary of the form {var[idx]: value}.
    :type value_dict:   dict[str, Number]
    :returns:           Dictionary of the form {idx: summed_value}.
    """
    # Raise an error if we're trying to sum something that isn't a number (incl. bool)
    value_types = [type(value) for value in value_dict.values()]
    if any([t == bool or not issubclass(t, Number) for t in value_types]):
        raise ValueError(
            f"Expected variable values to be numeric, but got the following types: {set(value_types)}"
        )

    result = defaultdict(float)
    for name, value in value_dict.items():
        index = safe_extract_variable_index(name)
        result[index] += value

    return dict(result)


def safe_mkdir(path: Path):
    """
    Creates a directory if it doesn't already exist. If the path exists but isn't a directory,
    raises a `FileExistsError`.

    :param path:    Path to new directory
    :type path:     pathlib.Path
    """
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"{path} exists and is not a directory")
    else:
        path.mkdir(parents=True)


def safe_write_dataframe_to_csv(
    dataframe: pd.DataFrame, directory: Path, filename: str
):
    """
    Writes a DataFrame to CSV, creating the directory if necessary.

    :param dataframe:   DataFrame to write.
    :type dataframe:    pandas.DataFrame
    :param directory:   Directory to write in.
    :type directory:    pathlib.Path
    :param filename:    Name of file to write to.
    :type filename:     str
    """
    safe_mkdir(directory)
    dataframe.to_csv((directory / filename).resolve(), index=False)


def populate_generators(
    data_input_path: Path, sol_object: ExpansionPlanningSolution, data_output_path: str
):
    """
    Takes a set of input generator data and updates it to reflect the solution object.

    :param data_input_path:     Path to folder with prescient generator data. Must contain a file named `"gen.csv"`.
    :param sol_object:          Solution object run with the prescient data at `data_input_path`
    :param data_output_path:    Path to write the generator data to
    :type data_input_path:      pathlib.Path
    :type sol_object:           gtep.gtep_solution.ExpansionPlanningSolution
    :type data_output_path:     str
    """

    # load existing and candidate generators from initial prescient data
    # note that -c in name indicates candidate
    input_df = pd.read_csv((data_input_path / "gen.csv").resolve())

    # get the sum by index of extended, operational, and installed variables for thermal gens during last investment period
    primals = extract_primals_last_investment_stage(sol_object)
    end_gen_idxs = sum_variable_values_by_index(extract_variable_values(primals, "gen"))
    end_gen_idx_list = [
        idx for idx, val in end_gen_idxs.items() if val > 0.5
    ]  # keep only gens which are "active"

    # get the sum by index of extended, operational, and installed variables for renewables during last investment period
    end_renew_idxs = sum_variable_values_by_index(
        extract_variable_values(primals, "renewable")
    )
    end_renew_idx_list = list(end_renew_idxs.keys())

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
    data_input_path: Path, sol_object: ExpansionPlanningSolution, data_output_path: Path
):
    """
    Takes a set of input transmission data and updates it to reflect the solution object.

    :param data_input_path:     Path to folder with prescient transmission data. Must contain a file named `"branch.csv"`.
    :param sol_object:          Solution object run with the prescient data at `data_input_path`
    :param data_output_path:    Path to write the transmission data to
    :type data_input_path:      pathlib.Path
    :type sol_object:           gtep.gtep_solution.ExpansionPlanningSolution
    :type data_output_path:     pathlib.Path
    """

    # load existing and candidate generators from initial prescient data
    # note that -c in name indicates candidate
    input_df = pd.read_csv((data_input_path / "branch.csv").resolve())

    # get the sum by index of extended, operational, and installed variables for branches during last investment period
    primals = extract_primals_last_investment_stage(sol_object)
    end_branch_idxs = sum_variable_values_by_index(
        extract_variable_values(primals, "branch")
    )
    end_branch_idx_list = [
        idx for idx, val in end_branch_idxs.items() if val > 0.5
    ]  # keep only branches which are active

    # define output dataframe
    output_df = input_df[input_df["UID"].isin(end_branch_idx_list)]

    safe_write_dataframe_to_csv(output_df, data_output_path, "branch.csv")


def filter_pointers(data_input_path: Path, data_output_path: Path):
    """
    Takes a set of input timeseries pointers and updates it to reflect the solution object.
    Must be run _after_ `populate_generators` and with the same `data_output_path`.

    :param data_input_path:     Path to folder with prescient timeseries pointers. Must contain a file named `"timeseries_pointers.csv"`.
    :param sol_object:          Solution object run with the prescient data at `data_input_path`
    :param data_output_path:    Path to write the timeseries pointers to. Must contain a file named `"gen.csv"`.
    :type data_input_path:      pathlib.Path
    :type sol_object:           gtep.gtep_solution.ExpansionPlanningSolution
    :type data_output_path:     pathlib.Path
    """

    # load initial timeseries pointers
    input_pointers_df = pd.read_csv(
        (data_input_path / "timeseries_pointers.csv").resolve()
    )

    # load final generators
    output_generators_df = pd.read_csv((data_output_path / "gen.csv").resolve())

    # keep generators that exist at the final investment stage and remove the rest
    # keep all non-generator timeseries pointers
    matching_gen_list = [gen for gen in output_generators_df["GEN UID"]]
    output_df = input_pointers_df[
        input_pointers_df["Object"].isin(matching_gen_list)
        | input_pointers_df["Category"]
        != "Generator"
    ]

    safe_write_dataframe_to_csv(output_df, data_output_path, "timeseries_pointers.csv")


def copy_prescient_inputs(data_input_path: Path, data_output_path: Path):
    """
    Copies all files at `data_input_path` to `data_output_path`, except for:
        - gen.csv
        - branch.csv
        - timeseries_pointers.csv

    These files are instead handled by other functions in this module (namely,
    `filter_pointers`, `populate_generators`, and `populate_transmission`).

    :param data_input_path:     Path to folder with files to copy.
    :param data_output_path:    Path to write the files to.
    :type data_input_path:      pathlib.Path
    :type data_output_path:     pathlib.Path
    """

    safe_mkdir(data_output_path)

    file_list = list(data_input_path.iterdir())

    # @jkskolf, I don't think I like this ...
    for from_fpath in file_list:
        if from_fpath.name in ["gen.csv", "branch.csv", "timeseries_pointers.csv"]:
            continue
        to_fpath = (data_output_path / from_fpath.name).resolve()
        shutil.copy(from_fpath, to_fpath)
