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

import pyomo.environ as pyo
import pyomo.gdp as gdp
from pyomo.environ import units as u
from pyomo.core.base.param import IndexedParam
from pyomo.core.base.expression import ScalarExpression, IndexedExpression

logger = logging.getLogger(__name__)


class ExpansionPlanningSolution:
    """This class stores the solution to the ExpansionPlanningModel
    class for writing and visualization.

    """

    def __init__(self, data_path):
        self.gen_df = pd.read_csv(f"{data_path}/gen.csv")
        self.gen_types = {
            gen_type: self.gen_df[self.gen_df["Unit Type"] == gen_type]["PMax MW"].sum()
            for gen_type in set(self.gen_df["Unit Type"])
        }

    def create_results_directory(self, dir_name="results"):
        """This function creates a directory to save model results.

        :param dir_name: Name or path of the directory where results
                         will be saved. Defaults to "results".
        :return: Directory name/path.

        """
        os.makedirs(dir_name, exist_ok=True)
        print(f"\nCreating the directory '{dir_name}' to save the results. ")

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

        print(f">>>Saved model configuration to: {config_csv_path}")

    def save_results_in_json_files(self, gtep_model, dir_name):
        """This functions saves model results to JSON files.

        Outputs include investments, load shed, costs, flows,
        generation, curtailment, loads, reserves, and storage
        charge/discharge. Creates the results directory if needed.

        :param gtep_model: Solved expansion planning model object.
        :param dir_name: Directory where JSON files are written.

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

        print(
            f"The following files have been created in the directory '{folder_name}':"
        )
        for name in output_files:
            print(f" - {folder_name}/{name}.json")
