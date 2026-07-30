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

from gtep.gtep_model import ExpansionPlanningModel

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

    def load_from_model(self, gtep_model):
        """This function loads the results from the solved model
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
        self.len_reps = gtep_model.len_reps
        self.num_commit = gtep_model.num_commit
        self.num_dispatch = gtep_model.num_dispatch

        # Store selected expression values for validation/reporting.
        self.expressions = {
            expr.name: pyo.value(expr)
            for expr in gtep_model.model.component_data_objects(pyo.Expression)
            if ("Commitment" in expr.name) or ("Investment" in expr.name)
        }

    def _to_dict(self) -> dict:
        """This function converts solution results into a
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

    def save_results_in_json_files(self, gtep_model, dir_name, value_threshold=1e-3):
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

        print(
            f"The following files have been created in the directory '{folder_name}':"
        )
        for name in output_files:
            print(f" - {folder_name}/{name}.json")
