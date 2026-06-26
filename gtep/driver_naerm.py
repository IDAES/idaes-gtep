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

import pyomo.environ as pyo
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

import xpress

import os


def _sanitize_path_part(name):
    return (
        str(name)
        .replace("[", "_")
        .replace("]", "")
        .replace(",", "_")
        .replace(" ", "_")
    )


def _block_path_components(block, root_block=None):
    """
    Convert a block name into hierarchical path components.

    Example:
        model.plant.unit[3].heater
    becomes:
        ["plant", "unit_3", "heater"]

    If root_block is provided and the block name starts with root_block.name,
    that prefix is removed.
    """
    full_name = block.name if block.name is not None else "anonymous_block"
    parts = full_name.split(".")

    if root_block is not None and root_block.name is not None:
        root_parts = root_block.name.split(".")
        if parts[:len(root_parts)] == root_parts:
            parts = parts[len(root_parts):]

    if not parts:
        parts = ["root"]

    return [_sanitize_path_part(p) for p in parts]


def export_block_constraints_hierarchical(
    root_block,
    output_folder="constraint_dump",
    individual_entries=False,
    recurse=True,
    active=True,
):
    """
    Export constraints from Pyomo blocks into a folder hierarchy matching
    the full block path.

    Parameters
    ----------
    root_block : pyo.Block
        Root Pyomo block/model.
    output_folder : str
        Top-level output folder.
    individual_entries : bool
        If False, use pprint() for each Constraint component.
        If True, write each indexed constraint entry separately.
    recurse : bool
        If True, include all sub-blocks.
    active : bool
        If True, only include active blocks/constraints.
    """
    os.makedirs(output_folder, exist_ok=True)

    if recurse:
        blocks = list(root_block.block_data_objects(active=active, descend_into=True))
    else:
        blocks = [root_block]

    for block in blocks:
        rel_parts = _block_path_components(block, root_block=root_block)
        block_dir = os.path.join(output_folder, *rel_parts)
        os.makedirs(block_dir, exist_ok=True)

        filename = os.path.join(block_dir, "constraints.txt")

        with open(filename, "w") as f:
            f.write(f"Block: {block.name}\n")
            f.write("=" * 80 + "\n\n")

            found_any = False

            for constr in block.component_objects(
                pyo.Constraint, active=active, descend_into=False
            ):
                found_any = True

                if not individual_entries:
                    constr.pprint(ostream=f)
                    f.write("\n")
                else:
                    f.write(f"Constraint component: {constr.name}\n")
                    f.write("-" * 80 + "\n")
                    if constr.is_indexed():
                        for index in constr:
                            f.write(f"{constr.name}[{index}] : {constr[index].expr}\n")
                    else:
                        f.write(f"{constr.name} : {constr.expr}\n")
                    f.write("\n")

            if not found_any:
                f.write("No constraints found on this block.\n")

# Add data
rep_days = [
    "2034-01-22 00:00",
    "2034-01-24 00:00",
    "2034-03-04 00:00",
    "2034-03-11 00:00",
    "2034-03-18 00:00",
    "2034-05-11 00:00",
    "2034-06-09 00:00",
    "2034-06-16 00:00",
    "2034-06-30 00:00",
    "2034-07-12 00:00",
    "2034-10-01 00:00",
    "2034-10-08 00:00",
    "2034-10-29 00:00",
    "2034-11-10 00:00",
    "2034-12-06 00:00",
]
rep_weights = [27, 32, 32, 37, 21, 29, 13, 25, 21, 21, 23, 26, 17, 23, 18]

data_date = "6-24-2026"
data_path = f"./gtep/data/WECC_ADS_PNNL_{data_date}"
data_object = ExpansionPlanningData(
    stages=1,
    num_reps=15,
    num_commit=24,
    num_dispatch=1,
    duration_representative_period=24,
    save_period_structure_file=False,
    period_structure_json_file=None,
    # period_structure_json_file="period_structure_from_gtep.json",
)

data_object.load_prescient(
    data_path, representative_dates=rep_days, representative_weights=rep_weights
)

# [ESR WIP: Add bus and cost data files to be used on the
# DataProcessing class. This class processes data for the following
# type of generators: Natural Gas, Solar, Wind, and Coal. Note that
# during this data processing stage, the generator type is matched to
# existent generators in the data. The data contains the following
# types: (a) Natural Gas: Combustion Turbine (CT) and Fuel Efficiency
# (FE) and (b) Solar: Utility PV and Concentrated Solar Power (CSP)

bus_data_path = "./data/costs/Bus_data_gen_weights_mappings.csv"
cost_data_path = "./data/costs/2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
candidate_gens = [
    "Natural Gas_CT",
    "Natural Gas_FE",
    "Solar - Utility PV",
    "Land-Based Wind",
]

data_processing_object = DataProcessing()
data_processing_object.load_gen_data(
    bus_data_path=bus_data_path,
    cost_data_path=cost_data_path,
    candidate_gens=candidate_gens,
    save_csv=False,
    candidate_gen_csv_path=f"{data_path}/gen.csv",
    candidate_storage_csv_path=f"{data_path}/storage.csv",
    candidate_branch_csv_path=f"{data_path}/branch.csv",
)

# Populate and create GTEP model
mod_object = ExpansionPlanningModel(
    data=data_object,
    cost_data=data_processing_object,
)

mod_object.config["include_investment"] = True
mod_object.config["include_commitment"] = False
mod_object.config["include_redispatch"] = True
mod_object.config["scale_loads"] = False
mod_object.config["transmission"] = True
mod_object.config["storage"] = False
mod_object.config["flow_model"] = "transport"
mod_object.config["advanced_hydro"] = True

mod_object.create_model()
print("model is created!")

mod_object.model.renewableQuota.pprint()
mod_object.model.investmentStage[1].renewableCurtailmentInvestment.pprint()

export_block_constraints_hierarchical(
    mod_object.model,
    output_folder="naerm_nonsense_invesstment",
    individual_entries=True,
    recurse=True,
)

# # print(mod_object.model.md.data['elements']['generator']['AESO_cc_gas'])
# from pyomo.core.expr.numvalue import as_numeric, is_numeric_data

# def get_variables_in_expr(expression):
#     # If it's a Pyomo variable, yield it
#     if hasattr(expression, "is_variable_type") and expression.is_variable_type():
#         yield expression
#     # If the expression has arguments (like a sum or product), recurse
#     elif hasattr(expression, "args"):
#         for arg in expression.args:
#             yield from get_variables_in_expr(arg)

# for var in get_variables_in_expr(mod_object.model.total_cost_objective_rule.expr):
#     if not is_numeric_data(var):
#         if not type(var) == bool:
#             print(f"Non-numeric term found: {type(var)} with value {var}")

# # Check specific terms during your constraint building logic
# # for term in mod_object.model.total_cost_objective_rule.expr.args:
# #     if not is_numeric_data(term):
# #         print(f"Non-numeric term found: {type(term)} with value {term}")
# raise SystemExit

# mod_object.model.total_cost_objective_rule.pprint()
# mod_object.model.investmentStage[1].genInstalled['AESO_cc_gas'].pprint()

pyo.TransformationFactory("gdp.bigm").apply_to(mod_object.model)
print("model is transformed!")

solver = "xpress"
opt = pyo.SolverFactory(solver)
xpress.init("/Users/jkskolf/naerm_xpauth.xpr")
if solver == "xpress":
    log_folder = "xpress_log_files"
    options_dict = {
        "logfile": log_folder + "/" + solver + ".log",
    }
    # print(dir(xpress.controls))
    xpress.controls.heurdivespeedup = 0
    xpress.controls.heursearchrootcutfreq = 1
    # xpress.controls.miprelstop = 0.1
    # xpress.controls.heurfreq = 0  # disable most heuristics
    # xpress.controls.threads = 0  # use all available threads
    # xpress.controls.presolve = 3
    # xpress.controls.cutstrategy = 3
    xpress.controls.scaling = 16
    xpress.controls.maxmiptasks = 16

    mod_object.results = opt.solve(
        mod_object.model,
        tee=True,
        # logfile=log_folder + "/" + solver + ".log",
    )
    mod_object.model.investmentStage[1].renewableCurtailmentInvestment.pprint()
else:
    options_dict = {"MIPGap": 0.05}
    mod_object.results = opt.solve(
        mod_object.model,
        tee=True,
        options=options_dict,
        # logfile=log_folder + "/" + solver + ".log",
    )
# mod_object.model.operatingCostTotal.display()
# mod_object.model.expansionCostTotal.display()
# mod_object.model.penaltyCostTotal.display()
# print(pyo.value(mod_object.model.total_cost_objective_rule))

# Save the results in .json files using the solution class
dir_name = f"NAERM_initial_testing_{data_date}"

# Define the plot type for the generationmix. The options are:
# stackplot, treemap or pie chart. If nothing is select, all the files
# will be produced, as the default
plot_type = "all"

sol_object = ExpansionPlanningSolution(data_path)
sol_object.save_results_in_json_files(mod_object, dir_name)

# Create stack plots, treemaps, and pie charts for gen mix for
# dispatchable and renewable generators
case_json = "dispatchables"
sol_object.create_plots(case_json, dir_name, data_path, plot_type)
case_json = "renewables"
sol_object.create_plots(case_json, dir_name, data_path, plot_type)
case_json = "combined"
sol_object.create_plots(case_json, dir_name, data_path, plot_type)

# Create stackgraph
day_hour_list = [("2034-07-12 00:00", 18), ("2034-07-12 00:00", 4)]
sol_object.create_stackgraph_and_metrics(dir_name, rep_days, day_hour_list)

# # Create report
# sol_object.create_html_report(dir_name, plot_type)
