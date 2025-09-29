from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from IPython import embed
import pyomo.environ as pyo
from forestlib.ph import stochastic_program
from forestlib.ph import ProgressiveHedgingSolver



def model_builder(scen, scen_args):
    model = pyo.ConcreteModel(scen["ID"])
    if scen["ID"] == "BelowAverageScenario":
        yields = [2, 2.4, 16]
    elif scen["ID"] == "AverageScenario":
        yields = [2.5, 3, 20]
    elif scen["ID"] == "AboveAverageScenario":
        yields = [3, 3.6, 24]
    # Variables
    model.X = pyo.Var(["WHEAT", "CORN", "BEETS"], within=pyo.NonNegativeReals)
    model.Y = pyo.Var(["WHEAT", "CORN"], within=pyo.NonNegativeReals)
    model.W = pyo.Var(
        ["WHEAT", "CORN", "BEETS_FAVORABLE", "BEETS_UNFAVORABLE"],
        within=pyo.NonNegativeReals,
    )

    # Objective function
    model.PLANTING_COST = (
        150 * model.X["WHEAT"] + 230 * model.X["CORN"] + 260 * model.X["BEETS"]
    )
    model.PURCHASE_COST = 238 * model.Y["WHEAT"] + 210 * model.Y["CORN"]
    model.SALES_REVENUE = (
        170 * model.W["WHEAT"]
        + 150 * model.W["CORN"]
        + 36 * model.W["BEETS_FAVORABLE"]
        + 10 * model.W["BEETS_UNFAVORABLE"]
    )
    model.OBJ = pyo.Objective(
        expr=model.PLANTING_COST + model.PURCHASE_COST - model.SALES_REVENUE,
        sense=pyo.minimize,
    )

    # Constraints
    model.CONSTR = pyo.ConstraintList()

    model.CONSTR.add(pyo.summation(model.X) <= 500)
    model.CONSTR.add(
        yields[0] * model.X["WHEAT"] + model.Y["WHEAT"] - model.W["WHEAT"] >= 200
    )
    model.CONSTR.add(
        yields[1] * model.X["CORN"] + model.Y["CORN"] - model.W["CORN"] >= 240
    )
    model.CONSTR.add(
        yields[2] * model.X["BEETS"]
        - model.W["BEETS_FAVORABLE"]
        - model.W["BEETS_UNFAVORABLE"]
        >= 0
    )
    model.W["BEETS_FAVORABLE"].setub(6000)

    return model

model_data = {
    "scenarios": [
        {
            "ID": "BelowAverageScenario",
            "Yield": {"WHEAT": 2.0, "CORN": 2.4, "SUGAR_BEETS": 16.0},
            "crops_multiplier": 1.0,
            "Probability": 0.3333333333333333,
        },
        {
            "ID": "AverageScenario",
            "Yield": {"WHEAT": 2.5, "CORN": 3.0, "SUGAR_BEETS": 20.0},
            "crops_multiplier": 1.0,
            "Probability": 0.3333333333333333,
        },
        {
            "ID": "AboveAverageScenario",
            "Yield": {"WHEAT": 3.0, "CORN": 3.6, "SUGAR_BEETS": 24.0},
            "crops_multiplier": 1.0,
            "Probability": 0.3333333333333333,
        },
    ]
}
model_data_gtep = {
    "scenarios": [
        {
            "ID": "NormalScenario",
            "Yield": {"WHEAT": 2.0, "CORN": 2.4, "SUGAR_BEETS": 16.0},
            "crops_multiplier": 1.0,
            "Probability": 1.0,
        }
    ]
}
def model_builder_wrapper(scen, scen_args):
    data_path = "./data/9bus"
    data_object = ExpansionPlanningData()
    data_object.load_prescient(data_path)

    mod_object = ExpansionPlanningModel(
        stages=2,
        data=data_object.md,
        num_reps=2,
        len_reps=1,
        num_commit=2,
        num_dispatch=2,
    )
    mod_object.create_model()
    TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model) 
    TransformationFactory("gdp.bigm").apply_to(mod_object.model)
    return mod_object.model


opt = Gurobi()
#opt = Highs()
GTEP_SP = stochastic_program(first_stage_variables=["investmentStage[1].genOperational[G3].binary_indicator_var"])
GTEP_SP.initialize_model(model_data=model_data_gtep, model_builder=model_builder_wrapper)
ph = ProgressiveHedgingSolver()
#embed()
results = ph.solve(GTEP_SP, max_iterations=10, solver="gurobi", rho=10)

#mod_object.results = opt.solve(mod_object.model) 
#sol_object = ExpansionPlanningSolution()
#sol_object.load_from_model(mod_object)
#sol_object.dump_json("./gtep_solution.json")

#sol_object.import_data_object(data_object)

#sol_object.plot_levels(save_dir="./plots/")
