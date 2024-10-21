import pyomo.common.unittest as unittest
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
import logging

from gtep.validation import clone_timeseries, filter_pointers, populate_generators, populate_transmission

input_data_source = "./gtep/data/5bus"
output_data_source = "./gtep/tests/data/5bus_out"

def test_solution():
    data_object = ExpansionPlanningData()
    data_object.load_prescient(input_data_source)

    mod_object = ExpansionPlanningModel(
        stages=2,
        data=data_object.md,
        num_reps=2,
        len_reps=1,
        num_commit=6,
        num_dispatch=4,
    )
    mod_object.create_model()
    TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model) 
    TransformationFactory("gdp.bigm").apply_to(mod_object.model) 
    # opt = SolverFactory("gurobi")
    # opt = Gurobi()
    opt = Highs()
    # # mod_object.results = opt.solve(mod_object.model, tee=True)
    mod_object.results = opt.solve(mod_object.model) 

    sol_object = ExpansionPlanningSolution()
    sol_object.load_from_model(mod_object)
    sol_object.dump_json("./gtep/tests/test_solution.json")
    return sol_object

solution = test_solution()

class TestValidation(unittest.TestCase):
    def test_populate_generators(self):
        populate_generators(input_data_source, solution, output_data_source)


    def test_populate_transmission(self):
        populate_transmission(input_data_source, solution, output_data_source)

    def test_filter_pointers(self):
        filter_pointers(input_data_source, output_data_source)

    def test_clone_timeseries(self):
        clone_timeseries(input_data_source, output_data_source)