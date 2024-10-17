import pyomo.common.unittest as unittest
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
import logging

input_data_source = "./gtep/data/5bus"

class TestValidation(unittest.TestCase):
    def test_populate_generators(self):
        pass

    def test_populate_transmission(self):
        pass

    def test_filter_pointers(self):
        pass

    def test_clone_timeseries(self):
        pass