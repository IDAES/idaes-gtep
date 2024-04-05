# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
import logging

logger = logging.getLogger(__name__)


class ExpansionPlanningSolution:
    def __init__(self):
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
        self.results = gtep_model.results
        self.stages = gtep_model.stages
        self.formulation = gtep_model.formulation
        self.data = gtep_model.data
        self.num_reps = gtep_model.num_reps
        self.len_reps = gtep_model.len_reps
        self.num_commit = gtep_model.num_commit
        self.num_dispatch = gtep_model.num_dispatch

    def dump_json(self):
        pass
