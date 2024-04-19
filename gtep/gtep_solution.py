# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
import logging

import json
from pathlib import Path

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
        self.results = gtep_model.results  # Highs results object
        self.stages = gtep_model.stages  # int
        self.formulation = gtep_model.formulation  # None (???)
        self.data = gtep_model.data  # ModelData object
        self.num_reps = gtep_model.num_reps  # int
        self.len_reps = gtep_model.len_reps  # int
        self.num_commit = gtep_model.num_commit  # int
        self.num_dispatch = gtep_model.num_dispatch  # int

        pass

    def dump_json(self, filename="./gtep_solution.json"):

        dump_filepath = Path(filename)
        with open(dump_filepath, "w") as fobj:
            json.dump(self._to_dict(), fobj)

    def _to_dict(self):

        results_dict = {
            "solution_loader": self.results.solution_loader,  # object
            "termination_condition": self.results.termination_condition,  # object
            "best_feasible_objective": self.results.best_feasible_objective,
            "best_objective_bound": self.results.best_objective_bound,
            "wallclock_time": self.results.wallclock_time,
        }

        # "best_feasible_objective", "best_objective_bound", and "wallclock_time" are all numbers, dont need subhandlers

        # subhandle "termination_condition"
        results_dict["termination_condition"] = {
            "value": self.results.termination_condition.value,
            "name": self.results.termination_condition.name,
        }

        # subhandle "solution_loader"
        results_dict["solution_loader"] = {"primals": {}}

        for key, val in self.results.solution_loader.get_primals()._dict.items():
            results_dict["solution_loader"]["primals"][key] = {
                "name": val[0].name,
                "upper": val[0].upper,
                "value": val[0].value,
                "bounds": val[0].bounds,
            }

        # renest "termination_condition" as a json-friendly dictionary
        # things are either vars (which have some sort of signifier in [] brackets) or are an attribute, which dont
        # the name variable will give it away
        results_dict["primals_tree"] = {}

        for key, val in self.results.solution_loader.get_primals()._dict.items():
            # split the name to figure out depth
            split_name = val[0].name.split(".")

            # start at the bottom and nest accordingly
            tmp_dict = {
                "name": val[0].name,
                "upper": val[0].upper,
                "value": val[0].value,
                "bounds": val[0].bounds,
            }

            # allocate the nested dictionary
            def nested_set(this_dict, key, val):
                if len(key) > 1:
                    this_dict.setdefault(key[0], {})
                    nested_set(this_dict[key[0]], key[1:], val)
                else:
                    this_dict[key[0]] = val

            nested_set(results_dict["primals_tree"], split_name, tmp_dict)
        out_dict = {"data": self.data.data, "results": results_dict}

        return out_dict
