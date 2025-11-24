#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pyomo.environ as pyo
import pyomo.gdp as gdp

valid_names = ["Inst", "Oper", "Disa", "Ext", "Ret"]
renewable_investments = {}
dispatchable_investments = {}
for var in mod_object.model.component_objects(pyo.Var, descend_into=True):
    for index in var:
        for name in valid_names:
            if name in var.name:
                print(var, index, pyo.value(var[index]))
                if pyo.value(var[index]) >= 0.001:
                    renewable_investments[var.name + str(index)] = pyo.value(var[index])
for var in mod_object.model.component_objects(gdp.Disjunct, descend_into=True):
    for index in var:
        for name in valid_names:
            if name in var.name:
                print(var.name)
                print(var, index, pyo.value(var[index].indicator_var))
                if pyo.value(var[index].indicator_var) == True:
                    dispatchable_investments[var.name + str(index)] = pyo.value(
                        var[index].indicator_var
                    )

import json

with open("renewable_investments.json", "w") as fil:
    json.dump(renewable_investments, fil)
with open("dispatchable_investments.json", "w") as fil:
    json.dump(dispatchable_investments, fil)
