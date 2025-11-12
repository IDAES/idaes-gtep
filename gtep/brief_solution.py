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
