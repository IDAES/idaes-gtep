import numpy as np
import matplotlib.pyplot as plt

try:
    import ujson as json
except ImportError:
    import json

import pyomo.environ as pyo

with open("generation.json", "r") as F:
    gen_data = json.load(F)


# CC	CC	CC	tab20	1
# CT	CT	CT	tab20	3
# COAL	COAL	CO	tab20	5
# NUCLEAR	NUCLEAR	NU	tab20	2
# PV	PV	PV	tab20	9
# WIND	WIND	WI	tab20	11
# THERMAL	THERMAL	TH	tab20	13
# HYDRO	HYDRO	HY	tab20	20
# BATT	BATT	BA	tab20	15
# ES4	ES4	ES	tab20	17
# PS	PS	PS	tab20	19
# Load Shed	Load Shed	SL	tab20	7


tab20 = plt.get_cmap("tab20")
GEN_TYPES = {
    "coal": tab20(5),
    "cc_gas": tab20(1),
    "ct_gas": tab20(3),
    "dr": tab20(19),
    "nuclear": tab20(2),
    "solar": tab20(9),
    "thermal_other": tab20(13),
    "wind": tab20(11),
    "gas_cc-c": tab20(0),
    "gas_ct-c": tab20(2),
    "pv-c": tab20(8),
    "wind-c": tab20(10),
}

generation = {}
for g, val in gen_data.items():
    c = list(pyo.ComponentUID(g)._cids)
    _, (stage,) = c.pop(0)
    if stage not in generation:
        generation[stage] = {}
    stage = generation[stage]
    _, (period,) = c.pop(0)
    if period not in stage:
        stage[period] = {}
    period = stage[period]
    _, (committment,) = c.pop(0)
    if committment not in period:
        period[committment] = {}
    committment = period[committment]
    _, (dispatch,) = c.pop(0)
    if dispatch not in committment:
        committment[dispatch] = dict.fromkeys(GEN_TYPES, 0)
    dispatch = committment[dispatch]
    gen_name = c[-1][0]
    _type = None
    for gt in GEN_TYPES:
        if gen_name.endswith(gt):
            _type = gt
    if _type is None:
        raise RuntimeError(f"Cannot map generator name '{gen_name}' to type")
    dispatch[_type] += val


time_periods = [
    (s, p, c, d)
    for s in generation
    for p in generation[s]
    for c in generation[s][p]
    for d in generation[s][p][c]
]
times = list(range(len(time_periods)))

fig, ax = plt.subplots()
bottom = np.zeros(len(time_periods))
width = 1
for name, color in GEN_TYPES.items():
    values = np.array([generation[s][p][c][d][name] for s, p, c, d in time_periods])
    p = ax.bar(
        times,
        values,
        width,
        label=name,
        color=color,
        bottom=bottom,
        edgecolor="#FFFFFF",
        linewidth=0.5,
    )
    bottom += values



# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# ax.set_title("Generation Mix")
# ax.set_ylabel('Nameplate Capacity [MW]')
# ax.set_xlabel('Investment Year')

plt.show()