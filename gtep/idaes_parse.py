from gtep.gtep_solution import ExpansionPlanningSolution
import pandas as pd
import json
import re
import matplotlib as mpl
from matplotlib import colormaps
import seaborn as sns


with open("./Sim Eng Viz/4Hourly/GTEP/idaes_solution.json") as fil:
    solutions = json.load(fil)

target_investment_stage = "2"
target_representative_period = "2"

filtered_solns = solutions["results"]["primals_tree"]["investmentStage[2]"][
    "representativePeriod[2]"
]

filtered_solns.pop("_logical_to_disjunctive")


def name_filter(gen_name):
    return re.search(r"\[.*\]", gen_name).group(0)[1:-1]


thermal_gen = {}
curtailment = {}
load_shed = {}
for key in filtered_solns.keys():
    for k, v in filtered_solns[key]["dispatchPeriod[1]"].items():
        if "Gen" in k:
            k = name_filter(k)
            if thermal_gen.get(k):
                thermal_gen[k].append(v["value"])
            else:
                thermal_gen[k] = [v["value"]]
        if "Curt" in k:
            if curtailment.get(k):
                curtailment[k].append(v["value"])
            else:
                curtailment[k] = [v["value"]]
        if "Shed" in k:
            if load_shed.get(k):
                load_shed[k].append(v["value"])
            else:
                load_shed[k] = [v["value"]]
    # only for the four hourly...
    if filtered_solns[key].get("dispatchPeriod[2]"):
        for k, v in filtered_solns[key]["dispatchPeriod[2]"].items():
            if "Gen" in k:
                k = name_filter(k)
                if thermal_gen.get(k):
                    thermal_gen[k].append(v["value"])
                else:
                    thermal_gen[k] = [v["value"]]
        for k, v in filtered_solns[key]["dispatchPeriod[3]"].items():
            if "Gen" in k:
                k = name_filter(k)
                if thermal_gen.get(k):
                    thermal_gen[k].append(v["value"])
                else:
                    thermal_gen[k] = [v["value"]]
        for k, v in filtered_solns[key]["dispatchPeriod[4]"].items():
            if "Gen" in k:
                k = name_filter(k)
                if thermal_gen.get(k):
                    thermal_gen[k].append(v["value"])
                else:
                    thermal_gen[k] = [v["value"]]


print(thermal_gen)
print(curtailment)
print(load_shed)

nonzero_gennames = ["10_STEAM", "3_CT", "4_STEAM"]
for gen in nonzero_gennames:
    thermal_gen.pop(gen)

gen_df = pd.DataFrame.from_dict(dict(sorted(thermal_gen.items())))
print(gen_df)
prescient_gen_df = pd.read_csv(
    "./Sim Eng Viz/4Hourly/Prescient/results/generation_out.csv",
    names=["gen", "val"],
)

print(prescient_gen_df)
for label, content in prescient_gen_df.items():
    print(f"label: {label}")
    print(f"content: {content}", sep="\n")

prescient_gen = {}
for row in prescient_gen_df.iterrows():
    gen = row[-1]["gen"]
    val = row[-1]["val"]
    if prescient_gen.get(gen):
        prescient_gen[gen].append(val)
    else:
        prescient_gen[gen] = [val]

for gen in nonzero_gennames:
    prescient_gen.pop(gen)
prescient_sorted_gen_df = pd.DataFrame.from_dict(dict(sorted(prescient_gen.items())))

my_cmap = sns.color_palette("muted")


ax1 = gen_df.plot.bar(stacked=True, color=my_cmap)
ax1.set_title("GTEP, 4-Hourly Commitment (04/01)")
ax1.set_xlabel("Hour")
ax1.set_ylabel("Generation (MW)")
ax1.figure.savefig("4blah.png")

ax2 = prescient_sorted_gen_df.plot.bar(stacked=True, color=my_cmap)
ax2.set_title("Prescient (04/01)")
ax2.set_xlabel("Hour")
ax2.set_ylabel("Generation (MW)")
ax2.figure.savefig("4blah2.png")

print(gen_df)
