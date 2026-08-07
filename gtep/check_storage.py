from pathlib import Path
import pandas as pd
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_model import ExpansionPlanningModel

data_dir = Path("gtep/data/9_bus_GTEP_dir")

data_object = ExpansionPlanningData(
    stages=1,
    num_reps=1,
    num_commit=1,
    num_dispatch=1,
)
data_object.load_prescient(data_dir)
mod_object = ExpansionPlanningModel(config={"storage": True}, data=data_object)
mod_object.create_model()

# This correctly prints the storage units
print("-" * 50)
print("storage Set on the model:")
print([stor for stor in mod_object.model.storage])

# This prints an empty list for every bus (!)
print("-" * 50)
print("storageByBus Param from the model:")
for bus, stors in mod_object.model.storageByBus.items():
    print(f"at bus {bus}:", [stor for stor in stors])

print("-" * 50)
print("First 5 storage ID and bus from the data object:")
n = 0
for stor, stor_data in data_object.md.data["elements"]["storage"].items():
    print(stor, "  \t", stor_data["bus"])
    n += 1
    if n >= 5:
        break

print("-" * 50)
print("Compare to first 5 generator ID and bus from the data object:")
n = 0
for gen, gen_data in data_object.md.data["elements"]["generator"].items():
    print(gen, "  \t", gen_data["bus"])
    n += 1
    if n >= 5:
        break

print("-" * 50)
print("Top of storage csv:")
print(pd.read_csv((data_dir / "storage.csv").resolve()).iloc[:5,:5])

print("-" * 50)
print("Compare to top of gen csv:")
print(pd.read_csv((data_dir / "gen.csv").resolve()).iloc[:5,:5])