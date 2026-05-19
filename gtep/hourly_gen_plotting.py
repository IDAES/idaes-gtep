import json
import re
from collections import defaultdict
import matplotlib.pyplot as plt

# ----------------------------
# Settings
# ----------------------------
INPUT_JSON = "generation.json"
OUTPUT_FIG = "hourly_generation_stackgraph_15_representative_days_broad.png"


# ----------------------------
# Helper: map detailed classes to broad categories
# ----------------------------
def map_to_broad_category(category, gen_name):
    """
    category: 'thermalGeneration' or 'renewableGeneration'
    gen_name: e.g. 'AESO_ct_gas', 'AESO_wind-c', 'BANC_thermal_other'
    """

    parts = gen_name.split("_", 1)
    tech = parts[1] if len(parts) == 2 else gen_name

    # Merge "-c" into equivalent class
    if tech.endswith("-c"):
        tech = tech[:-2]

    # Normalize some naming
    tech_lower = tech.lower()

    # Renewable side
    if category == "renewableGeneration":
        if "wind" in tech_lower:
            return "wind"
        if "solar" in tech_lower or "pv" in tech_lower:
            return "solar"
        return "other_renewable"

    # Thermal side
    if category == "thermalGeneration":
        if tech_lower == "cc_gas" or tech_lower == "gas_cc":
            return "gas_cc"
        if tech_lower == "ct_gas" or tech_lower == "gas_ct":
            return "gas_ct"

        # catch variants like gas_cc-c / gas_ct-c after -c removal
        if "cc" in tech_lower and "gas" in tech_lower:
            return "gas_cc"
        if "ct" in tech_lower and "gas" in tech_lower:
            return "gas_ct"

        if "coal" in tech_lower:
            return "coal"
        if tech_lower == "dr":
            return "dr"
        if "thermal_other" in tech_lower or tech_lower == "other":
            return "thermal_other"

        return "other_thermal"

    return "unknown"


# ----------------------------
# Load JSON
# ----------------------------
with open(INPUT_JSON, "r") as f:
    data = json.load(f)

# ----------------------------
# Parse and aggregate
# ----------------------------
# Build a continuous hour across representative periods:
# global_hour = (representativePeriod - 1) * 24 + commitmentPeriod

hourly_by_type = defaultdict(lambda: defaultdict(float))
global_hours_found = set()

for key, value in data.items():
    rp_match = re.search(r"representativePeriod\[(\d+)\]", key)
    cp_match = re.search(r"commitmentPeriod\[(\d+)\]", key)
    gen_match = re.search(r"(thermalGeneration|renewableGeneration)\.([^.]+)$", key)

    if not (rp_match and cp_match and gen_match):
        continue

    representative_day = int(rp_match.group(1))
    hour_in_day = int(cp_match.group(1))
    category = gen_match.group(1)
    gen_name = gen_match.group(2)

    global_hour = (representative_day - 1) * 24 + hour_in_day
    global_hours_found.add(global_hour)

    broad_type = map_to_broad_category(category, gen_name)
    hourly_by_type[global_hour][broad_type] += float(value)

# ----------------------------
# Build plotting arrays
# ----------------------------
global_hours = sorted(global_hours_found)

preferred_order = [
    "solar",
    "wind",
    "gas_ct",
    "gas_cc",
    "coal",
    "dr",
    "thermal_other",
    "other_thermal",
    "other_renewable",
    "unknown",
]

all_types_found = {t for h in global_hours for t in hourly_by_type[h].keys()}
all_types = [t for t in preferred_order if t in all_types_found] + sorted(
    all_types_found - set(preferred_order)
)

series = []
for t in all_types:
    series.append([hourly_by_type[h].get(t, 0.0) for h in global_hours])

# Optional color map
color_map = {
    "solar": "#FDB813",
    "wind": "#4C78A8",
    "gas_ct": "#E45756",
    "gas_cc": "#72B7B2",
    "coal": "#4D4D4D",
    "dr": "#B279A2",
    "thermal_other": "#8C564B",
    "other_thermal": "#9C755F",
    "other_renewable": "#54A24B",
    "unknown": "#BAB0AC",
}
colors = [color_map.get(t, None) for t in all_types]

# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots(figsize=(22, 9))
ax.stackplot(global_hours, series, labels=all_types, colors=colors)

ax.set_title("Hourly Generation by Broad Type Across 15 Representative Days")
ax.set_xlabel("Hour Across Representative Days")
ax.set_ylabel("Generation")

# Draw representative-day boundaries
num_days = max((h - 1) // 24 + 1 for h in global_hours)
for d in range(1, num_days):
    ax.axvline(d * 24 + 0.5, color="k", linestyle="--", linewidth=0.5, alpha=0.4)

# Label each representative period
day_centers = [d * 24 + 12.5 for d in range(num_days)]
day_labels = [f"RP{d+1}" for d in range(num_days)]
ax.set_xticks(day_centers)
ax.set_xticklabels(day_labels)

ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=9)
plt.tight_layout()

plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved figure to {OUTPUT_FIG}")
