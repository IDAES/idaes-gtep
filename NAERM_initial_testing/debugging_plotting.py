import json
import re
from collections import defaultdict

with open("generation.json", "r") as f:
    data = json.load(f)

rep_commit_dispatch = defaultdict(set)

for key in data:
    rp = re.search(r"representativePeriod\[(\d+)\]", key)
    cp = re.search(r"commitmentPeriod\[(\d+)\]", key)
    dp = re.search(r"dispatchPeriod\[(\d+)\]", key)

    if rp and cp and dp:
        rep = int(rp.group(1))
        com = int(cp.group(1))
        dis = int(dp.group(1))
        rep_commit_dispatch[rep].add((com, dis))

for rep in sorted(rep_commit_dispatch):
    cps = sorted({c for c, d in rep_commit_dispatch[rep]})
    dps = sorted({d for c, d in rep_commit_dispatch[rep]})
    print(f"RepresentativePeriod {rep}:")
    print(f"  commitmentPeriods: {cps[:10]} ... total={len(cps)}")
    print(f"  dispatchPeriods: {dps}")
