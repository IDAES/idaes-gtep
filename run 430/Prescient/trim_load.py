import pandas as pd

factor = 2

rt_load = pd.read_csv("/Users/jkskolf/idaes-gtep/run 430/Prescient/REAL_TIME_load.csv")
cols = [str(i) for i in range(1, 124)]
rt_load[cols] /= factor
rt_load.to_csv(
    "/Users/jkskolf/idaes-gtep/run 430/Prescient/REAL_TIME_load.csv", index=False
)

da_load = pd.read_csv("/Users/jkskolf/idaes-gtep/run 430/Prescient/DAY_AHEAD_load.csv")
cols = [str(i) for i in range(1, 124)]
da_load[cols] /= factor
da_load.to_csv(
    "/Users/jkskolf/idaes-gtep/run 430/Prescient/DAY_AHEAD_load.csv", index=False
)
