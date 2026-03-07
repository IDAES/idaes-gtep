import json
import pandas as pd

start_date = "12-30"
week = 52

with open(f"power_flows/power_flows_{start_date}.json", "r") as power_flows_file:
    power_flows = json.load(power_flows_file)
with open(f"max_flows/max_flows_{start_date}.json", "r") as max_flows_file:
    max_flows = json.load(max_flows_file)

branch_names = [
    "branch_2_3",
    "branch_1_2",
    "branch_1_4",
    "branch_4_10",
    "branch_1_10",
    "branch_3_4_0",
    "branch_3_4_1",
]

average_congestion = {}  # stores (per-line) average hourly congestion on each day
percent_over_99 = {}     # stores % of each day w/ congestion >= 0.99
percent_over_98 = {}     # stores % of each day w/ congestion >= 0.98
percent_over_95 = {}     # stores % of each day w/ congestion >= 0.95
percent_over_90 = {}     # stores % of each day w/ congestion >= 0.95
percent_over_70 = {}     # stores % of each day w/ congestion >= 0.70
percent_over_50 = {}     # stores % of each day w/ congestion >= 0.50

day = 0
while day < 1: #7:
    one_day = {}
    for i in range(24 * day + 1, 24 * day + 25):  # AVERAGE PER-LINE CONGESTION:
        one_day[i] = (
            sum(
                abs(
                    power_flows[
                        f"investmentStage[1].representativePeriod[1].commitmentPeriod[{i}].dispatchPeriod[1].powerFlow.{branch_name}"
                    ]
                )
                / max_flows[
                    f"investmentStage[1].representativePeriod[1].commitmentPeriod[{i}].dispatchPeriod[1].powerFlow.{branch_name}"
                ]
                for branch_name in branch_names
            )
        ) / len(branch_names)
    average_congestion[day + (7 * week)] = sum(one_day.values()) * (1 / 24)
    percent_over_99[day + (7 * week)] = sum(1 for c_val in one_day.values() if c_val >= 0.99)*(1/24)
    percent_over_98[day + (7 * week)] = sum(1 for c_val in one_day.values() if c_val >= 0.98)*(1/24)
    percent_over_95[day + (7 * week)] = sum(1 for c_val in one_day.values() if c_val >= 0.95)*(1/24)
    percent_over_90[day + (7 * week)] = sum(1 for c_val in one_day.values() if c_val >= 0.9)*(1/24)
    percent_over_70[day + (7 * week)] = sum(1 for c_val in one_day.values() if c_val >= 0.7)*(1/24)
    percent_over_50[day + (7 * week)] = sum(1 for c_val in one_day.values() if c_val >= 0.5)*(1/24)
    day += 1

#print(
#    f"highest_ave_congestion day is {max(average_congestion, key=average_congestion.get)}:",
#    max(average_congestion.values()),
#)
# print(f"highest_total_congestion day is {max(total_congestion, key=total_congestion.get)}", max(total_congestion.values()))

dict_to_dataframe = {
    'Average Per-Line Hourly Con.': average_congestion, 
    'Percent Con. >= 0.99': percent_over_99, 
    'Percent Con. >= 0.98': percent_over_98, 
    'Percent Con. >= 0.95': percent_over_95, 
    'Percent Con. >= 0.9': percent_over_90, 
    'Percent Con. >= 0.7': percent_over_70, 
    'Percent Con. >= 0.5': percent_over_50, 
}
df = pd.DataFrame.from_dict(dict_to_dataframe, orient='index')
#print(df)

if week == 0:
    df.to_csv('flows_ALLDAYS_5.csv')
else:
    main_df = pd.read_csv('flows_ALLDAYS_5.csv')
    main_df = main_df.reset_index(drop=True)
    df = df.reset_index(drop=True)
    final_df = pd.concat([main_df, df], axis=1)
    final_df.to_csv('flows_ALLDAYS_5.csv')
    print(final_df)

'''
    Q: how were the four default rep days in GTEP chosen???
    TODO:
    - what's happening w/ distribution of vals
    - visualize patters: are there lines that are constantly under-utilized/congested? 
'''
