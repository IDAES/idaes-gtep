import json
import pandas as pd

start_date = "12-30"
week = 52

with open(f'load_shed/load_shed_{start_date}.json', 'r') as file:
    load_shed = json.load(file)

average_load_shed = {} # stores average hourly load shed on each day
total_load_shed = {}   # stores total load shed on each day
percent_load_shed = {} # stores % of each day having non-zero load shed

day = 0
while day < 7: #2:
    one_day = {}
    for i in range(24*day + 1, 24*day + 25):
        one_day[i] = sum(load_shed[l] for l in load_shed.keys() if f'commitmentPeriod[{i}]' in l) # summing load shed from all buses
    nonzero_load_shed_count = sum(1 for ls_val in one_day.values() if ls_val > 0.001) # summing num. of hours with non-zero load shed
    total_load_shed[day + (7 * week)] = sum(one_day.values())
    average_load_shed[day + (7 * week)] = total_load_shed[day + (7 * week)]*(1/24)
    percent_load_shed[day + (7 * week)] = nonzero_load_shed_count*(1/24)
    day += 1

dict_to_dataframe = {
    'Average Load Shed': average_load_shed, 
    'Total Load Shed': total_load_shed, 
    'Percent Load Shed': percent_load_shed
}
df = pd.DataFrame.from_dict(dict_to_dataframe, orient='index')
#print(df)

#print(f"highest_ave_load_shed day is {max(average_load_shed, key=average_load_shed.get)}:", max(average_load_shed.values()))
#print(f"highest_total_load_shed day is {max(total_load_shed, key=total_load_shed.get)}:", max(total_load_shed.values()))
#print(f"largest_percent_hrs_load_shed day is {max(percent_load_shed, key=percent_load_shed.get)}:", max(percent_load_shed.values()))

if week == 0:
    df.to_csv('loadshed_ALLDAYS.csv')
else:
    main_df = pd.read_csv('loadshed_ALLDAYS.csv')
    main_df = main_df.reset_index(drop=True)
    df = df.reset_index(drop=True)
    final_df = pd.concat([main_df, df], axis=1)
    final_df.to_csv('loadshed_ALLDAYS.csv')
    print(final_df)

'''
    I think we want to use average hourly load shed as our metric here?? or do we want to use average per-bus load shed???
    - how many extreme days to select? 

    PROBLEM: choose day with largest % hours of load shed; almost every period has load shed

'''
