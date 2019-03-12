import pandas as pd
import copy

def create_dataframes(var_to_plot, month_list, Num_Realizations, Num_Years, TSID_key_list, Locations_to_Import_NEW,
                      TSID, simulation_title, how_type = 'mean'):
    units_conversion = 1.0
    master_column_list = []
    Scenario_Dictionary = {}
    new_index = {}
    for scenarios in TSID_key_list:
        for locations in Locations_to_Import_NEW[scenarios]:
            DF_SUM = TSID[scenarios][locations][var_to_plot].resample('M',how=how_type).copy(deep=True)*units_conversion
            DF_SUM.columns = ['Realization1' for i in range(Num_Realizations[scenarios])]  # Rename all columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.
            # In the sequence below, use DF_SUM to plot sums for months, or new_DF to plot all daily values for each month.
            column_list_arrays = [DF_SUM.ix[:,i] for i in range(1,Num_Realizations[scenarios])]  # have to skip zero-th element, because it's the first element included in the append call below.
            master_column_list.append(scenarios + ' ' + locations)
            Scenario_Dictionary[scenarios + ' ' + locations] = DF_SUM.ix[:,0].append(column_list_arrays)  # Take the first column of new_DF and append all the other columns to it.
        new_index[scenarios] = [i for i in range(int(Num_Realizations[scenarios]*Num_Years[scenarios]))]  # For now scenarios need to have same num years and num realizations or the cdf ranking wont work. so just take most recent value for num_realizations

    DF_SUM_STACKED = {}
    for scenario in master_column_list:
        DF_SUM_STACKED[scenario] = pd.DataFrame(Scenario_Dictionary[scenario])
        DF_SUM_STACKED[scenario].columns = [scenario]
    #Master_DF = DF_SUM_STACKED.copy(deep=True)
    Master_DF = copy.deepcopy(DF_SUM_STACKED)

    for scenario in master_column_list:
        DF_SUM_STACKED[scenario]['Month Number'] = DF_SUM_STACKED[scenario].index.month  # Create new column that stores
        Monthly_Energy_DF_GROUPED = DF_SUM_STACKED[scenario].groupby('Month Number')  # Create a group that slices the data by the Month Number column
        Monthly_Energy_DF_GROUPED.describe()  # Get statistics for each month
        # Generate a separate figure for every month, and output the monthly dataframe to the user.

    DF_Monthly_Columns = {}
    wse_traject_dict = {}
    for scenarios in TSID_key_list:
        locations = Locations_to_Import_NEW[scenarios][0]
        # For each scenario, create a dataframe that stores monthly info.
        DF_Monthly_Columns[scenarios] = pd.DataFrame(index=new_index[scenarios])  # Stores a column for each month with all values (Num_Realizations[sims]*(Sim_Dur/36500)) stored in each column, with generic numbered row index.
        i = 0
        wse_traject_dict[scenarios] = []
        # Create a subplot with 12 plots (one for each month). Save such a plot for each scenario.
        for mon in month_list:
            i+=1
            DF_Monthly_Columns[scenarios][mon] = DF_SUM_STACKED[scenarios + ' ' + locations][DF_SUM_STACKED[scenarios + ' ' + locations].index.month==i][[scenarios + ' ' + locations]].values
            wse_traject_dict[scenarios].append(DF_Monthly_Columns[scenarios][mon].mean())
    return new_index, DF_SUM_STACKED, wse_traject_dict