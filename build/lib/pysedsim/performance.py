# -*- coding: utf-8 -*-

'''
Purpose: This file evaluates the performance of a given simulation run for purposes of returning those values to an
external optimization model (e.g., Borg) as objective function values.
'''

# Import libraries
import pandas as pd
import copy

def Performance(output_dict, Sim_Dur, opt_dict=None, main_output_file_dir=None, sim_title=None):
    '''

    Purpose: Evaluates the performance of specified system elements in a particular PySedSim system configuration.

    The ideal use for this module is in an optimization setting, where the performance values returned by the
    function correspond to "objective function" values returned to the optimization model in each function
    evaluation.

    :param Sim_Dur: Simulation Duration in days
    :param opt_dict: optimization dictionary of preferences defined in optimization_pyedsim
    :param main_output_file_dir: string for directory in which output files are located (so can be imported/processed)
    :param sim_title: title of the simulation being run
    :param output_dict: output dictionary of data frames produced by PySedSim (stores variables user requests to export)
    :return: objs (a list of objective function values)
    '''

    # If this function is being used as part of an optimization routine, unpack optimization preferences/parameters.
    if opt_dict is not None:
        n_objs = opt_dict['Num Objectives']
        n_constrs = opt_dict['Num Constraints']
        simulation_title = opt_dict['Scenario Name']

        # Initialize arrays to store objective function values and constraints
        objs = [0.0]*n_objs
        if n_constrs > 0:
            constrs = [0.0]*n_constrs

    # Compute objective function value (performance). Borg is a minimizer, so negate values of all objectives that
    # are to be maximized, but don't negate values of objectives that are to be minimized.

    # Loop through objectives, computing relevant performance measure values according to user specifications.
    for i in range(n_objs):
        obj_value = 0
        obj_name = opt_dict['Objective Names Ordered List'][i]
        state_var_name = opt_dict['State Variable Ordered List'][i]
        locs = opt_dict['Objective Preferences'][obj_name]['Locations']  # List of location names
        # Loop through ordered objective variables and calculate objective function values to be returned to the
        # optimization model. Use objective preferences to guide calculation of objective values.

        for loc in locs:
            # If multiple locations apply to this objective, the objective values will be added.
            try:
                time_slice = opt_dict['Objective Preferences'][obj_name]['Time Slice']
                # time_slice is a number. Compute a single value sliced at a particular time.
                if time_slice is not None:
                    if time_slice > 1:
                        if time_slice not in [int, long]:
                            time_slice = int(time_slice)  # Force day value to be an integer
                        if time_slice > (Sim_Dur-1):
                            time_slice = int(Sim_Dur-1)  # Can't exceed duration of the simulation
                    else:
                        # Time slice is a fraction, so use it to multiply by simulation duration.
                        time_slice = int(time_slice*(Sim_Dur-1))

                # Compute statistic over time series.
                Locations_to_Import = {sim_title: [loc]}
                Input_Dict = {sim_title: output_dict}
                var_to_eval = state_var_name
                resample_freq = opt_dict['Objective Preferences'][obj_name]['Resample Frequency']
                resample_freq_orig = None
                if resample_freq == 'O':
                    resample_freq = str(time_slice) + 'D'  # Sample at a frequency of only the last day of simulation.
                    resample_freq_orig = opt_dict['Objective Preferences'][obj_name]['Resample Frequency']
                resample_stat = opt_dict['Objective Preferences'][obj_name]['Resample Stat']
                distribution_stat = opt_dict['Objective Preferences'][obj_name]['PM Distribution Stat']
                Num_Realizations = {sim_title: opt_dict['Num Realizations']}
                Num_Years = {sim_title: output_dict[loc][state_var_name].index[-1].year - output_dict[loc][state_var_name].index[0].year + 1}
                [temp_obj_value, Master_DF] = Performance_Measure_Calc(Input_Dict, Locations_to_Import, var_to_eval,
                                                                       resample_freq, resample_freq_orig, resample_stat,
                                                                       distribution_stat, Num_Realizations,
                                                                       Num_Years)
                obj_value += temp_obj_value
            except KeyError:
                raise KeyError('Variable %s is not stored in output dictionary. Either you did not specify that '
                               'this variable be exported in the "Export Preferences" worksheet of the input '
                               'file, or you incorrectly specified it as a variable for location %s '
                               'in the "Optimization" worksheet of the input file' % (state_var_name, loc))
        # Adjust objective value for use in Borg (Borg tries to minimize all objectives).
        if opt_dict['Objective Preferences'][obj_name]['Type'] == 'Max':
            objs[i] = -1*obj_value
        else:
            objs[i] = obj_value

    # Return final objective values
    if opt_dict is not None:
        if n_constrs > 0:
            return (objs, constrs)
        else:
            return objs


def Performance_Measure_Calc(Input_Dict, Locations_to_Import, var_to_eval, resample_freq,
                             resample_freq_orig, resample_stat, distribution_stat, Num_Realizations, Num_Years,
                             units_conversion=1, scenarios_list=None, after_n_years=None, combine_months=None,
                             months_to_plot=None, empirical_cdf = 'No'):
    '''
    Purpose: Use Pandas to compute and return performance measures (values summarizing simulation performance)

    Ideal use is in an optimization setting (function is called during each iteration of the evolutionary
    optimization algorithm, for example).

    :param Locations_to_Import: Dictionary, wherein keys are the names of simulation-optimizaiton scenarios being
    evaluated, such as ["Alternative with Flushing"], and the contents are a list of strings of locations for each of
    those scenarios for which data should be imported. Example: Locations_to_Import = {"Alternative with Flushing":
    ["Sambor Alternative"]}

    :param resample_freq: Frequency of resampling. Choices: 'A', 'M', or 'D' for Annual, Monthly, Daily, respectively.
    :param resample_stat: Method of manipulating resampled data. Choices: ['mean', 'variance', 'median', 'min', 'max',
    'sum']
    :param distribution_stat: Statistic to take over the resulting distribution of resampled statistics. Choices: [
    'mean', 'variance', 'median', 'min', 'max', 'sum'] as well as any quantile value in [0,1]
    :param units_conversion: Value (any number) to divide all values in DF by
    :param scenarios_list: A list of the scenario names + locations that is exact order in which you want them to be
    called from the dataframe when plotting. Warning: dataframe column names are set internally in function.
    :param after_n_years: After how many years to slice single values from every realization. List of year values (
    integers).
    :param months_to_plot: months to pull out of data and plot. Two choices: (1) Single list of month numbers (
    integers) to plot. e.g., [5, 12] or [1,2,3,4,5,6,7,8,9,10,11,12]; (2) Nested list, in which case each nested list
    will appear on the SAME figure, e.g. [[5, 12], [6,7,9], [1,2,3,4]]
    :param combine_months: whether to aggregate the specified month resampling into a single value (e.g., [5,6,7,8,
    9] into one value to represent dry season]. Choices: 1 or 0. If zero then things will produce n separate plots,
    e.g. 12 or 52 or 365.
    :param var_to_eval: name of variable to plot (this is the title/an axis of the panel). if you just feed in a
    dataframe itself, this doesn't need to be specified.
    :param empirical_CDF: 'Yes' or 'No', indicates whether data need to be sorted, rank-ordered, and assigned
    exceedance probabilities to create an empirical CDF.
    '''

    scenarios_to_import = Locations_to_Import.keys()

    # Step 1. Given a dataframe, resample the specified variable at the specified time frame. Store that resampled DF into a new DF.
    if scenarios_to_import is not None:
        Scenario_Dictionary = {}  # keys are scenario names, values are a stack of values. we then directly create a dataframe from this dict.
        new_index = {}  #  used to give numbered indices for purpose of weibull rank ordering in the CDF.
        for scenarios in scenarios_to_import:
            for locations in Locations_to_Import[scenarios]:
                if resample_freq != 'XXXXXXXXXXXXXXXXX':
                    # Create the "resampled" dataframe
                    DF_SUM = Input_Dict[scenarios][locations][var_to_eval].resample(resample_freq, how=resample_stat).copy(deep=True)/units_conversion
                    if resample_freq_orig == 'O':
                        DF_SUM = DF_SUM.tail(1)  # Only use last element if resample is supposed to be once.
                    # give the resampled DF column headings to reflect realization number (number of ensemble member)
                    DF_SUM.columns = ['Realization1' for i in range(int(Num_Realizations[scenarios]))]  # Rename all
                    # columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.
                    # In the sequence below, use DF_SUM to plot sums for months, or new_DF to plot all daily values for each month.
                    column_list_arrays = [DF_SUM.ix[:,i] for i in range(1,Num_Realizations[scenarios])]  # have to skip zero-th element, because it's the first element included in the append call below.
                    if column_list_arrays == []:
                        # Only one realization, so no need to append column_list_arrays to the DF_SUM DF.
                        Scenario_Dictionary[scenarios + ' ' + locations] = DF_SUM.ix[:,0]  # Take the first column of new_DF and append all the other columns to it.
                    else:
                        Scenario_Dictionary[scenarios + ' ' + locations] = DF_SUM.ix[:,0].append(column_list_arrays)  # Take the first column of new_DF and append all the other columns to it.
                elif resample_freq == 'XXXXXXXXXXXXXXXXX':
                    DF_SLICE = Input_Dict[scenarios][locations][var_to_eval].copy(deep=True)/units_conversion
                    DF_SLICE.columns = ['Realization1' for i in range(Num_Realizations[scenarios])]  # Rename all columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.
                    for years in after_n_years:
                        # Loop through years and create new "scenario" entiries in the scenario dictionary.
                        Scenario_Dictionary[locations + ' ' + scenarios + ' after ' + str(years) + ' years'] = DF_SLICE.ix[years*365].copy(deep=True)
            new_index[scenarios] = [i for i in range(int(Num_Realizations[scenarios]*Num_Years[scenarios]))]  # For now scenarios need to have same num years and num realizations or the cdf ranking wont work. so just take most recent value for num_realizations
        # Scenario dictionary has been created. Now create dataframe columns from those keys.
        DF_SUM_STACKED = pd.DataFrame(Scenario_Dictionary)  # Create a dataframe from the dictionary above
        Master_DF = DF_SUM_STACKED.copy(deep=True)
        if resample_freq == 'M':
            DF_SUM_STACKED['Month Number'] = DF_SUM_STACKED.index.month  # Create new column that stores the month value for every time series value (e.g., 1,2,...12)
            Monthly_Energy_DF_GROUPED = DF_SUM_STACKED.groupby('Month Number')  # Create a group that slices the data by the Month Number column
            Monthly_Energy_DF_GROUPED.describe()  # Get statistics for each month
            # Generate a separate figure for every month, and output the monthly dataframe to the user.
            Monthly_Energy_Figs = Monthly_Energy_DF_GROUPED.plot.hist(bins=20, histtype='step', cumulative=False, normed=True)  # Plots a separate histogram for each month
            # If user isn't specifying a particular month, then create a monthly data frame for every scenario.
            if months_to_plot == None:
                month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
                DF_Monthly_Columns = {}
                for scenarios in scenarios_to_import:
                    # For each scenario, create a dataframe that stores monthly info.
                    DF_Monthly_Columns[scenarios] = pd.DataFrame(index=new_index[scenarios])  # Stores a column for each month with all values (Num_Realizations[sims]*(Sim_Dur/36500)) stored in each column, with generic numbered row index.
                    i = 0
                    # Create a subplot with 12 plots (one for each month). Save such a plot for each scenario.
                    for mon in month_list:
                        i+=1
                        # Note: This just copies the VALUES from one data frame (without indices attached) and
                        # stores them in this dataframe, where dates dont really need to be stored.
                        DF_Monthly_Columns[scenarios][mon] = DF_SUM_STACKED[DF_SUM_STACKED.index.month==i][scenarios + ' ' + locations].values
            else:
                # A specific month (or subset of months) was specified.
                # Create a dataframe that stores the values for that month for each scenario.
                Scenario_Dictionary_Monthly = {}  # will create a dataframe from this dict. keys are scenarios.
                for scenarios in scenarios_to_import:
                    for locations in Locations_to_Import[scenarios]:
                        counter = 0
                        if type(months_to_plot[0]) is int:
                            # User fed in a list of month values, so aggregate those months and plot them on one figure.
                            for mon in months_to_plot:
                                if counter == 0:
                                    Scenario_Dictionary_Monthly[scenarios + ' ' + locations] = DF_SUM_STACKED[DF_SUM_STACKED.index.month == mon][scenarios + ' ' + locations]
                                else:
                                    Scenario_Dictionary_Monthly[scenarios + ' ' + locations] = Scenario_Dictionary_Monthly[scenarios + ' ' + locations].append(DF_SUM_STACKED[DF_SUM_STACKED.index.month == mon][scenarios + ' ' + locations])  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
                                    #Scenario_Dictionary_Monthly[scenarios + ' ' + locations] = DF_SUM_STACKED[DF_SUM_STACKED.index.month == mon[i]][scenarios + ' ' + locations].values  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
                                counter += 1
                        elif type(months_to_plot[0]) is list:
                            multi_month = 1
                            ct_3 = 0  # sub-list counter
                            for sub_list in range(len(months_to_plot)):
                                counter = 0
                                ct_3 += 1
                                for mon in months_to_plot[sub_list]:
                                    if counter == 0:
                                        Scenario_Dictionary_Monthly[scenarios + ' ' + locations + ' ' + str(ct_3)] = DF_SUM_STACKED[DF_SUM_STACKED.index.month == mon][scenarios + ' ' + locations]
                                    else:
                                        Scenario_Dictionary_Monthly[scenarios + ' ' + locations + ' ' + str(ct_3)] = Scenario_Dictionary_Monthly[scenarios + ' ' + locations + ' ' + str(ct_3)].append(DF_SUM_STACKED[DF_SUM_STACKED.index.month == mon][scenarios + ' ' + locations])  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
                                        #Scenario_Dictionary_Monthly[scenarios + ' ' + locations] = DF_SUM_STACKED[DF_SUM_STACKED.index.month == mon[i]][scenarios + ' ' + locations].values  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
                                    counter += 1

                                # User fed in a list of lists containing month values, so each month list needs to be a separate "scenario" that is placed on the same figure.
                if multi_month == 1:
                    DF_SUM_STACKED_MONTHLY = pd.DataFrame(index=[i for i in range(len(months_to_plot[0])*Num_Realizations[scenarios]*Num_Years[scenarios])])
                    for j in Scenario_Dictionary_Monthly.keys():
                        # CAUTION: CURRENTLY we can only plot cdfs for items that have the same number of values. So you cant do a 3-month and a 6-month on the same plot.
                        # Should be able to fix this. Just need to plot different x-values for the CDF, and cannot have everything in the same dataframe, which actually wasn't needed in the first place.
                        DF_SUM_STACKED_MONTHLY[j] = pd.DataFrame(Scenario_Dictionary_Monthly[j].values)
                else:
                    DF_SUM_STACKED_MONTHLY = pd.DataFrame(Scenario_Dictionary_Monthly)
                Master_DF = DF_SUM_STACKED_MONTHLY.copy(deep=True)

        if empirical_cdf == 'Yes':
            # Whether it has yearly values or monthly values (or whatever else) in it, now get the master_DF sorted and
            # numbered, and plot CDF using dataframe.
            cop = {}
            rank_ord_cum = {}
            new_index_2 = {}
            ct_2 = 0
            if scenarios_list is not None:
                plot_list = scenarios_list
            else:
                plot_list = copy.deepcopy(Master_DF.columns)
            for scenarios in plot_list:
                # Do the necessary calculations for a probability plot.
                cop[scenarios] = Master_DF[scenarios].copy(deep=True)
                cop[scenarios].sort_values(ascending = True, inplace=True)  # Sort the values.
                rank_ord_cum[scenarios] = pd.DataFrame(cop[scenarios].values)
                new_index_2[scenarios] = rank_ord_cum[scenarios].index + 1
                # Create DF that will store cumulative probability
                rank_ord_cum[scenarios] = rank_ord_cum[scenarios].set_index(new_index_2[scenarios])
                rank_ord_cum[scenarios]['cum. probability'] = rank_ord_cum[scenarios].index/len(rank_ord_cum[scenarios].index)
                rank_ord_cum[scenarios].columns = ['Data', 'cum. probability']

        if type(distribution_stat) in [str, unicode]:
            obj_value = getattr(Master_DF, distribution_stat)()[0]  # Don't need to specify column heading.
        else:
            # A quantile value has been specified. Use built-in pandas quantile function.
            obj_value = Master_DF.quantile(distribution_stat)[0]  # Don't need to specify column heading.

        if empirical_cdf == 'Yes':
            return obj_value, Master_DF, rank_ord_cum
        else:
            return obj_value, Master_DF