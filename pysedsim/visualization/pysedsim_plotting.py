# -*- coding: utf-8 -*-

'''
This module produces probability plots (PDFs or CDFs) of the outputs from monte carlo simulations.

The actual collection of data from monte carlo simulations into time series and/or statistics must be done externall
done by the user. This module simply takes data that have already been processed and creates plots (e.g.,
empirical CDFs) of those processed data.
'''

# Import relevant libraries
from __future__ import division
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from pylab import plot, show, hist, figure, title, suptitle, savefig
import os  # allows export of csv file by joining intended file name with the destination folder's path.
import pandas as pd
import copy
from pysedsim.data_processing.data_processing import Op_Sys_Folder_Operator

def Single_Time_Series_Value_Plotting(Sims_to_Import, Input_DF, after_n_years, save_type = None, save_image_as = None, cum_in=0, normed_in=0):
    # after_n_years = list of number of years after which to grab values. must be integer number of years <= simulation duration. Providing more than one value includes them on the same plot.
    # This function grabs a value at a single period in time for every realization for a given scenario/system element and creates either a histogram or cdf plot.
    # save_image_as: specify a name or nothing will be saved.

    # Establish histogram preferences
    if cum_in == 1:
        cum_in = True
    else:
        cum_in = False
    if normed_in == 1:
        normed_in = True
    else:
        normed_in = False

    Single_val_plot = plt.figure(1)
    for scenarios in Sims_to_Import:
        # for each scenario, loop through number of time slices that are desired and stick on plot.
        for year in after_n_years:
            # Loop adds more than one time slice onto plot if valid.
            Single_val_plot = Input_DF[scenarios].ix[year*365].hist(cumulative=cum_in, normed=normed_in)  # slice values for all realizations after n years (for the variable that the Input_DF represents)
            Single_val_plot.set_xlabel((save_image_as + ' after %s' + ' years') % (str(year)), fontsize=16)
            Single_val_plot.set_ylabel('Frequency', fontsize=16)
            rc('xtick', labelsize=14)
            rc('ytick', labelsize=14)
    if save_type is not None:
        # User wants to save an image
        save_image_as = save_image_as + ' After ' + (str(after_n_years)) + ' Years' + '.' + save_type
        plt.savefig(save_image_as, format=save_type, dpi=1200)


def Probability_Plot(Input_Dict, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=None,
                     units_conversion=1, var_to_plot=None, plot_title=None, after_n_days=None, resample_frequency=None,
                     resample_how=None, months_to_plot=None, combine_months=None, save_type=None, save_image_as=None,
                     plot_specs=None, file_loc=None, cum_in=None, normed_in=None, plot_type='CDF',
                     hist_preferences=None, x_label='', plt_colors=[None], plt_line_style=[None], legend_labels=[None],
                     plot_text=None, plot_text_position=None, scenarios_list=None, axis_range=None,
                     input_fig=None, fig_subplot_slot= None, log_x = None, log_y = None, tick_dict = None):

    # This function creates a CDF plot from a given dataframe and set of scenarios/locations. The purpose is to plot all scenarios/locations on the same plot. If you want separate plots, just call this function separately for each scenario.

    # resample_frequency = how often to resample. Choices: 'A', 'M', and 'Once'. 'A' and 'M' are for annual and monthly, respectively. 'Once' is for a one time slice (e.g., storage cap. after 40 years for all realizations).
    # resample_how = how to resample. Choices: 'sum', 'mean
    # after_n_days = after how many years to slice single values from every realization. List of year values (integers).
    # var_to_plot = name of variable to plot (this is the title/an axis of the panel). if you just feed in a dataframe itself, this doesn't need to be specified.
    # months_to_plot = months to pull out of data and plot. Two choices: (1) Single list of month numbers (integers) to plot. e.g., [5, 12] or [1,2,3,4,5,6,7,8,9,10,11,12]; (2) Nested list, in which case each nested list will appear on the SAME figure, e.g. [[5, 12], [6,7,9], [1,2,3,4]]
    # combine_months = whether to aggregate the specified month resampling into a single value (e.g., [5,6,7,8,9] into one value to represent dry season]. Choices: 1 or 0. If zero then things will produce n separate plots, e.g. 12 or 52 or 365.
    # save_type
    # save_image_as
    # Sims_to_Import
    # Locations_to_Import
    # Input_Dict
    # units_conversion = value to multiply all values in DF by
    # plot_specs = a dictionary of plot specs. Example: plot_params = {'axes.labelsize': 12, 'font.size': 20, 'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'text.usetex': False, 'figure.figsize': [4.5, 4.5]}
    # plot_type = plot type. Choices: 'CDF' or 'Histogram'. Default: 'CDF'.
    # legend_labels = labels to appear in legend. Choices: (1) list of strings that are labels for each scenario to be plotted. (2) "Off" or "off", which means now legend will appear. If provide nothing, default = None.
    # scenarios_list = A list of the scenario names + locations that is exact order in which you want them to be called from the dataframe when plotting. Warning: dataframe column names are set internally in function.

    os_fold = Op_Sys_Folder_Operator()

    # Establish histogram preferences
    if hist_preferences == None:
        cum_in = False
        normed_in = False
        histtype_in = None
        bins_in = None
        stacked = False
        fill = False
    else:
        try:
            cum_in = hist_preferences['cumulative']
        except KeyError:
            cum_in = False
        try:
            normed_in = hist_preferences['normed']
        except KeyError:
            normed_in = False
        try:
            histtype_in = hist_preferences['histtype']
        except KeyError:
            histtype_in = 'bar'
        try:
            bins_in = hist_preferences['bins']
        except KeyError:
            bins_in = 20
        try:
            stacked = hist_preferences['stacked']
        except KeyError:
            stacked = False
        try:
            fill = hist_preferences['fill']
        except KeyError:
            fill = False
    if plot_type == "Histogram":
        y_label = 'Frequency'
    else:
        y_label = 'CDF'
    # Step 1. Given a dataframe, resample the specified variable at the specified time frame. Store that resampled DF into a new DF.

    if input_fig is not None:
        main_fig, ax_arr = input_fig

    if Sims_to_Import is not None:
        Scenario_Dictionary = {}  # keys are scenario names, values are a stack of values. we then directly create a dataframe from this dict.
        new_index = {}  #  used to give numbered indices for purpose of weibull rank ordering in the CDF.
        master_column_list = []
        for scenarios in Sims_to_Import:
            for locations in Locations_to_Import[scenarios]:
                try:
                    if resample_frequency not in ['Once', 'once', 'o', 'O']:
                        DF_SUM = Input_Dict[scenarios][locations][var_to_plot].resample(resample_frequency,how=resample_how).copy(deep=True)*units_conversion
                        DF_SUM.columns = ['Realization1' for i in range(Num_Realizations[scenarios])]  # Rename all columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.
                        # In the sequence below, use DF_SUM to plot sums for months, or new_DF to plot all daily values for each month.
                        column_list_arrays = [DF_SUM.ix[:,i] for i in range(1,Num_Realizations[scenarios])]  # have to skip zero-th element, because it's the first element included in the append call below.
                        master_column_list.append(scenarios + ' ' + locations)
                        Scenario_Dictionary[scenarios + ' ' + locations] = DF_SUM.ix[:,0].append(column_list_arrays)  # Take the first column of new_DF and append all the other columns to it.
                    elif resample_frequency in ['Once', 'once', 'o', 'O']:
                        DF_SLICE = Input_Dict[scenarios][locations][var_to_plot].copy(deep=True)*units_conversion
                        DF_SLICE.columns = ['Realization1' for i in range(Num_Realizations[scenarios])]  # Rename all columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.
                        for days in after_n_days:
                            # Loop through years and create new "scenario" entries in the scenario dictionary.
                            Scenario_Dictionary[locations + ' ' + scenarios + ' after ' + str(days) + ' days'] = DF_SLICE.ix[days].copy(deep=True)
                            master_column_list.append(locations + ' ' + scenarios + ' after ' + str(days) + ' days')
                except KeyError:
                    pass
            new_index[scenarios] = [i for i in range(int(Num_Realizations[scenarios]*Num_Years[scenarios]))]  # For now scenarios need to have same num years and num realizations or the cdf ranking wont work. so just take most recent value for num_realizations
        # Scenario dictionary has been created. Now create dataframe columns from those keys.
        #DF_SUM_STACKED = pd.DataFrame(Scenario_Dictionary)  # Create a dataframe from the dictionary above
        #DF_SUM_STACKED = {scenario: pd.DataFrame(Scenario_Dictionary[scenario]) for scenario in master_column_list}
        DF_SUM_STACKED = {}
        for scenario in master_column_list:
            DF_SUM_STACKED[scenario] = pd.DataFrame(Scenario_Dictionary[scenario])
            DF_SUM_STACKED[scenario].columns = [scenario]
        #Master_DF = DF_SUM_STACKED.copy(deep=True)
        Master_DF = copy.deepcopy(DF_SUM_STACKED)
        if resample_frequency == 'M':
            for scenario in master_column_list:
                DF_SUM_STACKED[scenario]['Month Number'] = DF_SUM_STACKED[scenario].index.month  # Create new column that stores
                # the month value for every time series value (e.g., 1,2,...12)
                Monthly_Energy_DF_GROUPED = DF_SUM_STACKED[scenario].groupby('Month Number')  # Create a group that slices the data by the Month Number column
                Monthly_Energy_DF_GROUPED.describe()  # Get statistics for each month
                # Generate a separate figure for every month, and output the monthly dataframe to the user.
                Monthly_Energy_Figs = Monthly_Energy_DF_GROUPED.plot.hist(bins=20, histtype='step', cumulative=False, normed=True)  # Plots a separate histogram for each month
                plt.close("all")
            # If user isn't specifying a particular month, then create a monthly data frame for every scenario.
            if months_to_plot == None:
                test = 0
                if test == 1:
                    month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
                    DF_Monthly_Columns = {}
                    for scenarios in Sims_to_Import:
                        # For each scenario, create a dataframe that stores monthly info.
                        DF_Monthly_Columns[scenarios] = pd.DataFrame(index=new_index[scenarios])  # Stores a column for each month with all values (Num_Realizations[sims]*(Sim_Dur/36500)) stored in each column, with generic numbered row index.
                        i = 0
                        # Create a subplot with 12 plots (one for each month). Save such a plot for each scenario.
                        Mon_Energ_Fig = plt.figure(1)
                        suptitle(save_image_as)  # This gives your overall figure a title, whereas subplot figures are done below.
                        for mon in month_list:
                            i+=1
                            DF_Monthly_Columns[scenarios][mon] = DF_SUM_STACKED[DF_SUM_STACKED.index.month==i][scenarios + ' ' + locations].values  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
                            #New_Fig.add_subplot()
                            ax1 = plt.subplot(4,3,i)  # Could also make this P.subplot (using pylab; may work better for having a 4x3 subplot)
                            plt.hist(DF_Monthly_Columns[scenarios][mon], bins=30, cumulative=False, normed=True)  # Create monthly histogram.
                            plt.title(mon)  # If you want a title to appear in each subplot, use this.
                            plt.xlabel(x_label)
                            plt.ylabel(y_label)
                        plt.tight_layout()  # Makes everything fit more nicely into the figure area.
                        save_image_as = save_image_as + 'monthly subplot for scenario ' + scenarios + '.' + save_type
                        Mon_Energ_Fig.savefig(save_image_as, format=save_type, dpi=1200)
            else:
                # A specific month (or subset of months) was specified.
                # Create a dataframe that stores the values for that month for each scenario.
                Scenario_Dictionary_Monthly = {}  # will create a dataframe from this dict. keys are scenarios.
                for scenarios in Sims_to_Import:
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

                # Create a subplot with 12 plots (one for each month). Save such a plot for each scenario.

                # UNDER CONSTRUCTION. PRODUCE A HISTOGRAM, FOR ONLY THIS MONTH, THAT HAS ALL THE SCENARIOS ON IT.
                #Mon_Energ_Fig = plt.figure(1)
                #suptitle(save_image_as)  # This gives your overall figure a title, whereas subplot figures are done below.
                #for mon in months_to_plot:
                #    DF_Monthly_Columns[scenarios][mon] = DF_SUM_STACKED[DF_SUM_STACKED.index.month==i][scenario + ' ' + location].values  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
                #    #New_Fig.add_subplot()
                #    ax1 = plt.subplot(4,3,i)  # Could also make this P.subplot (using pylab; may work better for having a 4x3 subplot)
                #    plt.hist(DF_Monthly_Columns[scenarios][mon], bins=30, cumulative=False, normed=True)  # Create monthly histogram.
                #    plt.title(mon)  # If you want a title to appear in each subplot, use this.
                #    plt.xlabel(save_image_as)
                #    plt.ylabel('Frequency (#)')
                #plt.tight_layout()  # Makes everything fit more nicely into the figure area.
                #save_image_as = save_image_as + 'monthly subplot for scenario ' + scenarios + '.' + save_type
                #plt.savefig(save_image_as, format=save_type, dpi=1200)
        # Whether it has yearly values or monthly values (or whatever else) in it, now get the master_DF sorted and numbered, and plot CDF using dataframe.
        cop = {}
        new_df_3 = {}
        new_index_2 = {}
        if plot_specs is not None:
            rcParams.update(plot_specs)
        if input_fig is None:
            main_fig = plt.figure()
            fig_subplot_slot = (1,1)
            ax_arr = main_fig.add_subplot(111)
        ct_2 = 0
        if legend_labels[0] is None:
            # Initialize a legend label list if none was explicitly provided. In this case just use the scenario names as legend entries.
            legend_labels = [None for j in range(len(master_column_list))]
        elif legend_labels[0] in ['Off', 'off']:
            # Initialize a legend label list if none was explicitly provided. In this case just use the scenario names as legend entries.
            legend_labels = ['Off' for j in range(len(master_column_list))]
        else:
            # Use the user-specified legend labels
            legend_labels = legend_labels

        if scenarios_list is not None:
            plot_list = scenarios_list
        else:
            plot_list = copy.deepcopy(master_column_list)
        for scenarios in plot_list:
            # Handle plot line colors and line styles lists.
            ct_2 += 1
            if ct_2 == 1:
                if plt_line_style[ct_2 - 1] is None:
                    plt_line_style = ['-' for j in range(len(master_column_list))]  # Set default line type to solid line.
                if plt_colors[ct_2 - 1] is None:
                    plt_colors = [None for j in range(len(master_column_list))]  # make a list of "Nones", so that matplotlib will use default colors.
            # Do the necessary calculations for a probability plot.
            cop[scenarios] = Master_DF[scenarios][scenarios].copy(deep=True)
            cop[scenarios].sort_values(ascending = True, inplace=True)  # Sort the values.
            new_df_3[scenarios] = pd.DataFrame(cop[scenarios].values)
            new_index_2[scenarios] = new_df_3[scenarios].index + 1
            new_df_3[scenarios] = new_df_3[scenarios].set_index(new_index_2[scenarios])
            new_df_3[scenarios]['cum. probability'] = new_df_3[scenarios].index/len(new_df_3[scenarios].index)
            new_df_3[scenarios].columns = ['Data', 'cum. probability']
            if legend_labels[ct_2-1] is None:
                legend_labels[ct_2-1] = scenarios  # Use the scenario in the legend if no legend labels are specified.
            if plot_type == "CDF":
                if input_fig is not None:
                    if log_x in [None, 'no'] and log_y in [None, 'no']:
                        ax_arr[fig_subplot_slot].plot(new_df_3[scenarios]['Data'],
                                                      new_df_3[scenarios]['cum. probability'],
                                                      label=legend_labels[ct_2 - 1], color=plt_colors[ct_2 - 1],
                                                      linestyle=plt_line_style[ct_2 - 1])
                    elif log_x in ['yes'] and log_y in ['yes']:
                        ax_arr[fig_subplot_slot].loglog(new_df_3[scenarios]['Data'],
                                                      new_df_3[scenarios]['cum. probability'],
                                                      label=legend_labels[ct_2 - 1], color=plt_colors[ct_2 - 1],
                                                      linestyle=plt_line_style[ct_2 - 1])
                    elif log_x in ['yes']:
                        ax_arr[fig_subplot_slot].semilogx(new_df_3[scenarios]['Data'],
                                                          new_df_3[scenarios]['cum. probability'],
                                                          label=legend_labels[ct_2 - 1], color=plt_colors[ct_2 - 1],
                                                          linestyle=plt_line_style[ct_2 - 1])
                    elif log_y in ['yes']:
                        ax_arr[fig_subplot_slot].semilogy(new_df_3[scenarios]['Data'],
                                                          new_df_3[scenarios]['cum. probability'],
                                                          label=legend_labels[ct_2 - 1], color=plt_colors[ct_2 - 1],
                                                          linestyle=plt_line_style[ct_2 - 1])

                    #                else:
#                    ax.plot(new_df_3[scenarios]['Data'], new_df_3[scenarios]['cum. probability'],
#                            label=legend_labels[ct_2 - 1], color=plt_colors[ct_2 - 1],
#                            linestyle=plt_line_style[ct_2 - 1])
            elif plot_type == "Histogram":
                #plt.hist(new_df_3[scenarios]['Data'], label=scenarios, cumulative=cum_in, normed=normed_in,
                # bins=bins_in, histtype=histtype_in)
                new_df_3[scenarios]['Data'].hist(label=legend_labels[ct_2-1], cumulative=cum_in, normed=normed_in,
                                                 bins=bins_in, histtype=histtype_in, color = plt_colors[ct_2-1],
                                                 stacked=stacked, fill=fill)
#                new_df_3[scenarios]['Data'].hist()
                #ax.plot(new_df_3[scenarios]['Data'].hist(label=legend_labels[ct_2-1], cumulative=cum_in, normed=normed_in, bins=bins_in, histtype=histtype_in, color = plt_colors[ct_2-1], linestyle = plt_line_style[ct_2-1]))
        if 'off' in legend_labels or 'Off' in legend_labels:
            pass  # Don't turn on legend.
        else:
            # To create a legend box:

            # box = ax.get_position()
            #ax.set_position([box.x0, box.y0 - box.height * 0.2, box.width, box.height * 0.8])

            # To create a mini legend at the top:
            # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), fancybox=True, shadow=True, ncol=3)
            plt.legend(loc='best')  # may want this instead for histogram
        ax_arr[fig_subplot_slot].set_ylabel(y_label, fontsize = 8)
        ax_arr[fig_subplot_slot].set_xlabel(x_label, fontsize = 8)
        if axis_range is not None:
            ax_arr[fig_subplot_slot].set_xlim(axis_range[0])
            ax_arr[fig_subplot_slot].set_ylim(axis_range[1])
        # Set x tick locations and labels
        if tick_dict is not None:
            # User wants to adjust tick parameters
            try:
                ax_arr[fig_subplot_slot].set_xticks(tick_dict['x_tick_locs'])
            except KeyError:
                pass
            try:
                x_tick_vals = tick_dict['x_tick_vals']
                ax_arr[fig_subplot_slot].set_xticklabels(x_tick_vals, fontsize=10)
            except KeyError:
                pass
            # Set y tick locations and labels
            try:
                ax_arr[fig_subplot_slot].set_yticks(tick_dict['y_tick_locs'])
            except KeyError:
                pass
            try:
                y_tick_vals = tick_dict['y_tick_vals']
                ax_arr[fig_subplot_slot].set_yticklabels(y_tick_vals, fontsize=10)
            except KeyError:
                pass

        if plot_title is not None:
             ax_arr[fig_subplot_slot].set_title(plot_title, fontsize = 8)
        else:
            pass
        #ax_arr[fig_subplot_slot].set_xticks(np.arange(0,np.shape(table)[1],1))
        #ax_arr[fig_subplot_slot].set_xticklabels(fontsize=8)
        #ax_arr[fig_subplot_slot].tick_params(axis='both', which='major', labelsize=10)
        ax_arr[fig_subplot_slot].tick_params(axis='y', which='both', right='off', labelsize=8)
        ax_arr[fig_subplot_slot].tick_params(axis='x', which='both', top='off', labelsize=8)
        if plot_text is not None:
            ax_arr[fig_subplot_slot].text(plot_text_position[0], plot_text_position[1], plot_text,
                     transform=ax_arr[fig_subplot_slot].transAxes, color='k', fontsize=11,
                     horizontalalignment='left', verticalalignment='center')  # weight="bold",
        # save_image_as = save_image_as + ', for scenario ' + scenarios + '.' + save_type
        if file_loc is not None:
            save_image_as = file_loc + os_fold + save_image_as + '.' + save_type
            if not os.path.exists(file_loc):
                os.makedirs(file_loc)  # Create new folder if it doesn't exist.
            main_fig.savefig(save_image_as, format=save_type, dpi=200, bbox_inches='tight')  #bbox trims white space.
        else:
            save_image_as += '.' + save_type
            savefig(save_image_as, format=save_type, dpi=1200, bbox_inches='tight')  #bbox trims off white space.
            #savefig(os.path.join(file_loc, save_image_as), format=save_type, dpi=1200, bbox_inches='tight')
        plt.close("all")  # Close figures in case you're in a loop and need to restart.
    if input_fig is None:
        return Master_DF, new_df_3
    else:
        return Master_DF, new_df_3, main_fig, ax_arr