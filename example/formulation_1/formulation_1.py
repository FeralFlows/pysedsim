# This module contains the pysedsim function calls required to execute Formulation 1 in Wild et al. (in review).

from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer.
import sys
sys.path.append(r'\\essi12.umd.edu\documents\twild\Documents\Publications\2018\Wild et al. (2018) - EMS\PySedSim')
import pysedsim
import numpy as np
import processing_reference_set
from data_processing import Import_Simulation_Output
from data_processing import Total_Storage_Capacity
from pysedsim_plotting import Single_Time_Series_Value_Plotting
from pysedsim_plotting import Probability_Plot
from matplotlib import pyplot
import copy
import matplotlib.lines as mlines

# Import reference set files from formulation 2, so that the same color scheme and objective can be used for the
# colors of the rule curves from formulation 1.
def determine_rule_curve_colors(formulation_name, base_formulation_name, color_obj_num, color_obj_num_base,
                                policies=None, main_file_loc=None):
    '''

    Determines color array to use for plotting rule curves (formulation I) or policies (formulation II or III),
    based upon a colormap that uses a range from the reference set for a base case.

    :param formulation_name: indicates formulation being plotted. required, string, options: 'formulation_1',
    'formulation_2', 'formulation_3'
    :param base_formulation_name: indicates reference set to be used to define base colormap. required, string,
    :param policies: optional, indicates row numbers (policies) from formulation_name for which colors should be
    returned.
    :param color_obj_num_base: objective number (column number) in reference set file that has no decision variable
    values, for the formulation used to define the base color.
    :param color_obj_num: objective number (column number) in reference set file that has no decision variable
    values, for the formulation being plotted.
    options: 'formulation_1', 'formulation_2', 'formulation_3'.
    :return: color_values: array of values to be used as plotting colors
    '''

    # Load reference set that contains objective that will be used as the base color map
    reference_set_path_no_vars = r'E:\Publications in prep\Paper 2\Results\Output_Storage\Sambor_SA_' + \
                                 base_formulation_name + r'\sets\pysedsim_ref_set_no_vars.ref'
    ref_set_no_vars = np.loadtxt(reference_set_path_no_vars)# reference set for formulation II (use to obtain flow values for coloring)
    ref_set_no_vars = (-1)*ref_set_no_vars  # reformat reference set
    max_q = np.max(ref_set_no_vars[:,color_obj_num_base])  # Max objective (flow) values for objective to be colored, for formulation II
    min_q = np.min(ref_set_no_vars[:,color_obj_num_base])  # Min objective (flow) values for objective to be colored, for formulation II
    if formulation_name == 'formulation_1':
        # Find color values based on performance of three rule curve policies from formulation I
        flow_val_rule_curves = np.zeros(3)  # Array to store flow values
        [Time_Series_Import_Dictionary, RETURNED_master, df3] = formulation_I_histogram(produce_plot = 'No',
                                                                                        main_file_loc = main_file_loc)
        flow_val_rule_curves[0] = Time_Series_Import_Dictionary['rule_curve_1_sedmgmt_stochastic']['Bypass Channel 1'][
            'Q_out'].resample('A', how='mean').mean().mean()  # Mean annual flow value (m3/s), rule curve 1
        flow_val_rule_curves[1] = Time_Series_Import_Dictionary['rule_curve_2_sedmgmt_stochastic']['Bypass Channel 1'][
            'Q_out'].resample('A', how='mean').mean().mean()  # Mean annual flow value (m3/s), rule curve 1
        flow_val_rule_curves[2] = Time_Series_Import_Dictionary['rule_curve_3_sedmgmt_stochastic']['Bypass Channel 1'][
            'Q_out'].resample('A', how='mean').mean().mean()  # Mean annual flow value (m3/s), rule curve 1
        norm_array = [((flow_val_rule_curves[i]-min_q)/(max_q-min_q)) for i in range(len(flow_val_rule_curves))]  # Reset norm array values
    # Load reference set for formulation of interest, so the objective values for the specified policy (row) numbers
    # can be determined
    elif formulation_name in ['formulation_2', 'formulation_3']:
        # Find color values based on performance of provided policy numbers from formulation II or III
        reference_set_path_no_vars = r'E:\Publications in prep\Paper 2\Results\Output_Storage\Sambor_SA_' + \
                                     formulation_name + r'\sets\pysedsim_ref_set_no_vars.ref'
        ref_set_no_vars = np.loadtxt(reference_set_path_no_vars)
        ref_set_no_vars = (-1)*ref_set_no_vars
        norm_array = [((ref_set_no_vars[i,color_obj_num]-min_q)/(max_q-min_q)) for i in policies]  # Reset norm array values
        # Find color values based on performance of three rule curve policies from formulation I

    cmap = pyplot.cm.get_cmap("jet_r")  # Use same color map used in other cases
    color_values = cmap(norm_array)  # Array storing line color values for the three policies
    return color_values

def run_simulations_formulation_I():
    # Run deterministic and stochastic simulations of rule curves 1, 2, 3 and 4
    pysedsim.PySedSim(file_name = 'formulation_1.csv')  # Combines deterministic and stochastic simulations

def formulation_I_histogram(produce_plot = 'Yes', main_file_loc=None, var_sub_list=None):
    # As this is simulation only (no optimization), we directly call the Probability_Plot() function in
    # pysedsim_plotting, rather than using the Reevaluation() function in the processing_reference_set.py module.

    # Specify scenarios (i.e., names of simulations corresponding to folders for which simulation outputs exist)

    Sims_to_Import_1 = ["rule_curve_1_nosedmgmt_deterministic", "rule_curve_1_nosedmgmt_stochastic",
                        "rule_curve_2_nosedmgmt_deterministic", "rule_curve_2_nosedmgmt_stochastic",
                        "rule_curve_3_nosedmgmt_deterministic", "rule_curve_3_nosedmgmt_stochastic"]
    Sims_to_Import_1A = ["rule_curve_1_sedmgmt_deterministic", "rule_curve_1_sedmgmt_stochastic",
                        "rule_curve_2_sedmgmt_deterministic", "rule_curve_2_sedmgmt_stochastic",
                        "rule_curve_3_sedmgmt_deterministic", "rule_curve_3_sedmgmt_stochastic"]
    Sims_to_Import_2 = ["rule_curve_1_nosedmgmt_deterministic", "rule_curve_1_nosedmgmt_stochastic",
                        "rule_curve_2_nosedmgmt_deterministic", "rule_curve_2_nosedmgmt_stochastic",
                        "rule_curve_3_nosedmgmt_deterministic", "rule_curve_3_nosedmgmt_stochastic",
                        "rule_curve_1_sedmgmt_deterministic", "rule_curve_1_sedmgmt_stochastic",
                        "rule_curve_2_sedmgmt_deterministic", "rule_curve_2_sedmgmt_stochastic",
                        "rule_curve_3_sedmgmt_deterministic", "rule_curve_3_sedmgmt_stochastic"]

    #Sims_to_Import = [[Sims_to_Import_1, Sims_to_Import_1, Sims_to_Import_1],
    #                  [Sims_to_Import_1, Sims_to_Import_1A, Sims_to_Import_2]]
    Sims_to_Import = [[Sims_to_Import_1, Sims_to_Import_1, Sims_to_Import_1],
                      [Sims_to_Import_1, Sims_to_Import_1, Sims_to_Import_1],
                      [Sims_to_Import_1A, Sims_to_Import_1A, Sims_to_Import_1A]
                      ]



    # System elements (i.e., reservoir, channel and junction names) for which data should be imported for the specified
    # simulation scenarios.
    Locations_to_Import_1 = {
                             "rule_curve_1_nosedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_1_nosedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_2_nosedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_2_nosedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_3_nosedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_3_nosedmgmt_stochastic": ["Sambor SA"]
                             }

    Locations_to_Import_1A = {
                             "rule_curve_1_sedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_1_sedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_2_sedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_2_sedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_3_sedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_3_sedmgmt_stochastic": ["Sambor SA"]
                             }

    Locations_to_Import_2 = {
                             "rule_curve_1_nosedmgmt_deterministic": ["Bypass Channel 1"],
                             "rule_curve_1_nosedmgmt_stochastic": ["Bypass Channel 1"],
                             "rule_curve_2_nosedmgmt_deterministic": ["Bypass Channel 1"],
                             "rule_curve_2_nosedmgmt_stochastic": ["Bypass Channel 1"],
                             "rule_curve_3_nosedmgmt_deterministic": ["Bypass Channel 1"],
                             "rule_curve_3_nosedmgmt_stochastic": ["Bypass Channel 1"]
                             }

    Locations_to_Import_3 = {
                             "rule_curve_1_nosedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_1_nosedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_2_nosedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_2_nosedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_3_nosedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_3_nosedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_1_sedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_1_sedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_2_sedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_2_sedmgmt_stochastic": ["Sambor SA"],
                             "rule_curve_3_sedmgmt_deterministic": ["Sambor SA"],
                             "rule_curve_3_sedmgmt_stochastic": ["Sambor SA"]
                             }

    Locations_to_Import_4 = {
                             "rule_curve_1_nosedmgmt_deterministic": ["Junction 4"],
                             "rule_curve_1_nosedmgmt_stochastic": ["Junction 4"],
                             "rule_curve_2_nosedmgmt_deterministic": ["Junction 4"],
                             "rule_curve_2_nosedmgmt_stochastic": ["Junction 4"],
                             "rule_curve_3_nosedmgmt_deterministic": ["Junction 4"],
                             "rule_curve_3_nosedmgmt_stochastic": ["Junction 4"],
                             "rule_curve_1_sedmgmt_deterministic": ["Junction 4"],
                             "rule_curve_1_sedmgmt_stochastic": ["Junction 4"],
                             "rule_curve_2_sedmgmt_deterministic": ["Junction 4"],
                             "rule_curve_2_sedmgmt_stochastic": ["Junction 4"],
                             "rule_curve_3_sedmgmt_deterministic": ["Junction 4"],
                             "rule_curve_3_sedmgmt_stochastic": ["Junction 4"]
                             }

    Locations_to_Import_all = {
                             "rule_curve_1_nosedmgmt_deterministic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_1_nosedmgmt_stochastic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_2_nosedmgmt_deterministic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_2_nosedmgmt_stochastic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_3_nosedmgmt_deterministic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_3_nosedmgmt_stochastic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_1_sedmgmt_deterministic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_1_sedmgmt_stochastic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_2_sedmgmt_deterministic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_2_sedmgmt_stochastic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_3_sedmgmt_deterministic": ["Sambor SA", "Junction 4", 'Bypass Channel 1'],
                             "rule_curve_3_sedmgmt_stochastic": ["Sambor SA", "Junction 4", 'Bypass Channel 1']
                             }

    Locations_to_Import = [[Locations_to_Import_1, Locations_to_Import_1, Locations_to_Import_2],
                           [Locations_to_Import_1, Locations_to_Import_1, Locations_to_Import_1],
                           [Locations_to_Import_1A, Locations_to_Import_1A, Locations_to_Import_1A]
                           ]  # _4

    # Time series variables to be imported for the scenarios and system elements specified above.
    if var_sub_list is None:
        # Import default variable names
        var_sub_list = ['Hydropower_avg_MWH', 'larv_surv', 'larv_pass', 'Q_out', 'Settled_mass', 'Settled_mass',
                        'SS_W_out', 'res_flushed_load', 'larv_mass_out_surv_total']  # res_flushed_load
    else:
        var_sub_list = var_sub_list

    # Import all relevant data for all scenarios, system locations and variables listed above
    [Time_Series_Import_Dictionary, Num_Realizations, Num_Years, TSID_key_list] = Import_Simulation_Output(Sims_to_Import_2,
                                                                                                           Locations_to_Import_all,
                                                                                                           var_sub_list,
                                                                                                           main_file_loc)

    # Loop through scenarios, locations, and variables and replace the deterministic ones with a single value so they can
    #  be plotted this way.
    replace_deterministic_with_mean = 'Yes'
    scenario_sub_list = [
                         'rule_curve_1_nosedmgmt_deterministic',
                         'rule_curve_1_sedmgmt_deterministic',
                         'rule_curve_2_nosedmgmt_deterministic',
                         'rule_curve_2_sedmgmt_deterministic',
                         'rule_curve_3_nosedmgmt_deterministic',
                         'rule_curve_3_sedmgmt_deterministic'
                         ]
    # Create a copy of Time_Series_Import_Dictionary to manipulate
    tsid_copy = copy.deepcopy(Time_Series_Import_Dictionary)
    tsid_copy_2 = copy.deepcopy(Time_Series_Import_Dictionary)
    month_num_list = ['6'] #, '6', '7', '8'
    month_num_list_ints = [int(i) for i in month_num_list]
    if replace_deterministic_with_mean == 'Yes':
        for scenario in scenario_sub_list:
            for loc_key in tsid_copy[scenario].keys():
                for var_key in tsid_copy[scenario][loc_key]:
                    # Replace all time series values with the mean value
                    if var_key in ['SS_W_out']:
                        i=int(month_num_list[0])-1
                        for mon in month_num_list:
                            i+=1
                            tsid_copy_2[scenario][loc_key][var_key]['Month'] = tsid_copy_2[scenario][loc_key][var_key][
                                'Realization1'].index.month
                            tsid_copy_2[scenario][loc_key][var_key][str(i)] = tsid_copy_2[scenario][
                                loc_key][var_key]['Realization1'][tsid_copy_2[scenario][loc_key][var_key].index.month==i]
                        mean_value = tsid_copy_2[scenario][loc_key][var_key][month_num_list].mean(axis=1).mean()
                        tsid_copy[scenario][loc_key][var_key]['Realization1'] = mean_value  #put mean daily value in
    #                        tsid_copy_2[scenario][loc_key][var_key]['mon_val'][tsid_copy_2[scenario][loc_key][var_key][
    #                                                                               'mon_val'].index.month==i] = \
    #                            tsid_copy_2[scenario][loc_key][var_key]['Realization1'][tsid_copy_2[scenario][loc_key][
    #                                                                                        var_key][
    #
    # 'Realization1'].index.month==i].values
                            #    tsid_copy_2[scenario][loc_key][var_key].index.month==i]['Realization1']
                            #tsid_copy_2[scenario][loc_key][var_key]['Realization1'].loc[] = tsid_copy_2[scenario][loc_key][
                            #    var_key][
                            #    tsid_copy_2[scenario][loc_key][var_key].index.month==i]['Realization1']

                                #tsid_copy_2[scenario][loc_key][var_key][
                                #'Realization1'].index.month
                            #tsid_copy_2[scenario][loc_key][var_key][mon]  [tsid_copy_2.index.month==i].values
                        # Sediment load discharge is being computed monthly, so the average daily value being placed
                        # needs to be handled differently

                        #mean_5 =
                        #tsid_copy[scenario][loc_key][var_key].loc[
                        #    (tsid_copy[scenario][loc_key][var_key].index.month == 2) & (
                        #    tsid_copy[scenario][loc_key][var_key].index.day == 29)] = 0
                        #month_5_mean = tsid_copy[scenario][loc_key][var_key]['Realization1'].resample
                        #month_1_mean = tsid_copy[scenario][loc_key][var_key]['Realization1'].resample
                        #tsid_copy[scenario][loc_key][var_key]['Realization1'] = tsid_copy[scenario][loc_key][var_key][
                        #    'Realization1'].resample('M', how=sum)
                    else:
                        tsid_copy[scenario][loc_key][var_key]['Realization1'] = tsid_copy[scenario][loc_key][var_key][
                            'Realization1'].mean()
                    # Set any values on day 29 as == 0 for any variables that involve summing (e.g., energy)
                    if var_key in ['Hydropower_avg_MWH']:
                        tsid_copy[scenario][loc_key][var_key].loc[
                            (tsid_copy[scenario][loc_key][var_key].index.month == 2) & (
                            tsid_copy[scenario][loc_key][var_key].index.day == 29)] = 0

    if produce_plot == 'Yes':
        # Plot PDFs of performance of the four rule curve policies across six performance measures.
        # Top 3 PDFs will be un-related to sediment so will show rule curves 1 and 2; bottom 3 will show rule curves 1-4
        figure_storage_loc = r'E:\Publications in prep\Paper 2\Results\Figures and Tables\Figure Formulation 1 histogram'
        num_pm_plots = [3, 3, 3]  # number of performance measures plots to create. Each will be saved as a separate image
        create_superplot = 'Yes'  # 'Yes' if user wants to place all performance measure plots as subplots into a superplot.
        superplot_dimensions = [3, 3]
        superplot_name = 'Wild_et_al_EMS_Fig_7'
        var_to_plot = [['Hydropower_avg_MWH', 'larv_mass_out_surv_total', 'Q_out'],
                       ['Settled_mass', 'Settled_mass', 'SS_W_out'],
                       ['Settled_mass', 'Settled_mass', 'SS_W_out']
                      ]  # res_flushed_load

        save_image_as = [['Annual Energy Production', 'Larvae Flow Fraction', 'Bypass Flow'],
                         ['Annual Sediment Deposition-No Sed Mgmt', 'Annual Stor Cap Loss Perc-No Sed Mgmt', 'Flushed Load-No Sed Mgmt'],
                         ['Annual Sediment Deposition-Sed Mgmt', 'Annual Stor Cap Loss Perc-Sed Mgmt', 'Flushed Load-Sed Mgmt']]

        # First establish generic plotting preferences that will remain the same across the six probability plots.
        hist_preferences = {'bins': 10, 'cumulative': False, 'normed': False}  #'histtype': 'step',
        fig_num_position = [0.05, .9]  # [0.05, .875]
        #plot_params = {'axes.labelsize': 6, 'legend.fontsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
        #               'text.usetex': False, 'figure.figsize': [1.5, 1.5]}  # keys removed: 'font.size': 20,
        plot_params = None
        plot_type = 'CDF'  # 'Histogram'
        legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']

        # Performance Measure 1: Annual energy production PDF across two policies (determ/stochast)
        initial_total_stor_cap = 1280000000  # Sambor SA total reservoir storage capacity (m^3)
        perc_conv = 100
        mass_vol_conv = 1100 # density of compacted sediment

        stor_cap_loss_multiplier = perc_conv/(mass_vol_conv*initial_total_stor_cap)
        units_conversion = [
                           [0.000001, 1, 1],
                           [0.000000001, stor_cap_loss_multiplier, 0.000000001],
                           [0.000000001, stor_cap_loss_multiplier, 0.000000001]
                           ]  # Multipliers to convert units
        x_label = [['$J_{Energy}^{Y,M,A}$: Annual Energy (Terawatt-hours)', '$J_{Larvae}^{Y,M,A}$: Larvae Flow Fraction',
                    '$J_{Fish}^{Y,M,A}$: Bypass Flow Rate (m$\mathregular{^3}\!$/s)'],
                   ['$J_{Trap}^{Y,M,A}$: Annual Trapped Sediment (10$\mathregular{^6}\!$ tons)',
                    '$J_{Stor}^{Y,M,A}$: Annual Storage Capacity Loss (%)',
                    '$J_{Flush}^{Y,M}$: Monthly Sediment Load (10$\mathregular{^6}\!$ tons)'],
                   ['$J_{Sed}^{Y,M,A}$: Annual Trapped Sediment (10$\mathregular{^6}\!$ tons)',
                    '$J_{Stor}^{Y,M,A}$: Annual Storage Capacity Loss (%)',
                    '$J_{Flush}^{Y,M}$: Monthly Sediment Load (10$\mathregular{^6}\!$ tons)']
                   ]

        #color_values = determine_rule_curve_colors('formulation_1', 'formulation_2', 1, 1, main_file_loc =
        # main_file_loc)
        color_values = ['k', 'g', 'b']
        c1 = color_values[0]
        c2 = color_values[1]
        c3 = color_values[2]

        plt_colors = [[c1, c1, c2, c2, c3, c3],
                      [c1, c1, c2, c2, c3, c3],
                      [c1, c1, c2, c2, c3, c3]
                      ]
        #plt_colors = [['k', 'k', 'g', 'g', 'b', 'b'],
        #              ['k', 'k', 'g', 'g', 'b', 'b'],
        #              ['k', 'k', 'g', 'g', 'b', 'b']
        #              ]

        #plot_title = [['Annual Energy Production', 'Larvae Flow Fraction', 'Bypass Flow'],
        #              ['Sediment Title 1', 'Sediment Title 2', 'Sediment Title 3']]  # None
        plot_title = [[None, None, None],
                      [None, None, None],
                      [None, None, None]
                      ]
        #plt_line_style = [['-', '--', '-', '--', '-', '--'], ['-', '--', '-', '--','-', '--','-', '--', '-', '--', '-', '--']]
        plt_line_style = [['-', '--', '-', '--', '-', '--'],
                          ['-', '--', '-', '--', '-', '--'],
                          ['-', '--', '-', '--', '-', '--']
                          ]
        fig_num_label = [
                         ['(A)', '(B)', '(C)'],
                         ['(D)', '(E)', '(F)'],
                         ['(G)', '(H)', '(I)']
                         ]
        other_text_label = [[None, None, None],
                            ['No\nSediment\nManagement', 'No\nSediment\nManagement', 'No\nSediment\nManagement'],
                            ['Sluicing &\nFlushing', 'Sluicing &\nFlushing', 'Sluicing &\nFlushing']
                            ]
        other_text_position = [[None, None, None],
                            [[0.6,0.2], [0.6,0.2], [0.6,0.2]],
                            [[0.7,0.2], [0.7,0.2], [0.7,0.2]]
                            ]
        other_text_color = [[None, None, None],
                            ['k', 'k', 'k'],
                            ['darkorange', 'darkorange', 'darkorange']
                            ]
        resample_frequency = [['A', 'A', 'A'], ['A', 'A', 'M'], ['A', 'A', 'M']]
        resample_how = [[sum, 'mean', 'mean'], [sum, sum, sum], [sum, sum, sum]]
        log_x = [['no','no','no'], ['no','no','yes'], ['no','no','yes']]
        log_y = [['no','no','no'], ['no','no','no'], ['no','no','no']]
        label_y_axis = [['CDF',None,None], ['CDF',None,None], ['CDF',None,None]]
        tick_dict = {}
        tick_dict['x_tick_locs'] = [2000, 4000, 6000, 8000, 10000]
        tick_dict['x_tick_vals'] = [2000, 4000, 6000, 8000, 10000]
        tick_dict_list = [[None, None, tick_dict], [None, None, None], [None, None, None]]
        #months_to_plot = [[None, None, None], [None, None, [5, 6, 7, 8]]]
        months_to_plot = [[None, None, None], [None, None, month_num_list_ints], [None, None, month_num_list_ints]]
        ax_range = [[[[2.5, 5.5], [0, 1.02]], [[0.65, 1.00], [0, 1.02]], [[2000, 10000], [0, 1.02]]],
                    [[[0, 25], [0, 1.02]], [[0, 2], [0, 1.02]], [[0, 200], [0, 1.02]]],
                    [[[0, 25], [0, 1.02]], [[0, 2], [0, 1.02]], [[0, 200], [0, 1.02]]]
                    ]

        # Create main plot, onto which subplots will be placed.
        [main_fig, ax_arr] = pyplot.subplots(superplot_dimensions[0], superplot_dimensions[1], sharey='row')
        main_fig.tight_layout(h_pad=2, w_pad=0.3)

        # Loop through superplot rows and columns, storing each subplot onto the superplot if desired.
        fig_row_counter = 0
        fig_col_counter = 0
        for i in range(superplot_dimensions[0]):
            for j in range(num_pm_plots[i]):
                sp_slot_tuple = (fig_row_counter, fig_col_counter)
                [RETURNED_master, df3, main_fig, ax_arr] = Probability_Plot(tsid_copy,
                                                                            Locations_to_Import[i][j],
                                                                            Num_Realizations,
                                                                            Num_Years, Sims_to_Import=Sims_to_Import[i][j],
                                                                            units_conversion=units_conversion[i][j],
                                                                            var_to_plot=var_to_plot[i][j],
                                                                            resample_frequency=resample_frequency[i][j],
                                                                            resample_how=resample_how[i][j],
                                                                            save_type='png', plot_title=plot_title[i][j],
                                                                            save_image_as=save_image_as[i][j],
                                                                            plot_specs=plot_params, file_loc=figure_storage_loc,
                                                                            plot_type=plot_type,
                                                                            hist_preferences=hist_preferences,
                                                                            plt_colors=plt_colors[i],
                                                                            x_label=x_label[i][j],
                                                                            legend_labels=legend_labels,
                                                                            plot_text=fig_num_label[i][j],
                                                                            plot_text_position=fig_num_position,
                                                                            other_text_position = other_text_position[i][j],
                                                                            other_text_label = other_text_label[i][j],
                                                                            other_text_color = other_text_color[i][j],
                                                                            plt_line_style=plt_line_style[i],
                                                                            fig_subplot_slot=sp_slot_tuple,
                                                                            input_fig=[main_fig, ax_arr], log_x=log_x[i][j],
                                                                            log_y = log_y[i][j], tick_dict = tick_dict_list[i][j],
                                                                            months_to_plot=months_to_plot[i][j],
                                                                            axis_range=ax_range[i][j], label_y_axis =
                                                                            label_y_axis[i][j])
                #main_fig.tight_layout()
                fig_col_counter += 1
            fig_col_counter = 0  # Zero out each time move to new row
            fig_row_counter += 1

        # Create lines that will appear in legend
        line_rule_curve_1_det = mlines.Line2D([], [], color='k', label='Rule Curve 1 (Deterministic): Energy',
                                              linewidth=1, linestyle='-')
        line_rule_curve_1_stoch = mlines.Line2D([], [], color='k', label='Rule Curve 1 (Stochastic): Energy',
                                          linewidth = 1, linestyle='--')
        line_rule_curve_2_det = mlines.Line2D([], [], color='g', label='Rule Curve 2 (Deterministic): Larvae',
                                              linewidth=1, linestyle='-')
        line_rule_curve_2_stoch = mlines.Line2D([], [], color='g', label='Rule Curve 2 (Stochastic): Larvae',
                                          linewidth = 1, linestyle='--')
        line_rule_curve_3_det = mlines.Line2D([], [], color='b', label='Rule Curve 3 (Deterministic): Adult Fish',
                                              linewidth=1, linestyle='-')
        line_rule_curve_3_stoch = mlines.Line2D([], [], color='b', label='Rule Curve 3 (Stochastic): Adult Fish',
                                          linewidth = 1, linestyle='--')

        lgd = ax_arr[2,1].legend(
            handles=[line_rule_curve_1_det, line_rule_curve_1_stoch, line_rule_curve_2_det, line_rule_curve_2_stoch,
                     line_rule_curve_3_det, line_rule_curve_3_stoch], fontsize=8, frameon=True, ncol=3,
            loc='upper center', bbox_to_anchor=(0.5, -0.35))
        frame=lgd.get_frame()
        frame.set_facecolor('0.9')  # 0-1 is gray; larger number is lighter gray
        frame.set_edgecolor('k')

        #main_fig.tight_layout()
        #main_fig.subplots_adjust(top=0.90)
        #pyplot.show()
        #main_fig.tight_layout()
        #main_fig.set_size_inches(6, 6)
        #main_fig.tight_layout()
        #fig.figsize=[5,5]

        #pyplot.subplots_adjust(wspace = 0.8)

        # Save vectorized image (pdf)
        main_fig.savefig(figure_storage_loc + '\\' + superplot_name + '.pdf',bbox_extra_artists=(lgd,), bbox_inches='tight',
                         dpi=1200)
        # Save .png file
        main_fig.savefig(figure_storage_loc + '\\' + superplot_name + '.png',bbox_extra_artists=(lgd,), bbox_inches='tight',
                         dpi=1200)
    if produce_plot == 'Yes':
        return Time_Series_Import_Dictionary, RETURNED_master, df3
    else:
        return Time_Series_Import_Dictionary