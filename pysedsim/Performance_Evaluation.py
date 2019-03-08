# -*- coding: utf-8 -*-

'''
This is an example of a script that uses the Performance_Measure_Analysis method to generate plots of outputs from a
simulation(s).
'''

# Import relevant libraries:
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0).
from data_processing import Import_Simulation_Output
from data_processing import Total_Storage_Capacity
from pysedsim_plotting import Single_Time_Series_Value_Plotting
from pysedsim_plotting import Probability_Plot
import matplotlib.pyplot as plt

# Define information relevant for this particular performance measure evaluation

#Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]

Sims_to_Import = ["Alternative 7A with Flushing"]

#Locations_to_Import = {"Alternative 7A with Flushing": ["Sambor Alternative", 'Bypass Struct 1'],
#                       "Alternative No Flushing": ["Sambor Alternative", 'Bypass Struct 1'],
#                       "Sambor Original": ['Sambor']}

Locations_to_Import = {"Alternative 7A with Flushing": ["Sambor Alternative 7A"]}


#var_sub_list = ['BS_W', 'SS_W_in', 'Q_in', 'Q_out', 'SS_W_in', 'SS_W_out', 'water_surface_elevation', 'Egg_Passage_Velocity_Success',
# 'Bypass_Flow_Success', 'bypass_flow_resilience_tracker', 'theoretical_peaking_capacity', 'TE_avg', 'capacity_active_reservoir', 'capacity_dead_reservoir', 'Hydropower_avg_MWH', 'S']

var_sub_list = ['BS_W', 'water_surface_elevation', 'capacity_active_reservoir','capacity_dead_reservoir', 'theoretical_peaking_capacity',
                'Hydropower_avg_MWH', 'Egg_Passage_Velocity_Success', 'Bypass_Flow_Success', 'bypass_flow_resilience_tracker']

var_sub_list = ['Hydropower_avg_MWH', 'egg_pass', 'Q_in']
main_file_loc = r'E:\PySedSim\Model Files\PySedSim Files\PySedSim_Git_Repository\Output_Storage'  # For work computer  # r'C:\Users\Thomas Wild\Documents\PySedSim\12 9 15\Output Storage'  # For laptop

# Import all relevant data for scenarios, system locations and variables listed above
[Time_Series_Import_Dictionary, Num_Realizations, Num_Years] = Import_Simulation_Output(Sims_to_Import, Locations_to_Import, var_sub_list, main_file_loc)

# Determine total storage capacity given the dataframes for active and dead storage capacity
Time_Series_Import_Dictionary = Total_Storage_Capacity(Sims_to_Import, Locations_to_Import, Time_Series_Import_Dictionary)  # Modify the time series import dictionary so it includes metrics related to storage capacity loss for the specified scenarios/reservoirs.

# Establish generic plotting preferences. Establish them in the code for each individual plot if specific preferences vary.
hist_preferences = {'bins': 10, 'cumulative': False, 'normed': False}  #'histtype': 'step',

# 1. Sediment

# a. Total 100-year accumulation

# b. Storage capacity loss (%)

# Plot storage capacity at T = Concession Period time by grabbing last row. Note: this creates column, better for plotting using DF.plot. To create row, use aa colon in brackets: Storage_Capacity_DF.ix[-1:]
var_plotted = 'Total Storage Capacity Loss (Percent of Original)'
after_n_years = [40]  # must be an integer number of years
Single_Time_Series_Value_Plotting(Sims_to_Import, Storage_Capacity_Fract_Remaining, after_n_years, save_type = 'png', save_image_as = var_plotted, cum_in=0, normed_in=0)

after_n_years = [100]
Single_Time_Series_Value_Plotting(Sims_to_Import, Storage_Capacity_Fract_Remaining, after_n_years, save_type = 'png', save_image_as = var_plotted, cum_in=0, normed_in=0)

after_n_years = [40, 100]
Single_Time_Series_Value_Plotting(Sims_to_Import, Storage_Capacity_Fract_Remaining, after_n_years, save_type = 'png', save_image_as = var_plotted, cum_in=0, normed_in=0)

Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"], "Sambor Original": ['Sambor']}

# Figures capturing loss in total reservoir storage capacity:
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"], "Sambor Original": ['Sambor']}
scenarios_to_plot = 3  # Number of figures desired to be created.
var_to_plot = ['Storage_Capacity_Loss_Percent', 'Storage_Capacity_Loss', 'Storage_Capacity_Total_Remaining']
save_image_as = ['Relative Stor. Cap. Loss CDF', 'Total Stor. Cap. Loss CDF', 'Total Stor. Cap. Remaining (mil. m3) CDF']
x_label = ['Storage Capacity Loss (% of Initial)', 'Total Storage Capacity Loss ($\mathregular{10^9}$$\mathregular{m^3}$)', 'Total Storage Capacity Remaining ($\mathregular{10^9}$$\mathregular{m^3}$)']
plt_colors = [['g', 'g','darkorange', 'darkorange', 'k', 'k'], ['g', 'g', 'darkorange', 'darkorange', 'k', 'k'], ['g', 'g','darkorange',
                                                                                                                  'darkorange', 'k', 'k']]
plt_line_style = [['-', '--', '-', '--','-', '--',], ['-', '--', '-', '--','-', '--',], ['-', '--', '-', '--','-', '--',]]
legend_labels = ['Off']  # ['Alternative with Flushing: after 100 years', 'Alternative with Flushing: after 40 years','Alternative without Flushing: after 100 years', 'Alternative without Flushing: after 40 years','Sambor Original: after 100 years', 'Sambor Original: after 40 years']  # ['Off']
units_conversion = [1, 1E9, 1E9]
#plot_title = ['Relative (%) Storage Capacity Loss over 40 and 100 years', 'Total Storage Capacity Loss (mil. $m^3$) over 40 and 100 years', 'Total Storage Capacity Remaining (mil. m3) over 40 and 100 years']
plot_title = [None, None, None]
plot_text = ['(B)', '(C)', '(D)']
plot_text_position = [.9, 0.05]
after_n_years = [40, 100]
main_file_loc =  r'E:\PySedSim\Model Files\PySedSim Files\PySedSim_Git_Repository\Output_Storage'  # For work computer
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [6.5, 5]}  # keys removed: 'font.size': 20,
#plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False}  # REMOVE figure size if you want the entire legend to be spread across the top.
RETURNED_master = {}
df3 = {}
plot_type = 'CDF' # 'Histogram'
for i in range(scenarios_to_plot):
    [RETURNED_master[var_to_plot[i]], df3[var_to_plot[i]]] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, var_to_plot=var_to_plot[i], resample_frequency = 'Once', units_conversion=units_conversion[i], after_n_years = after_n_years, save_type = 'png', plot_title=plot_title[i], save_image_as = save_image_as[i], plot_specs = plot_params, file_loc = main_file_loc, plot_type = plot_type, hist_preferences = hist_preferences, x_label = x_label[i], plt_colors = plt_colors[i], plt_line_style = plt_line_style[i], legend_labels = legend_labels, plot_text = plot_text[i], plot_text_position=plot_text_position)

# Annual Sediment Load Trapped
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"], "Sambor Original": ['Sambor']}
var_to_plot = 'BS_W'
after_n_years = [100]
main_file_loc =  r'E:\PySedSim\Model Files\PySedSim Files\PySedSim_Git_Repository\Output_Storage'  # For work computer
scenarios_to_plot = 3  # Number of figures desired to be created.
x_label = r'Annual Trapped Sediment Load ($\mathregular{10^6}$$\mathregular{\frac{tonnes}{yr}}$)'
legend_labels = ['Off']
plot_title = None
plt_colors = ['g', 'darkorange', 'k']
save_image_as = 'Annual Sediment Deposition (Mt-yr)'
units_conversion = (Num_Years[Num_Years.keys()[0]])*1000*1000000  # Converts load into millions of tons/yr
plot_text = '(A)'
plot_text_position = [.9, 0.05]
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [6.5, 5]}  # keys removed: 'font.size': 20,
plot_type = 'CDF'
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'Once', save_type = 'png', after_n_years = after_n_years, plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc = main_file_loc, plot_type = plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position)

# Monthly Energy Production CDF--(1) All months in year; (2) dry season; (3) wet season; All separate plots.
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"], "Sambor Original": ['Sambor']}
scenarios_to_plot = 3  # Number of figures desired to be created.
var_to_plot = 'Hydropower_avg_MWH'
save_image_as = ['Monthly Energy Production (TWh)', 'Dry Season Monthly Energy Production (TWh)', 'Wet Season Monthly Energy Production (TWh)']
x_label = ['Monthly Energy Production (TWh)', 'Dry Season Monthly Energy Production (TWh)', 'Wet Season Monthly Energy Production (TWh)']
plt_colors = [['darkorange', 'g', 'k'], ['darkorange', 'g', 'k'], ['darkorange', 'g', 'k']]
legend_labels = ['Off']  # ['Alternative with Flushing: after 100 years', 'Alternative with Flushing: after 40 years','Alternative without Flushing: after 100 years', 'Alternative without Flushing: after 40 years','Sambor Original: after 100 years', 'Sambor Original: after 40 years']  # ['Off']
plot_title = [None, None, None]  # ['Monthly Energy Production (over the year)', 'Dry Season Monthly Energy Production', 'Wet Season Monthly Energy Production']
plot_text = ['(X)', '(B)', '(C)']
plot_text_position = [.9, 0.05]
units_conversion = 1000000  # Converts MWh to GWh
months_to_plot = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 12], [6, 7, 8, 9, 10, 11]]
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [6.5, 5]}  # keys removed: 'font.size': 20,
RETURNED_master = {}
df3 = {}
plot_type = 'CDF'
for i in range(scenarios_to_plot):
    [RETURNED_master[i], df3[i]] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = sum, months_to_plot = months_to_plot[i], save_type = 'png', plot_title=plot_title[i], save_image_as = save_image_as[i], plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors[i], x_label = x_label[i], legend_labels = legend_labels, plot_text = plot_text[i], plot_text_position=plot_text_position)

# Monthly Energy Production CDF--(1) Dry season; (2) wet season; All on the SAME plot.
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"], "Sambor Original": ['Sambor']}
scenarios_to_plot = 1  # Number of figures desired to be created.
scenarios_list = ['Alternative with Flushing Sambor Alternative 1', 'Alternative with Flushing Sambor Alternative 2', 'Alternative No Flushing Sambor Alternative 1', 'Alternative No Flushing Sambor Alternative 2', 'Sambor Original Sambor 1', 'Sambor Original Sambor 2']# used if you need to have things plotted in an exact order so the markers are correct.
var_to_plot = 'Hydropower_avg_MWH'
save_image_as = 'Monthly Energy Production'
x_label = 'Energy (Terawatt-hours)'
plt_colors = ['darkorange', 'darkorange', 'g', 'g', 'k', 'k']
plot_title = 'Monthly Energy Production'  # None  # ['Monthly Energy Production (over the year)', 'Dry Season Monthly Energy Production', 'Wet Season Monthly Energy Production']
plot_text = '(B)'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
plt_line_style = ['-', '--', '-', '--','-', '--']
plot_text_position = [.875, 0.05]
units_conversion = 1000000  # Converts MWh to GWh
months_to_plot = [[6, 7, 8, 9, 10, 11], [1, 2, 3, 4, 5, 12]]  # [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [5, 5]}  # keys removed: 'font.size': 20,
RETURNED_master = {}
df3 = {}
plot_type = 'CDF'
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, scenarios_list=scenarios_list, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = sum, save_type = 'png', months_to_plot = months_to_plot, plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position, plt_line_style=plt_line_style)

# Annual Energy Production CDF
#Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]
#Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor
# Alternative"], "Sambor Original": ['Sambor']}
var_to_plot = 'Hydropower_avg_MWH'
save_image_as = 'Annual Energy Production'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
units_conversion = 1000000  # Converts MWh to GWh
x_label = 'Energy (Terawatt-hours)'
plt_colors = ['darkorange', 'g', 'k']
plot_title = 'Annual Energy Production' # None
plot_text = '(A)'
plot_text_position = [.875, 0.05]
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [5, 5]}  # keys removed: 'font.size': 20,
plot_type = 'Histogram'  # 'CDF'
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'A', resample_how = sum, save_type = 'png', plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position)

# Seasonal Theoretical Peaking Capacity CDF (Seasons on DIFFERENT plots)
scenarios_to_plot = 3  # Number of figures desired to be created.
var_to_plot = 'theoretical_peaking_capacity'
save_image_as = ['Peaking Maximum Flow Duration - Annual', 'Dry Season Peaking Maximum Flow Duration', 'Wet Season  Peaking Maximum Flow Duration']
months_to_plot = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 12], [6, 7, 8, 9, 10, 11]]
plot_title = ['Peaking Maximum Flow Duration (hours)', 'Dry Season Peaking Maximum Flow Duration (hours)', 'Wet Season Peaking Maximum Flow Duration (hours)']
plot_params = {'axes.labelsize': 12, 'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'text.usetex': False}  # keys removed: 'font.size': 20, 'figure.figsize': [4.5, 4.5]
RETURNED_master = {}
df3 = {}
plot_type = 'CDF'
for i in range(scenarios_to_plot):
    [RETURNED_master[i], df3[i]] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', months_to_plot = months_to_plot[i], save_type = 'png', plot_title=plot_title[i], save_image_as = save_image_as[i], plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences)

# Theoretical Peaking Capacity CDF (Seasons ALL on SAME plot)
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing", "Sambor Original"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"], "Sambor Original": ['Sambor']}
scenarios_to_plot = 1  # Number of figures desired to be created.
scenarios_list = ['Alternative with Flushing Sambor Alternative 1', 'Alternative with Flushing Sambor Alternative 2', 'Alternative No Flushing Sambor Alternative 1', 'Alternative No Flushing Sambor Alternative 2', 'Sambor Original Sambor 1', 'Sambor Original Sambor 2']# used if you need to have things plotted in an exact order so the markers are correct.
var_to_plot = 'theoretical_peaking_capacity'
save_image_as = 'Maximum Peaking Flow Duration'
x_label = 'Duration (hours)'
months_to_plot = [[6, 7, 8, 9, 10, 11], [1, 2, 3, 4, 5, 12]]  # [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [5, 5]}  # keys removed: 'font.size': 20,
plot_title = 'Daily Max. Peaking Duration' # None
plt_colors = ['darkorange', 'darkorange', 'g', 'g', 'k', 'k']
plot_text = '(C)'
plt_line_style = ['-', '--', '-', '--','-', '--']
RETURNED_master = {}
df3 = {}
plot_type = 'CDF'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
plot_text_position = [.875, 0.05]
units_conversion = 1  # Converts MWh to GWh
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, scenarios_list=scenarios_list, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', save_type = 'png', months_to_plot = months_to_plot, plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position, plt_line_style=plt_line_style)


# Qin
var_to_plot = 'Q_in'
save_image_as = 'Daily Inflow'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
units_conversion = 1  # Converts MWh to GWh
x_label = 'Qin (m3/s)'
plot_title = 'Daily Inflow' # None
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
               'text.usetex': False, 'figure.figsize': [5, 5]}  # keys removed: 'font.size': 20,
plot_type = 'Histogram'  # 'CDF'
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations,
                                          Num_Years, Sims_to_Import=Sims_to_Import, units_conversion =
                                          units_conversion, var_to_plot=var_to_plot, resample_frequency = 'A',
                                          resample_how = 'mean', save_type = 'png', plot_title=plot_title,
                                          save_image_as = save_image_as, plot_specs = plot_params,
                                          file_loc=main_file_loc, plot_type=plot_type, hist_preferences =
                                          hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels
                                          = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position)


# Egg Passage

# a. EH1
var_to_plot = 'egg_pass'
save_image_as = 'Annual Larvae Passage'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
units_conversion = 1  # Converts MWh to GWh
x_label = 'Energy (Terawatt-hours)'
plt_colors = ['darkorange', 'g', 'k']
plot_title = 'Annual Larvae Passage (%)' # None
plot_text = '(A)'
plot_text_position = [.875, 0.05]
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [5, 5]}  # keys removed: 'font.size': 20,
plot_type = 'Histogram'  # 'CDF'
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations,
                                          Num_Years, Sims_to_Import=Sims_to_Import, units_conversion =
                                          units_conversion, var_to_plot=var_to_plot, resample_frequency = 'A',
                                          resample_how = 'mean', save_type = 'png', plot_title=plot_title,
                                          save_image_as = save_image_as, plot_specs = plot_params,
                                          file_loc=main_file_loc, plot_type=plot_type, hist_preferences =
                                          hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels
                                          = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position)

# a. Reliability of meeting velocity target--Seasons on DIFFERENT plots
scenarios_to_plot = 3  # Number of figures desired to be created.
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing"]  # Exclude Original dam, which has no egg passage.
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"], "Alternative No Flushing": ["Sambor Alternative"]}
var_to_plot = 'Egg_Passage_Velocity_Success'
save_image_as = ['Annual Egg Passage Minimum Velocity Reliability', 'Dry Season Egg Passage Minimum Velocity Reliability', 'Wet Season Egg Passage Minimum Velocity Reliability']
months_to_plot = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 12], [6, 7, 8, 9, 10, 11]]
plot_title = ['Annual Egg Passage Minimum Velocity Reliability', 'Dry Season Egg Passage Minimum Velocity Reliability', 'Wet Season Egg Passage Minimum Velocity Reliability']
plot_params = {'axes.labelsize': 12, 'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'text.usetex': False}  # keys removed: 'font.size': 20, 'figure.figsize': [4.5, 4.5]
RETURNED_master = {}
df3 = {}
plot_type = 'CDF'
for i in range(scenarios_to_plot):
    [RETURNED_master[i], df3[i]] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', months_to_plot = months_to_plot[i], save_type = 'png', plot_title=plot_title[i], save_image_as = save_image_as[i], plot_specs = plot_params, file_loc = main_file_loc, plot_type = plot_type, hist_preferences = hist_preferences)

# a. Reliability of meeting velocity target--Seasons on the SAME plots
Sims_to_Import = ["Alternative with Flushing"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"]}
scenarios_to_plot = 1  # Number of figures desired to be created.
scenarios_list = ['Alternative with Flushing Sambor Alternative 1', 'Alternative with Flushing Sambor Alternative 2']
var_to_plot = 'Egg_Passage_Velocity_Success'
save_image_as = 'Annual Egg Passage Velocity Reliability'
#x_label = 'Mean Monthly Reliability (%):\n Reservoir Velocity Target'
x_label = 'Mean Monthly Reliability (%)'
months_to_plot = [[6, 7, 8, 9, 10, 11], [1, 2, 3, 4, 5, 12]]
plot_title = 'Reliability: Reservoir Velocity Target' # None  # ['Annual Egg Passage Minimum Velocity Reliability', 'Dry Season Egg Passage Minimum Velocity Reliability', 'Wet Season Egg Passage Minimum Velocity Reliability']
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [6.5, 5]}  # keys removed: 'font.size': 20,
plt_colors = ['b', 'r']
axis_range = [-1,101,0,1]  # To make the vertical lines on the CDF visible.
plt_line_style = ['-', '-']
plot_text = '(B)'
plot_type = 'CDF'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
plot_text_position = [.875, 0.05]
units_conversion = .01
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, scenarios_list=scenarios_list, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', save_type = 'png', months_to_plot = months_to_plot, plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position, plt_line_style=plt_line_style, axis_range = axis_range)

# a. THIS IS JUST A TEST. Reliability of meeting velocity target.
Sims_to_Import = ["Alternative with Flushing"]
Locations_to_Import = {"Alternative with Flushing": ["Sambor Alternative"]}
scenarios_to_plot = 1  # Number of figures desired to be created.
scenarios_list = ['Alternative with Flushing Sambor Alternative 1', 'Alternative with Flushing Sambor Alternative 2']
var_to_plot = 'Q_in'  #'Egg_Passage_Velocity_Success'
save_image_as = 'Annual Egg Passage Velocity Reliability'
#x_label = 'Mean Monthly Reliability (%):\n Reservoir Velocity Target'
x_label = 'Mean Monthly Reliability (%)'
months_to_plot = [[6, 7, 8, 9, 10, 11], [1, 2, 3, 4, 5, 12]]
plot_title = 'Reliability: Reservoir Velocity Target' # None  # ['Annual Egg Passage Minimum Velocity Reliability', 'Dry Season Egg Passage Minimum Velocity Reliability', 'Wet Season Egg Passage Minimum Velocity Reliability']
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [6.5, 5]}  # keys removed: 'font.size': 20,
plt_colors = ['b', 'r']
axis_range = [-1,101,0,1]  # To make the vertical lines on the CDF visible.
plt_line_style = ['-', '-']
plot_text = '(B)'
plot_type = 'Histogram'
legend_labels = ['Off'] #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing: Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']
plot_text_position = [.875, 0.05]
units_conversion = 1
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, scenarios_list=scenarios_list, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', save_type = 'png', months_to_plot = months_to_plot, plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position, plt_line_style=plt_line_style, axis_range = axis_range)


# Fish Bypass -- different seasons on different plots

# a. Reliability of meeting bypass flow target
scenarios_to_plot = 3  # Number of figures desired to be created.
Sims_to_Import = ["Alternative with Flushing", "Alternative No Flushing"]  # Exclude Original dam, which has no fish bypass.
Locations_to_Import = {"Alternative with Flushing": ["Bypass Struct 1"], "Alternative No Flushing": ["Bypass Struct 1"]}
var_to_plot = 'Bypass_Flow_Success'
save_image_as = ['Annual Fish Pass Flow Target Reliability', 'Dry Season Fish Pass Flow Target Reliability', 'Wet Season Fish Pass Flow Target Reliability']
months_to_plot = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 12], [6, 7, 8, 9, 10, 11]]
plot_title = ['Annual Fish Pass Flow Target Reliability', 'Dry Season Fish Pass Flow Target Reliability', 'Wet Season Fish Pass Flow Target Reliability']
plot_params = {'axes.labelsize': 12, 'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'text.usetex': False}  # keys removed: 'font.size': 20, 'figure.figsize': [4.5, 4.5]
RETURNED_master = {}
df3 = {}
plot_type = 'CDF'
for i in range(scenarios_to_plot):
    [RETURNED_master[i], df3[i]] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', months_to_plot = months_to_plot[i], save_type = 'png', plot_title=plot_title[i], save_image_as = save_image_as[i], plot_specs = plot_params, file_loc = main_file_loc, plot_type = plot_type, hist_preferences = hist_preferences)

# a. Reliability of meeting bypass flow target -- different seasons on the SAME plot.
scenarios_to_plot = 1  # Number of figures desired to be created.
Sims_to_Import = ["Alternative with Flushing"]
Locations_to_Import = {"Alternative with Flushing": ["Bypass Struct 1"]}
scenarios_list = ['Alternative with Flushing Bypass Struct 1 1', 'Alternative with Flushing Bypass Struct 1 2']
var_to_plot = 'Bypass_Flow_Success'
save_image_as = 'Annual Fish Pass Flow Target Reliability'
#x_label = 'Mean Monthly Reliability (%):\nBypass Flow Target'
x_label = 'Mean Monthly Reliability (%)'
months_to_plot = [[6, 7, 8, 9, 10, 11], [1, 2, 3, 4, 5, 12]]
plot_title = 'Reliability: Bypass Flow Target' # None  # ['Annual Egg Passage Minimum Velocity Reliability', 'Dry Season Egg Passage Minimum Velocity Reliability', 'Wet Season Egg Passage Minimum Velocity Reliability']
plot_params = {'axes.labelsize': 18, 'legend.fontsize': 22, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'text.usetex': False, 'figure.figsize': [6.5, 5]}  # keys removed: 'font.size': 20,
plt_colors = ['b', 'r']
axis_range = [-1,101,0,1]  # To make the vertical lines on the CDF visible.
plt_line_style = ['-', '-']
plot_text = '(A)'
plot_type = 'CDF'
legend_labels = ['Off'] #  ['Alternative: Wet Season', 'Alternative: Dry Season']
plot_text_position = [.04, 0.95]
units_conversion = .01
[RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import, Num_Realizations, Num_Years, Sims_to_Import=Sims_to_Import, scenarios_list=scenarios_list, units_conversion = units_conversion, var_to_plot=var_to_plot, resample_frequency = 'M', resample_how = 'mean', save_type = 'png', months_to_plot = months_to_plot, plot_title=plot_title, save_image_as = save_image_as, plot_specs = plot_params, file_loc=main_file_loc, plot_type=plot_type, hist_preferences = hist_preferences, plt_colors=plt_colors, x_label = x_label, legend_labels = legend_labels, plot_text = plot_text, plot_text_position=plot_text_position, plt_line_style=plt_line_style, axis_range = axis_range)





# Object (System Element) to plot
Sys_Obj = ['Sambor Alternative']  # Plots will be produced for this system object/element.
var_name = 'SS_W_in'  # State variable to be plotted
Sim_name = Sims_to_Import

# 1. Sediment

# Daily sediment inflow - time series
Daily_BS_W_Fig = plt.figure(1)
Daily_BS_W_Fig = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]]['BS_W'].plot(legend=None)
Daily_BS_W_Fig.set_xlabel('Time', fontsize=18)
Daily_BS_W_Fig.set_ylabel('Deposited Sediment Mass (kg)', fontsize=18)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
plt.savefig('Bottom sediment time series.png', format='png', dpi=2400)

# Annual sediment inflow
# Annual sediment outflow

# Monthly sediment inflow

# Monthly sediment outflow

# Monthly sediment inflow and outflow COMBINED

########################################################################################################################
# 2. Energy production

# a. Annual production (distribution of Num_Realizations*100yrs values)
# i. Histogram
Ann_Sum_Series = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].resample('A', how=sum)
Ann_Sum_Series_Stack = Ann_Sum_Series.stack(0)  # Create single stack of all realizations to create one plot for all realizations.
plt.figure(1)
Ann_Sum_Series_Stack.hist()
plt.ylabel('Frequency')
plt.xlabel('Annual Energy Production (MWH)')
title('Annual Energy Production (MWH)')
plt.savefig('Annual Energy Prod Time Series.png', format='png', dpi=2400)
plt.figure(2)
plt.plot(Ann_Sum_Series_Stack)
plt.ylabel('Annual Energy Production (MWH)')
plt.xlabel('Year n of 10,000')
plt.savefig('Annual Energy Prod Histogram.png', format='png', dpi=2400)


# b. Monthly energy production (MWh)
new_DF = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].copy(deep=True)/1000  # use deep copy to ensure you are only creating a copy rather than permanently fixing these dataframes.
new_DF.columns = ['Realization1' for i in range(Num_Realizations[Sims_to_Import[0]])]  # Rename all columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.
new_DF_SUM = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].resample('M',how=sum).copy(deep=True)/1000
new_DF_SUM.columns = ['Realization1' for i in range(Num_Realizations[Sims_to_Import[0]])]  # Rename all columns as having the same name so that we can "stack" them up using .append without a multi-index appearing.

# In the sequence below, use new_DF_SUM to plot sums for months, or new_DF to plot all daily values for each month.
column_list_arrays = [new_DF_SUM.ix[:,i] for i in range(1,Num_Realizations[Sims_to_Import[0]])]  # have to skip zero-th element, because it's the first element included in the append call below.
Monthly_Energy_Dictionary = {'Monthly Energy Production': new_DF_SUM.ix[:,0].append(column_list_arrays)}  # Take the first column of new_DF and append all the other columns to it.
Monthly_Energy_DF = pd.DataFrame(Monthly_Energy_Dictionary)  # Create a dataframe from the dictionary above
Monthly_Energy_DF['Month Number'] = Monthly_Energy_DF.index.month  # Create new column that stores the month value for every time series value (e.g., 1,2,...12)
Monthly_Energy_DF_GROUPED = Monthly_Energy_DF.groupby('Month Number')  # Create a group that slices the data by the Month Number column
Monthly_Energy_DF_GROUPED.describe()  # Get statistics for each month
Monthly_Energy_Figs = Monthly_Energy_DF_GROUPED.plot.hist(bins=20, histtype='step', cumulative=False, normed=True)  # Plots a separate histogram for each month

# We may also want to plot a subplot of histograms.
new_index = [i for i in range(int(Num_Realizations[Sims_to_Import[0]]*Num_Years))]  # Monthly_Energy_DF.index #  simulation_dates
DF_Monthly_Columns = pd.DataFrame(index=new_index)  # Stores a column for each month with all values (Num_Realizations[sims]*(Sim_Dur/36500)) stored in each column, with generic numbered row index.
month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
i = 0
#New_Fig = plt.figure(1)
Mon_Energ_Fig = plt.figure(1)
suptitle("Monthly Inflow Volume")  # This gives your overall figure a title, whereas subplot figures are done below.
for mon in month_list:
    i+=1
    DF_Monthly_Columns[mon] = Monthly_Energy_DF[Monthly_Energy_DF.index.month==i]['Monthly Energy Production'].values  # Note: This just copies the VALUES from one data frame (without indices attached) and stores them in this dataframe, where dates dont really need to be stored.
    #New_Fig.add_subplot()
    ax1 = plt.subplot(4,3,i)  # Could also make this P.subplot (using pylab; may work better for having a 4x3 subplot)
    plt.hist(DF_Monthly_Columns[mon], bins=30, cumulative=False, normed=True)  # Create monthly histogram.
    plt.title(mon)  # If you want a title to appear in each subplot, use this.
    plt.xlabel("Monthly Inflow Volume")
    plt.ylabel('Frequency (#)')
plt.tight_layout()  # Makes everything fit more nicely into the figure area.
Mon_Energ_Fig.savefig('Monthly Flow Subplot.png', format='png', dpi=1200)

for dataframes try to use .plot.hist(legend=None


# Create monthly histogram just from dataframe using hist feature of pandas dataframe. Limited in axis labeling features.
Monthly_Energy_Subplot = DF_Monthly_Columns.hist()
# Create monthly histogram just from dataframe, but this time using plot feature of pyplot.
Monthly_Energy_Subplot = DF_Monthly_Columns.plot(kind='hist', stacked=False)  # Create monthly histogram.

# c. firm power


# d. Peaking capacity
# Plot an annual resampling of minimum values, to see how dry season peaking capacity is affected.
Peaking_Fig = plt.figure(1)
Peaking_Fig = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]]['theoretical_peaking_capacity'].resample("A",how=min).plot(legend=None)
Peaking_Fig.set_xlabel('Time', fontsize=18)
Peaking_Fig.set_ylabel('Number of hours peak discharge is possible', fontsize=18)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
plt.savefig('Peaking capability.png', format='png', dpi=1200)

# Plot a random year (1985) of Inflow, and then hours
Peaking_Fig_1985 = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]]['theoretical_peaking_capacity']['Realization1']['1985'].plot(legend=None)
Peaking_Fig_1985.set_xlabel('Time', fontsize=18)
Peaking_Fig_1985.set_ylabel('Number of hours peak discharge is possible', fontsize=18)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
plt.savefig('Peaking capability for 1 year.png', format='png', dpi=1200)

Peaking_Fig_Inflow = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]]['Q_in']['Realization1']['1985'].plot(legend=None)
Peaking_Fig_Inflow.set_xlabel('Time', fontsize=18)
Peaking_Fig_Inflow.set_ylabel('Inflow Rate (m^3/s)', fontsize=18)
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
plt.savefig('Daily inflow 1985.png', format='png', dpi=1200)
########################################################################################################################
# 3. Fish passage




# b. Vulnerability to providing less than target

# c. Resilience (bouncing back from target deficit)
########################################################################################################################
# 4. Egg passage

# b. Vulnerability to providing less than velocity target


# c. Resilience (bouncing back from velocity target deficit
Egg_Pass_Vel_DF = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]]['Egg_Passage_Velocity_Success']
Egg_Pass_Resil_DF = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]]['egg_pass_resilience_tracker']
Num_success_days = Egg_Pass_Vel_DF.sum()
Num_failed_days = len(Egg_Pass_Vel_DF.index) - Num_success_days
Egg_Passage_Resilience = 1 - Egg_Pass_Resil_DF.sum()/(Num_failed_days)
plt.figure(1)
Egg_Passage_Resilience.hist()  # Histogram plot. Use cumulative=True to plot a cumulative histogram.
plt.ylabel('Frequency')
plt.xlabel('Egg Passage Minimum Velocity Resilience')
plt.title('Egg Passage Velocity Resilience')
plt.savefig('Egg Passage Velocity Resilience.png', format='png', dpi=1200)


# a. Annual average daily inflow rate
# i. Histogram
var_name = 'Q_in'  # State variable to be plotted
Ann_Sum_Series = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].resample('A', how='mean')
Ann_Sum_Series_Stack = Ann_Sum_Series.stack(0)  # Create single stack of all realizations to create one plot for all realizations.
plt.figure(1)
Ann_Sum_Series_Stack.hist()
plt.ylabel('Frequency')
plt.xlabel('Inflow rate')
plt.title('Inflow rate')
plt.savefig('Inflow rate Time Series.png', format='png', dpi=2400)
plt.figure(2)
plt.plot(Ann_Sum_Series_Stack)
plt.ylabel('Inflow rate')
plt.xlabel('Inflow rate')
plt.savefig('Inflow rate Histogram.png', format='png', dpi=2400)















# CODE THAT WORKED BUT IS NO LONGER NEEDED.

# b. Monthly energy production. MAY ONLY WORK FOR A SINGLE REALIZATION, NOT FOR MULTIPLES.
Energy_DF = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name]
Energy_DF_STACKED = Energy_DF.stack(0)
Energy_DF_STACKED['mon-new-col'] = Energy_DF_STACKED.index.month  # Add new month column at end of dataframe.
Energy_DF_STACKED_GROUPED = Energy_DF_STACKED.groupby('mon-new-col')  # Add new "month number" Column to Energy_DF_STACKED
# Stack all Realizations on top of one another so one histogram

New_DF = Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].resample('M',how=sum).stack(0)
New_DF_2 = pd.DataFrame(index=simulation_dates)
month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
i = 0
for mon in month_list:
    i+=1
    New_DF_2[mon] = New_DF[New_DF.index[:][0].month==i]

New_DF_2.hist(sharey=True, sharex = True, layout=(4,3)) # Scale same on x-axis
suptitle("Monthly Energy Production")
New_DF_2.hist(sharey=True, sharex = False, layout=(4,3))  # Don't scale same on x-axis
suptitle("Monthly Energy Production 2")

New_DF_2.plot(kind='hist', stacked=False) # Scale same on x-axis

# Testing of distribution fitting features in scipy

# Distribution fitting
param = norm.fit(Ann_Sum_Series_Stack) # now, param[0] and param[1] are the mean and
# the standard deviation of the fitted distribution
x = linspace(0.6E07,1.1E07,100)
# fitted distribution
pdf_fitted = norm.pdf(Ann_Sum_Series_Stack, loc=param[0],scale=param[1])
# original distribution
plot(x,pdf_fitted,'r-')
#hist(samp,normed=1,alpha=.3)
show()

# Tips and tricks:

# Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].drop(Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].columns[[2]], axis=1, inplace=True) # To drop a column that you accidentally added to a DF
# Time_Series_Import_Dictionary[Sim_name[0]][Sys_Obj[0]][var_name].columns = ['Realization1', 'Realization2']  # Say you had accidentally named both columns by the same name. now you can rename the columns of the DF this way.
# new_DF.columns = ['Realization1', 'Realization3']  # If you have a new_DF with Columns named Realization1 and Realization2, here you can rename the second column.
# a = Ann_Sum_Series_Stack.as_matrix()  # What you'd do if you wanted to convert the long series to a frame.