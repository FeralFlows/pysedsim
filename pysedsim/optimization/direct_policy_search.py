import numpy as np
from pysedsim.data_processing.data_processing import Load_Input_File
from pysedsim.data_processing.data_processing import Excel_Data_Import
from pysedsim.visualization.processing_reference_set import *
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends import backend_agg as agg
import math
import logging

def DPS_Policy_Decision(DPS_inputs, DPS_dict):

    '''

    Purpose: Method to determine a reservoir's daily release/target water level according to a specified operating
    policy. Developed outside of the reservoir class so that the method can be accessed more efficiently externally.

    Before reservoir's mass balance is simulated for this time step, take externally specified reservoir
    operating policy (if user indicates such a RBF-based policy has been externally specified using DPS or other
    methods) and determine the required release/water level given system state variable values of interest. The
    method Set_DPS_Operations_Goal() is what actually determines the daily operational decision.

    Adapted from code originally written by Julianne Quinn, Cornell University (jdq2101@gmail.com)

    Args:
        DPS_inputs: list of daily inputs to DPS policy (e.g., [inflow[t], S[t]])
        DPS_dict: dictionary of DPS preferences created in Import_DPS_Preferences()
    '''

    # Unpack DPS preferences and parameters from dictionary
    M = DPS_dict['num_inputs']
    N = DPS_dict['RBF Specs']['Number']
    K = DPS_dict['num_reservoirs']
    DPS_input_ranges = DPS_dict['input_ranges']
    DPS_output_ranges = DPS_dict['output_ranges']
    C = DPS_dict['Policy Function']['Centers']
    R = DPS_dict['Policy Function']['Radii']
    W = DPS_dict['Policy Function']['Weights']

    # Load inputs to DPS (RBF) policy function into an array (inputs), then normalize that array in range [0,
    # 1] as normInputs
    DPS_normInputs = np.zeros(M)  # init normalized inputs array so it can store M inputs.
    for i in range(len(DPS_dict['ordered_input_names'])):
        name = DPS_dict['ordered_input_names'][i]
        DPS_normInputs[i] = (DPS_inputs[i] - DPS_input_ranges[name]['range'][0]) / (
        DPS_input_ranges[name]['range'][1] - DPS_input_ranges[name]['range'][0])

    # Calculate normalized releases (u) corresponding to each sample of normalized inputs
    # u is a 2-D matrix with K columns (1 for each output)
    # and as many rows as there are samples of inputs
    u = np.zeros((K,))
    use_sin_time_index = 'No'
    if DPS_dict['RBF Specs']['Type'] == 'Gaussian':
        for k in range(K):
            for n in range(N):
                BF = 0
                for m in range(M):
                    if R[m,n] > 10**-6:
                        if DPS_dict['ordered_input_names'][m] == 'current_date':
                            if use_sin_time_index == 'Yes':
                                BF = BF + ((math.sin(math.pi*DPS_normInputs[m])-C[m,n])/R[m,n])**2
                            else:
                                BF = BF + ((DPS_normInputs[m]-C[m,n])/R[m,n])**2
                        else:
                            BF = BF + ((DPS_normInputs[m]-C[m,n])/R[m,n])**2
                    else:
                        BF = BF + ((DPS_normInputs[m]-C[m,n])/(10**-6))**2

                u[k] = u[k] + W[k, n]*np.exp(-BF)

    # Convert normalized decision to actual decision for every reservoir in policy (e.g., release or water level)
    for i in range(len(u)):
        name = DPS_dict['ordered_output_names'][i]
        u[i] = DPS_output_ranges[name]['range'][0] + u[i] * (
        DPS_output_ranges[name]['range'][1] - DPS_output_ranges[name]['range'][0])
    return u

def Create_DPS_Policy(DPS_dict):

    # Unpack DPS preferences from DPS dictionary
    DPS_dict = DPS_dict
    M = DPS_dict['num_inputs']
    N = DPS_dict['RBF Specs']['Number']
    K = DPS_dict['num_reservoirs']
    DPS_Policy = DPS_dict['Policy Function']['Raw Parameters']

    # Calculate reservoir WSE corresponding to inputs at time t, K columns are the WSE targets for each reservoir
    # Re-organize decision variables into weight (W), center (C) and raddi (R) matrices
    # C and R are M x N, and W is K X N
    # Decision variables are currently arranged in 'DPS_policy' as
    # N consecutive sets of {M pairs of {C, R} followed by K Ws}, e.g. for N = 2, M = 3 and K = 4:
    # C, R, C, R, C, R, W, W, W, W, C, R, C, R, C, R, W, W, W, W
    C = np.zeros([M, N])
    R = np.zeros([M, N])
    W = np.zeros([K, N])
    for n in range(N):
        for m in range(M):
            C[m, n] = DPS_Policy[(2*M+K)*n + 2*m]
            R[m, n] = DPS_Policy[(2*M+K)*n + 2*m + 1]
        for k in range(K):
            W[k, n] = DPS_Policy[(2*M+K)*n + 2*M + k]

    # Normalize weights to sum to 1 across N RBFs (so each row of W should sum to 1)
    totals = np.sum(W,1)
    for k in range(K):
        if totals[k] > 10**-6:
            W[k, :] = W[k, :]/totals[k]

    # Having created the policy function, now store it in the DPS dict, and return that new dict.

    DPS_dict['Policy Function']['Centers'] = C
    DPS_dict['Policy Function']['Radii'] = R
    DPS_dict['Policy Function']['Weights'] = W

    return DPS_dict

def Import_DPS_Preferences(simulation_title=None, imported_specs=None, main_input_files_dir=None, input_data_file=None):
    '''

    Purpose: To import preferences related to Direct Policy Search (DPS).

    Module will often be used in an optimization setting, but can also be used if the user intends to simulate a
    reservoir operating policy that has been externally created using DPS.

    Args:
        simulation_title: Used to locate input data file. Must be specified if input_data_file is not specified.
        imported_specs: Used to locate input data file. Must be specified if input_data_file is not specified.
        main_input_files_dir: Used to locate input data file. Must be specified if input_data_file is not specified.
        input_data_file: The input data .xlsx file (openpyxl workbook object)
    Returns:
        DPS_dict: A dictionary storing all preferences relevant to DPS
    '''

    if input_data_file is None:
        Input_Data_File = Load_Input_File(simulation_title, main_input_files_dir, imported_specs)
    else:
        Input_Data_File = input_data_file
    opt_approach = 'DPS'

    # Gather optimization information from "DPS" worksheet if it exists.
    # Set the number of seeds for optimization
    try:
        num_RBF = Input_Data_File[opt_approach]['B3'].value
        if (num_RBF is None) or (type(num_RBF) not in [int, long]):
            # No value set by user. Set default:
            num_RBF = 2
        num_inputs = Input_Data_File[opt_approach]['B2'].value
        num_reservoirs = Input_Data_File[opt_approach]['B4'].value
        # Set the type (RBF = 'cubic', 'Gaussian' or 'TPS') and number (n) of RBFs
        RBF_type = 'Gaussian'
        # Set the number of decision variables, objectives and constraints
        n_vars = (2*num_inputs + num_reservoirs)*num_RBF

        DPS_input_ranges = {}
        var_names = []
        var_time_index = []
        for n in range(num_inputs):
            raw_current_input_name = Input_Data_File[opt_approach].cell(row=2 + n, column=3).value
            current_input_time_index = None
            if '(t-' in raw_current_input_name or '(t -' in raw_current_input_name:
                if '(t-' in raw_current_input_name:
                    raw_current_input_list = raw_current_input_name.split('(t-')
                else:
                    raw_current_input_list = raw_current_input_name.split('(t -')

                current_input_name = raw_current_input_list[0]  # Grab variable name
                current_input_time_index = raw_current_input_list[1]  # Get information relevant to time index.
                # Remove any parentheses that exist in the number.
                if ')' in current_input_time_index:
                    current_input_time_index = current_input_time_index.replace(')', '')
                current_input_time_index = int(current_input_time_index)
            else:
                current_input_name = Input_Data_File[opt_approach].cell(row=2 + n, column=3).value
            var_names.append(current_input_name)
            # Append a time index to use for the time series variable.
            if current_input_time_index is None:
                var_time_index.append(0)
            else:
                var_time_index.append(current_input_time_index)
            DPS_input_ranges[current_input_name] = {'range': np.zeros(2), 'rank': n}
            for col in range(2):
                DPS_input_ranges[current_input_name]['range'][col] = Input_Data_File[opt_approach].cell(row=2 + n,
                                                                                               column=4 + col).value

        # Loop through "DPS" worksheet of input file and store output variable ranges.
        DPS_output_ranges = {}
        output_names = []
        for n in range(num_reservoirs):
            current_output_name = Input_Data_File[opt_approach].cell(row=2 + n, column=6).value
            output_names.append(current_output_name)
            DPS_output_ranges[current_output_name] = {'range': np.zeros(2), 'rank': n}
            for col in range(2):
                DPS_output_ranges[current_output_name]['range'][col] = Input_Data_File[opt_approach].cell(row=2 + n,
                                                                                       column=7 + col).value

    except KeyError:
        raise KeyError(
            "ERROR: Required worksheet %s does not exist in input file for scenario %s" % (opt_approach,
                                                                                           simulation_title))

    # Read in preferences related to optimization of reservoir sediment flushing.
    element_name = 'Sambor SA'
    flushing_dict = {}
    n_flush_vars = 0
    if 'Flushing' in Input_Data_File.sheetnames:
        flushing_data = Excel_Data_Import(element_name, Input_Data_File, 'Flushing', 2, 16,
                                          max_distinct_data_types=None, data_name_offset=2)
        try:
            if flushing_data[12][0] in ['Yes', 'yes', 'y', 'Y']:
                # User wishes to include flushing parameters in the optimization/search.
                flushing_dict['Input variable'] = {}
                flushing_dict['Ordered input variable list'] = []
                if flushing_data[13] is not None:
                    flushing_dict['Input variable']['Frequency'] = [int(x) for x in flushing_data[13]]
                    n_flush_vars += 1
                if flushing_data[14] is not None:
                    flushing_dict['Input variable']['Day of year'] = [int(x) for x in flushing_data[14]]
                    n_flush_vars += 1
                if flushing_data[15] is not None:
                    flushing_dict['Input variable']['Inflow trigger'] = flushing_data[15]
                    n_flush_vars += 1
                flushing_dict['n_vars'] = n_flush_vars
                for key in flushing_dict['Input variable'].keys():
                    flushing_dict['Ordered input variable list'].append(key)
        except IndexError:
            pass  # Optimization is not being done, no value provided in that part of Flushing worksheet.
    # Determine total number of decision variables
    total_vars = 0
    try:
        total_vars += n_vars
    except KeyError:
        total_vars += 0
    try:
        total_vars += n_flush_vars
    except KeyError:
        total_vars += 0

    variable_type = [0 for v in range(n_vars + n_flush_vars)]
    for v1 in range(n_vars):
        variable_type[v1] = 'real valued'
    try:
        # If flushing included in the optimization
        for v2 in range(n_flush_vars):
            if flushing_dict['Ordered input variable list'][v2] in ['Frequency', 'Day of year']:
                variable_type[n_vars + v2] = 'integer'
            elif flushing_dict['Ordered input variable list'][v2] == 'Inflow trigger':
                variable_type[n_vars + v2] = 'real valued'
    except KeyError:
        pass

    # Define DPS dictionary to return:
    if opt_approach == 'DPS':
        DPS_dict = {
            'RBF Specs':
            {'Type': RBF_type, 'Number': num_RBF},
            'num_inputs': num_inputs,
            'num_reservoirs': num_reservoirs,
            'n_vars': n_vars,
            'input_ranges': DPS_input_ranges,
            'ordered_input_names': var_names,
            'ordered_input_time_index': var_time_index,
            'ordered_output_names': output_names,
            'output_ranges': DPS_output_ranges,
            'Policy Function': {},
            'Flushing Optimization': flushing_dict,
            'variable_type': variable_type,
            'total_vars': total_vars
        }
    else:
        DPS_dict = None

    return DPS_dict

def Import_Flushing_Preferences(simulation_title=None, imported_specs=None, main_input_files_dir=None,
                            input_data_file=None):
    '''

    Purpose: To import preferences related to Flushing.

    Module will only be used in an optimization setting, wherein optimal search is being performed on flushing
    parameters such as timing, trigger inflow, and frequency.

    Args:
        simulation_title: Used to locate input data file. Must be specified if input_data_file is not specified.
        imported_specs: Used to locate input data file. Must be specified if input_data_file is not specified.
        main_input_files_dir: Used to locate input data file. Must be specified if input_data_file is not specified.
        input_data_file: The input data .xlsx file (openpyxl workbook object)
    Returns:
        DPS_dict: A dictionary storing all preferences relevant to DPS
    '''

def DPS_Policy_Plotting(DPS_dict, dps_plt_dict, objs_to_plot, ref_set_pref_dict):
    '''

    Args:
        DPS_dict: Dictionary produced by calling Import_DPS_Preferences().
        dps_plt_dict: Dictionary that contains preferences for plots of policies the module will produce. An example
        appears like the following:

        dps_plt_dict = {'plot_x': 'Qin', 'plot_y': 'output_elev', 'plot_color': 'day of year', 'plot_const':
        'Water_surf_elev', 'xlabel': 'Qin', 'ylabel': 'WSE(t+1)', 'colorlabel': 'day of year', 'input_constant': {
        'Water_surf_elev': 1, 'day of year': 0, 'Qin': 0}, 'slice_num': 16, 'grid_size': 25}

        description of dps_plt_dict keys:
        input_constant: what factor you want to be sliced and have a different value for each subplot.

        Note that if xlabel, ylabel, and colorlabel keys aren't provided, axes just get labeled with DPS inputs
        names from the pysedsim input data sheet.

        objs_to_plot: A list of strings of objective value names, which will eventually serve as plot labels for any
        plots that plot by objective value, and will also be keys for a dictionary storing particular policies plotted.
        ref_set_pref: A dictionary containing the following keys:
        a. unit_conv: List of factors by which to multiply values in each array of objective values.
        b. perc_conv: List of 'Yes' or 'No' for each objective, for converting objective values to percent from
        fraction.
        c. invert: List of 'Yes' or 'No' for each objective, for reversing percentage (doing 100-value).
    Returns:
        various plots of reservoir release vs. DPS input variables.

    '''
    # If reference set preferences dictionary exists, then the user wishes to use it:
    # Unpack ref_set_pref_dict
    num_objs = ref_set_pref_dict['num_objs']
    n_inputs = DPS_dict['num_inputs']
    if n_inputs > 1:
        try:
            dps_plt_dict['colorlabel']
        except KeyError:
            # xlabel not provided, so fill in value with name of DPS input variable for color axis.
            dps_plt_dict['colorlabel'] = dps_plt_dict['plot_color']
    try:
        ref_set_file_name = ref_set_pref_dict['ref_set_file_name']
    except KeyError:
        ref_set_file_name = None
    try:
        unit_conv = ref_set_pref_dict['unit_conv']
    except KeyError:
        unit_conv = None
    try:
        perc_conv = ref_set_pref_dict['perc_conv']
    except KeyError:
        perc_conv = None
    try:
        invert = ref_set_pref_dict['invert']
    except KeyError:
        invert = None
    parse_objs = None  # Really only applicable for objective space plotting.

    # Unpack DPS plot dict
    try:
        dps_plt_dict['xlabel']
    except KeyError:
        # xlabel not provided, so fill in value with name of DPS input variable for x axis.
        dps_plt_dict['xlabel'] = dps_plt_dict['plot_x']
    try:
        dps_plt_dict['ylabel']
    except KeyError:
        # ylabel not provided, so fill in value with name of DPS input variable for y axis.
        dps_plt_dict['ylabel'] = dps_plt_dict['plot_y']
    try:
        if dps_plt_dict['single_policy_plot'] in ['Yes', 'yes', 'Y']:
            single_policy_plot = 'Yes'  # Each figure/image produced is for a single policy
        elif dps_plt_dict['single_policy_plot'] in ['No', 'no', 'N']:
            single_policy_plot = 'No'  # Each figure/image produced contains data for multiple policies.
        else:
            single_policy_plot = 'No'  # Each figure/image produced contains data for multiple policies.
    except KeyError:
        single_policy_plot = 'Yes'  # Default is to produce a single policy plot

    try:
        grid_size = dps_plt_dict['grid_size']  # number of increments for each variable to be plotted.
    except KeyError:
        grid_size = 25  # Default

    try:
        slice_num = dps_plt_dict['slice_num']  # number of increments for each variable to be plotted.
    except KeyError:
        if n_inputs == 1:
            slice_num = 1
        else:
            slice_num = 16  # # Default: num slices/cross-sections of "constant" input variable to take.

    try:
        fig_const_var_name = dps_plt_dict['figure_const']  # number of increments for each variable to be plotted.
        # Variable held constant in each figure, but is varied across multiple plots being created.
    except KeyError:
        fig_const_var_name = ''

    subplot_rows = int(math.sqrt(slice_num))  # number of subplot rows in figure
    subplot_cols = int(math.sqrt(slice_num))  # number of subplot columns in figure
    [ref_set_array, objective_values, dec_var_values] = Initial_Processing(num_objs,
                                                                                                    DPS_dict['n_vars'],
                                                                                                    ref_set_file_name,
                                                                                                    parse_objs,
                                                                                                    perc_conv, invert,
                                                                                                    unit_conv)

    # Code that plots Decision vs. WSE and Decision vs. Q_in, both with t as a color map. Snapshots with either Qin or
    # WSE over 12 values are taken and provided as subplots.

    input_val = np.zeros(n_inputs, dtype=object)
    const_val = np.zeros([n_inputs,slice_num])
    increment = np.zeros(n_inputs)

    if single_policy_plot == 'Yes':
        try:
            policies_to_plot = dps_plt_dict['policies_to_plot']  # User needs to provide a list of policies to make
            # plots for.
        except KeyError:
            logging.critical("User must supply row numbers of policies to plot from reference set file")
    else:
        try:
            policies_to_plot = dps_plt_dict['policies_to_plot']  # User needs to provide a list of policies to make
            # plots for.
        except KeyError:
            policies_to_plot = [pol for pol in range(len(ref_set_array))]  # List contains all policies from reference set.
        figure_const_value = dps_plt_dict['figure_const_value']

    fig = matplotlib.figure.Figure()
    agg.FigureCanvasAgg(fig)

    cmap = plt.cm.get_cmap('jet_r') #plt.cm.get_cmap('Spectral')
    inp_name = DPS_dict['ordered_input_names']

    try:
        title_subtitle = {'varied': [], 'not varied': []}
        for sub_keys in dps_plt_dict['input_constant']:
            if dps_plt_dict['input_constant'][sub_keys] == 0:
                title_subtitle['varied'].append(sub_keys)
            elif dps_plt_dict['input_constant'][sub_keys] == 1:
                title_subtitle['not varied'].append(sub_keys)
    except KeyError:
        pass

    if single_policy_plot == 'Yes':
        if n_inputs == 1:
            dps_plt_dict['main_plot_title'] = '%s varied' % (title_subtitle['varied'][0])
        elif n_inputs == 2:
            dps_plt_dict['main_plot_title'] = '%s and %s varied' % (title_subtitle['varied'][0], title_subtitle[
                'varied'][1])
        elif n_inputs == 3:
            dps_plt_dict['main_plot_title'] = '%s and %s varied, %s=snapshot' % (
            title_subtitle['varied'][0], title_subtitle['varied'][1], title_subtitle['not varied'][0])
    else:
        dps_plt_dict['main_plot_title'] = '%s varied, %s=snapshot, %s=constant=%s' % (
        title_subtitle['varied'][0], title_subtitle['not varied'][0], title_subtitle['not varied'][1],
        dps_plt_dict['figure_const_value'])

    # Determine how many values are being computed as inputs for each input variable, and size the abcis_ord array
    # accordingly.
    grid_size_var = np.zeros(n_inputs)  # Number of different values of each input variable to evaluate.
    len_abcis_ord = 1
    for i in range(n_inputs):
        if dps_plt_dict['input_constant'][inp_name[i]] == 0:
            grid_size_var[i] = grid_size  # Use regular grid size, not a constant snapshot input variable.
        else:
            grid_size_var[i] = 1
        len_abcis_ord = len_abcis_ord*grid_size_var[i]
    # Initialize size of abcis_ord so it reflects varying dimensions for each input variable.
    len_abcis_ord = int(len_abcis_ord)*len(policies_to_plot)  # Must get extended for multiple policies.
    abcis_ord = np.zeros([len_abcis_ord, n_inputs+1+num_objs])  # Initialize

    for z in range(1, slice_num + 1):
        for i in range(n_inputs):
            if dps_plt_dict['input_constant'][inp_name[i]] == 0:
                increment[i] = ((DPS_dict['input_ranges'][inp_name[i]]['range'][1] -
                                 DPS_dict['input_ranges'][inp_name[i]]['range'][0]) / grid_size_var[i])
                input_val[i] = np.arange(DPS_dict['input_ranges'][inp_name[i]]['range'][0],
                                         DPS_dict['input_ranges'][inp_name[i]]['range'][1], increment[i])
            else:
                if inp_name[i] == fig_const_var_name:
                    const_val[i, :] = figure_const_value
                else:
                    # Variable is regular constant, not being varied over multiple plots.
                    const_val[i, :] = np.arange(DPS_dict['input_ranges'][inp_name[i]]['range'][0],
                                                DPS_dict['input_ranges'][inp_name[i]]['range'][1]+1,
                                                ((DPS_dict['input_ranges'][inp_name[i]]['range'][1] - DPS_dict[
                                                    'input_ranges'][inp_name[i]]['range'][0]) / (slice_num-1)))
                input_val[i] = np.empty(grid_size_var[i])
                input_val[i].fill(const_val[i][z-1])

        ax = fig.add_subplot(subplot_rows, subplot_cols, z)
        counter = 0

        for policy in policies_to_plot:
            DPS_dict['Policy Function']['Raw Parameters'] = dec_var_values[policy]  # Add policy to dict.
            DPS_dict = Create_DPS_Policy(DPS_dict)
            if n_inputs == 1:
                for i in range(int(grid_size_var[0])):
                    DPS_inputs = [input_val[0][i]]
                    abcis_ord[counter, 0:n_inputs] = DPS_inputs
                    abcis_ord[counter, n_inputs] = DPS_Policy_Decision(DPS_inputs, DPS_dict)  # Release decision
                    for n in range(num_objs):
                        abcis_ord[counter, n_inputs + 1 + n] = objective_values[n][policy]
                    counter += 1
            elif n_inputs == 2:
                for i in range(int(grid_size_var[0])):
                    for j in range(int(grid_size_var[1])):
                        DPS_inputs = [input_val[0][i], input_val[1][j]]
                        abcis_ord[counter, 0:n_inputs] = DPS_inputs
                        abcis_ord[counter, n_inputs] = DPS_Policy_Decision(DPS_inputs, DPS_dict)  # Release decision
                        for n in range(num_objs):
                            abcis_ord[counter, n_inputs + 1 + n] = objective_values[n][policy]
                        counter += 1
            elif n_inputs == 3:
                for i in range(int(grid_size_var[0])):
                    for j in range(int(grid_size_var[1])):
                        for k in range(int(grid_size_var[2])):
                            DPS_inputs = [input_val[0][i], input_val[1][j], input_val[2][k]]
                            abcis_ord[counter, 0:n_inputs] = DPS_inputs
                            abcis_ord[counter, n_inputs] = DPS_Policy_Decision(DPS_inputs, DPS_dict)  # Release decision
                            for n in range(num_objs):
                                abcis_ord[counter, n_inputs + 1 + n] = objective_values[n][policy]
                            counter += 1
        xloc = DPS_dict['input_ranges'][dps_plt_dict['plot_x']]['rank']
        yloc = n_inputs  # + DPS_dict['output_ranges'][dps_plt_dict['plot_y']]['rank']
        if single_policy_plot == 'Yes':
            if n_inputs > 1:
                c_loc = DPS_dict['input_ranges'][dps_plt_dict['plot_color']]['rank']
                ax.scatter(abcis_ord[:,xloc],abcis_ord[:,yloc],c=abcis_ord[:,c_loc], cmap=cmap, s=10, linewidth=0)
            else:
                ax.scatter(abcis_ord[:,xloc],abcis_ord[:,yloc], s=10, linewidth=0)
        else:
            c_loc = n_inputs + 1 + objs_to_plot.index(dps_plt_dict['colorlabel'])
            try:
                s_loc = n_inputs + 1 + objs_to_plot.index(dps_plt_dict['sizelabel'])
            except KeyError:
                # No sizelabel specified
                ax.scatter(abcis_ord[:,xloc],abcis_ord[:,yloc],c=abcis_ord[:,c_loc], cmap=cmap, s=15, linewidth=0)

        ax.tick_params(axis='both',labelsize=6)
        ax.set_xlim([DPS_dict['input_ranges'][dps_plt_dict['plot_x']]['range'][0], DPS_dict['input_ranges'][
            dps_plt_dict['plot_x']]['range'][1]])
        ax.set_ylim([DPS_dict['output_ranges'][dps_plt_dict['plot_y']]['range'][0], DPS_dict['output_ranges'][
            dps_plt_dict['plot_y']]['range'][1]])
        # Give x axis scientific notation if valuable.
        if DPS_dict['input_ranges'][dps_plt_dict['plot_x']]['range'][1] > 1000:
            ax.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
        if DPS_dict['output_ranges'][dps_plt_dict['plot_y']]['range'][1] > 1000:
            ax.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
#            formatter = matplotlib.ticker.ScalarFormatter()
#            formatter.set_powerlimits((-3, 5))
#            ax.xaxis.set_major_formatter(formatter)
        ax.set_xlabel(dps_plt_dict['xlabel'],fontsize=8)
        ax.set_ylabel(dps_plt_dict['ylabel'],fontsize=8)
        if n_inputs > 1:
            current_const_val = int(const_val[inp_name.index(dps_plt_dict['plot_const'])][z-1])
            ax.set_title(dps_plt_dict['plot_const'] + '=' + '%s' % current_const_val, fontsize=4)
        #ax.plot(abcis_ord[:,2],abcis_ord[:,n_inputs],c=cmap(abcis_ord[:,1]))

    # Put finishing touches on plot (formatting of whole plot, not individual subplots)
    if n_inputs > 1:
        sm = matplotlib.cm.ScalarMappable(cmap=cmap)
        if single_policy_plot == 'Yes':
            sm.set_array([DPS_dict['input_ranges'][dps_plt_dict['plot_color']]['range'][0],
                          DPS_dict['input_ranges'][dps_plt_dict['plot_color']]['range'][1]])
        else:
            sm.set_array([min(objective_values[objs_to_plot.index(dps_plt_dict['colorlabel'])]),
                          max(objective_values[objs_to_plot.index(dps_plt_dict['colorlabel'])])])
    fig.tight_layout()
    if n_inputs > 1:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])  # OLD: cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.5])
        fig.colorbar(sm, cax=cbar_ax)
        fig.axes[-1].set_xlabel(dps_plt_dict['colorlabel'])  # Give colorbar a label/title
    #fig.suptitle(dps_plt_dict['main_plot_title'], fontsize=14)
    fig.suptitle('')
    plt.subplots_adjust(top=0.90)
    fig.savefig(dps_plt_dict['main_plot_title'] + '.png', dpi=1200)
    plt.show()
    #c=cmap(abcis_ord[:,1])

    return