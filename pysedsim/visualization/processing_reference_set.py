'''
The purpose of this module is to import an externally generated reference set and plot the objective function and
decision variable values.

The reference set that this file will import must be generated using a shell script file (e.g.,
pysedsim_reference_set.sh) that uses Pareto.py to sort through reference files produced for each seed by the
optimization algorithm and make a final reference set that is a .txt file called PySedSim_reference.txt. That text
file is imported and manipulated/plotted here. This module contains functions that process the data (extract the
objective values and the decision variable values), plots objective value space, plots decision variable space,
and produces text files that store particular operating policies of interest.

'''
from __future__ import division
import numpy as np
import matplotlib
import matplotlib.cm as cm
#matplotlib.use('Agg')
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
#import pylab as p
from matplotlib import rc, rcParams
#from movie_pareto_front import rotanimate
import subprocess  # for calling Pareto.py command-line style
import platform
import os
import shutil  # for copying Pareto.py over to directory where seed outputs from optimization are located.
from data_processing import Op_Sys_Folder_Operator
from data_processing import Determine_Num_Scenarios
import direct_policy_search
import itertools
from data_processing import Load_Input_File
import matplotlib.patheffects as pe

def Reference_Set(input_file_name = 'PySedSim_Input_Specifications.csv', ref_set_file_name = 'pysedsim_ref_set',
                  file_type = '.ref', objs_to_process = None, inline_optimization = None, eval_runtime_perf='No',
                  create_ref_set = 'Yes', thin_ref_set = 'No', thinned_ref_set_file_name='thinned_pysedsim_ref_set',
                  keep_min_max = 'Yes', Objective_Range=None):

    '''
        Purpose: To post-process simulation-optimization (S-O) seed results to produce a reference set. Creates
        reference set file named "ref_set_file_name.ref", where ref_set_file_name is defined above.

        This module would be called on the back end after all optimization "scenarios" are
        completed. It is not explicitly set up to be parallelized. The assumption is a single processor will be used
        to process reference sets for each scenario from N seeds. The seed files are assumed to be stored in the
        location for output data specified in your PySedSim input file (argument: input_file_name), in a directory
        called "seed". This is the file structure set up by PySedSim during the process of exporting
        S-O results, so if you are going to run this particular function on a separate computer,
        from where the S-O results were generated, it is suggested you replicate the same file structure, and modify
        your input_file_name.csv file to reflect the new working directory.

        Sorting is performed with Pareto.py over all objectives. This
        file must be in the present working directory. It can be downloaded here:
        https://github.com/matthewjwoodruff/pareto.py

        Args:
            :param input_file_name: Name of top level input file containing scenario (simulation) names,
            default: 'PySedSim_Input_Specifications.csv'. If no input file name is specified, then the user must supply
            num_objs and num_dec_vars, as the function does not otherwise know how to determine these values.

            :param ref_set_file_name: Desired name of reference set file to be produced (string). Each row will
            correspond to a different pareto solution, with the first num_dec_vars columns corresponding to decision
            variable values, and the remaining columns corresponding to num_objs columns. Default will be
            'pysedsim_ref_set.ref'

            :param  objs_to_process: Optional, List of strings of objective numbers in reference sets to process,
            if user doesn't want to create a reference set across all objectives present in the reference set files.
            Valid values: from zero to num_objs-1. Order of objectives is assumed to be same as appears in pysedsim
            optimization outputs.

            :param inline_optimization: Optional, used only when the Reference_Set() method is called inline after an
            optimization scenario has been completed. Contains the following list of inputs obtained directly from
            optimization_pysedsim: [simulation_title, Borg_dict, DPS_dict, main_output_file_dir]

            :param eval_runtime_perf: Optional, 'Yes' or 'No'. Used if user wants to evaluate the runtime performance
            of the MOEA used to generate solutions. Calls runtime_performance() function to evaluate performance
            across random seeds.

            :param thin_ref_set: Optional, 'Yes' or 'No'. Used if user wants to thin the reference set that has
            already been generated using epsilon dominance and the Pareto.py sorting algorithm. If using this option,
            the file name specified through the ref_set_file_name argument will be used as the reference set file. In
            this case, the user should specify create_ref_set = 'No', so the reference set is not generated again.

            :param thinned_ref_set_file_name: Optional, 'Yes' or 'No'. If user specifies thin_ref_set_argument='Yes',
            this file stores the new thinned reference set in this file name. Will be stored in the output file
            directory.

            :param keep_min_max: Optional, 'Yes' or 'No'. If user specfiies thin_ref_set_argument='Yes',
            this indicates whether the maximum and minimum values from the set should be kept in this file,
            to preserve the best/worst in various objectives.

            :param Objective_Range: Optional, this is a dictionary that specifies whether the minimum, maximum,
            or both values will be retained during the reference set thinning process. For example, to keep both the
            worst performing and best performing solutions for a six objective problem in the thinned set (they will
            be reinserted at the end of the set), specify the following dictionary as an argument (with preferences
            specified for every objective):


            Objective_Range = {'Objective 0': {'Range': ['max', 'min'], 'Type': 'Max'},
                               'Objective 1': {'Range': ['max', 'min'], 'Type': 'Max'},
                               'Objective 2': {'Range': ['max', 'min'], 'Type': 'Max'},
                               'Objective 3': {'Range': ['max', 'min'], 'Type': 'Max'},
                               'Objective 4': {'Range': ['max', 'min'], 'Type': 'Max'},
                               'Objective 5': {'Range': ['max', 'min'], 'Type': 'Min'},
                               }

            In the above example, 'max' indicates keep the best performing solution in the given objective,
            whereas 'min' indicates to keep the worst performing solution in the given objective. You can specify
            just 'max' or just 'min' if desired.

        Returns:
            :returns Borg_dict, DPS_dict
    '''

    from optimization_pysedsim import Import_Optimization_Preferences

    # Load basic information about the simulated scenarios, including input file directory and simulation names.
    os_fold = Op_Sys_Folder_Operator()
    ref_set_file_name_FINAL = ref_set_file_name + file_type
    ref_set_file_name_no_vars = ref_set_file_name + '_no_vars' + file_type

    if inline_optimization is None:
        # This file is not being called directly from within an optimization scenario (or loop of scenarios)
        file_name = input_file_name
        [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir,
         main_output_file_dir, Monte_Carlo_Parameters_File_List, external_mc_data] = Determine_Num_Scenarios(file_name)
    else:
        num_scenarios = 1
        main_output_file_dir = inline_optimization[3]

    for j in range(num_scenarios):
        if inline_optimization is None:
            # Import optimization details so reference set can be properly processed
            simulation_title = simulation_titles_list[j]
            Borg_dict = Import_Optimization_Preferences(simulation_title, imported_specs, main_input_files_dir)
            Borg_dict['opt_output_filename'] = ref_set_file_name_FINAL
            if Borg_dict['optimization approach'] == 'DPS':
                DPS_dict = direct_policy_search.Import_DPS_Preferences(simulation_title=simulation_title,
                                                                       imported_specs=imported_specs,
                                                                       main_input_files_dir=main_input_files_dir)
        else:
            simulation_title = inline_optimization[0]
            Borg_dict = inline_optimization[1]
            DPS_dict = inline_optimization[2]

        num_dec_vars = DPS_dict['total_vars']
        num_objs = Borg_dict['n_objs']
        epsilon_list = Borg_dict['epsilon_list']

        # Only create the Reference Set if user wants it to be created. This module may be called as a means of
        # generating plots from reference sets that have already been generated, as it produces dictionaries of
        # preferences related to optimization and operating policies.
        output_location = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets'

        if create_ref_set == 'Yes' or thin_ref_set == 'Yes':

            # Purpose is to generate reference set for all objectives from output of all seeds in optimization.
            # Can execute shell Pareto.py in shell if user operating in Linux. These steps must be done whether a
            # regular reference set is being generated or just thinning needs to occur.

            # Change to sets directory, and create the reference set.
            global sorting_file
            sorting_file = "Pareto.py"
            shutil.copy(sorting_file, output_location)  # Place a copy of Pareto.py temporarily in seed directory
            cwd = os.getcwd()  # Current working directory to return to after changing directories.
            os.chdir(output_location)

            # If reference set is being thinned, and user wishes to maintain the highest and lowest values from the
            # reference set in doing this, then store those values now, so they can be inserted later.
            if thin_ref_set == 'Yes' and keep_min_max == 'Yes':
                # If user did not specify the Objective_Range argument, then default is to assume user wants to keep
                # minimum and maximum values for all objectives.
                if Objective_Range is None:
                    Objective_Range = {}  # Initialize dictionary
                    for ob in range(num_objs):
                        Objective_Range['Objective %s' % ob] = {'Range': ['max', 'min'], 'Type': 'Max'}
                # Call separate module to actually do the thinning/brushing.
                [rowlist, refset_objs, refset] = brush_pareto_policies(Objective_Range,
                                                                      refset_file_path=ref_set_file_name_FINAL,
                                                                      objs_file_path=ref_set_file_name_no_vars,
                                                                      objective_union = 'Yes', one_best_worst='Yes')

            if objs_to_process is None:
                # User wants to create reference set across all objectives.
                obj_range = str(num_dec_vars) + '-' + str(num_dec_vars+num_objs-1)
                epsilon_strings = ""
                for i in range(num_objs):
                    epsilon_strings += str(epsilon_list[i])
                    epsilon_strings += " "
            else:
                # User has specified a subset of objectives for which to create a reference set.
                obj_range = ''
                epsilon_strings = ""
                for i in range(len(objs_to_process)):
                    obj_range += str(num_dec_vars + int(objs_to_process[i])) + " "
                    epsilon_strings += str(epsilon_list[int(objs_to_process[i])]) + " "

            if create_ref_set == 'Yes':
                file_to_sort = 'Borg_DPS_PySedSim*.set'
                file_to_produce_1 = ref_set_file_name_FINAL
                file_to_produce_2 = ref_set_file_name_no_vars
            elif thin_ref_set == 'Yes':
                file_to_sort = ref_set_file_name_FINAL
                file_to_produce_1 = thinned_ref_set_file_name + file_type
                file_to_produce_2 = thinned_ref_set_file_name + '_no_vars' + file_type

            if platform.system() == 'Linux':
                cmd_line_string_1 = "python Pareto.py " + file_to_sort + " -o " + obj_range + " -e " + epsilon_strings \
                                    + "--output " + file_to_produce_1 + " --delimiter=\" \" --comment=\"#\" --blank"
                cmd_line_string_2 = "python Pareto.py " + file_to_sort + " -o " + obj_range + " -e " + epsilon_strings \
                                    + "--output " + file_to_produce_2 + " --delimiter=\" \" --comment=\"#\" " \
                                                                                "--blank" + " --print-only-objectives"
                subprocess.call(cmd_line_string_1, shell=True)
                subprocess.call(cmd_line_string_2, shell=True)
                #subprocess.call("echo " + "\"#\" " + ">> " + ref_set_file_name_no_vars)  # Add # to EoL
                subprocess.call("rm " + sorting_file, shell=True)  # Remove the temporary Pareto.py
            elif platform.system() == 'Windows':
                # get list of files names
                file_list = []  # Initialize
                if create_ref_set == 'Yes':
                    # Only add seed file to the list of files to process if the file actually exists.
                    for i in range(Borg_dict['nSeeds']):
                        if os.path.isfile('Borg_DPS_PySedSim' + str(i+1) + '.set'):
                            file_list.append(' Borg_DPS_PySedSim' + str(i+1) + '.set')
                    file_list_string = ''.join(file_list)
                elif thin_ref_set == 'Yes':
                    file_list_string = ' ' + file_to_sort
                #file_list = [' Borg_DPS_PySedSim' + str(i+1) + '.set' for i in range(Borg_dict['nSeeds'])]

                cmd_line_string_1 = "python Pareto.py" + file_list_string + " -o " + obj_range + " -e " + epsilon_strings \
                                   + "--output " + file_to_produce_1 + " --delimiter=\" \" --comment=\"#\" --blank"
                subprocess.call(cmd_line_string_1, shell=True)
                cmd_line_string_2 = "python Pareto.py" + file_list_string + " -o " + obj_range + " -e " + epsilon_strings\
                                    + "--output " + file_to_produce_2 + " --delimiter=\" \" --comment=\"#\" " \
                                                                                "--blank" + " --print-only-objectives"
                subprocess.call(cmd_line_string_2, shell=True)
                subprocess.call("del " + sorting_file, shell=True)  # Remove the temporary Pareto.py

            # Open both of the created files (one with and without objectives) and place a # at the end of the file.
            # This signifies end of file, and is compatible with requirements for processing of outputs in
            # MOEAFramework.

            ref_file_list = [file_to_produce_1, file_to_produce_2]
            for x in range(len(ref_file_list)):
                with open(ref_file_list[x], "a") as ref_file:
                    # If user thinned set and kept min/max values to re-insert, then reinsert before completing file.
                    if thin_ref_set == 'Yes' and keep_min_max == 'Yes':
                        # Loop through policies and dump objs/vars into appropriate files. Begin w/ refset file.
                        for pol in range(len(rowlist)):
                            if x == 0:
                                ref_file.write(' '.join(map(str, list(refset[rowlist[pol]]))) + "\n")
                            elif x == 1:
                                ref_file.write(' '.join(map(str, list(refset_objs[rowlist[pol]]))) + "\n")
                        # Now move to file with no variable values (only objectives).
                    # Finish out file with hastag, and close it before moving to next file.
                    ref_file.write("#")  # Add hashtag
                    ref_file.close()  # Close file to free up memory

            os.chdir(cwd)  # Return to original directory

        # Evaluate runtime performance (across seeds) of reference set, if user desires.
        if eval_runtime_perf == 'Yes':
            runtime_performance(Borg_dict, DPS_dict, output_location)
    if thin_ref_set == 'Yes':
        return Borg_dict, DPS_dict, rowlist, refset_objs, refset
    else:
        return Borg_dict, DPS_dict

def Initial_Processing(num_objs, num_dec_vars, ref_set_file_name, parse_objs=None, perc_conv=None, invert=None,
                       unit_conv=None, reverse_sign_all_objs = None):
    '''

    Purpose: To do pre-processing of reference set data points (to organize the referenece set data for plotting),
    by organizing decision variables and objective function values into arrays. Will be called directly from
    Scatter_Plot_Objectives before plotting is done.

    Args:
        num_objs: Defined in Scatter_Plot_Objectives()
        num_dec_vars: Defined in Scatter_Plot_Objectives
        ref_set_file: Defined in Scatter_Plot_Objectives
        parse_objs: Defined in Scatter_Plot_Objectives
        perc_conv: Defined in Scatter_Plot_Objectives
        invert: Defined in Scatter_Plot_Objectives
        unit_conv: Defined in Scatter_Plot_Objectives

    Returns:
        1. objective_values: Objective function values of points in the reference set. Type: numpy array. Size:
        [num_objs, num_ref_points]
        2. ref_set_array: The reference set array imported from the externally generated text file. Type: numpy array
    '''

    ref_set_array = np.loadtxt(ref_set_file_name)
    if (ref_set_array.shape[1]-num_dec_vars) < num_objs:
        # Reference set in specified file contains performance for a different number of objectives than the user has
        #  specified in the optimization worksheet. This likely occured because the user wants to reevaluate
        # performance of reference set policies for new objectives and store the resulting input file. This code
        # below allows for this to happen when Initial_Processing() is called from the reevaluation process.
        new_cols = np.zeros([ref_set_array.shape[0],2])
        ref_set_array = np.append(ref_set_array, new_cols, 1)
    objective_values = [[0 for i in range(len(ref_set_array))] for j in range(num_objs)]
    dec_var_values = [[0 for i in range(num_dec_vars)] for j in range(len(ref_set_array))]
    if reverse_sign_all_objs is None:
        make_all_objs_pos = 1  # Make all objective values positive for plotting purposes
    else:
        make_all_objs_pos = 0
        reverse_sign_all_objs = 1
    if parse_objs is None:
        col_to_pull = [i for i in range(num_objs)]
    else:
        col_to_pull = parse_objs

    perc = [1 for i in range(num_objs)]
    if perc_conv is not None:
        for obj in range(num_objs):
            if perc_conv[col_to_pull[obj]] == 'Yes':
                perc[obj] = 100

    if unit_conv is None:
        unit_conv = [1 for i in range(num_objs)]

    if invert is None:
        invert = ['No' for i in range(num_objs)]

    for point in range(len(ref_set_array)):
        for dec in range(num_dec_vars):
            # Process decision variable values
            dec_var_values[point][dec] = ref_set_array[point][dec]
        for obj in range(num_objs):
            # Process objective function values
            objective_values[obj][point] = ref_set_array[point][num_dec_vars+col_to_pull[obj]]*unit_conv[col_to_pull[
                obj]]*perc[obj]
            if make_all_objs_pos == 1:
                if objective_values[obj][point] < 0:
                    objective_values[obj][point] = -objective_values[obj][point]
            elif reverse_sign_all_objs == 1:
                # Negate all values, regardless of sign. Used for parallel axis plots.
                objective_values[obj][point] = -objective_values[obj][point]
            if invert[col_to_pull[obj]] == 'Yes':
                objective_values[obj][point] = 100-objective_values[obj][point]

    objective_values = np.asarray(objective_values)  # Cast objective_values list as array
    dec_var_values = np.asarray(dec_var_values)  # Cast dec_var_values list as array

    return ref_set_array, objective_values, dec_var_values


def Policy_Reevaluation(ref_set_file='pysedsim_ref_set.ref', ref_set_file_no_vars = 'pysedsim_ref_set_no_vars.ref',
                        post_optimization='Yes', file_name='PySedSim_Input_Specifications.csv', dps_dict=None,
                        borg_dict=None, parallelize_policies='No', parallelize_realizations='No',
                        internal_opt_call=None, policies=None, reevaluation='Yes', plotting='No',
                        agg_sim_outputs='No', re_eval_opt_performance='No'):
    '''
    Purpose: Re-evaluates reservoir operating policies from the reference set (or from any provided set of policies).
    Reevaluations can be either stochastic or deterministic.

    This module can be called in a post-optimization setting. If you desire to run the policy re-evaluations in
    parallel in Linux using mpi4py, then you can call this method using mpirun. If instead you want to re-evaluate
    policies in serial on one processor, then call this method without the mpirun prefix.

    :param ref_set_file: A string representing the name of the reference set file or optimization result. This
    can be either (1) a file name only (e.g., 'file_name.ref'), in which case it must be located in the output storage
    directory specified in the input specifications file; or (2) a string containing the file path (e.g.,
    r'C:/dir1/dir2/file_name.ref'). The latter option should be used if post_optimization='No'.

    :param ref_set_file_no_vars: A string representing the name of the reference set file or optimization result,
    where no variable values are stored in the file. This can be either (1) a file name only (e.g., 'file_name.ref'),
    in which case it must be located in the output storage directory specified in the input specifications file; or (
    2) a string containing the file path (e.g., r'C:/dir1/dir2/file_name.ref'). The latter option should be used if
    post_optimization='No'.

    :param DPS_dict: DPS_dict as defined in direct_policy_search.Import_DPS_Preferences. Must be specified if this
    file is called internally from optimization.
    :param Borg_dict: Borg_dict as defined in direct_policy_search.Import_DPS_Preferences. Must be specified if this
    file is called internally from optimization.
    :param post_optimization: 'Yes' or 'No' (string), indicating whether this file is to be executed internally
    within an optimization run or batch of optimization runs, in which case all input data can be passed in directly,
    or of instead all relevant input data need to be imported.
    :param file_name: [OPTIONAL] is the name of the input file (as string), which must be a .csv file located in the
    same directory as the file that contains this pysedsim.py file. (Default: 'PySedSim_Input_Specifications.csv').
    :param parallelize_policies: [OPTIONAL] String, 'Yes' or 'No. Whether or not to parallelize in Linux the
    re-evaluations using  mpi4py. Using this approach will result in parallelization over the policies, not over the
    realizations.
    :param parallelize_realizations: [OPTIONAL] String, 'Yes' or 'No. Whether or not to parallelize in Linux the
    re-evaluations using  mpi4py. Using this approach will result in parallelization over the realizations,
    not over the policies.
    :param internal_opt_call: Optional. List containing one string that represents the name of the current optimization
    scenario from which this pysedsim function is being called. Example: ['Scenario 1']. Ensures that if
    pysedsim is called from within an optimzation loop (over multiple scenarios) for purposes of re-evaluating
    policies from a reference set, only the current optimization scenario gets re-evalutad here, rather than all
    scenarios listed in the input file. Will be fed into pysedsim call if not=None.
    Returns:
    :param policies: Optional. List of integers corresponding to row numbers of policies in the array of decision
    variable values imported from the reference set (dec_var_values, as shown below). Example: [19,3, 332]. List
    should be nested if multiple scenarios are being simulated. Example: [[4,232,1], [4, 23, 2323, 532]]. This list
    can also be specified in the "Reevaluation" worksheet of the input data file. If policies are not specified,
    all policies from the reference set will be reevaluated.
    :param reevaluation: Optional. 'Yes' or 'No'. Indicates whether or not to stochastically reevaluate specific
    policies.
    :param plotting: Optional. 'Yes' or 'No'. Indicates whether probability plots of reevaluated policies should be
    created.
    :param agg_sim_outputs: Optional. 'Yes' or 'No'. Aggregates outputs from simulations distributed across processors.
    :param re_eval_opt_performance: Optional. 'Yes' or 'No'. Indicates whether user wants to reevaluate performance of
    all policies in the pareto set. Most common use: A pareto set is generated from stochastic optimization, but actual
    performance is noisy because it reflects best performance for one particular stochastic realization,
    in circumstance where optimization is only for one realization. Or, used in cases where you are reporting a
    quantile from a distribution for a stochastic optimization as a representation of performance, but now want to
    report some other quantile for the policies in the reference set (i.e. the mean). To use this function,
    all policies from the pareto set must already have been reevaluated and results stored.

    '''

    # Imports
    from optimization_pysedsim import Import_Optimization_Preferences
    from pysedsim import PySedSim
    from data_processing import Import_Simulation_Output
    from pysedsim_plotting import Probability_Plot
    from cluster_output_parser import aggregate_sim_outputs
    from performance import Performance

    os_fold = Op_Sys_Folder_Operator()

    # Must load top level input file and relevant specifications if information cannot be passed directly through an
    # existing/ongoing simulation/optimization run.
    if post_optimization == 'Yes':
        [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
         Monte_Carlo_Parameters_File_List, external_mc_data] = Determine_Num_Scenarios(file_name)
    else:
        # User is calling this module directly from within an optimization call, so Borg/DPS data have already been
        # provided as inputs to this method.
        num_scenarios = 1

    # Loop over separate simulation/optimization "scenarios", re-evaluating each one.
    scenario_reeval_dict = {}  # Stores re_eval_dict for each scenario. Keys are scenario names.
    scenario_opt_dict_re_eval = {}  # Stores opt_dict_re_eval for each scenario. Keys are scenario names.
    policies_dict = {}
    comm = None  # Indicating MPI environment has not yet been initialized (in case of parallel processing)
    ref_set_dict = {'ref_set_array': None, 'objective_values': None, 'dec_var_values': None}
    for j in range(num_scenarios):
        re_eval_dict = None
        if post_optimization == 'Yes':
            simulation_title = simulation_titles_list[j]
            Input_Data_File = Load_Input_File(simulation_title, main_input_files_dir, imported_specs)
            # Read in preferences related to variable re-evaluation, from "Reevaluation" sheet.
            try:
                # User wishes to re-evaluate
                obj_row = 8
                opt_row = 16
                sim_type = Input_Data_File['Reevaluation']['B1'].value
                num_reevals = Input_Data_File['Reevaluation']['B2'].value
                num_reeval_vars = Input_Data_File['Reevaluation']['B3'].value
                policies = Input_Data_File['Reevaluation']['B4'].value
                figure_folder_name = Input_Data_File['Reevaluation']['B5'].value
                num_reeval_procs = Input_Data_File['Reevaluation']['B6'].value

                if figure_folder_name is None:
                    figure_folder_name = 'Policy Probability Plots'

                if type(policies) in [str, unicode]:
                    # User has specified a list of numbers, which will be imported as one long string with commas
                    policies = policies.split(', ')  # Create list of policies as strings
                    policies = [int(policies[i]) for i in range(len(policies))]  # Convert list elements to integers
                else:
                    if policies is not None:
                        policies = [int(policies)]  # Make it a list
                policies_dict[simulation_title] = policies
                re_eval_dict = {}
                for var in range(num_reeval_vars):
                    key = Input_Data_File['Reevaluation'].cell(row = obj_row+var, column=1).value
                    re_eval_dict[key] = {}
                    re_eval_dict[key]['State Variable Name'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var,
                                                                                           column=2).value
                    re_eval_dict[key]['Resample Frequency'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var,
                                                                                           column=3).value
                    re_eval_dict[key]['Resample Stat'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var, column=4).value
                    re_eval_dict[key]['Location'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var, column=7).value
                    re_eval_dict[key]['unit_conv'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var,
                                                                                          column=8).value
                    if re_eval_dict[key]['unit_conv'] is None:
                        re_eval_dict[key]['unit_conv'] = 1
                    re_eval_dict[key]['perc_conv'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var, column=9).value
                    re_eval_dict[key]['invert'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var, column=10).value
                    re_eval_dict[key]['plot_type'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var,column=11).value
                    re_eval_dict[key]['axis_range'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var,
                                                                                          column=12).value
                    if re_eval_dict[key]['axis_range'] is not None:
                        re_eval_dict[key]['axis_range'] = re_eval_dict[key]['axis_range'].split(', ')
                        re_eval_dict[key]['axis_range'] = [float(re_eval_dict[key]['axis_range'][i]) for i in range(len(
                            re_eval_dict[key]['axis_range']))]  # Convert list elements to integers
                    re_eval_dict[key]['plot_title'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var,
                                                                                          column=13).value

                    re_eval_dict[key]['num_reevals'] = num_reevals
                    re_eval_dict[key]['sim_type'] = sim_type
                    re_eval_dict[key]['num_reeval_procs'] = num_reeval_procs

                    re_eval_dict[key]['Time Slice'] = Input_Data_File['Reevaluation'].cell(row = obj_row+var, column=5).value
                    Sim_Dur = (Input_Data_File['Simulation Specifications']['B4'].value - Input_Data_File['Simulation Specifications']['B3'].value).days + 1
                    if type(re_eval_dict[key]['Time Slice']) not in [str, unicode, None]:
                        # Slice is a number. Will determine how to use it in performance.py
                        re_eval_dict[key]['Time Slice'] = [re_eval_dict[key]['Time Slice']]  # Make it a list
                    elif type(re_eval_dict[key]['Time Slice']) in [str, unicode]:
                        if re_eval_dict[key]['Time Slice'] in ['Last', 'last']:
                            re_eval_dict[key]['Time Slice'] = [Sim_Dur-1]  # Slice at simulation duration
                        else:
                            # No other strings currently supported. Set = 1 (last value).
                            re_eval_dict[key]['Time Slice'] = [Sim_Dur-1]
                    else:
                        re_eval_dict[key]['Time Slice'] = None  # No value specified, so don't do a time slice at all.
                scenario_reeval_dict[simulation_title] = re_eval_dict
            except KeyError:
                raise KeyError("ERROR: Required optimization worksheet 'Optimization' does not exist in input file for "
                       "scenario %s" % simulation_title)
            internal_opt_call = [simulation_title]  # Feed into pysedsim, so scenario looping doesnt also happen there.
            borg_dict = Import_Optimization_Preferences(simulation_title, imported_specs, main_input_files_dir)
            if borg_dict['optimization approach'] == 'DPS':
                dps_dict = direct_policy_search.Import_DPS_Preferences(simulation_title=simulation_title,
                                                                       imported_specs=imported_specs,
                                                                       main_input_files_dir=main_input_files_dir)
                borg_dict['total_vars'] = dps_dict['total_vars']  # Store num decision variables in Borg_dict.
            # Get reference set file name, but append file path so it includes output storage directory.
            ref_set_file_loc = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets'
            ref_set_file_name = ref_set_file_loc + os_fold + ref_set_file
            ref_set_file_name_no_vars = ref_set_file_loc + os_fold + ref_set_file_no_vars

            # Store optimization state variable names.
            opt_dict_re_eval = {}
            for var in range(borg_dict['n_objs']):
                try:
                    key = Input_Data_File['Optimization'].cell(row = opt_row + var, column=1).value
                    opt_dict_re_eval[key] = {}
                    opt_dict_re_eval[key]['State Variable Name'] = Input_Data_File['Optimization'].cell(
                        row=opt_row + var, column=2).value
                    opt_dict_re_eval[key]['Location'] = Input_Data_File['Optimization'].cell(row=opt_row + var,
                                                                                             column=9).value
                    scenario_opt_dict_re_eval[simulation_title] = opt_dict_re_eval
                except IndexError:
                    pass
        num_dec_vars = dps_dict['total_vars']
        num_objs = borg_dict['n_objs']

        # Extract decision variable values from reference set (optimization results)
        [ref_set_array, objective_values, dec_var_values] = Initial_Processing(num_objs, num_dec_vars,
                                                                               ref_set_file_name)
        if policies is None:
            # User did not specify particular policies to reevaluate. So reevaluate all policies from reference set.
            n_policies = len(dec_var_values)  # Number of pareto solutions (operating policies) to re-evaluate
            if parallelize_policies == 'Yes':
                # Call PySedSim model for (potentially multiple) policies on each processor.
                if comm is None:
                    # Need to initialized MPI environment within parallelize_pysedsim
                    [start, stop, rank, comm] = parallelize_pysedsim(n_policies)  # Distribute re-evaluation to processors.
                else:
                    # MPI environment has already been initialized. Distribute reevaluation to processors.
                    [start, stop, rank, comm] = parallelize_pysedsim(n_policies, comm=comm)
                policy_list = range(start, stop)
            else:
                # Each processor will evaluate every policy.
                policy_list = range(n_policies)
                if parallelize_realizations == 'Yes':
                    if comm is None:
                        # Need to initialized MPI environment within parallelize_pysedsim
                        [start, stop, rank, comm] = parallelize_pysedsim(num_reevals)  # Distribute re-evaluation
                    else:
                        # MPI environment has already been initialized. Distribute reevaluation to processors.
                        [start, stop, rank, comm] = parallelize_pysedsim(num_reevals, comm=comm)
        else:
            # Only reevaluate subset list of policies specified by user.
            if len(policies) > 1 and len(policies) == num_scenarios:
                # Nested list provided. Each sub-list is a list of policies for each scenario.
                policy_list = policies[j]  # Pull specified list of policies for this scenario
            else:
                policy_list = policies
            n_policies = len(policies)  # Number of pareto solutions (operating policies) to re-evaluate
            if parallelize_policies == 'Yes':
                # Call PySedSim model for (potentially multiple) policies on each processor.
                if comm is None:
                    # Need to initialized MPI environment within parallelize_pysedsim
                    [start, stop, rank, comm] = parallelize_pysedsim(n_policies)  # Distribute re-evaluation
                else:
                    # MPI environment has already been initialized. Distribute reevaluation to processors.
                    [start, stop, rank, comm] = parallelize_pysedsim(n_policies, comm=comm)
                policy_list = policy_list[start:stop]  # Slice off the parts of policy_list to be done by this processor
            else:
                # Each processor will evaluate every policy.
                if parallelize_realizations == 'Yes':
                    if comm is None:
                        # Need to initialized MPI environment within parallelize_pysedsim
                        [start, stop, rank, comm] = parallelize_pysedsim(num_reevals)  # Distribute re-evaluation
                    else:
                        # MPI environment has already been initialized. Distribute reevaluation to processors.
                        [start, stop, rank, comm] = parallelize_pysedsim(num_reevals, comm=comm)
        # Run PySedSim, and have it produce output for the re-evaluation variables (specified in the "Reevaluation"
        # sheet of the input file) by externally adding variables to the variable export list.
        if reevaluation == 'Yes':
            for i in policy_list:
                if parallelize_policies == 'Yes':
                    PySedSim(file_name = file_name, dps_dict = dps_dict, decision_vars=dec_var_values[i],
                             scenario_name=internal_opt_call, re_eval=re_eval_dict, policy_name=i)
                elif parallelize_realizations == 'Yes':
                    PySedSim(file_name = file_name, parallel_label = rank, start_stop = [start, stop], dps_dict = dps_dict,
                             decision_vars=dec_var_values[i], scenario_name=internal_opt_call, re_eval=re_eval_dict,
                             policy_name=i)
                else:
                    # No parallelization is taking place
                    PySedSim(file_name = file_name, dps_dict = dps_dict, decision_vars=dec_var_values[i],
                             scenario_name=internal_opt_call, re_eval=re_eval_dict, policy_name=i)

        ref_set_dict[simulation_title] = {'ref_set_array': ref_set_array, 'objective_values': objective_values,
                                          'dec_var_values': dec_var_values, 'num_dec_vars': num_dec_vars,
                                          'ref_set_file_name': ref_set_file_name, 'ref_set_file_name_no_vars':
                                          ref_set_file_name_no_vars, 'n_policies': n_policies, 'num_objs': num_objs,
                                          'Sim_Dur': Sim_Dur}

    # Realizations have been distributed across processors. Need to execute code to aggregate the
    # outputs from the different processors into a single file for each variable/location,
    # so it can be imported for probability plotting purposes. This should not be done in a loop of scenarios if
    # being done in parallel.
    if agg_sim_outputs == 'Yes':
        rep_scenario = simulation_titles_list[0]  # Representative scenario to grab input file of variable names.
        num_reeval_procs = scenario_reeval_dict[rep_scenario][scenario_reeval_dict[rep_scenario].keys()[0]]['num_reeval_procs']
        [var_sub_list, Locations_to_Import] = vars_locs(simulation_titles_list, scenario_reeval_dict)
        aggregate_sim_outputs(var_sub_list=var_sub_list, Locations_to_Import=Locations_to_Import,
                              policies=policy_list, num_reeval_procs=num_reeval_procs)

    if plotting == 'Yes':
        [var_sub_list, Locations_to_Import] = vars_locs(simulation_titles_list, scenario_reeval_dict)
        # Re-evaluations are completed. Import the outputs from these simulations.
        # Create a nested list in the event there are multiple scenarios with multiple policies in each scenario.
        if len(simulation_titles_list) > 1:
            policies_list = [[] for i in range(len(simulation_titles_list))]
            for sim in range(len(simulation_titles_list)):
                policies_list[sim] = policies_dict[simulation_titles_list[sim]]
        else:
            policies_list = policies

        # Import all relevant data for scenarios, system locations and variables listed above
        [Time_Series_Import_Dictionary, Num_Realizations, Num_Years, TSID_key_list] = Import_Simulation_Output(simulation_titles_list,
                                                                                                Locations_to_Import,
                                                                                                var_sub_list,
                                                                                                main_output_file_dir,
                                                                                                policies=policies_list)

        # Create new Locations_to_Import and Sims_to_Import lists to reflect the new sub-scenarios that have been
        # added that correspond to the selected policies within each scenario.
        Locations_to_Import_NEW = {}
        simulation_titles_list_NEW = TSID_key_list
        for sim in simulation_titles_list_NEW:
            Locations_to_Import_NEW[sim] = Time_Series_Import_Dictionary[sim].keys()

    if re_eval_opt_performance == 'Yes':
        for scen_num in range(len(simulation_titles_list)):
            simulation_titles_list_opt = [simulation_titles_list[scen_num]]
            [var_sub_list, Locations_to_Import] = vars_locs(simulation_titles_list_opt, scenario_opt_dict_re_eval)

            # Re-evaluations are completed. Import the outputs from these simulations.
            # Create a nested list in the event there are multiple scenarios with multiple policies in each scenario.
            if len(simulation_titles_list_opt) > 1:
                policies_list = [[] for i in range(len(simulation_titles_list_opt))]
                for sim in range(len(simulation_titles_list)):
                    policies_list[sim] = policies_dict[simulation_titles_list_opt[sim]]
            else:
                policies_list = policy_list

            # Import all relevant data for scenarios, system locations and variables listed above
            [Time_Series_Import_Dictionary, Num_Realizations, Num_Years, TSID_key_list] = Import_Simulation_Output(simulation_titles_list_opt,
                                                                                                    Locations_to_Import,
                                                                                                    var_sub_list,
                                                                                                    main_output_file_dir,
                                                                                                    policies=policies_list)

            # Create new Locations_to_Import and Sims_to_Import lists to reflect the new sub-scenarios that have been
            # added that correspond to the selected policies within each scenario.
            Locations_to_Import_NEW = {}
            simulation_titles_list_NEW = TSID_key_list
            for sim in simulation_titles_list_NEW:
                Locations_to_Import_NEW[sim] = Time_Series_Import_Dictionary[sim].keys()

            # Loop through policies, reevaluate performance across objectives for policy,
            policy_performance = np.zeros([ref_set_dict[simulation_titles_list_opt[0]]['n_policies'], ref_set_dict[
                simulation_titles_list_opt[0]]['num_objs']])

            for p in range(ref_set_dict[simulation_titles_list_opt[0]]['n_policies']):
                # Store performance in both reference set files.
                # Make call to performance.py
                # Store the performance, BY policy, back in BOTH the ref_set and ref_set_no_vars sheets.
                performance = Performance(Time_Series_Import_Dictionary[simulation_titles_list_NEW[p]],
                                          ref_set_dict[simulation_titles_list_opt[0]]['Sim_Dur'],
                                          opt_dict=borg_dict['opt_dict'],
                                          sim_title=simulation_titles_list_opt[0])
                policy_performance[p] = np.asarray(performance)
            # Store decision variable values
            new_ref_set_array = np.zeros([ref_set_dict[simulation_titles_list_opt[0]]['ref_set_array'].shape[0],
                                          ref_set_dict[simulation_titles_list_opt[0]]['num_dec_vars']+ref_set_dict[
                simulation_titles_list_opt[0]]['num_objs']])# create new
            # reference set array
            #new_ref_set_array = ref_set_dict[simulation_titles_list_opt[0]]['ref_set_array']
            new_ref_set_array[:,0:ref_set_dict[simulation_titles_list_opt[0]]['num_dec_vars']] = ref_set_dict[
                simulation_titles_list_opt[0]]['dec_var_values']
            # Write over previous objective values with new objective vals from performance.py call.
            new_ref_set_array[:, ref_set_dict[simulation_titles_list_opt[0]]['num_dec_vars']:] = policy_performance

            np.savetxt(ref_set_file_loc + os_fold + 'pysedsim_ref_set_reeval.ref', new_ref_set_array,
                       delimiter=' ')
            np.savetxt(ref_set_file_loc + os_fold + 'pysedsim_ref_set_no_vars_reeval.ref', policy_performance,
                       delimiter=' ')

    if plotting == 'Yes':
        # Plot performance histograms/CDFs. Loop through scenarios.

        # Establish generic plotting preferences. Establish them in the code for each individual plot if specific preferences vary.
        hist_preferences = {'bins': 10, 'cumulative': False, 'normed': True, 'histtype': 'step', 'stacked': True,
                            'fill': True}
        plot_params = {'axes.labelsize': 18, 'legend.fontsize': 10, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
                       'text.usetex': False}  # keys removed: 'figure.figsize': [5, 5], 'font.size': 20,

        rep_scenario = simulation_titles_list[0]  # Representative scenario to grab input file of variable names.
        # Settings for all plots
        # ['Off']  # [None]  #  ['Alternative with Flushing: Wet Season', 'Alternative with Flushing:
        # Dry Season','Alternative without Flushing: Wet Season', 'Alternative without Flushing: Dry Season','Sambor Original: Wet Season', 'Sambor Original: Dry Season']

        legend_labels = ['MC-1-Policy 1', 'MC-1-Policy 2', 'MC-2-Policy 1', 'MC-2-Policy 2', 'MC-5-Policy 1',
                         'MC-5-Policy 2', 'MC-10-Policy 1', 'MC-10-Policy 2']
        plt_colors = ['r', 'r', 'g', 'g', 'k', 'k', 'b', 'b']  # [None]
        plt_line_style = ['-', '-', '-', '-', '-', '-', '-', '-']

#        plt_colors = ['r', 'g', 'b', 'k']  # [None]
#        plt_line_style = ['-', '-', '-', '-']
#        legend_labels = ['Deterministic Opt.: Max Energy Policy', 'Deterministic Opt.: Compromise Policy',
#                         'Stochastic Opt.: Max 1% Energy Policy', 'Stochastic Opt.: Compromise Policy']

#        plt_colors = ['r', 'g', 'r', 'g']  # [None]
#        plt_line_style = ['-', '-', '--', '--']
#        legend_labels = ['Max Energy Policy: Deterministic Reeval.', 'Compromise Policy: Deterministic Reeval.',
#                         'Max Energy Policy: Stochastic Reeval.', 'Compromise Policy: Stochastic Reeval.']

        plot_text = None  # '(A)'
        plot_text_position = None  # [.875, 0.05]

        # Loop through variables and create a plot (histogram or CDF) across scenarios for each of those variables.
        for var_name in scenario_reeval_dict[rep_scenario]:
            save_image_as = var_name
            file_loc = main_output_file_dir + os_fold + figure_folder_name
            x_label = var_name
            if scenario_reeval_dict[rep_scenario][var_name]['plot_title'] is not None:
                plot_title = scenario_reeval_dict[rep_scenario][var_name]['plot_title']
            else:
                plot_title = var_name
            axis_range = scenario_reeval_dict[rep_scenario][var_name]['axis_range']
            var_to_plot = scenario_reeval_dict[rep_scenario][var_name]['State Variable Name']
            resample_frequency = scenario_reeval_dict[rep_scenario][var_name]['Resample Frequency']
            resample_how = scenario_reeval_dict[rep_scenario][var_name]['Resample Stat']
            units_conversion = scenario_reeval_dict[rep_scenario][var_name]['unit_conv']
            plot_type = scenario_reeval_dict[rep_scenario][var_name]['plot_type']
            after_n_days = scenario_reeval_dict[rep_scenario][var_name]['Time Slice']
            scenarios_list = []
            for scenario_num in range(len(simulation_titles_list_NEW)):
                scenarios_list.append(
                    simulation_titles_list_NEW[scenario_num] + " " + scenario_reeval_dict[rep_scenario][var_name][
                        'Location'])

            [RETURNED_master, df3] = Probability_Plot(Time_Series_Import_Dictionary, Locations_to_Import_NEW, Num_Realizations,
                                                      Num_Years, Sims_to_Import=simulation_titles_list_NEW,
                                                      units_conversion=units_conversion, var_to_plot=var_to_plot,
                                                      resample_frequency=resample_frequency, resample_how=resample_how, save_type='png',
                                                      plot_title=plot_title, save_image_as=save_image_as,
                                                      plot_specs=plot_params, file_loc=file_loc, plot_type=plot_type,
                                                      hist_preferences=hist_preferences, plt_colors=plt_colors, x_label=x_label,
                                                      legend_labels=legend_labels, plot_text=plot_text,
                                                      plot_text_position=plot_text_position,
                                                      after_n_days=after_n_days, scenarios_list=scenarios_list,
                                                      plt_line_style=plt_line_style, axis_range=axis_range)

    if plotting == 'Yes':
        return ref_set_array, objective_values, dec_var_values, RETURNED_master, df3
    else:
        return ref_set_array, objective_values, dec_var_values


def vars_locs(simulation_titles_list, scenario_reeval_dict):
    '''

    Purpose: Creates a list of pysedsim state variables, and dictionary of scenario names/locations to import.

    Args:
        simulation_titles_list: list of scenarios being simulated
        scenario_reeval_dict: dictionary of reevaluation specifications created in Policy_Reevaluation()

    Returns:
        var_sub_list: list of pysedsim state variables to import
        Locations_to_Import: dictionary, where keys are scenario names, and corresponding entry is a list of system
        locations for which data should be imported.
    '''
    var_sub_list = []
    Locations_to_Import = {}
    for scenario in simulation_titles_list:
        Locations_to_Import[scenario] = []
        for var_name in scenario_reeval_dict[scenario]:
            if scenario_reeval_dict[scenario][var_name]['State Variable Name'] not in var_sub_list:
                var_sub_list.append(scenario_reeval_dict[scenario][var_name]['State Variable Name'])
            if scenario_reeval_dict[scenario][var_name]['Location'] not in Locations_to_Import[scenario]:
                Locations_to_Import[scenario].append(scenario_reeval_dict[scenario][var_name]['Location'])
    return var_sub_list, Locations_to_Import

def ref_set_plot_loop(parallel = 'No', file_name = None):
    '''

    Purpose: Generates all possible permutations of 2-D tradeoff plots (with a third dimension being color,
    if desirable) from a reference set of many objectives. Only works for more than one objective.

    For example, if there are 3 objectives, 6 plots are produced, where 2 plots are possible for a given two
    objectives, because they can be plotted on different axes.

    :param parallel: String, either 'Yes' or 'No', indicates whether user will be calling ref_set_plot in parallel
    using mpi4py. The parallelization is over the number of plots to be produced, rather than over the potentially
    multiple scenarios of interest.
    :param file_name: Top-level pysedsim input file (.csv) that contains name of scenario(s) and input file(s)
    locations.
    :return:
    '''
    from optimization_pysedsim import Import_Optimization_Preferences
    from direct_policy_search import Import_DPS_Preferences

    [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
     os_fold] = import_specs(file_name=file_name)

    comm = None  # Indicating MPI environment has not yet been initialized (in case of parallel processing)

    # Loop through as many optimization scenarios as user has specified. Each processor gets a different scenario.
    for j in range(num_scenarios):
        simulation_title = simulation_titles_list[j]
        output_location = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets'
        Borg_dict = Import_Optimization_Preferences(simulation_title, imported_specs, main_input_files_dir)
        if Borg_dict['optimization approach'] == 'DPS':
            DPS_dict = Import_DPS_Preferences(simulation_title=simulation_title, imported_specs=imported_specs,
                                              main_input_files_dir=main_input_files_dir)
            Borg_dict['n_vars'] = DPS_dict['total_vars']  # Store num decision variables in Borg_dict, copied from
            # DPS_dict.

        # objs_to_plot may contain many objective variable strings. A subset of those strings of length 3 must now be
        # selected for 2D plotting, with color as a third dimension within the 2D space.
        objs_to_plot = Borg_dict['opt_dict']['Objective Names Ordered List']  # Plot all objectives in optimization.
        num_list = [i for i in range(Borg_dict['n_objs'])]
        num_objs_to_plot = min(3, Borg_dict['n_objs'])
        obj_num_list = list(itertools.permutations(num_list, num_objs_to_plot))
        ct=0
        plot_list = []
        for item in obj_num_list:
            if len(obj_num_list) == 3:
                plot_list.append([objs_to_plot[obj_num_list[ct][0]], objs_to_plot[obj_num_list[ct][1]], objs_to_plot[obj_num_list[ct][2]]])
            elif len(obj_num_list) == 2:
                plot_list.append([objs_to_plot[obj_num_list[ct][0]], objs_to_plot[obj_num_list[ct][1]]])
            ct+=1
        ref_set_pref_dict = {'num_objs': Borg_dict['n_objs'], 'num_dec_vars': DPS_dict['total_vars'],
                             'num_objs_to_plot': num_objs_to_plot}
        ref_set_pref_dict['invert'] = [i for i in range(Borg_dict['n_objs'])]
        ref_set_pref_dict['perc_conv'] = [i for i in range(Borg_dict['n_objs'])]
        ref_set_pref_dict['unit_conv'] = [i for i in range(Borg_dict['n_objs'])]
        ref_set_pref_dict['ref_set_file_name'] = output_location + os_fold + Borg_dict['opt_output_filename']

        # Loop through objective names, store preferences related to how the objective values are re-manipulated
        # after the reference set is processed.
        for x in range(len(Borg_dict['opt_dict']['Objective Names Ordered List'])):
            current_name = Borg_dict['opt_dict']['Objective Names Ordered List'][x]
            ref_set_pref_dict['invert'][x] = Borg_dict['opt_dict']['Objective Preferences'][current_name]['invert']
            ref_set_pref_dict['perc_conv'][x] = Borg_dict['opt_dict']['Objective Preferences'][current_name][
                'perc_conv']
            ref_set_pref_dict['unit_conv'][x] = Borg_dict['opt_dict']['Objective Preferences'][current_name][
                'unit_conv']

        if parallel == 'Yes':
            if comm is None:
                # Need to initialized MPI environment within parallelize_pysedsim
                [start, stop, rank, comm] = parallelize_pysedsim(len(obj_num_list))  # Distribute re-evaluation to processors.
            else:
                # MPI environment has already been initialized. Distribute reevaluation to processors.
                [start, stop, rank, comm] = parallelize_pysedsim(len(obj_num_list), comm=comm)
        else:
            [start, stop] = [0, len(obj_num_list)]

        for plot in range(start, stop):
            #plot_dict = {'Axis Range': [[0, 6000], [0, 100], [0, 8000]], 'plot_order': obj_num_list[plot],
            #             '3d_plot': 'No', 'compute_upper_ax_lim': 'Yes'}
            plot_dict = {'plot_order': obj_num_list[plot], '3d_plot': 'No', 'compute_upper_ax_lim': 'No'}

            # Produce/save plots.
            [ref_set_array, objective_values, dec_var_values] = Scatter_Plot_Objectives(ref_set_pref_dict, objs_to_plot,
                                                                            plot_dict=plot_dict,
                                                                            save_fig= output_location)
            pyplot.clf()  # Close all figures that are open, to avoid memory issues.

def parallelize_pysedsim(n_evaluations, comm=None):
    '''
    Purpose: Sets up parallelization of pysedsim simulation, or other processes (e.g., reference set generation) in
    parallel on a linux cluster using mpi4py. This can only be used currently in a Linux setting that has mpi4py
    installed. Windows performance is not yet verified.

    :param n_evaluations: Number of evaluations of pysedsim that need to be parallelized
    returns: start: integer, (first simulation out of total=n_evaluations to run)
    returns: stop: integer, (last simulation out of total=n_evaluations to run)

    :param comm: This is the initialized MPI environment object, called with mpi4py.MPI.COMM_WORLD. Specify this
    argument only if the MPI environment has already been initialized with a previous call.

    Note: if start, stop = 3,4, only one evaluation will actually be performed, not two.
    '''

    # Use MPI4PY to parallelize the policy re-evaluations
    from mpi4py import MPI

    if comm is None:
        # Initialize MPI environment.
        comm = MPI.COMM_WORLD

    # Get the number of processors and the rank of processors
    n_procs = comm.size  # Number of processors (e.g., 512)
    rank = comm.rank  # Rank of individual processor from 0 to (nprocs-1) (e.g., 0, 1, 2, ..., 511)

    # Use the processor rank to determine the chunk of work each processor will do. If a remainder exists,
    # some processors will be required to do extra simulations.
    count = int(n_evaluations/n_procs)  # Integer number of reevaluations to be executed per processor.
    remainder = n_evaluations % n_procs
    if rank < remainder:
        start = rank*(count+1)
        stop = start+count+1
    else:
        start = count*rank + remainder
        stop = start+count
    return start, stop, rank, comm

def Grab_Pareto_Min_Max(ref_set_array, objective_values, num_objs, num_dec_vars, objectives_names=[],
                        create_txt_file='No'):

    """
    Purposes: Identifies the operating policies producing the best and worst performance in each objective.

    Gets called automatically by processing_reference_set.Reference_Set()

    Required Args:
        1. ref_set_array: an array of P arrays, P=number of points in the reference set. Each of the P
        arrays contains N=num_objs+num_vars (number of objective values in optimization problem and number of
        decision variables). Decision variables come first, followed by objective values.
        2. objective_values: The objective value portion of the ref_set_array. It is also an array of P arrays,
        where each of the P arrays is of length num_objs.
        3. num_objs = integer number of objective values (e.g., 5)
        4. num_dec_vars: integer number of decision variable values (e.g., 30)
    Optional Args:
        5. objectives_names: (list of names of objectives for objective_values returned from Reference_Set,
        as defined above). Example: ['Sediment', 'Hydropower']. Used to name .txt file and provide output dictionary
        keys.
        6. create_text_file: String of 'Yes' or 'No'. Indicates whether users wants function to produce text files of
        operating policies.
    Returns:
        1. indices of ref_set_array that correspond to the points of the highest and lowest value in each of the
        objectives.
        2. Various text files that store the DPS parameters corresponding to the operating policy, if the user wishes
        to create such files.
    """

    # Find operating policy parameters corresponding to the largest objective
    # List of index of (1) highest and (2) lowest values (column in objective_values array)
    indices = [[0 for i in range(2)] for j in range(num_objs)]
    indices_dict = {}
    for obj in range(num_objs):
        indices[obj][0] = np.argmin(objective_values[obj])  # MIN for each objective
        indices[obj][1] = np.argmax(objective_values[obj])  # MAX for each objective
        if create_txt_file == 'Yes':
            # Save max and min policies so PySedSim can import them.
            np.savetxt('RBF_Parameters_max_' + objectives_names[obj] + '.txt',
                       ref_set_array[indices[1][obj]][0:num_dec_vars], newline=' ')
            np.savetxt('RBF_Parameters_min_' + objectives_names[obj] + '.txt',
                       ref_set_array[indices[0][obj]][0:num_dec_vars], newline=' ')
        indices = np.asarray(indices)  # cast list as array
        indices_dict[objectives_names[obj]] = {'Min': indices[obj][0], 'Max': indices[obj][1]}

    return indices_dict

def brush_pareto_policies(Objective_Range, refset_file_path = 'pysedsim_ref_set.ref',
                          objs_file_path='pysedsim_ref_set_no_vars.ref', objective_union = 'No', one_best_worst = 'No'):
    '''
    Purpose: Returns solutions from the pareto set (e.g., reference set) that are within the provided range.
    Solutions with objective values exactly equaling values in the provided range will be included.

    :param Objective_Range: Dictionary, where keys are the objective value range. Example: Objective_Range = {
    'Objective 1': {'Range': [13, 5500], 'Type': 'Max'}, 'Objective 2': {'Range': [5, 18], 'Type': 'Min'}. Instead
    of a list of values for the range, the user may instead specify strings including either 'max' or 'min', indicating
    whether to keep the best performing ('max') or worst performing ('min') solution. If no range is provided and no
    'max' or 'min' string is listed, the entire range will be included. The 'Type' sub-key is required.

    :param objective_union: Optional, 'Yes' or 'No'. 'Yes' indicates user wants to return any policies that meet
    performance critieria for individual objectives (the union), whereas 'No' indicates only policies meeting
    performance criteria for all objective will be returned (the intersection).

    :param one_best_worst: Optional, 'Yes' or 'No'. 'Yes' indicates that user is seeking the best and worst policies
    to be returned, but if multiple policies perform the best in a given objective, only one policy (the first to be
    located) will be kept. 'No' indicates all policies that have the best and/or worst performance will be kept. This
    is useful, for example, if multiple policies have zero values for a given objective.

    :return:
    rowlist: a list of the row numbers of solutions in the reference set
    refset: a numpy array that stores all soutions in the reference set (only objective values).

    '''

    rowlist = []
    refset_objs = np.loadtxt(objs_file_path)
    refset = np.loadtxt(refset_file_path)
    n_objs = np.size(refset_objs, 1)
    rng = [[] for i in range(n_objs)]
    pos_neg = [1 for i in range(n_objs)]  # Set default to 1. Used to handle assumed negativity of maximized
    # objective values.
    # Loop through objectives, and store provided range. If none provided, use min/max of column in refset_objs for that
    # variable.
    num_brushes = [1 for i in range(n_objs)]
    for obj in range(n_objs):
        if Objective_Range["Objective %s" % obj]['Type'] == 'Max':
            pos_neg[obj] = -1
        try:
            # Find max or min value and make that the range.
            # Account for the assumed negativity of objective values that are being maximized. Make both
            # values the same, since user wants to locate max or min only.
            if type(Objective_Range["Objective %s" % obj]['Range'][0]) in [str, unicode]:
                # User may be requesting both the max and min performing policies to be included.
                num_brushes[obj] = len(Objective_Range["Objective %s" % obj]['Range'])  # Reset value if min/max
                for brush_range_num in range(num_brushes[obj]):
                    if Objective_Range["Objective %s" % obj]['Range'][brush_range_num] in ['max', 'Max']:
                        # User wants to have the best performing objective value be retained. rng variable stores actual
                        # values (with negative sign removed)
                        if Objective_Range["Objective %s" % obj]['Type'] == 'Max':
                            # If maximization objective
                            rng[obj].append(-1*min(refset_objs[:, obj]))
                            rng[obj].append(-1*min(refset_objs[:, obj]))
                        else:
                            # If minimization objective
                            rng[obj].append(min(refset_objs[:, obj]))
                            rng[obj].append(min(refset_objs[:, obj]))
                    if Objective_Range["Objective %s" % obj]['Range'][brush_range_num] in ['min', 'Min']:
                        # User wants to retain the worst performing objective values
                        if Objective_Range["Objective %s" % obj]['Type'] == 'Max':
                            # If maximization objective
                            rng[obj].append(-1*max(refset_objs[:, obj]))
                            rng[obj].append(-1*max(refset_objs[:, obj]))
                        else:
                            # If minimization objective
                            rng[obj].append(max(refset_objs[:, obj]))
                            rng[obj].append(max(refset_objs[:, obj]))
            else:
                # User has provided an actual range of objective values. Replace default empty list with user list.
                rng[obj] = Objective_Range["Objective %s" % obj]['Range']
        except KeyError:
            # No range provided for this objective. Range is entire objective range.
            if Objective_Range["Objective %s" % obj]['Type'] == 'Max':
                # Account for the assumed negativity of objective values that are being maximized.
                rng[obj].append(-1*max(refset_objs[:, obj]))  # Lower value of brushed range
                rng[obj].append(-1*min(refset_objs[:, obj]))  # Upper value of brushed range
            else:
                # Objective is a minimization objective
                rng[obj].append(min(refset_objs[:, obj]))  # Lower value of brushed range
                rng[obj].append(max(refset_objs[:, obj]))  # Upper value of brushed range

    # Determine which policies (rows) meet individual objective criteria.
    if objective_union == 'Yes':
        for obj in range(n_objs):
            for n_brush in range(num_brushes[obj]):
                for row in range(len(refset_objs)):
                    if (pos_neg[obj] * refset_objs[row][obj] >= rng[obj][2*n_brush]) and (
                            pos_neg[obj] * refset_objs[row][obj] <= rng[obj][2*n_brush + 1]):
                        rowlist.append(row)  # Performance criterion satisfied. Add to list.
                        if one_best_worst == 'Yes':
                            # If user wants best/worst performing policies, and multiple have same performance,
                            # only keep one. Break for loop and move to next brush.
                            break
    else:
        # Return list of policies (row numbers) that meet the specified objective range criteria across all objectives.
        for row in range(len(refset_objs)):
            for obj in range(n_objs):
                if (pos_neg[obj] * refset_objs[row][obj] < rng[obj][0]) or (
                        pos_neg[obj] * refset_objs[row][obj] > rng[obj][1]):
                    break  # Condition for this variable not met. Exit and skip to next row (policy) in array.
                if obj == n_objs-1:
                    rowlist.append(row)  # All columns satisfied criteria. This policy meets all criteria. Add to list.
    return rowlist, refset_objs, refset

def runtime_performance(Borg_dict, DPS_dict, output_location, input_file_name = 'PySedSim_Input_Specifications.csv'):
    '''

    Purpose: Evaluates runtime performance (hypervolume, generational distance, etc.) across seeds using built-in
    functionality in MOEA Framework.

    Some code borrowed from Jazmin Zatarain-Salazar, Cornell University, blog post:
    https://waterprogramming.wordpress.com/2015/07/03/random-seed-analysis-for-the-borg-moea-using-dtlz2-3-objective
    -instance/

    :param Borg_dict:
    :param input_file_name:
    :return:
    '''

    # Import relevant preferences for the particular optimization run

    # 1. Produce runtime file that only contains objective values, with snapshots at each NFE separated by a hashtag.
    num_seeds = Borg_dict['nSeeds']
    runtime_interval = Borg_dict['runtime_preferences']['runtime_freq']
    os_fold = Op_Sys_Folder_Operator()
    seed_counter = 0  # used to make sure num_seeds worth of runtime and seed files actually exist.
    for seed in range(num_seeds):
        #print output_location
        #print output_location + os_fold + 'runtime_file_seed_%s' % str(seed) + '.runtime'
        try:
            f = open(output_location + os_fold + 'runtime_file_seed_%s' % str(seed+1) + '.runtime', "r")
            lines = f.readlines()
            f.close()

            for line in range(len(lines)):
                if lines[line][0] not in ["/", "#"]:
                    lines[line] = lines[line].split(' ')[DPS_dict['total_vars']:DPS_dict['total_vars']+Borg_dict['n_objs']]

            f = open(output_location + os_fold + 'runtime_file_seed_%s' % str(seed+1) + '.objs', "w")

            for line in lines:
                if line[0] != "/":
                    f.write(" ".join(line))
            f.close()
            seed_counter+=1
        except IOError:
            num_seeds = seed_counter  # Count number of runtime files, in case less seeds than spec'd in input file.
            break  # stop for loop that loops through seeds, as fewer seed files are available than user indicated.

    # Compute absolute runtime metrics (e.g., absolute hypervolume) for the reference set
	cwd = os.getcwd()  # Current working directory to return to after changing directories.
	os.chdir(output_location)  # Go to folder containing runtime file and reference set file.
	cmd_line_string_1 = "java -cp MOEAFramework-2.0-Executable.jar HypervolumeEval pysedsim_ref_set_no_vars.ref >> pysedsim_ref_set.metrics"
	subprocess.call(cmd_line_string_1, shell=True)
	os.chdir(cwd)  # Return to original directory

    # 2. Produce runtime metrics (for each seed, 6 metrics computed for every NFE value)
    for seed in range(num_seeds):
        cwd = os.getcwd()  # Current working directory to return to after changing directories.
        os.chdir(output_location)  # Go to folder containing runtime file and reference set file.
        cmd_line_string_1 = "java -cp MOEAFramework-2.0-Executable.jar " \
                            "org.moeaframework.analysis.sensitivity.ResultFileEvaluator -d " + str(Borg_dict['n_objs'])
        cmd_line_string_2 = " -i runtime_file_seed_%s.objs" % str(seed+1) + " -r pysedsim_ref_set_no_vars.ref"
        cmd_line_string_3 = " -o runtime_metrics_seed_%s.metrics" % str(seed+1)
        subprocess.call(cmd_line_string_1 + cmd_line_string_2 + cmd_line_string_3, shell=True)
        os.chdir(cwd)  # Return to original directory

    # 3. Produce average runtime metrics
    cwd = os.getcwd()  # Current working directory to return to after changing directories.
    os.chdir(output_location)  # Go to folder containing runtime file and reference set file.
    cmd_line_string_1 = "java -cp MOEAFramework-2.0-Executable.jar " \
                        "org.moeaframework.analysis.sensitivity.SimpleStatistics -m average -o pysedsim_run.average " + \
                        "runtime_metrics_seed_*.metrics"
    subprocess.call(cmd_line_string_1, shell=True)
    os.chdir(cwd)  # Return to original directory

    # 4. Produce plots of runtime metrics vs. NFE for each seed and for average.

    # Plot average average_metrics, as well as average_metrics by seed.
    average_metrics = np.loadtxt(output_location + os_fold + 'pysedsim_run.average', delimiter= ' ')
    runtime_metrics = {}
    for s in range(num_seeds):
        runtime_metrics['Seed' + str(s + 1)] = np.loadtxt(
            output_location + os_fold + 'runtime_metrics_seed_%s.metrics' % str(s + 1), delimiter=' ')

    fig, (ax1, ax2, ax3) = pyplot.subplots(3, sharex=True)
    NFE = [runtime_interval*(i+1) for i in range(len(average_metrics[:,0]))]
    with open(output_location + os_fold + 'pysedsim_ref_set.metrics','r') as f:
        hpv_ref = f.readline().rstrip()
        hpv_ref = float(hpv_ref)

    # Hypervolume subplot
    ax1.set_ylabel('Hypervolume')
    # Plot each seed separately in blue
    for s in range(num_seeds):
        ax1.plot(NFE, (runtime_metrics['Seed' + str(s+1)][:,0])/hpv_ref, color='b')
    #ax1.set_title('Runtime output')  # Set title for each subplot
    ax1.plot(NFE, (average_metrics[:,0])/hpv_ref, color='r')  # Average (across seeds) hypervolume plot

    # Generational Distance subplot
    ax2.set_ylabel('Generational Distance')
    # Plot each seed separately in blue
    for s in range(num_seeds):
        ax2.plot(NFE, runtime_metrics['Seed' + str(s+1)][:,1], color='b')
    ax2.plot(NFE, average_metrics[:,1], color='r')  # Average (across seeds) generational distance plot

    # Epsilon Indicator subplot
    ax3.set_ylabel('Additive $\epsilon$-indicator')
    # Plot each seed separately in blue
    for s in range(num_seeds):
        ax3.plot(NFE, runtime_metrics['Seed' + str(s+1)][:,4], color='b')
    ax3.plot(NFE, average_metrics[:,4], color='r')  # Average (across seeds) epsilon indicator plot
    ax3.set_xlabel('Number of Function Evaluations (NFE)')

    # Characteristics applying to entire plot (all of subplots)
    #p.suptitle('Runtime dynamics')  # Set title
    fig.subplots_adjust(hspace=0)
    pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    pyplot.tight_layout()
    pyplot.show()
    fig.savefig(output_location + os_fold + "runtime_plot.png", dpi=1200, bbox_inches='tight')

def import_specs(file_name=None):
    # Get operator for changing directory based on operating system, then import various assumptions about the simulation
    # from the top level input file.
    os_fold = Op_Sys_Folder_Operator()
    if file_name is not None:
        [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
         Monte_Carlo_Parameters_File_List, external_mc_data] = Determine_Num_Scenarios(file_name=file_name)
    else:
        [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
         Monte_Carlo_Parameters_File_List, external_mc_data] = Determine_Num_Scenarios()

    return num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir, os_fold

def Subplot_Scatter_Plot_Objectives(ref_set_pref_dict, objs_to_plot, parse_objs=None, plot_dict=None,
                            save_fig = None, gen_op_policy_file='No', num_objs_to_plot=None,
                            tradeoff_plotting_subplot=None, axis_test = None, sp_slot_tuple=None):

    '''

    Purpose: Produces a scatter plot of objective function values up to 4D (4th dimension is color).

    Args (if not defined here, it is defined in Reference_Set() method):
        3. objs_to_plot: list of objective names as they appear in the pysedsim input file and in the reference set
        file, e.g. ['Objective 1', 'Objective 2']. Only list the names of objectives that you want to plot.
        4. gen_op_policy_file: Optional, specify 'Yes' or 'No'. default = 'No'. Calls Grab_Pareto_Min_Max()
        function to generate text file of operating policy for minimum and maximum points in the reference set.
        5. parse_objs: List of objective numbers (column numbers in 'PySedSim_reference.txt' to plot from
        reference set
        6. ref_set_pref: A dictionary containing the following keys:
        a. num_objs: Number of objective function values
        b. num_dec_vars: Number of decision variable values
        c. ref_set_file_name: name of reference set file (e.g., 'ref_set_file.txt')
        d. unit_conv: List of factors by which to multiply values in each array of objective values.
        e. perc_conv: List of 'Yes' or 'No' for each objective, for converting objective values to percent from
        fraction.
        f. invert: List of 'Yes' or 'No' for each objective, for reversing percentage (doing 100-value).
        9. plot_dict: Dictionary of plotting preferences. Keys include:
        a. 'Axis Range': List of lists (sub-lists), containing axis limits for 2 or 3 primary axes: x, y,
        z. e.g., [[2,3],[3,5],[4,5]]
        b. 'x_labal': string to label x-axis.
        c. 'y_label': string to label y-axis.
        d. 'z_label': string to label z-axis.
        e. 'c_label': string to label "colored" axis (if there is a 4th dimension).
        f. 'Plot Specifications': dictionary, follows format of rcParams module from matplotlib. An example is:

            'Plot Specifications': {'axes.labelsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
                                 'text.usetex': False, 'figure.figsize': [9, 6], 'figure.figsize': [6.5, 5],
                                 'legend.fontsize': 22}
        g. '3d_plot': 'Yes' or 'No', indicates whether user wants to plot 3 objectives in 3D space, or in 2D space
        using color as the third objective.

        10. plt_order: Optional, this is a list that specifies in which order the objectives are plotted. Used if you
        want to be able to determine which of the objectives become a color/point size/orientation in the >3D case.
        Example: plt_order = [0,3,2,1] to have the third variable in the objective values array plotted in the 3D
        space, while the third item is plotted last.

        11. num_objs_to_plot: integer number of objectives you actually wish to plot. Will be computed as length of
        objs_to_plot list if nothing is specified.

        tradeoff_plotting_subplot: indicates whether call is being made from tradeoff_plotting module for a subplot,
        in which case the figure object should be returned so it can be included in a subplot.
    Returns:
        Nothing, but several figures may be saved
    '''

    # Load basic information about the simulated scenarios, including input file directory and simulation names.
    os_fold = Op_Sys_Folder_Operator()

    # Figure fontsize for all labels
    all_label_font = 8

    # Unpack ref_set_pref_dict
    try:
        ref_set_file_name = ref_set_pref_dict['ref_set_file_name']
        num_objs = ref_set_pref_dict['num_objs']
    except KeyError:
        ref_set_file_name = None
        num_objs = None
    try:
        num_dec_vars = ref_set_pref_dict['num_dec_vars']
    except KeyError:
        pass
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

    try:
        ideal_point_loc = plot_dict['ideal_point']['coordinates']
        plot_ideal_point = 'Yes'
        ideal_point_marker = plot_dict['ideal_point']['marker_type']
        ideal_point_color = plot_dict['ideal_point']['color']
        ideal_point_size = plot_dict['ideal_point']['size']
    except KeyError:
        plot_ideal_point = 'No'

    try:
        preference_arrow = plot_dict['preference_arrow']
    except KeyError:
        preference_arrow = 'No'

    # Unpack any preferences related to colormaps, which can be used in 2D and 3D plots to represent an additional
    # objective.
    try:
        cmap = plot_dict['cmap_name']
    except KeyError:
        try:
            # No colormap specified. See if user specified a file called 'custom_colormap.txt'. If so, use it to define
            # the colormap.
            new_cmap_matrix = np.loadtxt('custom_colormap.txt')
            cmap = matplotlib.colors.ListedColormap(new_cmap_matrix/255.0)  # Must be a 256 color map.
        except IOError:
            # Use default (jet)
            cmap = "jet_r"  # 'gray_4, 'jet', 'jet_r', cool_r, 'Reds_r','Blues_r','Greens_r',
    # If only a matplotlib default colormap name has been selected, then create a colormap from it.
    if type(cmap) in [str, unicode]:
        cmap = pyplot.cm.get_cmap(cmap)

    # Unpack any preferences related to selected policies to highlight in the figure
    try:
        pols_to_highlight = plot_dict['pols_to_highlight']
        transparency_main_plot = 0.5  # 0.35  # 1  # 0.3
        alpha = [1 for x in range(len(pols_to_highlight)+1)]
        alpha[0] = transparency_main_plot
        try:
            pol_names = plot_dict['pol_names']
        except KeyError:
            pol_names = ['Policy %s' %p for p in pols_to_highlight]
        try:
            label_policies = plot_dict['label_policies']
        except KeyError:
            label_policies = 'Yes'

        try:
            gray_policies_background = plot_dict['gray_policies_background']
        except KeyError:
            gray_policies_background = 'No'  # 'Yes'
        try:
            pol_name_locs = plot_dict['pol_name_locs']
        except KeyError:
            pol_name_locs = [
                            [
                            [(37, 5000), (75, 3750), (89, 2800)],
                            [(33, 5000), (60, 5000), (89, 400)],
                            ],
                            [
                            [(33, 700), (65, 550), (70, 50)],
                            [(12, 5000), (30, 2000), (15, 900)]
                            ]
                            ]
    except KeyError:
        alpha = [1]
        pols_to_highlight = None
        label_policies = 'No'
        gray_policies_background = 'No'  # Plot all in color as default

    if ref_set_file_name is not None:
        # User wishes to import reference set here.
        [ref_set_array, objective_values, dec_var_values] = Initial_Processing(num_objs, num_dec_vars, ref_set_file_name,
                                                                               parse_objs=parse_objs, perc_conv=perc_conv,
                                                                               invert=invert, unit_conv=unit_conv)
    else:
        # User has provided reference set information.
        ref_set_array = ref_set_pref_dict['ref_set_array']
        objective_values = ref_set_pref_dict['objective_values']
        dec_var_values = ref_set_pref_dict['dec_var_values']

    # Unpack plotting dictionary
    try:
        plt_order = plot_dict['plot_order']
    except KeyError:
        plt_order = None
    try:
        ax_lim = plot_dict['Axis Range']
        # Determine axis limits
        try:
            if plot_dict['compute_upper_ax_lim'] == 'Yes':
                for obj in range(len(objective_values)):
                    ax_lim[obj][1] = np.max(objective_values[obj]) + 0.1*np.max(objective_values[obj])
        except KeyError:
            pass
    except KeyError:
        ax_lim = None
    try:
        three_d_plot = plot_dict['3d_plot']
    except KeyError:
        three_d_plot = 'No'
    try:
        plot_specs = plot_dict['Plot Specifications']
    except KeyError:
        plot_specs = None
    try:
        invert_axis = plot_dict['Invert Axis']
    except KeyError:
        invert_axis = None

    try:
        num_objs_to_plot = ref_set_pref_dict['num_objs_to_plot']
    except KeyError:
        num_objs_to_plot = len(objs_to_plot)

    if num_objs_to_plot == 1:
        print("Only 1 objective; no plot will be produced. Result is a single value.")
    elif num_objs_to_plot == 2:
        fig = pyplot.figure()
        if plt_order is None:
            # Order of plotting not specified.
            plt_order = [0,1]

        if plot_specs is not None:
            rcParams.update(plot_specs)
        pyplot.scatter(objective_values[plt_order[0]],objective_values[plt_order[1]])
        try:
            pyplot.ylabel(plot_dict['y_label'])
            pyplot.xlabel(plot_dict['x_label'])
        except KeyError:
            pyplot.ylabel(objs_to_plot[plt_order[1]])
            pyplot.xlabel(objs_to_plot[plt_order[0]])
        if ax_lim is not None:
            pyplot.xlim(ax_lim[plt_order[0]][0], ax_lim[plt_order[0]][1])
            pyplot.ylim(ax_lim[plt_order[1]][0], ax_lim[plt_order[1]][1])
        if save_fig is not None:
            save_fig = save_fig + os_fold + objs_to_plot[plt_order[0]] + '-' + objs_to_plot[
                plt_order[1]] + ".png"
    elif num_objs_to_plot >= 3:
        if three_d_plot == 'Yes':
            fig = pyplot.figure()  # create the figure
        else:
            fig, ax_arr = axis_test
        if num_objs_to_plot == 3:
            # Create a 3D plot, with no color used to map a fourth objective.
            if plt_order is None:
                # Order of plotting not specified.
                plt_order = [0, 1, 2]

            if three_d_plot == 'Yes':
                ax3D = Axes3D(fig)
                pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[
                    plt_order[2]])
            else:
                # Update plot specifications
                if plot_specs is not None:
                    rcParams.update(plot_specs)
                if ax_lim[plt_order[2]] != []:
                    norm_array = np.array(
                        [(c - ax_lim[plt_order[2]][0]) / (ax_lim[plt_order[2]][1] - ax_lim[plt_order[2]][0]) for c in
                         objective_values[plt_order[2]]])
                else:
                    low = min(objective_values[plt_order[2]])
                    high = max(objective_values[plt_order[2]])
                    norm_array = np.array([(c-low)/(high-low) for c in objective_values[plt_order[2]]])
                # Plot main tradeoffs
                if label_policies == 'Yes':
                    if gray_policies_background == 'No':
                        ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                                                      color=cmap(norm_array), marker='o', s=10,
                                                      linewidth=0, alpha=alpha[0])
                    else:
                        ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                                                      color=(0.85,0.85,0.85), marker='o', s=10, linewidth=0,
                                                      alpha=alpha[0])
                else:
                    if gray_policies_background == 'No':
                        #ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                        #                              c=objective_values[plt_order[2]], cmap=cmap, marker='o', s=10,
                        #                              linewidth=0)
                        ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                                                      color=cmap(norm_array), marker='o', s=10,
                                                      linewidth=0)
                    #color=cm.ScalarMappable(norm=norm, cmap=cm.jet)
                    else:
                        ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                                                      color=(0.85,0.85,0.85), marker='o', s=10, linewidth=0)

                # Plot individual policies on same axis if user wishes to highlight particular policies in this way.
                if pols_to_highlight is not None:
                    for x in range(len(pols_to_highlight)):
                        low_value = int(min(objective_values[plt_order[2]]))
                        high_value = int(max(objective_values[plt_order[2]]))
                        diff = int(high_value-low_value)
                        color_value = pyplot.cm.jet_r((np.clip(objective_values[plt_order[2]][pols_to_highlight[x]],
                                                           low_value,high_value)-low_value)/diff)
                        if label_policies == 'Yes':
                            # Create thick marker linewidth for policies being highlighted, and use larger size,
                            # and don't have a transparent point.
                            ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]][pols_to_highlight[x]],
                                                          objective_values[plt_order[1]][pols_to_highlight[x]],
                                                          color=cmap(norm_array[pols_to_highlight[x]]), marker='o',
                                                          s=20, linewidth=1.0, alpha=alpha[x+1], edgecolor='k')
                        else:
                            # No points being highlighted. Create no marker linewidth, size=10, and no transparency.
                            ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]][pols_to_highlight[x]],
                                                          objective_values[plt_order[1]][pols_to_highlight[x]],
                                                          color=cmap(norm_array[pols_to_highlight[x]]), marker='o',
                                                          s=10, linewidth=0, alpha = 1.0)

                        xy_loc = (objective_values[plt_order[0]][pols_to_highlight[x]], objective_values[plt_order[
                            1]][pols_to_highlight[x]])
                        xy_text_loc = pol_name_locs[sp_slot_tuple[0]][sp_slot_tuple[1]][x]
                        if label_policies == 'Yes':
                            ax_arr[sp_slot_tuple].annotate(pol_names[x], fontsize=6, xy=xy_loc, xycoords='data',
                                        xytext=xy_text_loc, bbox=dict(boxstyle="round", fc='w'),
                                        arrowprops=dict(arrowstyle="->",connectionstyle="arc,angleA=0,angleB=90,"
                                                                                        "rad=10"))
                        #, fc="0.8"

                # Plot ideal point if user wants.
                if plot_ideal_point == 'Yes':
                    if plot_dict['ideal_point']['coordinates'] == 'best':
                        # Plot point that corresponds to best objective values across objectives
                        plot_dict['ideal_point']['coordinates'] = []
                        # X-axis coordinate for best point
                        if plot_dict['objective_type'][0] in ['max', 'Max']:
                            plot_dict['ideal_point']['coordinates'].append(max(objective_values[plt_order[0]]))
                        else:
                            plot_dict['ideal_point']['coordinates'].append(min(objective_values[plt_order[0]]))
                        # Y-axis coordinate for best point
                        if plot_dict['objective_type'][1] in ['max', 'Max']:
                            plot_dict['ideal_point']['coordinates'].append(max(objective_values[plt_order[1]]))
                        else:
                            plot_dict['ideal_point']['coordinates'].append(min(objective_values[plt_order[1]]))
                        ideal_point_loc = plot_dict['ideal_point']['coordinates']
                    else:
                        # Just plot axis limits in absence of other information.
                        if type(plot_dict['ideal_point']['coordinates']) not in [list] and ax_lim is not None:
                            plot_dict['ideal_point']['coordinates'] = [ax_lim[plt_order[0]][1],
                                                                       ax_lim[plt_order[1]][1]]
                            ideal_point_loc = plot_dict['ideal_point']['coordinates']

                    color_star = 'No'
                    # Uses the color scheme for the start (i.e., makes it blue if blue represents the best value for
                    # the objective corresponding to color.
                    if color_star == 'Yes':
                        ax_arr[sp_slot_tuple].scatter(ideal_point_loc[0], ideal_point_loc[1], marker=ideal_point_marker,
                                                      c=[cmap(max(norm_array))], s=ideal_point_size, linewidth=0)
                    else:
                        ax_arr[sp_slot_tuple].scatter(ideal_point_loc[0], ideal_point_loc[1], marker=ideal_point_marker,
                                                      c='k', s=ideal_point_size, linewidth=0)

                #cmap(max(norm_array))
                # Axis Labels
                try:
                    # Deal with formatting of x axis labels
                    x_ax_lab = ''
                    x_ax_lab_components = plot_dict['x_label'].split(' \\n ')
                    if len(x_ax_lab_components) > 1:
                        for k in range(len(x_ax_lab_components)):
                            if k < len(x_ax_lab_components) - 1:
                                x_ax_lab += x_ax_lab_components[k] + '\n'
                            else:
                                x_ax_lab += x_ax_lab_components[k]
                    else:
                        x_ax_lab = plot_dict['x_label']
                    # Adjust so powers ("^") show correctly in x-axis label
                    try:
                        power_index_x = x_ax_lab.index('^')
                        follow_val_x = x_ax_lab[power_index_x+1]
                    except ValueError:
                        power_index_x = None
                        follow_val_x = None
                    if power_index_x is not None:
                        x_ax_lab = x_ax_lab.replace('^%s' % follow_val_x, '$\mathregular{^%s}\!$' % follow_val_x)  #

                    # Deal with formatting of y axis labels
                    y_ax_lab = ''
                    y_ax_lab_components = plot_dict['y_label'].split(' \\n ')
                    if len(y_ax_lab_components) > 1:
                        for k in range(len(y_ax_lab_components)):
                            if k < len(y_ax_lab_components) - 1:
                                y_ax_lab += y_ax_lab_components[k] + '\n'
                            else:
                                y_ax_lab += y_ax_lab_components[k]
                    else:
                        y_ax_lab = plot_dict['y_label']
                    # Adjust so powers ("^") show correctly in y-axis label
                    try:
                        power_index_y = y_ax_lab.index('^')
                        follow_val_y = y_ax_lab[power_index_y+1]
                    except ValueError:
                        power_index_y = None
                        follow_val_y = None
                    if power_index_y is not None:
                        y_ax_lab = y_ax_lab.replace('^%s' % follow_val_y, '$\mathregular{^%s}\!$' % follow_val_y)  #

                    if preference_arrow == 'Yes':
                        if invert_axis is not None:
                            if invert_axis[plt_order[0]] == 'No':
                                # Arrows point to right/up
                                ax_arr[sp_slot_tuple].set_xlabel('$\\longrightarrow$' + '\n' + x_ax_lab,
                                                                 fontsize=all_label_font)
                            else:
                                ax_arr[sp_slot_tuple].set_xlabel('$\\longrightarrow$' + '\n' + x_ax_lab,
                                                                 fontsize=all_label_font)
                            if invert_axis[plt_order[1]] == 'No':
                                # Arrows point to right/up
                                ax_arr[sp_slot_tuple].set_ylabel(y_ax_lab + '\n' + '$\\longrightarrow$',
                                                                 fontsize=all_label_font)
                            else:
                                ax_arr[sp_slot_tuple].set_ylabel(y_ax_lab + '\n' + '$\\longrightarrow$',
                                                                 fontsize=all_label_font)
                        else:
                            # Arrows point to right/up
                            ax_arr[sp_slot_tuple].set_xlabel('$\\longrightarrow$' + '\n' + x_ax_lab,
                                     fontsize=all_label_font)
                            ax_arr[sp_slot_tuple].set_ylabel(y_ax_lab + '\n' + '$\\longrightarrow$',
                                                             fontsize=all_label_font)
                    else:
                        ax_arr[sp_slot_tuple].set_xlabel(x_ax_lab, fontsize=all_label_font)
                        ax_arr[sp_slot_tuple].set_ylabel(y_ax_lab, fontsize=all_label_font)
                except KeyError:
                    ax_arr[sp_slot_tuple].set_ylabel(objs_to_plot[plt_order[1]], fontsize=all_label_font)
                    ax_arr[sp_slot_tuple].set_xlabel(objs_to_plot[plt_order[0]], fontsize=all_label_font)
                # Axis Limits
                if ax_lim is not None:
                    ax_arr[sp_slot_tuple].set_xlim(ax_lim[plt_order[0]][0], ax_lim[plt_order[0]][1])
                    ax_arr[sp_slot_tuple].set_ylim(ax_lim[plt_order[1]][0], ax_lim[plt_order[1]][1])

                # Set fontsize for y tick labels
                for tick in ax_arr[sp_slot_tuple].xaxis.get_major_ticks():
                    tick.label.set_fontsize(all_label_font)
                # Set fontsize for x tick labels
                for tick in ax_arr[sp_slot_tuple].yaxis.get_major_ticks():
                    tick.label.set_fontsize(all_label_font)

                if invert_axis is not None:
                    if invert_axis[plt_order[0]] == 'Yes':
                        ax_arr[sp_slot_tuple].invert_xaxis()
                    if invert_axis[plt_order[1]] == 'Yes':
                        ax_arr[sp_slot_tuple].invert_yaxis()

                # Turn top-most and right-most ticks off

                ax_arr[sp_slot_tuple].tick_params(axis='x', which='both', top='off')
                ax_arr[sp_slot_tuple].tick_params(axis='y', which='both', right='off')

                # Turn top and right axes off (no line appears, so axes are not a box)
                turn_off_top_right_axes = 'No'
                if turn_off_top_right_axes == 'Yes':
                    ax_arr[sp_slot_tuple].spines["top"].set_visible(False)
                    ax_arr[sp_slot_tuple].spines["right"].set_visible(False)

                # Add Axis arrows on x-axis as annotated text, indicating direction of increasing preference.
                place_arrow_set_position = 'No'  # Default is to place arrows under axis labels
                if preference_arrow == 'Yes':
                    if place_arrow_set_position == 'Yes':
                        ax_arr[sp_slot_tuple].annotate('', xy=(0.25,-1.15), xycoords = 'axes fraction', xytext=(0.75,-1.15),
                                                       arrowprops=dict(arrowstyle='<-, head_width=0.15', color='k'))

                        # Add Axis arrows on y-axis as annotated text, indicating direction of increasing preference.
                        ax_arr[sp_slot_tuple].annotate('', xy=(-1.25,0.25), xycoords = 'axes fraction', xytext=(-1.25,0.75),
                                                       arrowprops=dict(arrowstyle='<-, head_width=0.15', color='k'))

                try:
                    ax_arr[sp_slot_tuple].set_title(plot_dict['subplot_title'], fontweight='bold')
                except KeyError:
                    pass

                #fig.colorbar(ax_arr[sp_slot_tuple], shrink=0.75)
                #fig.axes[-1].set_ylabel(objs_to_plot[plt_order[2]], fontsize=18)
            #p.show()
            if save_fig is not None:
                save_fig = save_fig + os_fold + objs_to_plot[plt_order[0]] + '-' + objs_to_plot[
                    plt_order[1]] + '-' + objs_to_plot[plt_order[2]] + ".png"
        elif num_objs_to_plot == 4:
            if three_d_plot == 'Yes':
                fig = pyplot.figure()  # create the figure
            else:
                fig, ax_arr = axis_test
            cmap = pyplot.cm.get_cmap("jet_r")
            if three_d_plot == 'Yes':
                ax3D = Axes3D(fig)
                pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[plt_order[2]])
                if plt_order is None:
                    plt_order = [0, 1, 2, 3]
                    pts = ax3D.scatter(objective_values[0], objective_values[1], objective_values[2], c= objective_values[3],
                    cmap=cmap, linewidth=0)
                    fig.colorbar(pts,ax=ax3D,shrink=0.75)
                    fig.axes[-1].set_ylabel(objs_to_plot[3])
                else:
                    # Order of plotting specified.
                    pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[
                        plt_order[2]], c= objective_values[plt_order[3]], cmap=cmap, linewidth=0)
                    fig.colorbar(pts, ax=ax3D, shrink=0.75)
                    try:
                        fig.axes[-1].set_ylabel(plot_dict['c_label'])
                    except KeyError:
                        fig.axes[-1].set_ylabel(objs_to_plot[plt_order[3]])
            else:
                if plot_specs is not None:
                    rcParams.update(plot_specs)
                ax_arr[sp_slot_tuple].scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                                              c=objective_values[plt_order[2]], cmap=cmap, marker='o',
                                              s = objective_values[plt_order[3]]/100, linewidth=0)
                # Axis labels.
                try:
                    x_ax_lab = ''
                    x_ax_lab_components = plot_dict['x_label'].split(' \\n ')
                    y_ax_lab = ''
                    y_ax_lab_components = plot_dict['y_label'].split(' \\n ')
                    if len(x_ax_lab_components) > 1:
                        for k in range(len(x_ax_lab_components)):
                            if k < len(x_ax_lab_components) - 1:
                                x_ax_lab += x_ax_lab_components[k] + '\n'
                            else:
                                x_ax_lab += x_ax_lab_components[k]
                    else:
                        x_ax_lab = plot_dict['x_label']
                    if len(y_ax_lab_components) > 1:
                        for k in range(len(y_ax_lab_components)):
                            if k < len(y_ax_lab_components) - 1:
                                y_ax_lab += y_ax_lab_components[k] + '\n'
                            else:
                                y_ax_lab += y_ax_lab_components[k]
                    else:
                        y_ax_lab = plot_dict['y_label']
                    ax_arr[sp_slot_tuple].set_ylabel(y_ax_lab)
                    ax_arr[sp_slot_tuple].set_xlabel(x_ax_lab)
                except KeyError:
                    ax_arr[sp_slot_tuple].set_ylabel(objs_to_plot[plt_order[1]])
                    ax_arr[sp_slot_tuple].set_xlabel(objs_to_plot[plt_order[0]])
                # Axis limits
                if ax_lim is not None:
                    ax_arr[sp_slot_tuple].set_xlim(ax_lim[plt_order[0]][0],ax_lim[plt_order[0]][1])
                    ax_arr[sp_slot_tuple].set_ylim(ax_lim[plt_order[1]][0],ax_lim[plt_order[1]][1])
                # Plot title
                try:
                    ax_arr[sp_slot_tuple].set_title(plot_dict['subplot_title'], fontweight='bold')
                except KeyError:
                    pass
                # Set fontsize for y tick labels
                for tick in ax_arr[sp_slot_tuple].xaxis.get_major_ticks():
                    tick.label.set_fontsize(all_label_font)
                # Set fontsize for x tick labels
                for tick in ax_arr[sp_slot_tuple].yaxis.get_major_ticks():
                    tick.label.set_fontsize(all_label_font)
        elif num_objs_to_plot == 5:
            ax3D = Axes3D(fig)
            cmap = pyplot.cm.get_cmap("jet_r")
            if plt_order is None:
                plt_order = [0, 1, 2, 3, 4]
                pts = ax3D.scatter(objective_values[0], objective_values[1], objective_values[2], c= objective_values[3],
                cmap=cmap, marker='o', s=objective_values[4]/10, linewidth=0)
                fig.colorbar(pts,ax=ax3D,shrink=0.75)
                fig.axes[-1].set_ylabel(objs_to_plot[3])
            else:
                # Order of plotting specified.
                pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[
                    plt_order[2]], c = objective_values[plt_order[3]], cmap=cmap, marker='o', s=objective_values[
                    plt_order[4]]/10, linewidth=0)
                fig.colorbar(pts, ax=ax3D, shrink=0.75)
                try:
                    fig.axes[-1].set_ylabel(plot_dict['c_label'])
                except KeyError:
                    fig.axes[-1].set_ylabel(objs_to_plot[plt_order[3]])

        if (three_d_plot == 'Yes') or num_objs_to_plot == 5:
            # Label axes, title plot, and save resulting plot.
            try:
                ax3D.set_xlabel(plot_dict['x_label'], fontsize=16, labelpad=10)
                ax3D.set_ylabel(plot_dict['y_label'], fontsize=16, labelpad=10)
                ax3D.set_zlabel(plot_dict['z_label'], fontsize=16, labelpad=10)
            except KeyError:
                ax3D.set_xlabel(objs_to_plot[plt_order[0]], fontsize=16, labelpad=10)
                ax3D.set_ylabel(objs_to_plot[plt_order[1]], fontsize=16, labelpad=10)
                ax3D.set_zlabel(objs_to_plot[plt_order[2]], fontsize=16, labelpad=10)
            if ax_lim is not None:
                # User set axis limits.
                ax3D.set_xlim3d(ax_lim[plt_order[0]][0], ax_lim[plt_order[0]][1])
                ax3D.set_ylim3d(ax_lim[plt_order[1]][0], ax_lim[plt_order[1]][1])
                ax3D.set_zlim3d(ax_lim[plt_order[2]][0], ax_lim[plt_order[2]][1])

    #p.show()
    if save_fig is not None:
        if type(save_fig) not in [str, unicode]:
            # Appropriate figure name not given. Give figure a title and save in pysedsim current working directory.
            fig.savefig("scatter_plot.png", dpi=1200, bbox_inches='tight')
        else:
            fig.savefig(save_fig, dpi=1200, bbox_inches='tight')

    # Create text files of operating policy parameters for particular points on pareto front, if user desires.
    if gen_op_policy_file == 'Yes':
        Grab_Pareto_Min_Max(ref_set_array, objective_values, num_objs, num_dec_vars)
    if ref_set_file_name is not None and tradeoff_plotting_subplot is None:
        # Return reference set information if user is internally processing reference set here.
        return ref_set_array, objective_values, dec_var_values
    if ref_set_file_name is not None and tradeoff_plotting_subplot is not None:
        # Return reference set information if user is internally processing reference set here.
        return ref_set_array, objective_values, dec_var_values, fig, ax_arr


def Scatter_Plot_Objectives(ref_set_pref_dict, objs_to_plot, parse_objs=None, plot_dict=None,
                            save_fig = None, gen_op_policy_file='No', num_objs_to_plot=None,
                            tradeoff_plotting_subplot=None, axis_test = None, provided_3D_fig_axis=None):

    '''

    Purpose: Produces a scatter plot of objective function values up to 4D (4th dimension is color).

    Args (if not defined here, it is defined in Reference_Set() method):
        3. objs_to_plot: list of objective names as they appear in the pysedsim input file and in the reference set
        file, e.g. ['Objective 1', 'Objective 2']. Only list the names of objectives that you want to plot.
        4. gen_op_policy_file: Optional, specify 'Yes' or 'No'. default = 'No'. Calls Grab_Pareto_Min_Max()
        function to generate text file of operating policy for minimum and maximum points in the reference set.
        5. parse_objs: List of objective numbers (column numbers in 'PySedSim_reference.txt' to plot from
        reference set
        6. ref_set_pref: A dictionary containing the following keys:
        a. num_objs: Number of objective function values
        b. num_dec_vars: Number of decision variable values
        c. ref_set_file_name: name of reference set file (e.g., 'ref_set_file.txt')
        d. unit_conv: List of factors by which to multiply values in each array of objective values.
        e. perc_conv: List of 'Yes' or 'No' for each objective, for converting objective values to percent from
        fraction.
        f. invert: List of 'Yes' or 'No' for each objective, for reversing percentage (doing 100-value).
        9. plot_dict: Dictionary of plotting preferences. Keys include:
        a. 'Axis Range': List of lists (sub-lists), containing axis limits for 2 or 3 primary axes: x, y,
        z. e.g., [[2,3],[3,5],[4,5]]
        b. 'x_labal': string to label x-axis.
        c. 'y_label': string to label y-axis.
        d. 'z_label': string to label z-axis.
        e. 'c_label': string to label "colored" axis (if there is a 4th dimension).
        f. 'Plot Specifications': dictionary, follows format of rcParams module from matplotlib. An example is:

            'Plot Specifications': {'axes.labelsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18,
                                 'text.usetex': False, 'figure.figsize': [9, 6], 'figure.figsize': [6.5, 5],
                                 'legend.fontsize': 22}
        g. '3d_plot': 'Yes' or 'No', indicates whether user wants to plot 3 objectives in 3D space, or in 2D space
        using color as the third objective.

        10. plt_order: Optional, this is a list that specifies in which order the objectives are plotted. Used if you
        want to be able to determine which of the objectives become a color/point size/orientation in the >3D case.
        Example: plt_order = [0,3,2,1] to have the third variable in the objective values array plotted in the 3D
        space, while the third item is plotted last.

        11. num_objs_to_plot: integer number of objectives you actually wish to plot. Will be computed as length of
        objs_to_plot list if nothing is specified.

        tradeoff_plotting_subplot: indicates whether call is being made from tradeoff_plotting module for a subplot,
        in which case the figure object should be returned so it can be included in a subplot.
    Returns:
        Nothing, but several figures may be saved
    '''

    # Load basic information about the simulated scenarios, including input file directory and simulation names.
    os_fold = Op_Sys_Folder_Operator()

    # Unpack ref_set_pref_dict
    try:
        ref_set_file_name = ref_set_pref_dict['ref_set_file_name']
        num_objs = ref_set_pref_dict['num_objs']
    except KeyError:
        ref_set_file_name = None
        num_objs = None
    try:
        num_dec_vars = ref_set_pref_dict['num_dec_vars']
    except KeyError:
        pass
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

    if ref_set_file_name is not None:
        # User wishes to import reference set here.
        [ref_set_array, objective_values, dec_var_values] = Initial_Processing(num_objs, num_dec_vars, ref_set_file_name,
                                                                               parse_objs=parse_objs, perc_conv=perc_conv,
                                                                               invert=invert, unit_conv=unit_conv)
    else:
        # User has provided reference set information.
        ref_set_array = ref_set_pref_dict['ref_set_array']
        objective_values = ref_set_pref_dict['objective_values']
        dec_var_values = ref_set_pref_dict['dec_var_values']

    # Unpack plotting dictionary
    try:
        plt_order = plot_dict['plot_order']
    except KeyError:
        plt_order = None
    try:
        ax_lim = plot_dict['Axis Range']
        # Determine axis limits
        try:
            if plot_dict['compute_upper_ax_lim'] == 'Yes':
                for obj in range(len(objective_values)):
                    ax_lim[obj][1] = np.max(objective_values[obj]) + 0.1*np.max(objective_values[obj])
        except KeyError:
            pass
    except KeyError:
        ax_lim = None
    try:
        pols_to_highlight = plot_dict['pols_to_highlight']
    except KeyError:
        pols_to_highlight = None
    try:
        policy_labels = plot_dict['policy_labels']
    except KeyError:
        policy_labels = None
    try:
        three_d_plot = plot_dict['3d_plot']
    except KeyError:
        three_d_plot = 'No'
    try:
        plot_specs = plot_dict['Plot Specifications']
    except KeyError:
        plot_specs = None
    try:
        include_colorbar = plot_dict['include_colorbar']
    except KeyError:
        include_colorbar = 'Yes'
    try:
        snapshot_distance = plot_dict['snapshot_distance']
    except KeyError:
        snapshot_distance = None
    try:
        viewpoint = plot_dict['viewpoint']
    except KeyError:
        viewpoint = None
    try:
        invert_x_axis = plot_dict['invert_x_axis']
    except KeyError:
        invert_x_axis = None
    try:
        invert_y_axis = plot_dict['invert_y_axis']
    except KeyError:
        invert_y_axis = None
    try:
        invert_z_axis = plot_dict['invert_z_axis']
    except KeyError:
        invert_z_axis = None
    try:
        plot_ideal_point = plot_dict['plot_ideal_point']
    except KeyError:
        plot_ideal_point = None
    try:
        num_objs_to_plot = ref_set_pref_dict['num_objs_to_plot']
    except KeyError:
        num_objs_to_plot = len(objs_to_plot)

    try:
        pol_name_locs = plot_dict['pol_name_locs']
    except KeyError:
        pass

    if num_objs_to_plot == 1:
        print("Only 1 objective; no plot will be produced. Result is a single value.")
    elif num_objs_to_plot == 2:
        if plot_specs is not None:
            rcParams.update(plot_specs)
        fig = pyplot.figure()
        if plt_order is None:
            # Order of plotting not specified.
            plt_order = [0,1]
        pyplot.scatter(objective_values[plt_order[0]],objective_values[plt_order[1]])
        try:
            pyplot.ylabel(plot_dict['y_label'])
            pyplot.xlabel(plot_dict['x_label'])
        except KeyError:
            pyplot.ylabel(objs_to_plot[plt_order[1]])
            pyplot.xlabel(objs_to_plot[plt_order[0]])
        if ax_lim is not None:
            pyplot.xlim(ax_lim[plt_order[0]][0],ax_lim[plt_order[0]][1])
            pyplot.ylim(ax_lim[plt_order[1]][0],ax_lim[plt_order[1]][1])
        if save_fig is not None:
            save_fig = save_fig + os_fold + objs_to_plot[plt_order[0]] + '-' + objs_to_plot[
                plt_order[1]] + ".png"
    elif num_objs_to_plot >= 3:
        fig = pyplot.figure()  # create the figure
        if num_objs_to_plot == 3:
            # Create a 3D plot, with no color used to map a fourth objective.
            if plt_order is None:
                # Order of plotting not specified.
                plt_order = [0,1,2]
            if plot_specs is not None:
                rcParams.update(plot_specs)
            if three_d_plot == 'Yes':
                ax3D = Axes3D(fig)
                pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[plt_order[2]])
            else:
                # Create a 2D plot, with color used to map the third objective.
                cmap = pyplot.cm.get_cmap("jet_r")
                new_plot = pyplot.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]],
                               c=objective_values[plt_order[2]], cmap=cmap, marker='o', s=30)
                try:
                    pyplot.ylabel(plot_dict['y_label'])
                    pyplot.xlabel(plot_dict['x_label'])
                except KeyError:
                    pyplot.ylabel(objs_to_plot[plt_order[1]])
                    pyplot.xlabel(objs_to_plot[plt_order[0]])
                if ax_lim is not None:
                    pyplot.xlim(ax_lim[plt_order[0]][0],ax_lim[plt_order[0]][1])
                    pyplot.ylim(ax_lim[plt_order[1]][0],ax_lim[plt_order[1]][1])
                fig.colorbar(new_plot, shrink=0.75)
                fig.axes[-1].set_ylabel(objs_to_plot[plt_order[2]], fontsize=18)
            pyplot.show()
            if save_fig is not None:
                save_fig = save_fig + os_fold + objs_to_plot[plt_order[0]] + '-' + objs_to_plot[
                    plt_order[1]] + '-' + objs_to_plot[plt_order[2]] + ".png"
        elif num_objs_to_plot == 4:
            if plot_specs is not None:
                rcParams.update(plot_specs)
            if provided_3D_fig_axis is None:
                ax3D = Axes3D(fig)
            else:
                ax3D = provided_3D_fig_axis
            cmap = pyplot.cm.get_cmap("jet_r")
            if plt_order is None:
                plt_order = [0,1,2,3]
                pts = ax3D.scatter(objective_values[0], objective_values[1], objective_values[2], c= objective_values[3],
                cmap=cmap, linewidth=0)
                if include_colorbar == 'Yes':
                    fig.colorbar(pts,ax=ax3D,shrink=0.75)
                    fig.axes[-1].set_ylabel(objs_to_plot[3])
            else:
                # Order of plotting specified.
                pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[
                    plt_order[2]], c= objective_values[plt_order[3]], cmap=cmap, linewidth=0)
                if include_colorbar == 'Yes':
                    fig.colorbar(pts, ax=ax3D, shrink=0.75)
                    try:
                        fig.axes[-1].set_ylabel(plot_dict['c_label'])
                    except KeyError:
                        fig.axes[-1].set_ylabel(objs_to_plot[plt_order[3]])
        elif num_objs_to_plot == 5:
            ax3D = Axes3D(fig)
            cmap = pyplot.cm.get_cmap("jet_r")
            if plt_order is None:
                plt_order = [0, 1, 2, 3, 4]
                pts = ax3D.scatter(objective_values[0], objective_values[1], objective_values[2], c= objective_values[3],
                cmap=cmap, marker='o', s=objective_values[4]/10, linewidth=0)
                fig.colorbar(pts,ax=ax3D,shrink=0.75)
                fig.axes[-1].set_ylabel(objs_to_plot[3])
            else:
                # Order of plotting specified.
                pts = ax3D.scatter(objective_values[plt_order[0]], objective_values[plt_order[1]], objective_values[
                    plt_order[2]], c = objective_values[plt_order[3]], cmap=cmap, marker='o', s=objective_values[
                    plt_order[4]]/10, linewidth=0)
                fig.colorbar(pts, ax=ax3D, shrink=0.75)
                try:
                    fig.axes[-1].set_ylabel(plot_dict['c_label'])
                except KeyError:
                    fig.axes[-1].set_ylabel(objs_to_plot[plt_order[3]])

        if (three_d_plot == 'Yes') or num_objs_to_plot == 4 or num_objs_to_plot == 5:
            # Label axes, title plot, and save resulting plot.
            try:
                ax3D.set_xlabel(plot_dict['x_label'], labelpad=10)
                ax3D.set_ylabel(plot_dict['y_label'], labelpad=50)
                tmp_planes = ax3D.zaxis._PLANES
                ax3D.zaxis._PLANES = (
                tmp_planes[2], tmp_planes[3], tmp_planes[0], tmp_planes[1], tmp_planes[4], tmp_planes[5])
                ax3D.set_zlabel(plot_dict['z_label'], labelpad=10, rotation=90)
            except KeyError:
                ax3D.set_xlabel(objs_to_plot[plt_order[0]], labelpad=10)
                ax3D.set_ylabel(objs_to_plot[plt_order[1]], labelpad=10)
                ax3D.set_zlabel(objs_to_plot[plt_order[2]], labelpad=10)
            if ax_lim is not None:
                # User set axis limits.
                ax3D.set_xlim3d(ax_lim[plt_order[0]][0], ax_lim[plt_order[0]][1])
                ax3D.set_ylim3d(ax_lim[plt_order[1]][0], ax_lim[plt_order[1]][1])
                ax3D.set_zlim3d(ax_lim[plt_order[2]][0], ax_lim[plt_order[2]][1])
            else:
                # Set axis limits based on best/worst performing values
                best_obj_list = []
                for obj in range(num_objs):
                    if plot_dict['objectives_type'][plt_order[obj]] == 'max':
                        if obj == 0:
                            ax3D.set_xlim3d(min(objective_values[plt_order[obj]]), max(objective_values[plt_order[
                                obj]]))
                        if obj == 1:
                            ax3D.set_ylim3d(min(objective_values[plt_order[obj]]), max(objective_values[plt_order[
                                obj]]))
                        if obj == 2:
                            ax3D.set_zlim3d(min(objective_values[plt_order[obj]]), max(objective_values[plt_order[
                                obj]]))
                        best_obj_list.append(max(objective_values[plt_order[obj]]))
                    else:
                        if obj == 0:
                            ax3D.set_xlim3d(max(objective_values[plt_order[obj]]), min(objective_values[plt_order[
                                obj]]))
                        if obj == 1:
                            ax3D.set_ylim3d(max(objective_values[plt_order[obj]]), min(objective_values[plt_order[
                                obj]]))
                        if obj == 2:
                            ax3D.set_zlim3d(max(objective_values[plt_order[obj]]), min(objective_values[plt_order[
                                obj]]))
                        best_obj_list.append(min(objective_values[plt_order[obj]]))
            if pols_to_highlight is not None:
                low = min(objective_values[plt_order[3]])
                high = max(objective_values[plt_order[3]])
                norm_array = np.array([(c-low)/(high-low) for c in objective_values[plt_order[3]]])
                for x in range(len(pols_to_highlight)):
                    low_value = int(min(objective_values[plt_order[3]]))
                    high_value = int(max(objective_values[plt_order[3]]))
                    diff = int(high_value-low_value)
                    color_value = pyplot.cm.jet_r((np.clip(objective_values[plt_order[3]][pols_to_highlight[x]],
                                                       low_value,high_value)-low_value)/diff)
                    if policy_labels is not None:
                        # Create thick marker linewidth for policies being highlighted, and use larger size,
                        # and don't have a transparent point.
                        ax3D.scatter(objective_values[plt_order[0]][pols_to_highlight[x]],
                                     objective_values[plt_order[1]][pols_to_highlight[x]],
                                     objective_values[plt_order[2]][pols_to_highlight[x]],
                                     color=cmap(norm_array[pols_to_highlight[x]]), marker='o',
                                     s=50, linewidth=1.0, edgecolor='k')

                    #xy_loc = (objective_values[plt_order[0]][pols_to_highlight[x]], objective_values[plt_order[
                    #    1]][pols_to_highlight[x]], objective_values[plt_order[1]][pols_to_highlight[x]])
                    #xy_text_loc = pol_name_locs[x]
                    #if policy_labels is not None:
                    #    ax3D.annotate(policy_labels[x], fontsize=6, xy=xy_loc, xycoords='data',
                    #                xytext=xy_text_loc, bbox=dict(boxstyle="round", fc='w'),
                    #                arrowprops=dict(arrowstyle="->",connectionstyle="arc,angleA=0,angleB=90,"
                    #                                                                "rad=10"))

            if invert_x_axis not in [None, 'No']:
                ax3D.invert_xaxis()
            if invert_y_axis not in [None, 'No']:
                ax3D.invert_yaxis()
            if invert_z_axis not in [None, 'No']:
                ax3D.invert_zaxis()
    #fig.tight_layout()
    if snapshot_distance is not None:
        ax3D.dist=snapshot_distance
    if viewpoint is not None:
        ax3D.view_init(elev=viewpoint[0], azim=viewpoint[1])

    if plot_ideal_point == 'Yes':
        ax3D.scatter(best_obj_list[0], best_obj_list[1], best_obj_list[2], marker="*", s=50, color='darkblue',
                     edgecolor='darkblue')
    ax3D.grid(False)
    pyplot.show()

    fig.set_size_inches(4,4)
    #fig.tight_layout()
    if save_fig is not None:
        if type(save_fig) not in [str, unicode]:
            # Appropriate figure name not given. Give figure a title and save in pysedsim current working directory.
            fig.savefig("scatter_plot.png", dpi=1200, bbox_inches='tight')
        else:
            fig.savefig(save_fig, dpi=1200, bbox_inches='tight')

    # Package figure and axis for return
    if ax3D is not None:
        fig_axis = [fig, ax3D]
    else:
        fig_axis = []

    # Create text files of operating policy parameters for particular points on pareto front, if user desires.
    if gen_op_policy_file == 'Yes':
        Grab_Pareto_Min_Max(ref_set_array, objective_values, num_objs, num_dec_vars)
    if ref_set_file_name is not None and tradeoff_plotting_subplot is None:
        # Return reference set information if user is internally processing reference set here.
        return ref_set_array, objective_values, dec_var_values, fig_axis
    if ref_set_file_name is not None and tradeoff_plotting_subplot is not None:
        # Return reference set information if user is internally processing reference set here.
        return ref_set_array, objective_values, dec_var_values, fig


def tradeoff_plotting(input_file_name = 'PySedSim_Input_Specifications.csv', ref_set_file_name=None):
    '''

    Purpose: Generates plots of reference set objectives, according to preferences in 'Tradeoff Plotting' sheet of
    input file.

    :param input_file_name: Top-level pysedsim input file (.csv) that contains name of scenario(s) and input file(s)
    locations. Can be a file path, but must include actual file name.

    :param ref_set_file_name: reference set file (can be file path). Default is to assume it is located in
    directory Output_Storage/Scenario_Name/sets/pysedsim_ref_set.ref

    :return:
    '''

    [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
     os_fold] = import_specs(file_name=input_file_name)

    # Loop through as many optimization scenarios as user has specified.
    tradeoff_plot_pref_dict = {}
    axis_offset_list = [13, 21, 29, 38, 46]  # Position of data in 'Tradeoff Plotting' sheet for each axis.
    plot_list = []  # names of all plots (subplots do not count; subplots make up a plot). Goes across scenarios.
    opt_dicts_dict = {}
    plot_list_dict = {}
    num_plots_dict = {}
    axis_list = ['x_axis', 'y_axis', 'color_axis', 'z_axis', 'size_axis']
    axis_loc_dict = {'x_axis': 0, 'y_axis': 1, 'color_axis': 2, 'z_axis': 4, 'size_axis': 3}
    sp_list = {}
    input_ref_set_file_path = ref_set_file_name

    for j in range(num_scenarios):
        simulation_title = simulation_titles_list[j]
        output_location = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets'
        plot_list_dict[simulation_title] = []
        # 1. Create reference set file if not already created. Run this code before running any of the figure generation code
        #  below.
        [Borg_dict, DPS_dict] = Reference_Set(input_file_name=input_file_name, create_ref_set='No')
        if Borg_dict['optimization approach'] == 'DPS':
            Borg_dict['n_vars'] = DPS_dict['total_vars']  # Store num decision variables in Borg_dict, copied from
            # DPS_dict.

        sp_list[simulation_title] = []

        # Import plot preferences
        opt_dicts_dict[simulation_title] = {'Borg_dict': Borg_dict, 'DPS_dict': DPS_dict}
        Input_Data_File = Load_Input_File(simulation_title, main_input_files_dir, imported_specs)

        # Read in preferences related to tradeoff plotting from the "Tradeoff Plotting" worksheet in input file.
        try:
            tradeoff_plot_pref_dict[simulation_title] = {}
            num_plots_dict[simulation_title] = Input_Data_File['Tradeoff Plotting']['B1'].value
            num_plots = num_plots_dict[simulation_title]
            for p in range(num_plots):
                plot_name = Input_Data_File['Tradeoff Plotting'].cell(row = 2, column = 2 + p).value
                pols_to_highlight = Input_Data_File['Tradeoff Plotting'].cell(row = 5, column = 2 + p).value
                pol_names = Input_Data_File['Tradeoff Plotting'].cell(row = 6, column = 2 + p).value
                if pols_to_highlight is not None:
                    # Handle policy numbers
                    if type(pols_to_highlight) in [str, unicode]:
                        # List of policies has been specified.
                        pols_to_highlight = pols_to_highlight.split(', ')  # Create list of policies as strings
                        # Convert list elements to integers
                        pols_to_highlight = [int(pols_to_highlight[i]) for i in range(len(pols_to_highlight))]
                    else:
                        # Only one policy specified. Put that value into a list.
                        pols_to_highlight = [pols_to_highlight]
                    # Handle policy names
                    if type(pol_names) in [str, unicode]:
                        # List of policies has been specified.
                        pol_names = pol_names.split(', ')  # Create list of policies as strings
                        # Convert list elements to integers
                        pol_names = [str(pol_names[i]) for i in range(len(pol_names))]
                    else:
                        # Only one policy specified. Put that value into a list.
                        pol_names = [pol_names]

                if plot_name not in plot_list:
                    plot_list.append(plot_name)
                if plot_name not in plot_list_dict[simulation_title]:
                    plot_list_dict[simulation_title].append(plot_name)
                    tradeoff_plot_pref_dict[simulation_title][plot_name] = {}
                    if Input_Data_File['Tradeoff Plotting'].cell(row = 7, column = 2 + p).value is not None:
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['main_plot_figsize_inches'] = \
                        [float(c) for c in Input_Data_File['Tradeoff Plotting'].cell(row=7, column=2 + p).value.split(', ')]
                    if Input_Data_File['Tradeoff Plotting'].cell(row = 8, column = 2 + p).value is not None:
                        # Choices: 'Yes' or 'No'
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['include_colorbar'] = Input_Data_File[
                            'Tradeoff Plotting'].cell(row=8, column=2 + p).value
                    else:
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['include_colorbar'] = 'No'
                    if Input_Data_File['Tradeoff Plotting'].cell(row = 9, column = 2 + p).value is not None:
                        # Choices: 'Yes' or 'No'
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['ideal_point'] = Input_Data_File[
                            'Tradeoff Plotting'].cell(row=9, column=2 + p).value
                    else:
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['ideal_point'] = 'No'
                    if Input_Data_File['Tradeoff Plotting'].cell(row = 10, column = 2 + p).value is not None:
                        # Choices: 'Yes' or 'No'
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['preference_arrow'] = Input_Data_File[
                            'Tradeoff Plotting'].cell(row=10, column=2 + p).value
                    else:
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['preference_arrow'] = 'No'
                    if Input_Data_File['Tradeoff Plotting'].cell(row = 11, column = 2 + p).value is not None:
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'] = Input_Data_File[
                            'Tradeoff Plotting'].cell(row = 11, column = 2 + p).value
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'] = [int(h) for h in
                                                                                                    tradeoff_plot_pref_dict[
                                                                                                        simulation_title][
                                                                                                        plot_name][
                                                                                                        'figure_rows_cols'].split(
                                                                                                        ', ')]
                    else:
                        tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'] = None
                if Input_Data_File['Tradeoff Plotting'].cell(row = 3, column = 2 + p).value is not None:
                    subplot_number = str(Input_Data_File['Tradeoff Plotting'].cell(row = 3, column = 2 + p).value)
                    sp_list[simulation_title].append(subplot_number)
                else:
                    subplot_number = None

                tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number] = {}
                if subplot_number is not None:
                    # For this subplot, load preferences for each axis.
                    for axis in axis_list:
                        ax_ind = axis_list.index(axis)
                        row_ind = axis_offset_list[ax_ind]
                        if Input_Data_File['Tradeoff Plotting'].cell(row = row_ind, column = 2 + p).value is not None:
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis] = {}
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['objective_name'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['unit_conv'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+1, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['perc_conv'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+2, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['invert'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+3, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['axis_range'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+4, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['title'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+5, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['invert_axis']\
                                = Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+6, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['cmap_name'] = \
                                Input_Data_File['Tradeoff Plotting'].cell(row = row_ind+7, column = 2 + p).value
                            # Subplot title: Same for every axis
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis][
                                'subplot_title'] = Input_Data_File['Tradeoff Plotting'].cell(row=4, column=2 + p).value
        except KeyError:
            print('Proper Tradeoff Plotting worksheet does not exist to generate plots')

    ref_set_pref_dict = {}
    for plot_name in plot_list:
        subplot_loc_list = []
        if num_scenarios > 1:
            # Leave subplot axes as an array
            [main_fig, ax_arr] = pyplot.subplots(num_plots, num_scenarios)
        else:
            # Unpack axes immediately for a horizontal plot.
            if tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'] is None:
                # No preferences about plot orientation have been specified. Default is horizontal row of subplots.
                [main_fig, ax_arr] = pyplot.subplots(1, num_plots)
                for a in range(num_plots):
                    subplot_loc_list.append(tuple([0,a]))
            else:
                # Create list of tuples that indicate location of subplot
                for a in range(tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'][0]):
                    for b in range(tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'][1]):
                        subplot_loc_list.append(tuple([a,b]))

                [main_fig, ax_arr] = pyplot.subplots(
                    tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'][0],
                    tradeoff_plot_pref_dict[simulation_title][plot_name]['figure_rows_cols'][1])

            #main_fig, ()
        # Create all of plots user specifies
        scenario_counter = 0
        for j in range(num_scenarios):
            simulation_title = simulation_titles_list[j]
            sub_plt_counter = 0
            if input_ref_set_file_path is None:
                # Use default file path
                ref_set_file_name = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets' + os_fold + \
                                  'pysedsim_ref_set.ref'
            objs_to_plot = opt_dicts_dict[simulation_title]['Borg_dict']['opt_dict']['Objective Names Ordered List']
            ref_set_pref_dict['ref_set_file_name'] = ref_set_file_name
            ref_set_pref_dict['num_objs'] = opt_dicts_dict[simulation_title]['Borg_dict']['n_objs']
            ref_set_pref_dict['num_dec_vars'] = opt_dicts_dict[simulation_title]['DPS_dict']['total_vars']
            for subplot_number in sp_list[simulation_title]:
                ref_set_pref_dict['unit_conv'] = [1 for i in range(len(objs_to_plot))]
                ref_set_pref_dict['perc_conv'] = ['No' for i in range(len(objs_to_plot))]
                ref_set_pref_dict['invert'] = ['No' for i in range(len(objs_to_plot))]
                ref_set_pref_dict['invert_axis'] = ['No' for i in range(len(objs_to_plot))]
                plot_dict = {}
                spn_str = subplot_number
                ref_set_pref_dict['num_objs_to_plot'] = len(tradeoff_plot_pref_dict[simulation_title][plot_name][
                    spn_str].keys())
                plot_dict['Axis Range'] = [[] for i in range(ref_set_pref_dict['num_objs'])]
                plot_dict['Invert Axis'] = ['No' for i in range(ref_set_pref_dict['num_objs'])]
                plot_dict['plot_order'] = [i for i in range(ref_set_pref_dict['num_objs_to_plot'])]
                plot_dict['objective_type'] = ['Max' for i in range(ref_set_pref_dict['num_objs_to_plot'])]
                if pols_to_highlight is not None:
                    plot_dict['pols_to_highlight'] = pols_to_highlight
                if pol_names is not None:
                    plot_dict['pol_names'] = pol_names
                plot_dict['Plot Specifications'] = {'axes.labelsize': 8, 'xtick.labelsize':
                    8, 'ytick.labelsize': 8, 'text.usetex': False}#, 'figure.figsize': [6.5, 5]}  # keys removed:
                # 'font.size': 20, 'legend.fontsize': 22,
                if tradeoff_plot_pref_dict[simulation_title][plot_name]['ideal_point'] == 'Yes':
                    plot_dict['ideal_point'] = {}
                    plot_dict['ideal_point']['coordinates'] = 'best'
                    plot_dict['ideal_point']['marker_type'] = "*"
                    plot_dict['ideal_point']['color'] = 'bl'
                    plot_dict['ideal_point']['size'] = 50
                if tradeoff_plot_pref_dict[simulation_title][plot_name]['preference_arrow'] == 'Yes':
                    plot_dict['preference_arrow'] = 'Yes'
                try:
                    if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['x_axis']['subplot_title'] is \
                            not None:
                        plot_dict['subplot_title'] = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['x_axis']['subplot_title']
                except KeyError:
                    pass
                for obj in objs_to_plot:
                    for axis in axis_list:
                        try:
                            if obj == tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['objective_name']:
                                try:
                                    ref_set_pref_dict['unit_conv'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['unit_conv']
                                except KeyError:
                                    pass
                                try:
                                    ref_set_pref_dict['perc_conv'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['perc_conv']
                                except KeyError:
                                    pass
                                try:
                                    ref_set_pref_dict['invert'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['invert']
                                except KeyError:
                                    pass

                                try:
                                    ref_set_pref_dict['invert_axis'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['invert_axis']
                                except KeyError:
                                    pass

                                # Order matters
                                if axis == 'x_axis':
                                    plot_dict['x_label'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                        spn_str][axis]['title']
                                elif axis == 'y_axis':
                                    plot_dict['y_label'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                        spn_str][axis]['title']
                                elif axis == 'color_axis':
                                    plot_dict['c_label'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                        spn_str][axis]['title']
                                    if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis][
                                        'cmap_name'] is not None:
                                        plot_dict['cmap_name'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                            spn_str][axis]['cmap_name']
                                elif axis == 'z_axis':
                                    plot_dict['z_label'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                        spn_str][axis]['title']
                                elif axis == 'size_axis':
                                    plot_dict['s_label'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                        spn_str][axis]['title']

                                plot_dict['plot_order'][axis_loc_dict[axis]] = objs_to_plot.index(obj)
                                plot_dict['objective_type'][axis_loc_dict[axis]] = \
                                Borg_dict['opt_dict']['Objective Preferences'][obj]['Type']

                                if type(tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['axis_range']) in [str, unicode]:
                                    # User has specified a list of numbers, which will be imported as one long string with commas
                                    ax_rng = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis][
                                        'axis_range'].split(', ')  # Create list of policies as strings
                                    # Convert to list of integers
                                    plot_dict['Axis Range'][objs_to_plot.index(obj)] = [float(ax_rng[z]) for z in
                                                                                        range(len(ax_rng))]
                                plot_dict['Invert Axis'][objs_to_plot.index(obj)] = \
                                tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['invert_axis']
                        except KeyError:
                            pass
                # Create the plot
                if num_scenarios > 1:
                    # Can access ax_arr like a true numpy array of axes.
                    sp_slot_tuple = (sub_plt_counter, scenario_counter)
                else:
                    sp_slot_tuple = sub_plt_counter
                    sp_slot_tuple = subplot_loc_list[sub_plt_counter]
                [ref_set_array, objective_values, dec_var_values, main_fig, ax_arr] = \
                    Subplot_Scatter_Plot_Objectives(ref_set_pref_dict,
                                                    objs_to_plot,
                                                    plot_dict=plot_dict,
                                                    gen_op_policy_file='No',
                                                    tradeoff_plotting_subplot='Yes',
                                                    axis_test = [main_fig, ax_arr],
                                                    sp_slot_tuple=sp_slot_tuple)
                # To keep track of where in main figure to place each subplot.
                sub_plt_counter += 1
            # To keep track of where in main figure to place each subplot.
            scenario_counter += 1

                #save_fig=output_location
                #ax_arr[0, 0] = fig
                #main_fig.show()

        # Put finishing touches on plot (formatting of whole plot, not individual subplots)
        if ref_set_pref_dict['num_objs_to_plot'] >= 3:
            cmap = pyplot.cm.get_cmap("jet_r")
            if tradeoff_plot_pref_dict[simulation_title][plot_name]['include_colorbar'] == 'Yes':
                sm = matplotlib.cm.ScalarMappable(cmap=cmap)
                if ref_set_pref_dict['num_objs_to_plot'] == 3:
                    if plot_dict['Axis Range'][2] != []:
                        # User specified a range for the color axis
                        sm.set_array(plot_dict['Axis Range'][2])
                    else:
                        sm.set_array()
                elif ref_set_pref_dict['num_objs_to_plot'] == 4:
                    if plot_dict['Axis Range'][3] != []:
                        # User specified a range for the color axis
                        sm.set_array(plot_dict['Axis Range'][3])
                    else:
                        sm.set_array()

            main_fig.tight_layout()  # main_fig.tight_layout()
            if tradeoff_plot_pref_dict[simulation_title][plot_name]['include_colorbar'] == 'Yes':
                main_fig.subplots_adjust(right=0.8)
                cbar_ax = main_fig.add_axes([0.85, 0.15, 0.05, 0.7])  # OLD: cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.5])
                main_fig.colorbar(sm, cax=cbar_ax)
                color_label_components = plot_dict['c_label'].split(' ')
                # Give colorbar a label/title
                col_lab = ''
                for k in range(len(color_label_components)):
                    col_lab += color_label_components[k] + '\n'
                main_fig.axes[-1].set_xlabel(col_lab)
        #main_fig.suptitle(dps_plt_dict['main_plot_title'], fontsize=14)
        #main_fig.suptitle('')

        #if tradeoff_plot_pref_dict[simulation_title][plot_name]['include_colorbar'] == 'Yes':
            # Only adjust figure size if there is a color bar.
            #main_fig.subplots_adjust(top=0.90)
        #main_fig.tight_layout()

        if plot_dict['Plot Specifications'] is not None:
            rcParams.update(plot_dict['Plot Specifications'])
        pyplot.show()
        try:
            main_fig.set_size_inches(tradeoff_plot_pref_dict[simulation_title][plot_name]['main_plot_figsize_inches'])
        except KeyError:
            pass
        main_fig.tight_layout()
        main_fig.savefig(plot_name + '.pdf', bbox_inches='tight', dpi=1200)  # '.png'

def tradeoff_plotting_parallel_axis(input_file_name = 'PySedSim_Input_Specifications.csv', ref_set_file_name=None):
    '''

    Purpose: Generates parallel axis plots of reference set objectives, according to preferences in 'Parallel Axis
    Plotting' sheet of input file.

    :param input_file_name: Top-level pysedsim input file (.csv) that contains name of scenario(s) and input file(s)
    locations. Can be a file path, but must include actual file name.

    :param ref_set_file_name: reference set file path. If you specify nothing, the default is to assume the file is
    located as follows: Output_Storage\Scenario_Name\sets\pysedsim_ref_set.ref, where output storage is named
    appropriately in the input_file_name file.

    :return:
    '''

    [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
     os_fold] = import_specs(file_name=input_file_name)

    # Loop through as many optimization scenarios as user has specified.
    tradeoff_plot_pref_dict = {}
    plot_list = []  # names of all plots (subplots do not count; subplots make up a plot). Goes across scenarios.
    opt_dicts_dict = {}
    plot_list_dict = {}
    num_plots_dict = {}
    num_axes_dict = {}
    sp_list = {}

    start_row_par_axis_data = 11

    for j in range(num_scenarios):
        simulation_title = simulation_titles_list[j]
        output_location = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets'
        plot_list_dict[simulation_title] = []
        # 1. Create reference set file if not already created. Run this code before running any of the figure generation code
        #  below.
        [Borg_dict, DPS_dict] = Reference_Set(input_file_name=input_file_name, create_ref_set='No')
        if Borg_dict['optimization approach'] == 'DPS':
            Borg_dict['n_vars'] = DPS_dict['total_vars']  # Store num decision variables in Borg_dict, copied from
            # DPS_dict.

        sp_list[simulation_title] = []

        # Import plot preferences
        opt_dicts_dict[simulation_title] = {'Borg_dict': Borg_dict, 'DPS_dict': DPS_dict}
        Input_Data_File = Load_Input_File(simulation_title, main_input_files_dir, imported_specs)

        # Read in preferences related to tradeoff plotting from the "Parallel Axis Plotting" worksheet in input file.
        try:
            tradeoff_plot_pref_dict[simulation_title] = {}
            num_plots_dict[simulation_title] = Input_Data_File['Parallel Axis Plotting']['B1'].value
            num_plots = num_plots_dict[simulation_title]
            num_axes_dict[simulation_title] = Input_Data_File['Parallel Axis Plotting']['B3'].value
            num_axes = num_axes_dict[simulation_title]
            # Position of data in 'Parallel Axis Plotting' sheet for each axis
            axis_offset_list = [int(start_row_par_axis_data+7*i) for i in range(num_axes)]

            axis_list = ['Axis ' + str(i) for i in range(num_axes)]
            axis_loc_dict = {}
            for i in range(len(axis_list)):
                axis_loc_dict[axis_list[i]] = i

            for p in range(num_plots):
                plot_name = Input_Data_File['Parallel Axis Plotting'].cell(row = 2, column = 2 + p).value
                if plot_name is None:
                    plot_name = 'parallel_axis_plot'  # default figure name
                if plot_name not in plot_list:
                    plot_list.append(plot_name)
                if plot_name not in plot_list_dict[simulation_title]:
                    plot_list_dict[simulation_title].append(plot_name)
                    tradeoff_plot_pref_dict[simulation_title][plot_name] = {}
                #subplot_number = str(Input_Data_File['Parallel Axis Plotting'].cell(row = 3, column = 2 + p).value)
                subplot_number = 0  # For now, no subplots. Just single parallel axis plots.
                sp_list[simulation_title] = [subplot_number]

                tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number] = {}
                if subplot_number is not None:
                    tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number]['Color Axis Num'] = \
                        Input_Data_File['Parallel Axis Plotting'].cell(row = 5, column = 2 + p).value
                    tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number]['Colormap Name'] = \
                        Input_Data_File['Parallel Axis Plotting'].cell(row = 6, column = 2 + p).value
                    tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number]['Colorbar Title'] = \
                        Input_Data_File['Parallel Axis Plotting'].cell(row = 7, column = 2 + p).value
                    tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number]['Policies to Brush'] = \
                        Input_Data_File['Parallel Axis Plotting'].cell(row = 8, column = 2 + p).value
                    tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number]['Policy Labels'] = \
                        Input_Data_File['Parallel Axis Plotting'].cell(row = 9, column = 2 + p).value

                    # For this subplot, load preferences for each axis.
                    for axis in axis_list:
                        ax_ind = axis_list.index(axis)
                        row_ind = axis_offset_list[ax_ind]
                        if Input_Data_File['Parallel Axis Plotting'].cell(row = row_ind, column = 2 + p).value is not None:
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis] = {}
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['objective_name'] = \
                                Input_Data_File['Parallel Axis Plotting'].cell(row = row_ind, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['unit_conv'] = \
                                Input_Data_File['Parallel Axis Plotting'].cell(row = row_ind+1, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['perc_conv'] = \
                                Input_Data_File['Parallel Axis Plotting'].cell(row = row_ind+2, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['invert'] = \
                                Input_Data_File['Parallel Axis Plotting'].cell(row = row_ind+3, column = 2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['brush_range']\
                                = Input_Data_File['Parallel Axis Plotting'].cell(row=row_ind + 4, column=2 + p).value
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis]['title'] = \
                                Input_Data_File['Parallel Axis Plotting'].cell(row = row_ind+5, column = 2 + p).value
                            # Subplot title: Same for every axis
                            tradeoff_plot_pref_dict[simulation_title][plot_name][subplot_number][axis][
                                'subplot_title'] = Input_Data_File['Parallel Axis Plotting'].cell(row=4, column=2 + p).value
        except KeyError:
            print('Proper Parallel Axis Plotting worksheet does not exist to generate plots')

    ref_set_pref_dict = {}
    for plot_name in plot_list:
        # Create all of plots user specifies
        scenario_counter = 0
        for j in range(num_scenarios):
            simulation_title = simulation_titles_list[j]
            sub_plt_counter = 0
            if ref_set_file_name is None:
                ref_set_file_name = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets' + os_fold + \
                                  'pysedsim_ref_set.ref'

            objs_to_plot = opt_dicts_dict[simulation_title]['Borg_dict']['opt_dict']['Objective Names Ordered List']
            ref_set_pref_dict['ref_set_file_name'] = ref_set_file_name
            ref_set_pref_dict['num_objs'] = opt_dicts_dict[simulation_title]['Borg_dict']['n_objs']
            ref_set_pref_dict['num_dec_vars'] = opt_dicts_dict[simulation_title]['DPS_dict']['total_vars']
            ref_set_pref_dict['Borg_dict'] = opt_dicts_dict[simulation_title]['Borg_dict']
            for subplot_number in sp_list[simulation_title]:
                ref_set_pref_dict['unit_conv'] = [1 for i in range(len(objs_to_plot))]
                ref_set_pref_dict['perc_conv'] = ['No' for i in range(len(objs_to_plot))]
                ref_set_pref_dict['invert'] = ['No' for i in range(len(objs_to_plot))]
                plot_dict = {}
                spn_str = subplot_number
                plot_dict['num_axes'] = num_axes
                ref_set_pref_dict['num_objs_to_plot'] = num_axes
                plot_dict['brush_range'] = [[] for i in range(ref_set_pref_dict['num_objs'])]
                plot_dict['plot_order'] = [i for i in range(ref_set_pref_dict['num_objs_to_plot'])]

                # Import preferences for colorbar if user provided an axis number for colorbar.
                if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['Color Axis Num'] is not None:
                    plot_dict['Color Axis Num'] = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][
                        'Color Axis Num']
                if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['Colormap Name'] is not None:
                    plot_dict['Colormap Name'] = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][
                        'Colormap Name']
                if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['Colorbar Title'] is not None:
                    plot_dict['Colorbar Title'] = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][
                        'Colorbar Title']
                if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['Policies to Brush'] is not None:
                    plot_dict['Policies to Brush'] = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][
                        'Policies to Brush']
                if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['Policy Labels'] is not None:
                    plot_dict['Policy Labels'] = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][
                        'Policy Labels']
                try:
                    if tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['x_axis']['subplot_title'] is \
                            not None:
                        plot_dict['subplot_title'] = \
                        tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str]['x_axis']['subplot_title']
                except KeyError:
                    pass
                for obj in objs_to_plot:
                    for axis in axis_list:
                        try:
                            if obj == tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['objective_name']:
                                try:
                                    ref_set_pref_dict['unit_conv'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['unit_conv']
                                except KeyError:
                                    pass
                                try:
                                    ref_set_pref_dict['perc_conv'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['perc_conv']
                                except KeyError:
                                    pass
                                try:
                                    ref_set_pref_dict['invert'][objs_to_plot.index(obj)] = \
                                    tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis]['invert']
                                except KeyError:
                                    pass

                                # Order matters
                                plot_dict[axis + ' label'] = tradeoff_plot_pref_dict[simulation_title][plot_name][
                                    spn_str][axis]['title']

                                plot_dict['plot_order'][axis_loc_dict[axis]] = objs_to_plot.index(obj)

                                if type(tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis][
                                            'brush_range']) in [str, unicode]:
                                    # User has specified a list of numbers for the brushing range, which will be
                                    # imported as one long string, including commas.
                                    ax_rng = tradeoff_plot_pref_dict[simulation_title][plot_name][spn_str][axis][
                                        'brush_range'].split(', ')  # Create list of policies as strings
                                    # Convert to list
                                    plot_dict['brush_range'][objs_to_plot.index(obj)] = []
                                    for z in range(len(ax_rng)):
                                        try:
                                            plot_dict['brush_range'][objs_to_plot.index(obj)].append(float(ax_rng[z]))
                                        except ValueError:
                                            # User provided a string (either "min, min" or "max, max"). Can't do float.
                                            plot_dict['brush_range'][objs_to_plot.index(obj)].append(str(ax_rng[z]))
                        except KeyError:
                            pass
                # Create the parallel axis plot
                parallel_axis_plot(ref_set_pref_dict, plot_name=plot_name, plot_dict=plot_dict)
                # To keep track of where in main figure to place each subplot.
                sub_plt_counter += 1
            # To keep track of where in main figure to place each subplot.
            scenario_counter += 1


def parallel_axis_plot(ref_set_pref_dict, plot_dict = None, parse_objs = None, objs_to_plot = None, plot_name=None):
    #import seaborn.apionly as sns
    from copy import deepcopy

    # Load basic information about the simulated scenarios, including input file directory and simulation names.
    os_fold = Op_Sys_Folder_Operator()

    # Unpack ref_set_pref_dict
    try:
        ref_set_file_name = ref_set_pref_dict['ref_set_file_name']
        num_objs = ref_set_pref_dict['num_objs']
    except KeyError:
        ref_set_file_name = None
        num_objs = None
    try:
        num_dec_vars = ref_set_pref_dict['num_dec_vars']
    except KeyError:
        pass
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

    try:
        num_axes = plot_dict['num_axes']
    except KeyError:
        print("User did not specify number of axes")

    if ref_set_file_name is not None:
        # User wishes to import reference set here.
        [ref_set_array, objective_values, dec_var_values] = Initial_Processing(num_objs, num_dec_vars, ref_set_file_name,
                                                                               parse_objs=parse_objs, perc_conv=perc_conv,
                                                                               invert=invert, unit_conv=unit_conv,
                                                                               reverse_sign_all_objs = 'Yes')

        # Get positive values of objective values as well
        obj_vals_pos = Initial_Processing(num_objs, num_dec_vars, ref_set_file_name,
                                                                               parse_objs=parse_objs, perc_conv=perc_conv,
                                                                               invert=invert, unit_conv=unit_conv, reverse_sign_all_objs = 'Yes')[1]

    else:
        # User has provided reference set information.
        ref_set_array = ref_set_pref_dict['ref_set_array']
        objective_values = ref_set_pref_dict['objective_values']
        dec_var_values = ref_set_pref_dict['dec_var_values']

    # In case user specified a policy number to be brushed instead of a brushing range, specify the objective values
    policies_to_brush = None
    try:
        if plot_dict['Policies to Brush'] is not None:
            # Create list of policy number strings
            if type(plot_dict['Policies to Brush']) in [str, unicode]:
                policies_to_brush = plot_dict['Policies to Brush'].split(', ')
            else:
                policies_to_brush = [plot_dict['Policies to Brush']]
            # Create list of policy number integers
            policies_to_brush = [int(policies_to_brush[pol]) for pol in range(len(policies_to_brush))]
            for pol in policies_to_brush:
                for obj in range(num_objs):
                    # Append objective value for specified policy twice, so range is limited to exactly this value.
                    plot_dict['brush_range'][obj].append(obj_vals_pos[obj][pol])
                    plot_dict['brush_range'][obj].append(obj_vals_pos[obj][pol])
            # Replace list of policy strings with policies_to_brush integer list
            plot_dict['Policies to Brush'] = policies_to_brush
    except KeyError:
        plot_dict['label_policies'] = 'No'  # No policies to label of none are being brushed

    try:
        if plot_dict['Policy Labels'] is not None:
            # Create list of policy labels for those policies being highlighted
            plot_dict['label_policies'] = 'Yes'
            # Create list of policy number strings
            if type(plot_dict['Policy Labels']) in [str, unicode]:
                pol_names = plot_dict['Policy Labels'].split(', ')
            else:
                pol_names = None
            if pol_names is not None:
                plot_dict['Policy Labels'] = pol_names
            plot_dict['Label Locations'] = [(0.45, 0.05), (0.20, 0.95), (0.3, 0.8)]
            #(0.2, 0.95), (0.6, 0.95), (0.8, 0.95)
            plot_dict['Data Point Location'] = [(0.69, 0.33), (0.94, 0.71), (0.5,0.65)]
            # (0.22, 0.83), (0.62, 0.63), (0.82,0.83)
        else:
            plot_dict['label_policies'] = 'No'
    except KeyError:
        plot_dict['label_policies'] = 'No'

    # Unpack plotting dictionary
    try:
        plt_order = plot_dict['plot_order']
    except KeyError:
        plt_order = None
    try:
        three_d_plot = plot_dict['3d_plot']
    except KeyError:
        three_d_plot = 'No'
    try:
        plot_specs = plot_dict['Plot Specifications']
    except KeyError:
        plot_specs = None

    try:
        num_objs_to_plot = ref_set_pref_dict['num_objs_to_plot']
    except KeyError:
        num_objs_to_plot = len(objs_to_plot)

    # Set "Style"
    # sns.set_style("dark")  # Use seaborn to set gray background

    # set plotting characteristics
    #re_formatted_obj_vals = ref_set_array[:, num_dec_vars:]
    obj_vals_transposed = np.transpose(obj_vals_pos)  # Transpose to get num_columns = num_objs
    re_formatted_obj_values = np.copy(obj_vals_transposed)  # Initially populate to get correct dimensions,
    # then re-order cols.
    for x in range(num_objs):
        # Reformat so objective values are in the right order.
        re_formatted_obj_values[:, x] = obj_vals_transposed[:, plt_order[x]]
    formulations = [re_formatted_obj_values]
    labels = [[plot_dict['Axis ' + str(i) + ' label'] for i in range(num_axes)]]
    labels_orig = deepcopy(labels)
    try:
        cmaps = [plot_dict['Colormap Name']]
    except KeyError:
        try:
            # No colormap specified. See if user specified a file called 'custom_colormap.txt'. If so, use it to define
            # the colormap.
            new_cmap_matrix = np.loadtxt('custom_colormap.txt')
            cmaps = matplotlib.colors.ListedColormap(new_cmap_matrix/255.0)  # Must be a 256 color map.
        except IOError:
            # Use default (jet)
            cmaps = ["jet_r"]  # 'gray_4, 'jet', 'jet_r', cool_r, ['Reds_r','Blues_r','Greens_r','Purples_r']

    precision = [[0, 0, 0, 0, 0, 0]]
    # titles = ['Worst Case (WC)', 'Worst 1st Percentile (WP1)']
    # make 2 x 2 subplot with parallel axes for each problem formulation

    create_par_axis_plot(formulations, labels, labels_orig, precision, cmaps, plot_name, plot_dict, ref_set_pref_dict,
                         num_axes)  # titles

def create_par_axis_plot(formulations, labels, labels_orig, precision, cmaps, plot_name, plot_dict, ref_set_pref_dict,
                         num_axes, titles=None):

    # Imports
    from matplotlib.backends import backend_agg as agg  # raster backend
    import pandas as pd

    fig = matplotlib.figure.Figure()  # create the figure
    agg.FigureCanvasAgg(fig)         # attach the rasterizer
    ax = fig.add_subplot(1, 1, 1)    # make axes to plot on
    fontsize_all = 9
    top_bottom_label_font = 9

    shadeIndex = [0]
    for i in range(1):
        # used to be range(4)
        table = pd.DataFrame(formulations[i], columns=labels[i])
        mins = np.min(formulations[i], 0)
        maxs = np.max(formulations[i], 0)

        # Create bottom labels, which include both a numerical value and an objective label
        # round number of significant digits shown on objective labels
        for j in range(len(labels[i])):

            # First, ensure that /n specified by user are properly handled.
            x_ax_lab = ''
            x_ax_lab_components = labels[i][j].split(' \\n ')
            if len(x_ax_lab_components) > 1:
                for k in range(len(x_ax_lab_components)):
                    if k < len(x_ax_lab_components) - 1:
                        x_ax_lab += x_ax_lab_components[k] + '\n'
                    else:
                        x_ax_lab += x_ax_lab_components[k]
            else:
                x_ax_lab = labels[i][j]
            labels[i][j] = x_ax_lab
            # Handle labels that raise any number to a power (get rid of carrot symbol).
            try:
                power_index = labels[i][j].index('^')
                follow_val = labels[i][j][power_index+1]
            except ValueError:
                power_index = None
                follow_val = None
            if power_index is not None:
                labels[i][j] = labels[i][j].replace('^%s' % follow_val, '$\mathregular{^%s}\!$' % follow_val)  #

            if precision[i][j] != 0:
                labels[i][j] = str(np.round(mins[j], precision[i][j])) + '\n' + labels[i][j]
            else:
                labels[i][j] = str(int(mins[j])) + '\n' + labels[i][j]
            # don't show negative sign on maximization objectives
            if mins[j] < 0:
                labels[i][j] = labels[i][j][1:]

        # Create top labels, which don't include text (only numbers)
        toplabels = []
        # round number of significant digits shown on objective labels
        remove_minus_sign = 0
        for x in range(len(labels[i])):
            if precision[i][x] != 0:
                toplabels.append(str(np.round(maxs[x], precision[i][x])))
                if np.round(maxs[x], precision[i][x])<0:
                    remove_minus_sign = 1
            else:
                toplabels.append(str(int(maxs[x])))
                if int(maxs[x]) < 0:
                    remove_minus_sign = 1
            #if maxs[x] < 0:
            if remove_minus_sign == 1:
                # don't show negative sign on maximization objectives
                toplabels[x] = toplabels[x][1:]
            # Place a star (for theoretical ideal solution) at end of axis label
            toplabels[x] += ''  # '\n' # r'$\star$'
            #toplabels[-1] += r'$textcolor{red}{Today}$'

        table.columns = labels[i]

        if type(cmaps) in [list]:
            # User specified matplotlib standard color map name, or a default name was used. Input file of custom
            # colors was not specified.
            cmaps = pyplot.cm.get_cmap(cmaps[i])
        scaled = table.copy()

        index = 0
        try:
            color_index = plot_dict['Color Axis Num']  # Use specified column values for color bar.
        except KeyError:
            color_index = 0  # Use zero-th axis as default value for colorbar if none specified by user.

        for column in table.columns:
            scaled[column] = (table[column] - mins[index]) / (maxs[index] - mins[index])
            index = index + 1

        it_scaled = scaled.copy()  # Table to be iteratively paired down to reflect brushing preferences.

        # Handle any brushing the user wants to occur.
        num_brushes = 0
        for z in range(num_axes):
            if plot_dict['brush_range'][plot_dict['plot_order'][z]] != []:
                num_brushes = int(len(plot_dict['brush_range'][plot_dict['plot_order'][z]])/2)
                break  # Know how many brushes based on first axis for which brushed range is listed.

        # Brushed policy subset 1
        brush_col_names = [[] for b in range(num_brushes)]
        brush_col_range = [[] for b in range(num_brushes)]
        for brush_num in range(num_brushes):
            brush_col_names[brush_num] = list(scaled.columns)
            for z in range(num_axes):
                if plot_dict['brush_range'][plot_dict['plot_order'][z]] == []:
                    value = [0, 1]
                    brush_col_range[brush_num].append(value)
                else:
                    value = plot_dict['brush_range'][plot_dict['plot_order'][z]][2*brush_num:2*brush_num+2]
                    # User may have specified either "min, min" or "max, max" to get highest line.
                    for m in range(len(value)):
                        if type(value[m]) in [str, unicode]:
                            if value[m] in ['Min', 'min']:
                                value[m] = mins[z]
                            elif value[m] in ['Max', 'max']:
                                value[m] = maxs[z]
                    value = [float((value[m] - mins[z])/(maxs[z] - mins[z])) for m in range(len(value))]
                    brush_col_range[brush_num].append(value)

        #brush_col_name_1 = [scaled.columns[1]]  # Max annual energy policy
        #brush_col_name_2 = [scaled.columns[0], scaled.columns[4]]  # Max wet/dry larvae passage policy
#        brush3_cols = [scaled.columns[0], scaled.columns[4]]  # Compromise policy

        #brushed_set = [it_scaled.copy() for b in range(num_brushes)]
        brushed_set = [b for b in range(num_brushes)]
        for brush_num in range(num_brushes):
            for a in range(num_axes):
                if num_axes == 6:
                    brushed_set[brush_num] = \
                        (scaled[brush_col_names[brush_num][0]] >= brush_col_range[brush_num][0][0]) & \
                        (scaled[brush_col_names[brush_num][0]] <= brush_col_range[brush_num][0][1]) & \
                        (scaled[brush_col_names[brush_num][1]] >= brush_col_range[brush_num][1][0]) & \
                        (scaled[brush_col_names[brush_num][1]] <= brush_col_range[brush_num][1][1]) & \
                        (scaled[brush_col_names[brush_num][2]] >= brush_col_range[brush_num][2][0]) & \
                        (scaled[brush_col_names[brush_num][2]] <= brush_col_range[brush_num][2][1]) & \
                        (scaled[brush_col_names[brush_num][3]] >= brush_col_range[brush_num][3][0]) & \
                        (scaled[brush_col_names[brush_num][3]] <= brush_col_range[brush_num][3][1]) & \
                        (scaled[brush_col_names[brush_num][4]] >= brush_col_range[brush_num][4][0]) & \
                        (scaled[brush_col_names[brush_num][4]] <= brush_col_range[brush_num][4][1]) & \
                        (scaled[brush_col_names[brush_num][5]] >= brush_col_range[brush_num][5][0]) & \
                        (scaled[brush_col_names[brush_num][5]] <= brush_col_range[brush_num][5][1])
                elif num_axes == 5:
                    brushed_set[brush_num] = \
                        (scaled[brush_col_names[brush_num][0]] >= brush_col_range[brush_num][0][0]) & \
                        (scaled[brush_col_names[brush_num][0]] <= brush_col_range[brush_num][0][1]) & \
                        (scaled[brush_col_names[brush_num][1]] >= brush_col_range[brush_num][1][0]) & \
                        (scaled[brush_col_names[brush_num][1]] <= brush_col_range[brush_num][1][1]) & \
                        (scaled[brush_col_names[brush_num][2]] >= brush_col_range[brush_num][2][0]) & \
                        (scaled[brush_col_names[brush_num][2]] <= brush_col_range[brush_num][2][1]) & \
                        (scaled[brush_col_names[brush_num][3]] >= brush_col_range[brush_num][3][0]) & \
                        (scaled[brush_col_names[brush_num][3]] <= brush_col_range[brush_num][3][1]) & \
                        (scaled[brush_col_names[brush_num][4]] >= brush_col_range[brush_num][4][0]) & \
                        (scaled[brush_col_names[brush_num][4]] <= brush_col_range[brush_num][4][1])
                elif num_axes == 4:
                    brushed_set[brush_num] = \
                        (scaled[brush_col_names[brush_num][0]] >= brush_col_range[brush_num][0][0]) & \
                        (scaled[brush_col_names[brush_num][0]] <= brush_col_range[brush_num][0][1]) & \
                        (scaled[brush_col_names[brush_num][1]] >= brush_col_range[brush_num][1][0]) & \
                        (scaled[brush_col_names[brush_num][1]] <= brush_col_range[brush_num][1][1]) & \
                        (scaled[brush_col_names[brush_num][2]] >= brush_col_range[brush_num][2][0]) & \
                        (scaled[brush_col_names[brush_num][2]] <= brush_col_range[brush_num][2][1]) & \
                        (scaled[brush_col_names[brush_num][3]] >= brush_col_range[brush_num][3][0]) & \
                        (scaled[brush_col_names[brush_num][3]] <= brush_col_range[brush_num][3][1])
                elif num_axes == 3:
                    brushed_set[brush_num] = \
                        (scaled[brush_col_names[brush_num][0]] >= brush_col_range[brush_num][0][0]) & \
                        (scaled[brush_col_names[brush_num][0]] <= brush_col_range[brush_num][0][1]) & \
                        (scaled[brush_col_names[brush_num][1]] >= brush_col_range[brush_num][1][0]) & \
                        (scaled[brush_col_names[brush_num][1]] <= brush_col_range[brush_num][1][1]) & \
                        (scaled[brush_col_names[brush_num][2]] >= brush_col_range[brush_num][2][0]) & \
                        (scaled[brush_col_names[brush_num][2]] <= brush_col_range[brush_num][2][1])
                elif num_axes == 2:
                    brushed_set[brush_num] = \
                        (scaled[brush_col_names[brush_num][0]] >= brush_col_range[brush_num][0][0]) & \
                        (scaled[brush_col_names[brush_num][0]] <= brush_col_range[brush_num][0][1]) & \
                        (scaled[brush_col_names[brush_num][1]] >= brush_col_range[brush_num][1][0]) & \
                        (scaled[brush_col_names[brush_num][1]] <= brush_col_range[brush_num][1][1])

    #                (brushed_set[brush_num][brush_col_names[brush_num][a]] >= brush_col_range[brush_num][a][0]) & \
    #                (brushed_set[brush_num][brush_col_names[brush_num][a]] <= brush_col_range[brush_num][a][1])

        # Brushed policy subset 2
        #brushed_set[1] = (scaled[brush_col_name_2[0]] > 0.98) & (scaled[brush_col_name_2[1]] > 0.98)

        # Brushed policy subset 3
        #brushing_3 = (scaled[brush3_cols[0]] > 0.9) & (scaled[brush3_cols[1]] < 0.5) & (scaled[brush3_cols[2]] < 0.5) \
        #             & (scaled[brush3_cols[3]] < 0.5) & (scaled[brush3_cols[4]] < 0.5) & (scaled[brush3_cols[5]] < 0.5)

        # Plot in gray ALL policies in set, but only if brushing of some sort is going to occur. Use color axis.
        # See if brushing is actually producing lines. If not, don't plot anything, not even gray lines in
        # background. Used for purposes of producing a blank par axis plot for figure build.

        solution_counter = 0
        plot_all_policies_in_background = 'Yes'
        plot_all_gray_background = 'No'
        plot_background_lighter = 'Yes'

        if num_brushes == 0:
            # No brushing will occur.
            for solution in scaled.iterrows():
                # for solution in scaled[brushing_1 ^ True].iterrows():
                ys = solution[1]
                xs = range(len(ys))
                #ax.plot(xs, ys, c=cmaps(ys[shadeIndex[i]]), linewidth=0.5)
                ax.plot(xs, ys, c=cmaps(ys[color_index]), alpha = 1.0, linewidth=0.5)  # Colorbar, set to column 1
        else:
            # Some brushing is going to occur. So, bring all solutions on in gray first, then brush out some to color.
            # User must indicate some form of brushing in input worksheet to get all gray lines in background.
            for solution in scaled.iterrows():
                # for solution in scaled[brushing_1 ^ True].iterrows():
                for set_mem in range(num_brushes):
                    for solution2 in scaled[brushed_set[set_mem]].iterrows():
                        solution_counter += 1  # At least one brushed line exists
                if solution_counter >= 1 or plot_all_policies_in_background == 'Yes':
                    ys = solution[1]
                    xs = range(len(ys))
                    #ax.plot(xs, ys, c=cmaps(ys[shadeIndex[i]]), linewidth=0.5)
                    if plot_all_gray_background == 'Yes':
                        ax.plot(xs, ys, color=(0.85, 0.85, 0.85), linewidth=0.5)  # Gray lines
                    else:
                        if plot_background_lighter == 'Yes':
                            ax.plot(xs, ys, c=cmaps(ys[color_index]), alpha=0.25, linewidth=0.5)  # 0.5, 0.1
                        else:
                            ax.plot(xs, ys, c=cmaps(ys[color_index]), alpha=1, linewidth=0.5)  # 0.5, 0.1
            # Plot policies in brushing groups
            for set_mem in range(num_brushes):
                for solution in scaled[brushed_set[set_mem]].iterrows():
                    ys = solution[1]
                    xs = range(len(ys))
                    #ax.plot(xs, ys, c='k', linewidth=1.5)  # Colorbar, set to column 1
                    ax.plot(xs, ys, c=cmaps(ys[color_index]), linewidth=2,
                            path_effects=[pe.Stroke(linewidth=3.0, foreground='k'), pe.Normal()])  # 1.5,  '--',

        # Having plotted solutions, now label them if user desires
        # Get location
        if plot_dict['label_policies'] == 'Yes':
            for x in range(len(plot_dict['Policies to Brush'])):
                ax.annotate(plot_dict['Policy Labels'][x], fontsize=8, xy=plot_dict['Data Point Location'][x],
                            xycoords='data', xytext=plot_dict['Label Locations'][x],
                            bbox=dict(boxstyle="round", fc='w'),
                            arrowprops=dict(arrowstyle="->", connectionstyle="arc,angleA=0,angleB=90,"
                                                                             "rad=10"))

        # Create and label colorbar
        # Create and label colorbar
        include_colorbar = 'Yes'
        if include_colorbar == 'Yes':
            # User wants a colorbar. Create a colorbar label that is formatted well, taking into account user
            # preferences regarding line returns and superscripts.
            sm = matplotlib.cm.ScalarMappable(cmap=cmaps)
            #sm.set_array([0.0,1.0])
            sm.set_array([mins[color_index], maxs[color_index]])
            cbar = fig.colorbar(sm)  # , shrink=0.75
            # If colorbar title has '\n' values in it, reformat them.
            try:
                cbar_lab = ''
                cbar_lab_components = plot_dict['Colorbar Title'].split(' \\n ')
                if len(cbar_lab_components) > 1:
                    for k in range(len(cbar_lab_components)):
                        if k < len(cbar_lab_components) - 1:
                            cbar_lab += cbar_lab_components[k] + '\n'
                        else:
                            cbar_lab += cbar_lab_components[k]
                else:
                    cbar_lab = plot_dict['Colorbar Title']

                # Replace ^ with sup-s
                try:
                    power_index = cbar_lab.index('^')
                    follow_val = cbar_lab[power_index+1]
                except ValueError:
                    power_index = None
                    follow_val = None
                if power_index is not None:
                    cbar_lab = cbar_lab.replace('^%s' % follow_val, '$\mathregular{^%s}\!$' % follow_val)  #
                actual_colorbar_label = cbar_lab
            except KeyError:
                pass
            # Create the colorbar label
            cbar.ax.set_ylabel(actual_colorbar_label, fontsize=fontsize_all)
            cbar.ax.tick_params(labelsize=fontsize_all)  # Set size of tick labels

        #ax.set_title(title[i],fontsize=fontsize_all,y=1.1)
        ax.set_xticks(np.arange(0,np.shape(table)[1],1))
        ax.set_xticklabels(labels[i], fontsize=top_bottom_label_font)
        ax.tick_params(axis='y', which='both', labelleft='off', left='off', right='off')
        ax.tick_params(axis='x', which='both', top='off', bottom='off')

        ax.set_ylabel('Direction of Preference $\\longrightarrow$', fontsize=fontsize_all)

        # make subplot frames invisible
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["right"].set_visible(True)

        # draw in axes
        for k in np.arange(0,np.shape(table)[1],1):
            ax.plot([k, k],[0,1],c='k')

        # create twin y axis to put x tick labels on top
        ax2 = ax.twiny()
        ax2.set_xticks(np.arange(0,np.shape(table)[1],1))
        ax2.set_xticklabels(toplabels, fontsize=top_bottom_label_font)

        # Attempt to use scientific notation on axis
        #formatter = matplotlib.ticker.ScalarFormatter()
        #formatter.set_powerlimits((-3, 5))
        #ax2.yaxis.set_major_formatter(formatter)
        #ax2.xaxis.set_major_formatter(formatter)

        ax2.tick_params(axis='y',which='both',labelleft='off',left='off',right='off')
        ax2.tick_params(axis='x',which='both',top='off',bottom='off')

        # make subplot frames invisible
        ax2.spines["top"].set_visible(False)
        ax2.spines["bottom"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        #ax.figsize=[4, 4]

    fig.set_size_inches(18.5/2.54, 3.5)  # 4
    fig.tight_layout()
    #fig.figsize=[5,5]
    fig.savefig(plot_name + '.pdf', bbox_inches='tight', dpi=1200, pad_inches=0.01)