# -*- coding: utf-8 -*-

'''
Purpose: To call the Borg MOEA (through the Borg wrapper) for purposes of identifying a pareto set of non-dominated
points in the objective space, each of which corresponds to a different reservoir operating policy.

For now, this particular module should not be called in parallel. Simulations are done in parallel (distributed
over the processors) for each seed (e.g., 1-50) in the optimization. The "function evaluation" is the PySedSim model,
which requires no inputs when being called in deterministic mode.
'''

# Imports
import numpy as np
from pysedsim.data_processing.data_processing import Op_Sys_Folder_Operator
from pysedsim.data_processing.data_processing import Determine_Num_Scenarios
from pysedsim.data_processing.data_processing import Load_Input_File
import pysedsim.optimization.direct_policy_search
import os
import platform
import logging

# Load in user-specified optimization assumptions from input file, which first must be located using its path that is
# stored in a top level input file. This should be done only once for the entire optimization.

def PySedSim_Caller(vars, additional_inputs):
    '''
    Purpose: Calls the PySedSim simulation model. Designed to be called by Borg wrapper in each function evaluation
    as part of a simulation-optimization experiment.

    Decision variable values from Borg can either be
    sent directly into PySedSim through the function interface (ideal for parallel function evaluations being
    performed by Master-Slave Borg, but works for Serial Borg as well), or can be dumped into a text file (works only
    for Serial Borg evaluations).

    Args:
        vars: A list of decision variable values from Borg
        additional_inputs: A list of two dictionaries and the input file name. The first dictionary is a dictionary of
        Borg preferences, and the second is a dictionary of Direct Policy Search optimization preferences.

    Returns:
        performance: objective value performance. A list of objective values, one for each of the objectives.
    '''
    Borg_dict = additional_inputs[0]
    DPS_Dict = additional_inputs[1]
    top_level_input_file_name = additional_inputs[2]
    borg_vars = vars
    # If any variables are discrete valued (e.g., integer), re-format borg variable values to reflect.
    for v in range(len(DPS_Dict['variable_type'])):
        if DPS_Dict['variable_type'][v] == 'integer':
            borg_vars[v] = int(borg_vars[v])

    # Call PySedSim according to user preferences for interface.
    if Borg_dict['link_to_borg'] == 0:
        # Use for deterministic simulation
        np.savetxt('RBF_Parameters.txt', borg_vars, newline = ' ')  # Save policy so PySedSim can import it.
        performance = pysedsim.PySedSim(borg_dict=Borg_dict, file_name=top_level_input_file_name)
    else:
        # Use for stochastic simulation/optimization, particularly if being run in parallel.
        op_policy_params = np.asarray(borg_vars)  # Cast borg output parameters as array for use in simulation.
        performance = pysedsim.PySedSim(borg_dict=Borg_dict, dps_dict = DPS_Dict, decision_vars = op_policy_params,
                                        file_name = top_level_input_file_name)  #
    return performance

def Optimization_DPS_Borg(file_name = 'PySedSim_Input_Specifications.csv'):

    '''

    Purpose: Perform main Borg optimization loop through scenarios and seeds. Calls PySedSim for simulation.

    Args:
        file_name: Name of the top-level input file that indicates where input files are located, etc. Default is
        'PySedSim_Input_Specifications.csv'. Must be the same file used for the simulation that is to be connected to
        the optimization.

    Returns:
        Produces pareto outputs for each seed, exported to separate files in a "sets" folder created in the model
        output location specified in the input file.
    '''

    import borg as bg  # Import borg wrapper
    from processing_reference_set import Reference_Set
    from processing_reference_set import Policy_Reevaluation
    from processing_reference_set import Scatter_Plot_Objectives

    # Get operator for changing directory based on operating system, then import various assumptions about the simulation
    # from the top level input file.
    os_fold = Op_Sys_Folder_Operator()
    [num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir,
     main_output_file_dir, Monte_Carlo_Parameters_File_List, external_mc_data] = Determine_Num_Scenarios(file_name)

    # Loop through as many optimization scenarios as user has specified
    for j in range(num_scenarios):
        simulation_title = simulation_titles_list[j]
        Borg_dict = Import_Optimization_Preferences(simulation_title, imported_specs, main_input_files_dir)
        if Borg_dict['optimization approach'] == 'DPS':
            DPS_dict = direct_policy_search.Import_DPS_Preferences(simulation_title=simulation_title,
                                                                   imported_specs=imported_specs,
                                                                   main_input_files_dir=main_input_files_dir)
            Borg_dict['n_vars'] = DPS_dict['n_vars']  # Store num decision variables in Borg_dict, copied from DPS_dict.

        # Populate DPS input variable ranges to specify for borg
        borg_variable_range = []
        # Direct Policy Search Radial Basis Function Parameters:
        for r in range(DPS_dict['RBF Specs']['Number']):
            for m in range(DPS_dict['num_inputs']):
                borg_variable_range.append([-1, 1])  # for normalized center
                borg_variable_range.append([0, 1])  # for normalized radius
            for k in range(DPS_dict['num_reservoirs']):
                borg_variable_range.append([0, 1])  # for weights, which are normalized later on as they must sum to 1.

        # Populate Flushing-related input variable ranges to specify for borg
        try:
            for flush_var in DPS_dict['Flushing Optimization']['Ordered input variable list']:
                borg_variable_range.append(DPS_dict['Flushing Optimization']['Input variable'][flush_var])
        except KeyError:
            pass

        # Compute total number of variables, which is a combination of DPS policy decision vars and flushing decision
        #  vars.
        Borg_dict['total_vars'] = DPS_dict['total_vars']

        # Create directory in which to store optimization output
        output_location = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets'
        if not os.path.exists(output_location):
            os.makedirs(output_location)

        # Set up interface with Borg for parallel function evaluation, if Borg Maser-Slave is being used, and if this
        #  is the first optimization scenario being considered. Can't start MPI twice.
        if Borg_dict['borg_type'] == 1 and j == 0:
            bg.Configuration.startMPI()  # start parallelization with MPI

        # Set up interface with Borg for optimization
        for j in range(Borg_dict['nSeeds']):
            borg = bg.Borg(Borg_dict['total_vars'], Borg_dict['n_objs'], Borg_dict['n_constrs'], PySedSim_Caller,
                        add_pysedsim_inputs = [Borg_dict, DPS_dict, file_name])  # Create instance of Borg class
            borg.setBounds(*borg_variable_range)  # Set decision variable bounds
            borg.setEpsilons(*Borg_dict['epsilon_list'])  # Set epsilon values
            runtime_filename = main_output_file_dir + os_fold + simulation_title + os_fold + 'sets' + os_fold + \
                               'runtime_file_seed_' + str(j+1) + '.runtime'
            if Borg_dict['borg_type'] == 0:
                # Run serial Borg
                if Borg_dict['runtime_preferences']['runtime_choice'] == 'Yes':
                    result = borg.solve({"maxEvaluations": Borg_dict['num_func_evals'], "runtimeformat": 'borg',
                                         "frequency": Borg_dict['runtime_preferences']['runtime_freq'],
                                         "runtimefile": runtime_filename})
                else:
                    result = borg.solve({"maxEvaluations": Borg_dict['num_func_evals']})

                #result = borg.solve({"maxEvaluations":Borg_dict['num_func_evals']})  # Minimum is 100 in Borg
            elif Borg_dict['borg_type'] == 1:
                # Run parallel Borg
                if Borg_dict['runtime_preferences']['runtime_choice'] == 'Yes':
                    result = borg.solveMPI(maxEvaluations=Borg_dict['num_func_evals'], runtime=runtime_filename,
                                           frequency = Borg_dict['runtime_preferences']['runtime_freq'])
                else:
                    result = borg.solveMPI(maxEvaluations=Borg_dict['num_func_evals'])
            if result:
                # This particular seed is now finished being run in parallel. The result will only be returned from
                # one node in case running Master-Slave Borg.
                result.display()

                # Create/write objective values and decision variable values to files in folder "sets", 1 file per seed.
                f = open(output_location + os_fold + 'Borg_DPS_PySedSim' + str(j+1) + '.set', 'w')
                f.write('#Borg Optimization Results\n')
                f.write('#First ' + str(Borg_dict['total_vars']) + ' are the decision variables, ' + 'last ' + str(
                    Borg_dict['n_objs']) + ' are the ' + 'objective values\n')
                f.write('#List of objective names (in order of appearance): ' + ', '.join(
                    Borg_dict['opt_dict']['Objective Names Ordered List']) + '\n')
                for solution in result:
                    line = ''
                    for i in range(len(solution.getVariables())):
                        line = line + (str(solution.getVariables()[i])) + ' '

                    for i in range(len(solution.getObjectives())):
                        line = line + (str(solution.getObjectives()[i])) + ' '

                    f.write(line[0:-1]+'\n')
                f.write("#")
                f.close()

                # Create/write objective values only to files in folder "sets", 1 file per seed. Purpose is so that
                # the file can be processed in MOEAFramework, where performance metrics may be evaluated across seeds.
                f2 = open(output_location + os_fold + 'Borg_DPS_PySedSim_no_vars' + str(j+1) + '.set', 'w')
                for solution in result:
                    line = ''
                    for i in range(len(solution.getObjectives())):
                        line = line + (str(solution.getObjectives()[i])) + ' '

                    f2.write(line[0:-1]+'\n')
                f2.write("#")
                f2.close()

                logging.info("Seed {0} complete".format(j))

        # All seeds for this scenario are complete. Now generate reference set if user desires.
        if Borg_dict['ref_set_yes_no'] == 'Yes':
            # Ensure that this reference set processing is done only once in event this is being run in parallel.
            #if platform.system() == 'Linux':
                #[start, stop, rank] = parallelize_pysedsim(1)  # Distribute re-evaluation work to processors.
                #if rank == 0:
            [Borg_dict, DPS_dict] = Reference_Set(inline_optimization=[simulation_title, Borg_dict, DPS_dict,
                                                                   main_output_file_dir])
            #else:
                #[Borg_dict, DPS_dict] = Reference_Set(inline_optimization=[simulation_title, Borg_dict, DPS_dict,main_output_file_dir])

        # Re-evaluate (simulate) all of the policies from the reference set.
        if Borg_dict['re_eval'] == 'Yes':
            if Borg_dict['ref_set_yes_no'] == 'Yes':
                ref_set_file_name = output_location + os_fold + Borg_dict['opt_output_filename']
                [ref_set_array, objective_values, dec_var_values] = Policy_Reevaluation(
                    ref_set_file_name=ref_set_file_name, post_optimization='No', file_name=file_name, dps_dict=DPS_dict,
                    borg_dict=Borg_dict, parallel=Borg_dict['parallel_choice'],
                    internal_opt_call=[Borg_dict['opt_dict']['Scenario Name']])
            else:
                # Do not specify a reference set file name, as user did not generate a reference set with
                # optimization. Instead, just use the default reference set file name from the Policy_Reevaluation
                # method.
                [ref_set_array, objective_values, dec_var_values, dps_dict, borg_dict] = Policy_Reevaluation(post_optimization='No',
                                                                                        file_name=file_name,
                                                                                        dps_dict=DPS_dict,
                                                                                        borg_dict=Borg_dict,
                                                                                        parallel=Borg_dict['parallel_choice'],
                                                                                        internal_opt_call=[Borg_dict['opt_dict']['Scenario Name']])

        if Borg_dict['ref_set_plot'] == 'Yes':
            ref_set_pref_dict = {'num_objs': Borg_dict['n_objs'], 'num_dec_vars': DPS_dict['n_vars'], 'invert': [
                'No', 'No', 'No'], 'perc_conv': ['No', 'Yes', 'No'], 'unit_conv': [.365, 1, 1]}
            movie_dict = {'Create Movie': 'No'}
            plot_dict = {'Axis Range': [[0, 6000], [0, 100], [0, 8000]], 'plot_order': [1, 0, 2], '3d_plot': 'No'}
            objs_to_plot = Borg_dict['opt_dict']['Objective Names Ordered List']  # Plot all objectives in optimization.
            if Borg_dict['re_eval'] == 'Yes':
                # Reference set has already been processed. Store those values in the reference set preferences dict.
                ref_set_pref_dict['ref_set_array'] = ref_set_array
                ref_set_pref_dict['objective_values'] = objective_values
                ref_set_pref_dict['dec_var_values'] = dec_var_values
                # Produce/save plots.
                Scatter_Plot_Objectives(ref_set_pref_dict, objs_to_plot, plot_dict=plot_dict, save_fig=output_location,
                                        movie_dict=movie_dict)
            else:
                ref_set_pref_dict['ref_set_file_name'] = output_location + os_fold + Borg_dict['opt_output_filename']
                ref_set_pref_dict['num_dec_vars'] = DPS_dict['n_vars']
                ref_set_pref_dict['num_objs'] = Borg_dict['n_objs']
                # Produce/save plots.
                [ref_set_array, objective_values, dec_var_values] = Scatter_Plot_Objectives(ref_set_pref_dict, objs_to_plot,
                                                                                plot_dict=plot_dict,
                                                                                save_fig= output_location,
                                                                                movie_dict=movie_dict)

    if Borg_dict['borg_type'] == 1:
        bg.Configuration.stopMPI()  # stop parallel function evaluation process

    return

def Import_Optimization_Preferences(simulation_title, imported_specs, main_input_files_dir):

    '''

    Purpose: To import preferences related to PySedSim reservoir operating policy optimization. Specify these
    preferences in the "Optimization" worksheet in the input data (.xlsx) file.

    Args:
        simulation_title: A list of strings of runs to conduct ("scenarios", each of which has a separate input file)
        imported_specs: A list of details regarding file locations for the scenario run(s) to be conducted,
        which represents the data contained in the top-level PySedSim input file (default file name:
        "PySedSim_Input_Specifications.csv")
        main_input_files_dir: String that contains the file path to the directory that contains the input file(s) for
        a given scenario.

    Returns:
        Borg_dict: A dictionary of Borg preferences, specified by the user in the "Optimization" worksheet of the
        pysedsim input file for the current scenario.

    '''

    Input_Data_File = Load_Input_File(simulation_title, main_input_files_dir, imported_specs)

    # Gather user-specified optimization preferences from "Optimization" worksheet if it exists.
    try:
        opt_approach = Input_Data_File['Optimization']['B2'].value
        if (opt_approach is None) or (type(opt_approach) not in [str, unicode]):
            opt_approach = 'DPS'  # Set default. Currently, only DPS is supported.

        nSeeds = Input_Data_File['Optimization']['B3'].value
        if (nSeeds is None) or (type(nSeeds) not in [int, long]):
            nSeeds = 1  # Default is 1 seeds if user did not specify anything, or didn't specify a number

        link_to_borg = Input_Data_File['Optimization']['B4'].value
        if link_to_borg in ['Yes', 'yes']:
            link_to_borg = 1
        elif link_to_borg in ['No', 'no']:
            link_to_borg = 0
        else:
            link_to_borg = 0  # Default if no user specified value
        borg_type = Input_Data_File['Optimization']['B5'].value
        if borg_type in ['Serial', 'serial']:
            borg_type = 0
        elif borg_type in ['Master-Slave', 'master-slave', 'Master Slave', 'master slave']:
            borg_type = 1
        else:
            borg_type = 0  # Default is Serial if no user specified value

        # Function evaluations
        num_func_evals = Input_Data_File['Optimization']['B6'].value
        if (num_func_evals is None) or (type(num_func_evals) not in [int, long]):
            num_func_evals = 100  # Set minimum number of function evaluations if user did not specify an integer.

        # Number of objectives
        n_objs = Input_Data_File['Optimization']['B7'].value
        if (n_objs is None) or (type(n_objs) not in [int, long]):
            n_objs = 1  # Default is a 1 objective optimization problem.
        else:
            n_objs = int(n_objs)

        n_constrs = Input_Data_File['Optimization']['B8'].value
        if (n_constrs is None) or (type(n_constrs) not in [int, long]):
            n_constrs = 0  # Default is a zero constraints in objective space
        else:
            n_constrs = int(n_constrs)

        maxtime = Input_Data_File['Optimization']['B9'].value
        if (maxtime is None) or (type(maxtime) not in [int, long]):
            maxtime = None  # Default is a zero constraints in objective space

        # Whether or not to produce reference set using Pareto.py
        ref_set = Input_Data_File['Optimization']['B10'].value
        if (ref_set is None) or (type(ref_set) not in [str, unicode]):
            ref_set = 'No'  # Default is a 1 objective optimization problem.
        else:
            if ref_set in ['yes', 'Yes', 'y', 'Y']:
                ref_set = 'Yes'
            else:
                ref_set = 'No'

        # Whether or not to produce reference set using Pareto.py, and if yes, whether or not to re-evaluate policies
        # in parallel after each optimization scenario is completed.
        re_eval = Input_Data_File['Optimization']['B11'].value
        if (re_eval is None) or (type(re_eval) not in [str, unicode]):
            re_eval = 'No'  # Default is not to automatically re-evaluate (simulate) all points in reference set.
            parallel_choice = 'No'
        else:
            if re_eval in ['yes', 'Yes', 'y', 'Y']:
                re_eval = 'Yes'
                # Whether or not to re-evaluate reference set in parallel.
                parallel_choice = Input_Data_File['Optimization']['B12'].value
                if (parallel_choice is None) or (type(parallel_choice) not in [str, unicode]):
                    parallel_choice = 'No'  # Default is not to automatically re-evaluate (simulate) all points in reference set.
                else:
                    if parallel_choice in ['yes', 'Yes', 'y', 'Y']:
                        parallel_choice = 'Yes'
                    else:
                        parallel_choice = 'No'
            else:
                re_eval = 'No'
                parallel_choice = 'No'

        # Establish user's preferences regarding plotting 2D tradeoffs in objective space, and saving plots?
        ref_set_plot = Input_Data_File['Optimization']['B13'].value
        if (ref_set_plot is None) or (type(ref_set_plot) not in [str, unicode]):
            ref_set_plot = 'No'  # Default is a 1 objective optimization problem.
        else:
            if ref_set_plot in ['yes', 'Yes', 'y', 'Y']:
                ref_set_plot = 'Yes'
            else:
                ref_set_plot = 'No'

        # Establish user's preferences regarding printing Borg runtime dynamics
        runtime_freq = Input_Data_File['Optimization']['B14'].value
        if (runtime_freq is None) or (type(runtime_freq) not in [int, long, float]):
            runtime_choice = 'No'  # Default is a 1 objective optimization problem.
        else:
            runtime_choice = 'Yes'
            runtime_freq = int(runtime_freq)
        runtime_preferences = {'runtime_choice': runtime_choice, 'runtime_freq': runtime_freq}
        # Establish preference for name of reference set file produced.
        # Currently, user cannot choose file name, but this section can be modified later, with a place added in the
        # input file for this string to be specified.

        opt_output_filename = 'pysedsim_ref_set.ref'
        obj_row = 16
        # Read in objective variable names and preferences.
        obj_var_pref = {}  # Objective variable preferences. Dictionary of variable names.
        obj_names_ordered = [0 for i in range(n_objs)]
        for i in range(n_objs):
            obj_name = Input_Data_File['Optimization'].cell(row = obj_row+i, column=1).value
            obj_names_ordered[i] = obj_name
            obj_var_pref[obj_name] = {}  # Each variable name key stores a sub dictionary of related preferences.
            obj_var_pref[obj_name]['State Variable'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=2).value

            # Determine whether variable is minimization or maximization
            obj_var_pref[obj_name]['Type'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=3).value
            # Handle various ways of indicating whether variable is minimization or maximization.
            if obj_var_pref[obj_name]['Type'] not in ['Max', 'Min']:
                if obj_var_pref[obj_name]['Type'] in ['max', 'MAX', 'maximize', 'MAXIMIZE']:
                    obj_var_pref[obj_name]['Type'] = 'Max'
                elif obj_var_pref[obj_name]['Type'] in ['min', 'MIN', 'minimize', 'MINIMIZE']:
                    obj_var_pref[obj_name]['Type'] = 'Min'
                elif obj_var_pref[obj_name]['Type'] in [None]:
                    obj_var_pref[obj_name]['Type'] = 'Min'  # Assume minimization if nothing specified

            # Handle various ways of indicating what the sampling frequency in pandas should be. Options: Daily,
            # Monthly, Annual. Default: Daily.
            # For example, in a 100-year simulation, a monthly resampling would mean you draw 100*12 different
            # samples for the corresponding state variable. For each of those 1200 samples, some statistic is
            # calculated (e.g., a sum, a mean, etc.). That statistic preference is in the next section ('Resample Stat')
            obj_var_pref[obj_name]['Resample Frequency'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=4).value
            if obj_var_pref[obj_name]['Resample Frequency'] not in ['A', 'M', 'D', 'O']:
                if obj_var_pref[obj_name]['Resample Frequency'] in ['Annual', 'annual', 'Annually', 'annually', 'a',
                                                                    'ann', 'Ann']:
                    obj_var_pref[obj_name]['Resample Frequency'] = 'A'
                elif obj_var_pref[obj_name]['Resample Frequency'] in ['Monthly', 'month', 'Month', 'monthly', 'm',
                                                                      'mon', 'Mon']:
                    obj_var_pref[obj_name]['Resample Frequency'] = 'M'
                elif obj_var_pref[obj_name]['Resample Frequency'] in ['Daily', 'day', 'Day', 'daily', 'd']:
                    obj_var_pref[obj_name]['Resample Frequency'] = 'D'
                elif obj_var_pref[obj_name]['Resample Frequency'] in ['Once', 'once', 'o']:
                    obj_var_pref[obj_name]['Resample Frequency'] = 'O'
                elif obj_var_pref[obj_name]['Resample Frequency'] in [None]:
                    obj_var_pref[obj_name]['Resample Frequency'] = 'D'  # Assume user is interested in mean of this
                    # variable over Monte Carlo runs.

            # Handle various ways of indicating what statistic is taken within the resampling period. For example,
            # if the resampling is Annual, and the Resample stat is Sum, then your PM is the annual sum,
            # for every year and every realization.
            # Options: Mean, Variance, Median, Max, Min, Sum. Default = Mean.
            obj_var_pref[obj_name]['Resample Stat'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=5).value
            if obj_var_pref[obj_name]['Resample Stat'] not in ['mean', 'variance', 'median', 'min', 'max', 'sum']:
                if obj_var_pref[obj_name]['Resample Stat'] in ['Mean', 'MEAN', 'avg', 'AVG', 'average', 'AVERAGE']:
                    obj_var_pref[obj_name]['Resample Stat'] = 'mean'
                elif obj_var_pref[obj_name]['Resample Stat'] in ['VARIANCE', 'Variance']:
                    obj_var_pref[obj_name]['Resample Stat'] = 'Variance'
                elif obj_var_pref[obj_name]['Resample Stat'] in ['MEDIAN', 'Med']:
                    obj_var_pref[obj_name]['Resample Stat'] = 'median'
                elif obj_var_pref[obj_name]['Resample Stat'] in ['Min', 'MIN', 'minimize', 'MINIMIZE']:
                    obj_var_pref[obj_name]['Resample Stat'] = 'min'
                elif obj_var_pref[obj_name]['Resample Stat'] in ['Max', 'MAX', 'maximize', 'MAXIMIZE']:
                    obj_var_pref[obj_name]['Resample Stat'] = 'max'
                elif obj_var_pref[obj_name]['Resample Stat'] in ['SUM', 'Sum']:
                    obj_var_pref[obj_name]['Resample Stat'] = 'sum'
                elif obj_var_pref[obj_name]['Resample Stat'] in [None]:
                    obj_var_pref[obj_name]['Resample Stat'] = 'mean'  # Default

            obj_var_pref[obj_name]['Time Slice'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=6).value
            if type(obj_var_pref[obj_name]['Time Slice']) not in [str, unicode, None]:
                pass  # Slice is a number. Will determine how to use it in performance.py
            elif type(obj_var_pref[obj_name]['Time Slice']) in [str, unicode]:
                if obj_var_pref[obj_name]['Time Slice'] in ['Last', 'last']:
                    obj_var_pref[obj_name]['Time Slice'] = 1.0  # Slice at simulation duration
                else:
                    # No other strings currently supported. Set = 1 (last value).
                    obj_var_pref[obj_name]['Time Slice'] = 1.0
            else:
                obj_var_pref[obj_name]['Time Slice'] = None  # No value specified, so don't do a time slice at all.

            # Handle various ways of indicating what statistic is taken over all of the samples that are collected
            # given the specified resampling frequency and statistic. Example: 10th percentile of all monthly sums.
            # Options: Mean, Variance, Median, Max, Min, Sum, or ANY QUANTILE FRACTION [0,1]. Default = Mean.
            obj_var_pref[obj_name]['PM Distribution Stat'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=7).value
            if obj_var_pref[obj_name]['PM Distribution Stat'] not in ['mean', 'variance', 'median', 'min', 'max', 'sum']:
                if type(obj_var_pref[obj_name]['PM Distribution Stat']) in [float, int, long]:
                    # User is specifying a quantile or percentile. If > 1 it will be assumed to be a percentage
                    # value, and will be normalized to 0-1.
                    if obj_var_pref[obj_name]['PM Distribution Stat'] <= 0 or obj_var_pref[obj_name]['PM Distribution Stat'] >= 1:
                        # User has not specified a valid fraciton. Use the mean as a default.
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'mean'
                else:
                    if obj_var_pref[obj_name]['PM Distribution Stat'] in ['Mean', 'MEAN', 'avg', 'AVG', 'average',
                                                                           'AVERAGE']:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'mean'
                    elif obj_var_pref[obj_name]['PM Distribution Stat'] in ['VARIANCE', 'Variance']:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'Variance'
                    elif obj_var_pref[obj_name]['PM Distribution Stat'] in ['MEDIAN', 'Med']:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'median'
                    elif obj_var_pref[obj_name]['PM Distribution Stat'] in ['Min', 'MIN', 'minimize', 'MINIMIZE']:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'min'
                    elif obj_var_pref[obj_name]['PM Distribution Stat'] in ['Max', 'MAX', 'maximize', 'MAXIMIZE']:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'max'
                    elif obj_var_pref[obj_name]['PM Distribution Stat'] in ['SUM', 'Sum']:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'sum'
                    elif obj_var_pref[obj_name]['PM Distribution Stat'] in [None]:
                        obj_var_pref[obj_name]['PM Distribution Stat'] = 'mean'  # Default

            obj_var_pref[obj_name]['Epsilon'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=8).value
            if type(obj_var_pref[obj_name]['Epsilon']) not in [int, float, long]:
                obj_var_pref[obj_name]['Epsilon'] = 1  # Default epsilon value == 1

            # Read in system locations for which objective applies.
            obj_var_pref[obj_name]['Locations'] = Input_Data_File['Optimization'].cell(row=obj_row + i,
                                                                                       column=9).value.split(', ')

            obj_var_pref[obj_name]['unit_conv'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=10).value
            if type(obj_var_pref[obj_name]['unit_conv']) not in [int, float, long]:
                obj_var_pref[obj_name]['unit_conv'] = 1  # Default epsilon value == 1

            obj_var_pref[obj_name]['perc_conv'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=11).value
            if (obj_var_pref[obj_name]['perc_conv'] is None) or (type(obj_var_pref[obj_name]['perc_conv']) not in [str, unicode]):
                obj_var_pref[obj_name]['perc_conv'] = 'No'  # Default is a 1 objective optimization problem.
            else:
                if obj_var_pref[obj_name]['perc_conv'] in ['yes', 'Yes', 'y', 'Y']:
                    obj_var_pref[obj_name]['perc_conv'] = 'Yes'
                else:
                    obj_var_pref[obj_name]['perc_conv'] = 'No'

            obj_var_pref[obj_name]['invert'] = Input_Data_File['Optimization'].cell(row = obj_row+i, column=12).value
            if (obj_var_pref[obj_name]['invert'] is None) or (type(obj_var_pref[obj_name]['invert']) not in [str, unicode]):
                obj_var_pref[obj_name]['invert'] = 'No'  # Default is a 1 objective optimization problem.
            else:
                if obj_var_pref[obj_name]['invert'] in ['yes', 'Yes', 'y', 'Y']:
                    obj_var_pref[obj_name]['invert'] = 'Yes'
                else:
                    obj_var_pref[obj_name]['invert'] = 'No'

        # Now store object_order in the obj_var_pref dictionary
        opt_dict = {}  # Stores all data relevant to the optimization that PySedSim will need to know.
        opt_dict['Objective Names Ordered List'] = obj_names_ordered
        # Store epsilons in same order as objective names
        epsilon_list = [obj_var_pref[names]['Epsilon'] for names in opt_dict['Objective Names Ordered List']]
        opt_dict['State Variable Ordered List'] = [obj_var_pref[names]['State Variable'] for names in
                                                                  opt_dict['Objective Names Ordered List']]
        opt_dict['Objective Preferences'] = obj_var_pref
        opt_dict['Num Objectives'] = n_objs
        opt_dict['Num Constraints'] = n_constrs
        opt_dict['Scenario Name'] = simulation_title
    except KeyError:
        raise KeyError("ERROR: Required optimization worksheet 'Optimization' does not exist in input file for "
                       "scenario %s" % simulation_title)

    # Read in general simulation preferences that are relevant to the optimization
    try:
        if Input_Data_File['Simulation Specifications']['B5'].value == 'Stochastic':
            Num_Realizations = Input_Data_File['Simulation Specifications']['B6'].value
        else:
            Num_Realizations = 1
    except KeyError:
        Num_Realizations = 1
    opt_dict['Num Realizations'] = Num_Realizations

    # Define Borg dictionary to return:
    Borg_dict = {
        'opt_dict': opt_dict,
        'link_to_borg': link_to_borg,
        'borg_type': borg_type,
        'num_func_evals': num_func_evals,
        'nSeeds': nSeeds,
        'n_objs': n_objs,
        'n_constrs': n_constrs,
        'maxtime': maxtime,
        'epsilon_list': epsilon_list,
        'optimization approach': opt_approach,
        'ref_set_yes_no': ref_set,
        're_eval': re_eval,
        'parallel_choice': parallel_choice,
        'opt_output_filename': opt_output_filename,
        'ref_set_plot': ref_set_plot,
        'runtime_preferences': runtime_preferences
    }

    return Borg_dict