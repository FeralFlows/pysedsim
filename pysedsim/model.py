"""
@author   Thomas B. Wild
@email:   twild@umd.edu
@Project: pysedsim 1.0.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Purpose: This file houses a top-level PySedSim function that conducts all processes necessary to run a simulation, from
data import to network and system element creation, to simulation and export of results. Most of the work is done
with imported functions. This file contains no function definitions.

"""

# -*- coding: utf-8 -*-

# Imports
from __future__ import division
from pysedsim.data_processing.data_processing import *
from pysedsim.river_basin_elements.system_element_creation import Simulation_Network_Creation
from pysedsim.simulation.pysedsim_main_simulation import SedSim_Main_Simulation
from pysedsim.data_processing.performance import Performance
from pysedsim.optimization.direct_policy_search import *
import logging
import os
import sys
import numpy as np
import time  # for keeping track of model execution time

class PySedSim():
	'''

	An Open Source Reservoir and Sediment Simulation Screening Framework for Identifying and Evaluating Dam Siting,
	Design, and Operation Alternatives.

    '''

	def __init__(self, file_name = 'PySedSim_Input_Specifications.csv', Simulation_mode = 'regular', start_stop =None,
					 parallel_label = None, scenario_name=None, re_eval = None, policy_name=None):
		"""
		Initialize pysedsim

		Purpose: Conducts all processes necessary to run a simulation, including basic data import, monte carlo data
		import, network and system element creation, simulation, and export of results. All of the work in this function
		is done by imported functions/methods.

		Note: The user is required to provide a .csv file titled "PySedSim_Input_Specifications" that contains
		assumptions related to input and output file paths, simulation title, and Monte Carlo parameter generation
		preferences. This file must be located in the same directory as this top-level file.

		Inputs:
		(1) file_name [OPTIONAL] is the name of the input file, which must be a .csv file (as string) located in the same
		directory as the file that contains this pysedsim.py file. (Default: 'PySedSim_Input_Specifications.csv').

		(2) simulation_mode [OPTIONAL] (Options: 'regular' or 'debug'. 'debug' stores all variable output temporarily and
		uses more RAM. Default: 'regular'.).

		(3) start_stop [OPTIONAL] is a list containing the start and stop integers defining the range of simulations to
		be conducted on a given processor (either through distributed computing over processors on a single computer,
		or through many processors on a cluster). Example of list: [8, 15]. This would mean 7 simulations are to be
		conducted. This 7 simulations may represent a subset of simulations being conducted as part of a monte carlo (
		i.e., this processor is being asked to conduct 7 simulations out of 1000 simulations).

		(4) parallel_label is a label to append to the end of any output files produced from a given simulation or batch
		of simulations, to reflect the output produced by different processors. For example, in a parallellized monte
		carlo simulation, this could be the rank of the processor on which the pysedsim_wrapper calling this function is
		running, so that the distributed output from the processors can be aggregated into a single file.

		(7) scenario_name: Optional. List containing one string that represents the name of the current scenario (
		internal to optimization or post-optimization) from which this pysedsim function is being called for purposes of
		policy reevaluation. Example: ['Scenario 1']. Ensures that if pysedsim is called from within an optimzation loop
		(over multiple scenarios) for purposes of re-evaluating policies from a reference set, only the current
		optimization scenario gets re-evalutad here, rather than all scenarios listed in the input file.

		(8) re_eval: List of preferences related to reevaluation of operating policies from the reference set,
		such as what variables/locations should be exported, whether simulation is deterministic/stochastic,
		and if stochastic how many realizations.

		(9) policy_name: Name of the policy being evaluated, if pysedsim is being run in a reevaluation setting (from
		processing_reference_set.Policy_Reevaluation()). Policy number is used as the name, and sub-folder is created.

		Outputs:
		1. Data are written to files and saved
		2. If user has connected PySedSim to an external optimization model, wherein PySedSim() is the "function
		evaluation", then PySedSim returns objs, the objective function values.

		Example instantiation:
		(1) PySedSim(file_name = 'PySedSim_Input_File.csv', simulation_mode = 'debug')
		(2) PySedSim()
		"""

		self.start_time = time.time()  # Keep track of PySedSim model execution time
		self.file_name = file_name
		self.Simulation_mode = Simulation_mode
		self.start_stop = start_stop
		self.parallel_label = parallel_label
		self.scenario_name = scenario_name
		self.re_eval = re_eval
		self.policy_name = policy_name
		# Import entire file "PySedSim_Input_Specification.csv", then load specific information to variables
		[self.num_scenarios, self.simulation_titles_list, self.imported_specs, self.main_input_files_dir,
		 self.main_output_file_dir, self.monte_carlo_file_list, self.external_mc_data] = Determine_Num_Scenarios(
			self.file_name)
		# create output directory
		self.make_dir(self.main_output_file_dir)
		# Initialize logging file
		self.init_log(self.main_output_file_dir)
		self.os_fold = Op_Sys_Folder_Operator()		


	@staticmethod
	def make_dir(pth):
		"""Create dir if not exists."""
		if not os.path.exists(pth):
			os.makedirs(pth)


	def execute_simulation(self, borg_dict=None, dps_dict=None, decision_vars=None, sim_opt_log_detail = 0):
		'''

		:param borg_dict: an optimization dictionary. Defined automatically in the optimization_pysedsim.py module.
		:param decision_vars: decision_vars (optional) is a list of decision variable values specified by an external
		optimization model (e.g., Borg) in the form of a numpy array. Currently, the only thing this input can be used
		for is to populate the parameters for a RBF-type reservoir operating policy.
		:param dps_dict: dictionary defining direct policy search preferences
		:param sim_opt_log_detail: level of detail in simulation logging during simulation=optimization run. 0=no
		output, 1 = all output (for every function evaluation/simulation).
		:return:
		'''

		if borg_dict is None:
			logging.info("Executing a simulation experiment")
		else:
			if sim_opt_log_detail == 1:
				logging.info("Executing a simulation embedded in a broader simulation-optimization experiment")

		self.borg_dict = borg_dict
		self.dps_dict = dps_dict
		self.decision_vars = decision_vars

		if self.borg_dict is not None:
			self.optimization = 1  # An optimization is being performed.
			# is being performed.
			opt_dict = self.borg_dict['opt_dict']
			self.num_scenarios = 1
			self.simulation_titles_list = [opt_dict['Scenario Name']]
			export_data = 'No'  # Do not dump output to .csv file in optimization setting.
			export_data_scenario = 'No'  # Ensures export_data is 'No' for all scenarios being looped through.
		else:
			self.optimization = 0  # Optimization will not take place.
			opt_dict = None
			if self.scenario_name is not None:
				self.num_scenarios = 1
				self.simulation_titles_list = self.scenario_name
			export_data_scenario = None  # Data export preferences for particular scenario in chain, to prevent carryover.
		if self.dps_dict is not None:
			self.dps_dict['Policy Function'] = {'Raw Parameters': self.decision_vars[0:self.dps_dict['n_vars']]}
			try:
				self.dps_dict['Flushing Optimization']['Flushing Parameters'] = self.decision_vars[self.dps_dict['n_vars']:self.dps_dict['total_vars']]
			except KeyError:
				pass

		if self.start_stop is not None:
			parallelize = 1  # Processing in parallel. Has implications for how to read in data if stochastic simulation.
		else:
			parallelize = None

		for i in range(self.num_scenarios):
			# Initialize logging file
			if (self.optimization == 0) or (sim_opt_log_detail == 1):
				logging.info('Beginning scenario: ' + self.simulation_titles_list[i])
			# In case looping through scenarios, re-initiate start_stop
			if parallelize is None:
				self.start_stop = None
			if export_data_scenario is None:
				export_data = None
			if self.optimization == 0:
				self.scen_start_time = time.time()  # Keep track of PySedSim model execution time
			simulation_title = self.simulation_titles_list[i]  # Feed in a string as the simulation title
			monte_carlo_file = self.monte_carlo_file_list[i]
			# Import simulation preferences from user-provided input file.
			[export_data, export_file_type, var_sub_list, element_export_list,
			 Input_Data_File, T, Sim_Dur, Stochastic_Sim, Num_Realizations, simulation_dates,
			 simulation_dates_no_leap, Col_Names, Monte_Carlo_Parameters_File_Name, self.external_mc_data, simulation_title,
			 start_stop, Sampled_Parameter_Dict, Synthetic_Inflow_dataframe_name_LIST,
			 Synthetic_Inflows_dictionary] = Import_Simulation_Preferences(self.imported_specs, simulation_title,
																		   self.main_input_files_dir, re_eval=self.re_eval,
																		   start_stop=self.start_stop,
																		   export = export_data,
																		   Monte_Carlo_Parameters_File_Name=monte_carlo_file,
																		   external_mc_data = self.external_mc_data)

			# Create an initial network of reservoirs, channels, juntions, etc. that will be simulated.
			[SystemObjects, Time_Series_Output_Dictionary, Output_Object_Dict, Element_Dict, Flushing_Group_Dict, element_stochastic_components,
			 Parameter_Input_Dictionary] = Simulation_Network_Creation(Input_Data_File)

			if Stochastic_Sim == 1:
				[Parameter_Input_Dictionary, element_stochastic_components, Sampled_Parameter_Dict, Synthetic_Inflow_dataframe_name_LIST,
				 Synthetic_Inflows_dictionary] = Monte_Carlo_Import(self.main_input_files_dir, Col_Names, element_stochastic_components,
																	SystemObjects["Ordered Simulation List"], self.external_mc_data,
																	Monte_Carlo_Parameters_File_Name, simulation_dates,
																	simulation_dates_no_leap, Num_Realizations, Sampled_Parameter_Dict,
																	Synthetic_Inflow_dataframe_name_LIST, Synthetic_Inflows_dictionary,
																	start_stop=start_stop, parallelize_import=parallelize)

			# Run main simulation
			[state_list_excel, Time_Series_Output_Dictionary, Output_Object_Dict] = SedSim_Main_Simulation(Num_Realizations, T, Input_Data_File,
																										   element_stochastic_components,
																										   SystemObjects, Element_Dict,
																										   Flushing_Group_Dict, Stochastic_Sim,
																										   Parameter_Input_Dictionary,
																										   simulation_dates_no_leap,
																										   Time_Series_Output_Dictionary,
																										   Sim_Dur,
																										   simulation_dates,
																										   Col_Names, self.Simulation_mode,
																										   Output_Object_Dict, var_sub_list, element_export_list,
																										   Sampled_Parameter_Dict,
																										   Synthetic_Inflow_dataframe_name_LIST,
																										   Synthetic_Inflows_dictionary,
																										   op_policy_params=self.dps_dict,
																										   optimization = self.optimization)

			# Export output data if user indicates this should be conducted, or if pysedsim is being called as part of a
			# re-evaluation process, in which case output is (assumed to be) desired to be stored.
			if export_data == 'Yes' or self.scenario_name is not None:
				Export_Simulation_Output(export_file_type, Time_Series_Output_Dictionary, state_list_excel,
										 self.main_output_file_dir, simulation_title, var_sub_list, rank =
										 self.parallel_label, policy_name=self.policy_name)

			if self.optimization == 0:
				logging.info("--- Simulation Experiment: {1} Complete in {0} minutes ---".format(
					(time.time() - self.scen_start_time) / 60, self.simulation_titles_list[i]))

		# If PySedSim() is connected to an external optimization model, evaluate performance over the realizations.
		if opt_dict is not None:
			objs = Performance(Time_Series_Output_Dictionary, Sim_Dur, opt_dict=opt_dict, sim_title=simulation_title,
							   optimization=self.optimization)
			return objs

		if opt_dict is None:
			self.cleanup()  # uncoupled simulation is over; close out log file

	def PySedSim_Caller(self, vars, additional_inputs):
		'''
        Purpose: Calls the PySedSim simulation model. Designed for use in an optimization setting. Designed to be called by
        Borg optimization model wrapper in each function evaluation as part of a simulation-optimization experiment.

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
			np.savetxt('RBF_Parameters.txt', borg_vars, newline=' ')  # Save policy so PySedSim can import it.
			performance = self.execute_simulation(borg_dict=Borg_dict)
		else:
			# Use for stochastic simulation/optimization, particularly if being run in parallel.
			op_policy_params = np.asarray(borg_vars)  # Cast borg output parameters as array for use in simulation.
			performance = self.execute_simulation(borg_dict=Borg_dict, dps_dict=DPS_Dict,
												  decision_vars=op_policy_params)
		return performance


	def execute_optimization(self, file_name='PySedSim_Input_Specifications.csv', borg_path = None,
							 sim_opt_log_detail = 0):

		'''

        Purpose: Perform main Borg optimization loop through scenarios and seeds. Calls PySedSim for simulation.

        Args:
            file_name: Name of the top-level input file that indicates where input files are located, etc. Default is
            'PySedSim_Input_Specifications.csv'. Must be the same file used for the simulation that is to be connected to
            the optimization.
            :param sim_opt_log_detail: level of detail to display in "log.txt" during execution. = minimal detail,
            1=maximum detail. Note: Detailed output is always printed in the .runtime file for each
            optimization seed.

        Returns:
            Produces pareto outputs for each seed, exported to separate files in a "sets" folder created in the model
            output location specified in the input file.
        '''

		logging.info("Executing a coupled simulation-optimization experiment")

		# Import borg wrapper given optimization is to be perfored
		from pysedsim.optimization import borg as bg

		# Loop through as many optimization scenarios as user has specified
		for j in range(self.num_scenarios):
			self.scen_start_time = time.time()
			simulation_title = self.simulation_titles_list[j]
			logging.info("Simulation-optimization experiment name: {0}".format(simulation_title))
			Borg_dict = Import_Optimization_Preferences(simulation_title, self.imported_specs, self.main_input_files_dir)
			if Borg_dict['optimization approach'] == 'DPS':
				DPS_dict = Import_DPS_Preferences(simulation_title=simulation_title,
												  imported_specs=self.imported_specs,
												  main_input_files_dir=self.main_input_files_dir)
				Borg_dict['n_vars'] = DPS_dict[
					'n_vars']  # Store num decision variables in Borg_dict, copied from DPS_dict.

			# Populate DPS input variable ranges to specify for borg
			borg_variable_range = []
			# Direct Policy Search Radial Basis Function Parameters:
			for r in range(DPS_dict['RBF Specs']['Number']):
				for m in range(DPS_dict['num_inputs']):
					borg_variable_range.append([-1, 1])  # for normalized center
					borg_variable_range.append([0, 1])  # for normalized radius
				for k in range(DPS_dict['num_reservoirs']):
					borg_variable_range.append(
						[0, 1])  # for weights, which are normalized later on as they must sum to 1.

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
			output_location = self.main_output_file_dir + self.os_fold + simulation_title + self.os_fold + 'sets'
			if not os.path.exists(output_location):
				os.makedirs(output_location)

			# Initialize Borg configuration, for purposes of specifying borg a path to Borg DLL/SO
			bg.Configuration.initialize(borg_path=borg_path)  # 

			# Set up interface with Borg for parallel function evaluation, if Borg Maser-Slave is being used, and if this
			#  is the first optimization scenario being considered. Can't start MPI twice.
			if Borg_dict['borg_type'] == 1 and j == 0:
				bg.Configuration.startMPI()  # start parallelization with MPI

			# Set up interface with Borg for optimization
			for sd in range(Borg_dict['nSeeds']):
				logging.info("Running random seed: {0}".format(sd))
				self.seed_start_time = time.time()
				# Create instance of Borg class
				borg = bg.Borg(Borg_dict['total_vars'], Borg_dict['n_objs'], Borg_dict['n_constrs'], self.PySedSim_Caller,
							   add_pysedsim_inputs=[Borg_dict, DPS_dict, file_name])
				borg.setBounds(*borg_variable_range)  # Set decision variable bounds
				borg.setEpsilons(*Borg_dict['epsilon_list'])  # Set epsilon values
				runtime_filename = self.main_output_file_dir + self.os_fold + simulation_title + self.os_fold + 'sets' + self.os_fold + \
								   'runtime_file_seed_' + str(sd + 1) + '.runtime'
				if Borg_dict['borg_type'] == 0:
					# Run serial Borg
					if Borg_dict['runtime_preferences']['runtime_choice'] == 'Yes':
						result = borg.solve({"maxEvaluations": Borg_dict['num_func_evals'], "runtimeformat": 'borg',
											 "frequency": Borg_dict['runtime_preferences']['runtime_freq'],
											 "runtimefile": runtime_filename})
					else:
						result = borg.solve({"maxEvaluations": Borg_dict['num_func_evals']})

				# result = borg.solve({"maxEvaluations":Borg_dict['num_func_evals']})  # Minimum is 100 in Borg
				elif Borg_dict['borg_type'] == 1:
					# Run parallel Borg
					if Borg_dict['runtime_preferences']['runtime_choice'] == 'Yes':
						result = borg.solveMPI(maxEvaluations=Borg_dict['num_func_evals'], runtime=runtime_filename,
											   frequency=Borg_dict['runtime_preferences']['runtime_freq'])
					else:
						result = borg.solveMPI(maxEvaluations=Borg_dict['num_func_evals'])
				if result:
					# This particular seed is now finished being run in parallel. The result will only be returned from
					# one node in case running Master-Slave Borg.
					if sim_opt_log_detail == 1:
						result.display()

					# Create/write objective values and decision variable values to files in folder "sets", 1 file per seed.
					f = open(output_location + self.os_fold + 'Borg_DPS_PySedSim' + str(sd + 1) + '.set', 'w')
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

						f.write(line[0:-1] + '\n')
					f.write("#")
					f.close()

					# Create/write objective values only to files in folder "sets", 1 file per seed. Purpose is so that
					# the file can be processed in MOEAFramework, where performance metrics may be evaluated across seeds.
					f2 = open(output_location + self.os_fold + 'Borg_DPS_PySedSim_no_vars' + str(sd + 1) + '.set', 'w')
					for solution in result:
						line = ''
						for i in range(len(solution.getObjectives())):
							line = line + (str(solution.getObjectives()[i])) + ' '

						f2.write(line[0:-1] + '\n')
					f2.write("#")
					f2.close()

					logging.info(
						"Seed {0} complete in {1} minutes".format(sd, ((time.time() - self.seed_start_time) / 60)))

			# All seeds for this scenario are complete. Now generate reference set if user desires.
			if Borg_dict['ref_set_yes_no'] == 'Yes':
				# Ensure that this reference set processing is done only once in event this is being run in parallel.
				# if platform.system() == 'Linux':
				# [start, stop, rank] = parallelize_pysedsim(1)  # Distribute re-evaluation work to processors.
				# if rank == 0:
				[Borg_dict, DPS_dict] = Reference_Set(inline_optimization=[simulation_title, Borg_dict, DPS_dict,
																		   self.main_output_file_dir])
			# else:
			# [Borg_dict, DPS_dict] = Reference_Set(inline_optimization=[simulation_title, Borg_dict, DPS_dict,main_output_file_dir])

			# Re-evaluate (simulate) all of the policies from the reference set.
			if Borg_dict['re_eval'] == 'Yes':
				if Borg_dict['ref_set_yes_no'] == 'Yes':
					ref_set_file_name = output_location + self.os_fold + Borg_dict['opt_output_filename']
					[ref_set_array, objective_values, dec_var_values] = Policy_Reevaluation(
						ref_set_file_name=ref_set_file_name, post_optimization='No', file_name=file_name,
						dps_dict=DPS_dict,
						borg_dict=Borg_dict, parallel=Borg_dict['parallel_choice'],
						internal_opt_call=[Borg_dict['opt_dict']['Scenario Name']])
				else:
					# Do not specify a reference set file name, as user did not generate a reference set with
					# optimization. Instead, just use the default reference set file name from the Policy_Reevaluation
					# method.
					[ref_set_array, objective_values, dec_var_values, dps_dict, borg_dict] = Policy_Reevaluation(
						post_optimization='No',
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
				objs_to_plot = Borg_dict['opt_dict'][
					'Objective Names Ordered List']  # Plot all objectives in optimization.
				if Borg_dict['re_eval'] == 'Yes':
					# Reference set has already been processed. Store those values in the reference set preferences dict.
					ref_set_pref_dict['ref_set_array'] = ref_set_array
					ref_set_pref_dict['objective_values'] = objective_values
					ref_set_pref_dict['dec_var_values'] = dec_var_values
					# Produce/save plots.
					Scatter_Plot_Objectives(ref_set_pref_dict, objs_to_plot, plot_dict=plot_dict,
											save_fig=output_location,
											movie_dict=movie_dict)
				else:
					ref_set_pref_dict['ref_set_file_name'] = output_location + self.os_fold + Borg_dict[
						'opt_output_filename']
					ref_set_pref_dict['num_dec_vars'] = DPS_dict['n_vars']
					ref_set_pref_dict['num_objs'] = Borg_dict['n_objs']
					# Produce/save plots.
					[ref_set_array, objective_values, dec_var_values] = Scatter_Plot_Objectives(ref_set_pref_dict,
																								objs_to_plot,
																								plot_dict=plot_dict,
																								save_fig=output_location,
																								movie_dict=movie_dict)
			# Optimization Scenario complete.
			logging.info("--- Simulation-Optimization Experiment: {1} Complete in {0} minutes ---".format((
				(time.time() - self.scen_start_time) / 60), simulation_title))

		if Borg_dict['borg_type'] == 1:
			bg.Configuration.stopMPI()  # stop parallel function evaluation process
		self.cleanup()  # close out log file
		return


	def init_log(self, outdir):
		"""
        Initialize project-wide logger.
        The logger outputs to both stdout and a file
        """
		log_format = logging.Formatter('%(levelname)s: %(message)s')
		log_level = logging.DEBUG
		log_file = os.path.join(outdir, 'logfile.log')
		# Delete previous logs; avoids merging outputs from multiple sessions  into a single log
		if os.path.exists(log_file):
			os.remove(log_file)

		logger = logging.getLogger()
		logger.setLevel(log_level)

		# logger console handler
		c_handler = logging.StreamHandler(sys.stdout)
		c_handler.setLevel(log_level)
		c_handler.setFormatter(log_format)
		logger.addHandler(c_handler)

		# logger file handler
		f_handler = logging.FileHandler(log_file)
		c_handler.setLevel(log_level)
		c_handler.setFormatter(log_format)
		logger.addHandler(f_handler)


	def cleanup(self):
		"""Close log files."""
		if self.optimization == 0:
			logging.info("--- All Simulation Experiments Complete in {0} minutes ---".format((time.time() -
																						 self.start_time)/60))
		else:
			logging.info("--- All Simulation-Optimization Experiments Complete in {0} minutes ---".format((time.time() -
																						 self.start_time)/60))


		# Remove logging handlers - they are initialized at the module level, so this prevents duplicate logs from
		# being created if pysedsim is run multiple times.
		logger = logging.getLogger()
		for handler in logger.handlers[:]:
			logger.removeHandler(handler)