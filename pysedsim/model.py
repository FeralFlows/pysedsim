"""
@author   Thomas B. Wild
@email:   twild@umd.edu
@Project: pysedsim 1.0.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Purpose: This file houses a top-level PySedSim function that conducts all processes necessary to run a simulation, from data import,
to network and system element creation, to simulation and export of results. Most of the work is done with imported functions. This
file contains no function definitions.

"""

# -*- coding: utf-8 -*-

# Imports
from __future__ import division
from pysedsim.data_processing.data_processing import Import_Simulation_Preferences
from pysedsim.data_processing.data_processing import Determine_Num_Scenarios
from pysedsim.data_processing.data_processing import Monte_Carlo_Import
from pysedsim.data_processing.data_processing import Export_Simulation_Output
from pysedsim.river_basin_elements.system_element_creation import Simulation_Network_Creation
from pysedsim.simulation.pysedsim_main_simulation import SedSim_Main_Simulation
import time  # for keeping track of model execution time
from pysedsim.data_processing.performance import Performance
import logging
import os
import sys

class PySedSim():
	'''

	An Open Source Reservoir and Sediment Simulation Screening Framework for Identifying and Evaluating Dam Siting, Design, and Operation Alternatives.

    '''

	def __init__(self, file_name = 'PySedSim_Input_Specifications.csv', Simulation_mode = 'regular', start_stop =None,
					 parallel_label = None, borg_dict = None, dps_dict = None, decision_vars = None, scenario_name=None,
					 re_eval = None, policy_name=None):
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

		(3) start_stop [OPTIONAL] is a list containing the start and stop integers defining the range of simulations to be
		conducted on a given processor (either through distributed computing over processors on a single computer,
		or through many processors on a cluster). Example of list: [8, 15]. This would mean 7 simulations are to be
		conducted. This 7 simulations may represent a subset of simulations being conducted as part of a monte carlo (
		i.e., this processor is being asked to conduct 7 simulations out of 1000 simulations).

		(4) parallel_label is a label to append to the end of any output files produced from a given simulation or batch
		of simulations, to reflect the output produced by different processors. For example, in a parallellized monte
		carlo simulation, this could be the rank of the processor on which the pysedsim_wrapper calling this function is
		running, so that the distributed output from the processors can be aggregated into a single file.

		(5) Borg_dict is an optimization dictionary. This is defined automatically in the optimization_pysedsim.py module.

		(6) decision_vars (optional) is a list of decision variable values specified by an external optimization model (
		e.g., Borg) in the form of a numpy array. Currently, the only thing this input can be used for is to populate the
		parameters for a RBF-type reservoir operating policy.

		(7) scenario_name: Optional. List containing one string that represents the name of the current scenario (
		internal to optimization or post-optimization) from which this pysedsim function is being called for purposes of
		policy reevaluation. Example: ['Scenario 1']. Ensures that if pysedsim is called from within an optimzation loop (
		over multiple scenarios) for purposes of re-evaluating policies from a reference set, only the current
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

		Example function call:
		(1) PySedSim(file_name = 'PySedSim_Input_File.csv', simulation_mode = 'debug')
		(2) PySedSim()
		"""

		self.file_name = file_name
		self.Simulation_mode = Simulation_mode
		self.start_stop = start_stop
		self.parallel_label = parallel_label
		self.borg_dict = borg_dict
		self.dps_dict = dps_dict
		self.decision_vars = decision_vars
		self.scenario_name = scenario_name
		self.re_eval = re_eval
		self.policy_name = policy_name


	@staticmethod
	def make_dir(pth):
		"""Create dir if not exists."""
		if not os.path.exists(pth):
			os.makedirs(pth)


	def execute(self):
		# Import entire file "PySedSim_Input_Specification.csv", then load specific information to variables
		[num_scenarios, simulation_titles_list, imported_specs, main_input_files_dir, main_output_file_dir,
		 monte_carlo_file_list, external_mc_data] = Determine_Num_Scenarios(self.file_name)

		# create output directory
		self.make_dir(main_output_file_dir)		 
		 
		# Initialize logging file
		self.init_log(main_output_file_dir)

		if self.borg_dict is not None:
			# An optimization is being performed. Scenarios can be looped through in optimization_pysedsim(), but not here given optimization
			# is being performed.
			opt_dict = self.borg_dict['opt_dict']
			num_scenarios = 1
			simulation_titles_list = [opt_dict['Scenario Name']]
			export_data = 'No'  # Do not dump output to .csv file in optimization setting.
			export_data_scenario = 'No'  # Ensures export_data is 'No' for all scenarios being looped through.
		else:
			# Optimization will not take place.
			opt_dict = None
			if self.scenario_name is not None:
				num_scenarios = 1
				simulation_titles_list = self.scenario_name
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

		for i in range(num_scenarios):
			# Initialize logging file
			logging.info('Beginning simulation scenario: ' + simulation_titles_list[i])
			# In case looping through scenarios, re-initiate start_stop
			if parallelize is None:
				self.start_stop = None
			if export_data_scenario is None:
				export_data = None
			start_time = time.time()  # Keep track of PySedSim model execution time
			simulation_title = simulation_titles_list[i]  # Feed in a string as the simulation title
			monte_carlo_file = monte_carlo_file_list[i]
			# Import simulation preferences from user-provided input file.
			[export_data, export_file_type, var_sub_list, element_export_list,
			 Input_Data_File, T, Sim_Dur, Stochastic_Sim, Num_Realizations, simulation_dates,
			 simulation_dates_no_leap, Col_Names, Monte_Carlo_Parameters_File_Name, external_mc_data, simulation_title,
			 start_stop, Sampled_Parameter_Dict, Synthetic_Inflow_dataframe_name_LIST,
			 Synthetic_Inflows_dictionary] = Import_Simulation_Preferences(imported_specs, simulation_title,
																		   main_input_files_dir, re_eval=self.re_eval,
																		   start_stop=self.start_stop,
																		   export = export_data,
																		   Monte_Carlo_Parameters_File_Name=monte_carlo_file,
																		   external_mc_data = external_mc_data)

			# Create an initial network of reservoirs, channels, juntions, etc. that will be simulated.
			[SystemObjects, Time_Series_Output_Dictionary, Output_Object_Dict, Element_Dict, Flushing_Group_Dict, element_stochastic_components,
			 Parameter_Input_Dictionary] = Simulation_Network_Creation(Input_Data_File)

			if Stochastic_Sim == 1:
				[Parameter_Input_Dictionary, element_stochastic_components, Sampled_Parameter_Dict, Synthetic_Inflow_dataframe_name_LIST,
				 Synthetic_Inflows_dictionary] = Monte_Carlo_Import(main_input_files_dir, Col_Names, element_stochastic_components,
																	SystemObjects["Ordered Simulation List"], external_mc_data,
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
																										   Synthetic_Inflows_dictionary, op_policy_params=self.dps_dict)

			print("--- Simulation(s) Complete in %s seconds ---" % (time.time() - start_time))
			logging.info("--- Simulation(s) Complete in %s seconds ---")

			# Export output data if user indicates this should be conducted, or if pysedsim is being called as part of a
			# re-evaluation process, in which case output is (assumed to be) desired to be stored.
			if export_data == 'Yes' or self.scenario_name is not None:
				Export_Simulation_Output(export_file_type, Time_Series_Output_Dictionary, state_list_excel, main_output_file_dir,
										 simulation_title, var_sub_list, rank = self.parallel_label,
										 policy_name=self.policy_name)

		self.cleanup()

		# If PySedSim() is connected to an external optimization model, evaluate performance over the realizations.
		if opt_dict is not None:
			objs = Performance(Time_Series_Output_Dictionary, Sim_Dur, opt_dict=opt_dict, sim_title = simulation_title)
			return objs


	def init_log(self, outdir):
		"""
        Initialize project-wide logger.
        The logger outputs to both stdout and a file
        """
		log_format = logging.Formatter('%(levelname)s: %(message)s')
		log_level = logging.DEBUG
		log_file = os.path.join(outdir, 'logfile.log')

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
		logging.info("End of simulation")

		# Remove logging handlers - they are initialized at the module level, so this prevents duplicate logs from
		# being created if pysedsim is run multiple times.
		logger = logging.getLogger()
		for handler in logger.handlers[:]:
			logger.removeHandler(handler)