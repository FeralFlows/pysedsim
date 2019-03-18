# -*- coding: utf-8 -*-

'''

This module defines the Reservoir Class and methods.

Reservoirs are large pools of stored water from which the releases of water and sediment are controlled with various
outlet structures.

'''

# import relevant modules
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0)
import numpy as np
from pysedsim.river_basin_elements.storage_element import Storage_Element
from outlet import Outlet
from pysedsim.data_processing.data_processing import *
from pysedsim.data_processing.matrix_interpolation import Matrix_Interpolation
from pysedsim.sediment_management.flushing import Flushing
from pysedsim.sediment_management.sluicing import Sluicing
from pysedsim.sediment_management.dredging import Dredging
from pysedsim.sediment_management.density_current_venting import Density_Current_Venting
# from bypassing import Bypassing
from datetime import datetime  # Used to work with date objects (especially date arithmetic)
from datetime import timedelta  # Used to add days/months/years to a datetime object
import pysedsim.optimization.direct_policy_search
import logging

class Reservoir(Storage_Element):
    def __init__(self, name, T, Input_Data_File, Element_Sub_Dict, stochastic_components=None, op_policy_params=None):
        # Note: The order in which method calls appear below is important. They should not be rearranged.
        if hasattr(Storage_Element, '__init__'):
            Storage_Element.__init__(self, name, T, Element_Sub_Dict)    # Call parent constructor.
        Reservoir.Array_Initialization(self, T)  # Initialize arrays reservoirs will have.
        self.Orifices = Outlet(self.name, T, Input_Data_File)  # Every reservoir must have outlets (orifices) of some sort
        Reservoir.Import_Data(self, T, Input_Data_File, stochastic_components,
                              self.Element_Sub_Dict["Simulation Start Date"], op_policy_params)
        Storage_Element.Time_Zero_Initialization(self)

    def Master_Method_Caller(self, t, Flow_in, Sed_in, flushing_dict):
        # This method calls all other relevant methods in this class in the correct order.
        Storage_Element.Master_Method_Caller(self, t)  # Increments date
        Reservoir.Import_External_Element_State(self, t)  # Imports relevant real-time data from other elements.
        Storage_Element.Element_Inflows(self, t, Flow_in, Sed_in)
        self.Flushing_Dictionary = flushing_dict  # Update dictionary so you know exactly where flushing is happening in system.
        if self.Reservoir_Operations_Goal_ID == 5:
            Reservoir.DPS_Set_WSE_Target(self, t)  # Set daily operations target if op. policy externally specified.
        Reservoir.Sediment_Management_Daily_Import(self, t)
        Reservoir.Trapped_Load(self, t)
        Reservoir.Reservoir_Volume_Reduction(self, t, self.E_Sed[1])
        # If flushing can occur, execute flushing related method:
        try:
            self.Operating_Policy["Flushing"].Available_Flushing_Load(t, self.Settled_volume[t],
                                                                      self.capacity_active_reservoir[t + 1],
                                                                      self.capacity_dead_reservoir[t + 1],
                                                                      self.capacity_active_reservoir[0],
                                                                      self.capacity_dead_reservoir[0],
                                                                      self.Original_River_Bed_Elevation)
        except KeyError:
            pass
        Reservoir.Elev_Stor_Target(self, t)
        # If water level needs to be capped due to upstream uncontrolled outlet, execute related method:
        if self.Cap_Res_WL == 1:
            Reservoir.Terminal_WL_Calc(self,t)
        self.Evap[t] = self.Orifices.Outlet_Euler_Capacity_Estimation(self.name, t, self.current_date,
                                                                      self.Monthly_Evap_Data, self.S[t],
                                                                      self.water_surface_elevation[t],
                                                                      self.water_surface_area[t],
                                                                      self.Environmental_Release_Goal[t],
                                                                      self.flushing_today,
                                                                      self.flushing_group_drawn_down_today,
                                                                      self.sluicing_today,
                                                                      self.density_current_venting_today,
                                                                      self.Reservoir_Type, self.Settled_volume[t],
                                                                      self.Q_in[t], self.E_V_A)
        Reservoir.Set_total_release_goal(self, t)
        Reservoir.Water_Storage_Mass_Balance(self, t)
        Reservoir.Distribute_Outflow_Among_Outlets(self, t)
        Reservoir.Outflow_Storage_Check(self, t)
        Reservoir.SedMgmt_Main(self, t)
        Reservoir.Sediment_Mass_Balance(self, t)
        Reservoir.Hydropower_Calculations(self, t)
        Storage_Element.Outflow_Element_Allocation(self, t)  # This is a reservoir, so flow will simply be allocated to single downstream junction.
        Storage_Element.Post_Mass_Balance_Calcs(self, t)
        Reservoir.Element_Data_Transfer(self, t)
        Reservoir.Performance_Measure_Pre_Calc(self, t)

    def Array_Initialization(self, T):
        # Initialize Arrays, dictionaries and variables that are required for simulation. Many other arrays are
        # initialized in the Import_Data module, as they only need to be conditionally created if the user is
        # simulating particular processes.
        self.water_surface_elevation = np.zeros(T+1)
        self.water_surface_area = np.zeros(T+1)
        self.S_a = np.zeros(T+1)
        self.S_d = np.zeros(T+1)
        self.capacity_active_reservoir = np.zeros(T+1)
        self.capacity_dead_reservoir = np.zeros(T+1)
        self.capacity_total_reservoir = np.zeros(T+1)
        self.Evap = np.zeros(T)
        self.TE_avg = np.zeros(T)
        self.Residence_Time = np.zeros(T)  # Average residence time of water in the reservoir (days). Running average
        #  can be taken, as given by Brune sediment trapping preferences.
        self.STORAGE_TARGET_elevation = np.zeros(T+1)
        self.ELEVATION_TARGET = np.zeros(T+1)
        self.STORAGE_TARGET = np.zeros(T+1)
        self.reservoir_is_full = 0  # Variable denoting whether reservoir is currently full of sediment (= 1) or not (=0, default)
        self.Settled_mass = np.zeros(T)
        self.Settled_volume = np.zeros(T)
        self.Operating_Policy = {}  # Dictionary that stores different operating policy components, including Sediment management components.
        self.Flushing_Dictionary = 0  # This will be populated during reservoir creation if Flushing groups of any type exist. DO NOT move this into a flushing module or conditionally import it, it must be kept here.
        self.Environmental_Release_Goal = np.zeros(T)
        self.Storage_Target_Release_Goal = np.zeros(T)
        self.S_delta = np.zeros(T)
        self.Dam_Tailwater_Elevation = np.zeros(T+1)  # stores tailwater elevation for reservoir.
        self.DCV_Flow = 0 # For use in preventing overestimate of sediment released during DCV. Must exist regardless of whether DCV exists.
        self.Hydropower_avg_MW = np.zeros(T)  # Tracks daily power production. Is initialized/used regardless of whether dam produces power.
        self.Hydropower_avg_MWH = np.zeros(T)  # Tracks daily energy production. Is initialized/used regardless of whether dam produces power.
        self.Curtailment_limited = np.zeros(T)  # Energy not produced because powerhouse operating at capacity
        self.Curtailment_unlimited = np.zeros(T)  # Energy not produced because powerhouse operating at capacity
        self.SS_C_adjusted_for_evap = np.zeros(T)  # Reservoir suspended sediment concentration (adjusted for evaporation)
        self.flushing_today = 0  # Initialize to zero. Equals reservoir's .Current value for sediment management type.
        self.res_flushed_load = np.zeros(T)  # daily sediment load released during a flushing event
        self.res_flushed_load[:] = 'NaN'  # Set all values to NaN initially, so zeros dont appear on non-flushing days.
        self.flushing_yesterday = 0  # Initialize to zero. Equals reservoir's yesterday .Current value for sediment management type.
        self.sluicing_today = 0  # Initialize to zero. Equals reservoir's .Current value for sediment management type.
        self.sluicing_yesterday = 0  # Initialize to zero. Equals reservoir's yesterday .Current value for sediment management type.
        self.density_current_venting_today = 0  # Initialize to zero. Equals reservoir's .Current value for sediment management type.
        self.density_current_venting_yesterday = 0  # Initialize to zero. Equals reservoir's yesterday .Current value for sediment management type.
        #self.E_V_A_TRACK = [[] for varx in range(T+1)]  # This tracks the evolution of the E_V_A double list over
        # time. It's a 3-D array.
        self.stoch_trap_adjust = 0  # Adds noise to Churchll trapping (by adding +/- % to % passing value if stochastic simulation). Initialize to zero, reset externally if it has a value.
        self.theoretical_peaking_capacity = np.zeros(T)  # Number of hours for which reservoir could theoretically operate at powerhouse design flow.
        self.Q_Downstream_Junction = np.zeros(T)  # Flow at junction downstream of dam to calculate tailwater elev.
        self.re_calc_energy = 0  # Re-calc energy after other elements simulated in pysedsim_main_simulation.py? Default is = 0 (No) outlet
        self.fraction_Q_jct = np.zeros(T)  # Fraction of flow present at upstream junction that enters reservoir.
        self.Q_jct = np.zeros(T)  # Upstream junction inflow, before any flow is distributed among outflow elements

    def Import_Data(self, T, Input_Data_File, stochastic_components, start_date, op_policy_params):
        self.op_policy_params = op_policy_params  # Will be reset if DPS. Set here to capture flushing preferences.
        if 'E-V-A-S' in Input_Data_File.sheetnames:
            # Import Elevation-Volume-Area-Sediment data
            self.E_V_A = Excel_Data_Import(self.name, Input_Data_File, 'E-V-A-S', 2, 4, max_distinct_data_types = None,
                                           data_name_offset = 3)
            #self.E_V_A_TRACK[0] = deepcopy(self.E_V_A)  # Time zero value of E_V_A is first curve before sedimentation.
            self.Original_River_Bed_Elevation = self.E_V_A[0][0]  # Initialize flushing channel elevation to original river bed elevation at each reservoir
            self.Settled_Volume_to_Distribute = 0  # Initialize variable used to distribute deposited sediment into the reservoir's storage volume/elevation profile.
            # Store small volumes in each reservoir for checking later whether reservoir is full of sediment or not. Miniscule volumes of water are left remaining at each elevation for model execution purposes.
            self.storage_sum_full_res = 0  # Initialize to zero
        else:
            logging.critical("Required worksheet 'E-V-A-S' is not provided in input data file")
        if 'Evaporation Data' in Input_Data_File.sheetnames:
            self.Monthly_Evap_Data = Excel_Data_Import(self.name, Input_Data_File, 'Evaporation Data',
                                                       1, 12, max_distinct_data_types=None,
                                                       data_name_offset=None)
        self.Environmental_Flow_Data = [None for i in range(12)]  # Default in case user does not specify worksheet.
        if 'Environmental Flow Data' in Input_Data_File.sheetnames:
            self.Environmental_Flow_Data = Excel_Data_Import(self.name, Input_Data_File,
                                                             'Environmental Flow Data', 1, 12,
                                                             max_distinct_data_types=None,
                                                             data_name_offset=None)
        # Set any "None" values to zero in environmental flow data array
        for i in range(len(self.Environmental_Flow_Data)):
            if self.Environmental_Flow_Data[i] is None:
                self.Environmental_Flow_Data[i] = 0
        if 'Reservoir Specifications' in Input_Data_File.sheetnames:
            [self.Act_Stor_Elev, self.Dead_Stor_Elev, self.Sed_Trapping_Curve_Spec, self.Brune_Time_Scale,
             self.Initial_Sediment_Mass, self.Reservoir_Capabilities, self.Reservoir_Operations_Goal,
             self.Initial_Storage, self.Dam_Tailwater_Elevation_1, self.MW_Capacity, self.Turbine_Efficiency,
             self.Reservoir_Length, self.Perform_Flushing, self.Perform_Dredging, self.Perform_Bypassing,
             self.Perform_Density_Current_Venting, self.Perform_Sluicing, self.Group_ID, self.Flushing_Group_ID,
             self.Sed_Group_ID, self.Non_Flushing_Group_ID, self.density_SS, self.DPS_Policy_Name,
             self.Min_Net_Head] = Excel_Data_Import(self.name, Input_Data_File,
                                                    'Reservoir Specifications', 1, 24,
                                                    max_distinct_data_types=None, data_name_offset=None)
            # Set default values for parameters read in from Reservoir Specifications.
            # Create entry in operating policy dictionary for each type of sediment management the user requests to have performed.
            # The order in which these imports occur often matters, so it is suggested not to reorganize code in this section.
            if self.Initial_Storage is not None:
                if type(self.Initial_Storage) in [int, float, long]:
                    pass  # If user specifies a valid number, proceed.
                else:
                    logging.critical("Initial reservoir water storage value {0} for Reservoir {1} is not "
                                     "valid".format(self.Initial_Storage, self.name))
            else:
                if self.Act_Stor_Elev is not None:
                    # Set time zero value. Default initial (time = 0) water storage (m^3) in reservoir is equal to
                    # volume associated with active storage.
                    self.Initial_Storage = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                self.Act_Stor_Elev)
                elif self.Dead_Stor_Elev is not None:
                    # Set time zero value. Default initial (time = 0) water storage (m^3) in reservoir is equal to
                    # volume associated with active storage
                    self.Initial_Storage = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                self.Dead_Stor_Elev)
                else:
                    self.Initial_Storage = 0  # If no active storage elevation is provided
            # Set time zero value for water surface elev and area, corresponding to initial storage.
            self.water_surface_elevation[0] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation",
                                                                   self.Initial_Storage)
            self.water_surface_area[0] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "area",
                                                              self.Initial_Storage)
            if self.Initial_Sediment_Mass is not None:
                if type(self.Initial_Sediment_Mass) in [int, float, long]:
                    pass  # If user specifies a valid number, proceed.
                else:
                    logging.critical("Initial reservoir deposited sediment value {0} for Reservoir {1} is not "
                                  "valid".format(self.Initial_Sediment_Mass, self.name))
            else:
                self.Initial_Sediment_Mass = 0
            if self.density_SS > 0:
                pass
            else:
                self.density_SS = 1000  # default (if no value provided by user): 1000 kg/m^3
            if self.Sed_Trapping_Curve_Spec is None:
                self.Sed_Trapping_Curve_Spec = 0  # if no value specified, set default trapping efficiency = 0
            elif type(self.Sed_Trapping_Curve_Spec) not in [str, unicode]:
                # A number is specified, but make sure it is <=1.
                if self.Sed_Trapping_Curve_Spec > 1 or self.Sed_Trapping_Curve_Spec < 0:
                    self.Sed_Trapping_Curve_Spec = 0  # if invalid value specified, set default trapping efficiency = 0
                    logging.critical("Invalid constant trapping efficiency fraction specified for reservoir {0}".format(
                        self.name))
            elif self.Sed_Trapping_Curve_Spec in ["L", "M", "H"]:
                if self.Sed_Trapping_Curve_Spec == "L":
                    self.m_br = 3
                    self.a_br = 0.00000102
                    self.b_br = -0.00013
                    self.c_br = 0.0262
                    self.d_br = 1.0266
                elif self.Sed_Trapping_Curve_Spec == "M":
                    self.m_br = 1
                    self.a_br = 0.012
                    self.b_br = 1.02
                    self.c_br = 0
                    self.d_br = 0
                elif self.Sed_Trapping_Curve_Spec == "H":
                    self.m_br = 2
                    self.a_br = 0.000003
                    self.b_br = 0.0063
                    self.c_br = 0.9947
                    self.d_br = 0
                if self.Brune_Time_Scale == "A":
                    self.days_in_sum = 365
                elif self.Brune_Time_Scale == "M":
                    self.days_in_sum = 30
                elif self.Brune_Time_Scale == "D":
                    self.days_in_sum = 1
                else:
                    self.Brune_Time_Scale = "A"  # Default Brune (1953) Curve calculation time scale is annual.
                    self.days_in_sum = 365
            elif self.Sed_Trapping_Curve_Spec == "C":
                self.days_in_sum = 1
            if self.days_in_sum > 0:
                self.S_over_last_period = np.zeros(self.days_in_sum)
                self.V_out_over_last_period = np.zeros(self.days_in_sum)
            if self.Act_Stor_Elev is not None:
                # Determine active and dead storage values
                self.Upper_Storage = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                          self.Act_Stor_Elev)  # Find water storage corresponding to active storage elevation.
                self.Lower_Storage = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                          self.Dead_Stor_Elev)  # Find water storage corresponding to dead storage elevation.
                self.Res_Avg_Area = self.Upper_Storage / (self.Act_Stor_Elev - self.Original_River_Bed_Elevation)  # Determine average surface area over which sediment settles (this is a simple average, used to track elevation of sediment accumulate outside of flushing channel)
                self.Active_Storage = self.Upper_Storage - self.Lower_Storage
                self.Dead_Storage = self.Lower_Storage - self.Initial_Sediment_Mass / self.density_SS
                if self.Dead_Storage < 0:
                    self.Active_Storage = max(0, self.Active_Storage - abs(self.Dead_Storage), 0)
                    self.Dead_Storage = 0
                else:
                    pass
                # Store information about initial active and dead storage in additional matrices
                self.capacity_active_reservoir[0] = self.Active_Storage
                self.capacity_dead_reservoir[0] = self.Dead_Storage
            else:
                logging.critical("No active storage capacity value provided for reservoir {0} in Reservoir " \
                               "Specifications worksheet for".format(self.name))

            # Assign Number to each Reservoir Operations Goal Type, and set defaults
            if self.Reservoir_Operations_Goal == "Meet specified daily water storage targets (m^3)":
                self.Reservoir_Operations_Goal_ID = 1
            elif self.Reservoir_Operations_Goal == "Meet specified daily water elevation (mamsl) targets":
                self.Reservoir_Operations_Goal_ID = 2
            elif self.Reservoir_Operations_Goal == "Meet daily energy (MWH) targets":
                self.Reservoir_Operations_Goal_ID = 3
            elif self.Reservoir_Operations_Goal == "Meet daily water release (m^3/s) targets":
                self.Reservoir_Operations_Goal_ID = 4
            elif self.Reservoir_Operations_Goal == "Specify DPS policy":
                # User can either provide a file in which DPS policy parameters are contained. Or,
                # the Borg optimization model can call PySedSim automatically and create this DPS policy parameter
                # file automatically in each optimization iteration.
                self.Reservoir_Operations_Goal_ID = 5

                if op_policy_params is None:
                    # User has not directly linked an optimization model to PySedSim, and/or wishes for the
                    # operating policy parameters to be read in from a text file rather than fed in directly.
                    try:
                        # Load all relevant user specifications from DPS worksheet
                        self.op_policy_params = direct_policy_search.Import_DPS_Preferences(input_data_file=Input_Data_File)
                        self.op_policy_params['Policy Function'] = {'Raw Parameters': np.loadtxt('RBF_Parameters.txt')}
                    except KeyError:
                        logging.critical("Required Worksheet 'DPS' does not exist in input file.")
                else:
                    self.op_policy_params = op_policy_params
                # Initialize DPS-relevant arrays
                # Update DPS dictionary to store C, R, and W for policy function.
                self.op_policy_params = direct_policy_search.Create_DPS_Policy(self.op_policy_params)
                self.DPS_inputs = np.zeros(self.op_policy_params['num_inputs'])  # Stores daily inputs to DPS policy
                self.DPS_Decision = np.zeros(T+1)
            elif self.Reservoir_Operations_Goal == "Meet specified daily water elevation (mamsl) targets - annually recurring":
                self.Reservoir_Operations_Goal_ID = 2  # value is still 2, because it's still an elevation target. We just must modify its form so it is a longer matrix.
                self.recurring_daily_elevation_target = 1  # If this is 1, the model will not load daily elevation targets for every simulation day for this reservoir
            elif self.Reservoir_Operations_Goal == "Water surface elevation targets as a function of inflow":
                self.Reservoir_Operations_Goal_ID = 6
            else:
                self.Reservoir_Operations_Goal_ID = 2  # Default is to meet daily water elevation targets

            if self.Reservoir_Operations_Goal == "Meet specified daily water elevation (mamsl) targets - annually recurring":
                if 'Elevation Target Recurring' in Input_Data_File.sheetnames:
                    self.Recurring_elevation_target = Excel_Data_Import(self.name, Input_Data_File,
                                                                        'Elevation Target Recurring', 0,
                                                                        366,
                                                                        max_distinct_data_types=None,
                                                                        data_name_offset=None,
                                                                        start_date=start_date)
                    # Assume these dates run from Jan 1 to Dec 31, and include a Feb 29 value. (Hence, there are assumed to be 366 values).
                    self.Recurring_elev_target_dict = {}  # Annually recurring reservoir elevation target dictionary that stores targets by month/day.

                    for z in range(1, 12+1):
                        self.Recurring_elev_target_dict[z] = {}  # Initialize month dictionary, so can store daily dict in it later.

                    i = datetime(1988, 1, 1)  # Pick a random leap year that will let us store month/day values for 366 days.
                    for j in range(366):
                        # Store 366 values for the common year (guide curve type operating policy).
                        self.Recurring_elev_target_dict[i.month][i.day] = self.Recurring_elevation_target[j]
                        i += timedelta(days = 1)  # increment by one day

                    i = self.Element_Sub_Dict["Simulation Start Date"]
                    for j in range(self.T):
                        self.ELEVATION_TARGET[j+1] = self.Recurring_elev_target_dict[i.month][i.day]  # Stores elevation targets for end of every day (so only Sim_Dur = T values are stored, as we don't set t = 0 target).
                        i += timedelta(days = 1)  # increment by one day
                    # Store a time 0 storage target corresponding to specified time zero elevation target.
                    self.STORAGE_TARGET[0] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                  self.ELEVATION_TARGET[0])
                else:
                    logging.critical("Worksheet 'Elevation Target Recurring' is not provided in input data file, "
                                     "but is required given user specifications for reservoir operations related to  "
                                     "Reservoir {0}".format(self.name))
            if self.Reservoir_Operations_Goal_ID == 1:
                if 'Storage Volume Target' in Input_Data_File.sheetnames:
                    self.STORAGE_TARGET = Excel_Data_Import(self.name, Input_Data_File,
                                                            'Storage Volume Target', 0, T + 1,
                                                            max_distinct_data_types=None,
                                                            data_name_offset=None,
                                                            start_date=start_date)
                    # Store a time 0 elevation target corresponding to specified time zero storage target.
                    self.ELEVATION_TARGET[0] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation",
                                                                    self.STORAGE_TARGET[0])
                else:
                    logging.critical("worksheet 'Storage Volume Target' is not provided in input data file, " \
                          "but is required given user specifications for reservoir operations related to Reservoir " \
                          "{0}".format(self.name))
            elif self.Reservoir_Operations_Goal == "Meet specified daily water elevation (mamsl) targets":
                if 'Storage Volume Elevation Target' in Input_Data_File.sheetnames:
                    self.ELEVATION_TARGET = Excel_Data_Import(self.name, Input_Data_File,
                                                              'Storage Volume Elevation Target', 0,
                                                              T + 1, max_distinct_data_types=None,
                                                              data_name_offset=None,
                                                              start_date=start_date)
                    # Store a time 0 storage target corresponding to specified time zero elevation target.
                    self.STORAGE_TARGET[0] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                  self.ELEVATION_TARGET[0])
                else:
                    logging.critical("worksheet 'Storage Volume Elevation Target' is not provided in input data file, " \
                          "but is required given user specifications for reservoir operations related to Reservoir " \
                          "{0}".format(self.name))
            elif self.Reservoir_Operations_Goal_ID == 3:
                if 'Daily Energy Targets' in Input_Data_File.sheetnames:
                    pass
                    # SET THESE LATER WHEN NEEDED:
                    # self.Daily_Energy_Targets = Excel_Data_Import(self.name, Input_Data_File, 'Daily Energy Targets' , 0, T,
                    # max_distinct_data_types = None, data_name_offset = 4, start_date = start_date)
                    # self.Energy_Release_Goal = np.zeros(T)
                else:
                    logging.critical("worksheet 'Daily Energy Targets' is not provided in input data file, " \
                          "but is required given user specifications for reservoir operations related to Reservoir " \
                          "{0}".format(self.name))
            elif self.Reservoir_Operations_Goal_ID == 4:
                if 'Daily Release Targets' in Input_Data_File.sheetnames:
                    pass  # self.Daily_Release_Target = data_processing.Excel_Data_Import(self.name, Input_Data_File, 'Daily Release Targets', 0, T,
                    # max_distinct_data_types = None, data_name_offset = 4, start_date = start_date)
                else:
                    logging.critical("worksheet 'Daily Release Targets' is not provided in input data file, " \
                          "but is required given user specifications for reservoir operations related to Reservoir " \
                          "{0}".format(self.name))
            elif self.Reservoir_Operations_Goal_ID == 6:
                if 'WSE-Inflow Policy' in Input_Data_File.sheetnames:
                    self.wse_inflow_target = Excel_Data_Import(self.name, Input_Data_File,
                                                               'WSE-Inflow Policy', 2, 2,
                                                               max_distinct_data_types=None,
                                                               data_name_offset=3)
                    # max_distinct_data_types = None, data_name_offset = 4, start_date = start_date)
                else:
                    logging.critical("worksheet 'WSE-Inflow Policy' is not provided in input data file, but is required " \
                          "given user specifications for reservoir operations related to Reservoir {0}".format(self.name))

            # Sediment management
            if self.Perform_Flushing == 'Yes':
                # Instantiate flushing object, and initialize some additional attribute values.
                self.Operating_Policy["Flushing"] = Flushing(self.name, T, Input_Data_File, self.density_SS,
                                                             self.Flushing_Group_ID, self.Non_Flushing_Group_ID,
                                                             self.Act_Stor_Elev,
                                                             self.Orifices.Low_Level_Outlet_Invert_Elevation,
                                                             self.Res_Avg_Area,
                                                             self.Original_River_Bed_Elevation, start_date,
                                                             opt_params = self.op_policy_params)
                self.Operating_Policy["Flushing"].Flushing_Channel_Elevation[0:2] = self.E_V_A[0][0]  # Flushing channel is assumed to built from the ground up
                self.Operating_Policy["Flushing"].Original_River_Bed_Elevation = self.Original_River_Bed_Elevation
                self.Operating_Policy["Flushing"].water_surface_elevation_pre_flush[0] = Matrix_Interpolation(self.name,
                                                                                                              self.E_V_A,
                                                                                                              "storage",
                                                                                                              "elevation",
                                                                                                              self.S[0])
            else:
                self.Perform_Flushing = 'No'  # Set default
            if self.Perform_Sluicing == 'Yes':
                self.Operating_Policy["Sluicing"] = Sluicing(self.name, T, Input_Data_File, self.Element_Sub_Dict["Simulation Start Date"])  # Instantiate sluicing class
            else:
                self.Perform_Sluicing = 'No'  # Set default
            if self.Perform_Dredging == 'Yes':
                self.Operating_Policy["Dredging"] = Dredging(self.name, T, Input_Data_File)  # Instantiate dredging object
                self.Element_Sub_Dict["Dredging Outflow Element"] = self.Operating_Policy["Dredging"]["Sediment Destination"]  # Set destination element for dredged sediment.
            else:
                self.Perform_Dredging = 'No'  # Set default
            if self.Perform_Density_Current_Venting == 'Yes':
                self.Operating_Policy["Density Current Venting"] = Density_Current_Venting(self.name, T, Input_Data_File, self.Reservoir_Length)  # Instantiate DCV object
            else:
                self.Perform_Density_Current_Venting = 'No'  # Set default

            # Set reservoir type as one of these possible types: "Storage", "Power", "Diversion", "Diversion and Power"
            try:
                if self.Element_Sub_Dict["Reservoir Type"] == "Diversion":
                    if self.MW_Capacity is not None:
                        self.Reservoir_Type = "Diversion and Power"
                    else:
                        self.Reservoir_Type = "Diversion"
            except KeyError:
                if self.MW_Capacity is not None:
                    self.Reservoir_Type = "Power"
                else:
                    self.Reservoir_Type = "Storage"

            # If reservoir has power, initialize relevant matrices.
            if self.Reservoir_Type in ["Power", "Diversion and Power"]:
                self.hydropower_head_average = np.zeros(T)  # Computes average hydropower head. Only need if dam produces hydropower.
                self.Max_MW_Turbine_Flow = np.zeros(T)  # The maximum flow (m^3/s) allowed for a certain water level elevation such that the MW capacity of plant is not exceeded. Only need if dam produces hydropower.
        else:
            logging.critical("required worksheet 'Reservoir Specifications' is not provided in input data file")

        if 'Tailwater Rating Curve' in Input_Data_File.sheetnames:
            self.Tailwater_Rating_Curve = Excel_Data_Import(self.name, Input_Data_File, 'Tailwater Rating Curve', 2, 2,
                                                            max_distinct_data_types=None, data_name_offset=3)
        if 'Egg Passage Targets' in Input_Data_File.sheetnames:
            # Initialize relevant arrays
            self.larv_pass = np.zeros(T)  # Binary time series variable that records success or failure to meet particular elevation targets given daily inflow conditions.
            self.Egg_Passage_Elevation_Target = np.zeros(T+1)  # End of day elevation target that will lead to egg passage velocity success given inflow conditions and an Egg Passage elevation-inflow success table.
            self.egg_pass_resilience = np.zeros(T)  # Tracks days when a failure to achieve adequate egg passage conditions occurs after the very same failure condition has occurred on a previous day.
            self.larv_stor = np.zeros(T+1)  # Tracks alive larvae "mass" contained in reservoir.
            self.larv_mass_out = np.zeros(T)  # Tracks larvae mass release
            self.larv_mass_out_surv = np.zeros(T)  # Tracks whether released larvae mass is alive or dead
            self.larv_mass_out_surv_total = np.zeros(T)  # Total mass released from site, accounting for spillage into
            # bypass if one exists.
            # reservoir every day (=1 unit of larvae).
            # Import data from input sheet
            self.Egg_Passage_Rating_Curve = Excel_Data_Import(self.name, Input_Data_File,
                                                              'Egg Passage Targets', 2, 4,
                                                              max_distinct_data_types=None,
                                                              data_name_offset=3)
            if len(self.Egg_Passage_Rating_Curve[0]) > 0:
                self.Evaluate_Egg_Passage = 1
                # Larvae Passage
                self.larv_pass_over_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin > provided flow threshold.
                self.larv_pass_over_thold[:] = 'NaN'  # Fill array with 'NaN'
                self.larv_pass_under_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.larv_pass_under_thold[:] = 'NaN'  # Fill array with 'NaN'

                # Larvae Survival
                self.larv_surv_over_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.larv_surv_over_thold[:] = 'NaN'  # Fill array with 'NaN'
                self.larv_surv_under_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.larv_surv_under_thold[:] = 'NaN'  # Fill array with 'NaN'

                # Residence Time
                self.res_time_over_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.res_time_over_thold[:] = 'NaN'  # Fill array with 'NaN'
                self.res_time_under_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.res_time_under_thold[:] = 'NaN'  # Fill array with 'NaN'

                # Residence Time Reliability
                self.res_time_reliab_over_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.res_time_reliab_over_thold[:] = 'NaN'  # Fill array with 'NaN'
                self.res_time_reliab_under_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.res_time_reliab_under_thold[:] = 'NaN'  # Fill array with 'NaN'

                # Energy production above/below threshold
                self.Hydropower_avg_MWH_over_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.Hydropower_avg_MWH_over_thold[:] = 'NaN'  # Fill array with 'NaN'
                self.Hydropower_avg_MWH_under_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
                self.Hydropower_avg_MWH_under_thold[:] = 'NaN'  # Fill array with 'NaN'

                # Larvae passage discharge (flow) threshold (m^3/s), above and below which to track passage success
                if self.Egg_Passage_Rating_Curve[2] != []:
                    self.larv_pass_flow_thold = self.Egg_Passage_Rating_Curve[2][0]
                else:
                    self.larv_pass_flow_thold = 0  # Track passage for all inflows

                # Residence time threshold (days), above and below which to track passage success
                if self.Egg_Passage_Rating_Curve[3] != []:
                    self.res_time_thold = (self.Egg_Passage_Rating_Curve[3][0])/365.25
                else:
                    self.res_time_thold = 0  # Track passage for all inflows
            else:
                self.Evaluate_Egg_Passage = 0
        else:
            self.Evaluate_Egg_Passage = 0

        # Import data related to Maximum Water Level, if reservoir has an unregulated outflow point at the upstream
        # end of the reservoir. Not a common scenario; use this feature with caution, as it assumes there is not
        # complete control of the reservoir's water levels.
        if 'Max Reservoir WSE' in Input_Data_File.sheetnames:
            self.TERMINAL_WL_Rating_Curve = Excel_Data_Import(self.name, Input_Data_File,
                                                              'Max Reservoir WSE', 2, 2,
                                                              max_distinct_data_types=None,
                                                              data_name_offset=3)
            if len(self.TERMINAL_WL_Rating_Curve[0]) > 0:
                self.Cap_Res_WL = 1
                self.Terminal_WSE = np.zeros(T+1)  # Maximum WSE of reservoir if uncontrolled upstream outlet/spillway
            else:
                self.Cap_Res_WL = 0
        else:
            self.Cap_Res_WL = 0

        # Import/create data related to distribution of deposited sediment within the reservoir's elevation/storage profile.
        for i in range(1, len(self.E_V_A[0])):
            self.storage_sum_full_res += i
        # If user has provided trapped sediment reservoir elevation distribution data in 4th column, import it. Otherwise set default values.
        self.E_Sed = [[], [], []]
        self.E_Sed[0] = self.E_V_A[0]  # E_Sed first column is identical to E_V_A first column (=elevations)
        if '9' not in stochastic_components:
            # Distribution of sediment within E-V-A curve is assumed to be deterministic. Don't randomly generate E_Sed curve.
            if not self.E_V_A[3]:
                # If sediment fraction list in E-V-A is empty (user did not specify), specify defaults here.
                if self.Act_Stor_Elev is not None:
                    elevation_val = self.Act_Stor_Elev  # User user-specified active storage elevation
                else:
                    elevation_val = self.E_V_A[0][-1]  # Just use last value in the E_V_A elevation table (using python's "-1" index)
                self.E_Sed[1].append(0)  # Set first sediment deposition fraction value = 0, since loop below begins at 1st element
                self.E_Sed[2].append(0)  # Set first sediment removal fraction value = 0, since loop below begins at 1st element
                # Load fractions based on linear relationship using elevation_val
                for i in range(1, len(self.E_V_A[0])):
                    self.elev_counter = i
                    if self.E_Sed[0][i] <= elevation_val:
                        self.E_Sed[1].append((self.E_Sed[0][i] - self.Original_River_Bed_Elevation)/(elevation_val - self.Original_River_Bed_Elevation))  # For sediment deposition
                        self.E_Sed[2].append((self.E_Sed[0][i] - self.Original_River_Bed_Elevation)/(elevation_val - self.Original_River_Bed_Elevation))  # For sediment removal
                        if self.E_Sed[0][i] == elevation_val:
                            break  # break for loop
                    else:
                        # Elevation value in E_Sed (during looping) exceeds active storage elevation (or max elevation value in EVA table), so store this as last value you'll loop to for distributing sediment.
                        self.E_Sed[1].append((self.E_Sed[0][i] - self.Original_River_Bed_Elevation)/(elevation_val - self.Original_River_Bed_Elevation))  # For sediment deposition
                        self.E_Sed[2].append((self.E_Sed[0][i] - self.Original_River_Bed_Elevation)/(elevation_val - self.Original_River_Bed_Elevation))  # For sediment removal
                        break  # break for loop
            else:
                self.elev_counter = len(self.E_V_A[3]) - 1 # Only loop to end of this list when assigning deposited sediment to different elevations.
                self.E_Sed[1] = self.E_V_A[3]  # If sediment fraction list is not empty (user did specify), set values here. First E_Sed list is the E-V-A fraction list.
                self.E_Sed[2] = self.E_V_A[3] # Assume second list in E_Sed (which is for distribution of sediment removal during flushing and general sediment removal) for now is an equal distribution across the dead and active storage elevation range.
        else:
            pass  # Will externally randomly generate entire E_Sed curve.

    def Import_External_Element_State(self, t):
        '''
        Imports data from other elements being simulated (e.g., junctions, reaches and reservoirs) for use in the
        simulation of this element.

        Input data are stored in the elem_xfer_input_dict (which originally came from other elements'
        elem_xfer_output_dict). This module simply takes this information and converts it into a specific use for
        this element. For this reason, each element will have their own unique implementation of this method. This
        method is executed at the beginning of every time step, before any processes are simulated. The module must be run a second time
        for reservoir elements that have a junction downstream with tailwater impacted by multiple elements flowing into junction.
        '''

        # Loop through elements sending data, and unpack it to provide use for this element.
        for keys in self.elem_xfer_input_dict:
            if keys in self.Element_Sub_Dict['Outflow Elements']:
                # Located data to be transferred from downstream junction regarding water level.
                if self.elem_xfer_input_dict[keys]["Num Inflow Elements"] > 1:
                    # Downstream junction has multiple inflow points. Will need to re-calculate tailwater after all
                    # elements have been completed in pysedsim_main_simulation.py since tailwater not properly
                    # computed.
                    self.re_calc_energy = 1
                for sub_keys in self.elem_xfer_input_dict[keys]:
                    if sub_keys == "Downstream Junction Flow":
                        self.Q_Downstream_Junction[t] = self.elem_xfer_input_dict[keys][sub_keys]
            elif keys in self.Element_Sub_Dict['Inflow Elements']:
                for sub_keys in self.elem_xfer_input_dict[keys]:
                    if sub_keys == "Upstream Junction Inflow Fraction":
                        self.fraction_Q_jct[t] = self.elem_xfer_input_dict[keys][sub_keys]
                    elif sub_keys == "Upstream Junction Inflow":
                        self.Q_jct[t] = self.elem_xfer_input_dict[keys][sub_keys]

    def Water_Storage_Mass_Balance(self, t):
        # Determine final water storage (t+1) in reservoir n for internal hydrologic simulation case
        # If sluicing exists for this reservoir, this method also checks to see if post_sluicing is occurring, and if so, if it can be ended based on storage/elevation conditions.
        self.S[t+1] = self.S[t] + self.V_in[t] - (self.Q_out[t] + self.Evap[t]) * self.dt * 86400 - self.Settled_volume[t]
        if self.S[t + 1] < 0:
            self.S[t + 1] = 0
        try:
            self.Operating_Policy["Sluicing"].Post_Sluicing_Check(self.STORAGE_TARGET[t+1], self.S[t+1])
        except KeyError:
            pass  # Sluicing does not exist for this reservoir.

    def Trapped_Load(self, t):
        # This method computes reservoir sediment trapping efficiency (TE), then determines how much sediment is trapped, adjusts mass BS_W, BS_V, and SS_C, and determines new reservoir active and dead storage capacities.
        exit_loop = 0
        while exit_loop == 0:
            if t != 0:
                # Compute residence time, separately from sediment trapping efficiency calculations.
                if self.Brune_Time_Scale == "A":
                    # TE computed over annual time scale
                    if t <= 365:
                        self.days_in_sum = t
                    else:
                        self.days_in_sum = 365
                elif self.Brune_Time_Scale == "M":
                    # TE computed over monthly time scale
                    if t <= 30:
                        self.days_in_sum = t
                    else:
                        self.days_in_sum = 30
                elif self.Brune_Time_Scale == "D":
                    self.days_in_sum = 1
                for x in range(1, self.days_in_sum + 1):
                    self.S_over_last_period[x-1] = self.S[t - x]
                    self.V_out_over_last_period[x-1] = self.V_out[t - x]
                self.Average_S_over_last_period = sum(self.S_over_last_period) / self.days_in_sum
                self.Sum_V_out_over_last_period = sum(self.V_out_over_last_period) / (self.days_in_sum / 365)
                if self.Sum_V_out_over_last_period > 0:
                    self.Residence_Time[t] = self.Average_S_over_last_period / self.Sum_V_out_over_last_period
                else:
                    self.Residence_Time[t] = 5E+19  # Really should be inf (and not counted in sum and average).

            # Order of calls below must be maintained
            try:
                if self.Operating_Policy["Flushing"].Current == 1 or self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down == 1:
                    self.TE_avg[t] = 0
                    break
            except KeyError:
                pass
            try:
                if self.Operating_Policy["Density Current Venting"].Current == 1 and self.Operating_Policy["Sluicing"].Current == 1:
                    # This section must come before the separate sections for individual density currents and sluicing are done, otherwise those others will be picked up first.
                    # MODEL DOES NOT CURRENTLY ACCOUNT FOR THIS EVENT.
                    break
            except KeyError:
                pass
            try:
                if self.Operating_Policy["Density Current Venting"].Current == 1:
                    self.TE_avg[t] = 0  # Actual TE associated with DCV will be set later in DCV mass balance routine. Be careful adding a direct connection to this routine in DCV file, as you are already adding sed mass if dreding occurred as well.
                    break
            except KeyError:
                pass
            try:
                if self.Operating_Policy["Sluicing"].Current == 1 or self.Sed_Trapping_Curve_Spec == "C":
                    # Apply Churchill method (either because sluicing is occurring, or because user has externally specified Churchill always be used).
                    # First, determine current reservoir length
                    self.Reservoir_Length_Current = (self.water_surface_elevation[t] - self.Original_River_Bed_Elevation) / ((self.Act_Stor_Elev - self.Original_River_Bed_Elevation) / self.Reservoir_Length)
                    if self.Reservoir_Length_Current > 0:
                        if self.Q_in[t] > 0:
                            Sedimentation_Index = (((self.S[t]/self.Q_in[t]) ** 2) / self.Reservoir_Length_Current)
                            if self.S[t] > 0:
                                if 800 * ((Sedimentation_Index / 3.28) ** (-0.2)) - 12 + self.stoch_trap_adjust > 100:
                                    self.TE_avg[t] = 0
                                elif 800 * ((Sedimentation_Index / 3.28) ** (-0.2)) - 12 + self.stoch_trap_adjust < 0:
                                    self.TE_avg[t] = 1
                                else:
                                    self.TE_avg[t] = (100 - (800 * ((Sedimentation_Index / 3.28) ** (-0.2)) - 12 + self.stoch_trap_adjust)) / 100
                            else:
                                self.TE_avg[t] = 0
                            if self.TE_avg[t] < 0:
                                self.TE_avg[t] = 0
                        else:
                            self.TE_avg[t] = 0
                    else:
                        self.TE_avg[t] = 0
                    break
            except KeyError:
                pass
            try:
                if t != 0:
                    if type(self.Sed_Trapping_Curve_Spec) not in [str, unicode]:
                        self.TE_avg[t] = self.Sed_Trapping_Curve_Spec
                    else:
                        if self.Residence_Time[t] == 0:
                            self.TE_avg[t] = 0
                        elif self.Residence_Time[t] == 5E+19:
                            self.TE_avg[t] = 1
                        else:
                            # Use piecewise functionalized Brune curve to determine TE.
                            if self.Sed_Trapping_Curve_Spec in ["L", "M", "H"]:
                                # Try Gill approximation of Brune instead
                                self.TE_avg[t] = (self.Residence_Time[t] ** self.m_br) / (self.a_br + self.b_br * self.Residence_Time[t] + self.c_br * (self.Residence_Time[t] ** 2) + self.d_br * (self.Residence_Time[t] ** 3))
                        if self.TE_avg[t] < 0:
                            self.TE_avg[t] = 0
                        if self.TE_avg[t] > 1:
                            self.TE_avg[t] = 1
                break
            except KeyError:
                pass
        # Compute impact of sedimentation on reservoir storage volume, and adjust total storage volume (and active/dead storages) accordingly.
        self.Settled_mass[t] = self.TE_avg[t] * self.SS_W_in[t]
        self.Settled_volume[t] = self.Settled_mass[t] / self.density_SS
        Reservoir.Sediment_Storage_Capacity_Check(self, t)
        self.BS_W[t + 1] += self.BS_W[t] + self.Settled_mass[t]  # IMPORTANT: the "+=" is used here because BS_W[t+1] could have a value in it, placed there as a result of dredging, if the element where dredging occurred was (upstream and) simulated before this element is simulated.
        self.BS_V[t + 1] = self.BS_W[t + 1] / self.density_SS

        if self.S[t] + self.V_in[t] - self.Settled_volume[t] <= 0:
            self.SS_C[t] = 0
        else:
            self.SS_C[t] = ((self.SS_W[t] + self.SS_W_in[t]) - self.Settled_mass[t]) / (self.S[t] + self.V_in[t])
        Reservoir.Active_Dead_Storage_Capacity_Calculation(self, t)  # Compute new active storage capacity since settled mass has been computed.

    def Sediment_Storage_Capacity_Check(self, t):
        # This method ensures that if the reservoir is nearly full to capacity (dead + active), no more sediment can settle than the storage capacity available.
        if self.Settled_volume[t] >= self.capacity_active_reservoir[t] + self.capacity_dead_reservoir[t]:
            self.reservoir_is_full = 1
            self.Settled_volume[t] = self.capacity_active_reservoir[t] + self.capacity_dead_reservoir[t]
            self.Settled_mass[t] = self.Settled_volume[t] * self.density_SS
            self.Q_out[t] = self.Q_in[t]  # Reservoir is full of sediment, so dam now becomes man-made waterfall with no storage (run-of-river, but with no routing).
        elif self.capacity_active_reservoir[t] + self.capacity_dead_reservoir[t] <= self.storage_sum_full_res:
            # Some small storage may remain for purposes of maintaining some semblance of an E-V-A curve, in which case we will still consider the reservoir to be full despite this minuscule storage capacity we assume exists.
            self.reservoir_is_full = 1
            self.Q_out[t] = self.Q_in[t]  # Reservoir is full of sediment, so dam now becomes man-made waterfall with no storage (run-of-river, but with no routing).

    def Active_Dead_Storage_Capacity_Calculation(self, t):
        # This method computes the active and dead reservoir storage capacity after: trapping is computed, or various forms of sediment REMOVAL take place (thereby affecting E-V-A), such as flushing, dredging, and DCV.
        if self.reservoir_is_full == 0:
            # Reservoir is not full of sediment. Compute new active/dead reservoir storage capacity normally.
            self.capacity_dead_reservoir[t + 1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                       self.Dead_Stor_Elev)
            self.capacity_active_reservoir[t + 1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                         self.Act_Stor_Elev) - self.capacity_dead_reservoir[t + 1]
        else:
            self.capacity_active_reservoir[t+1] = 0
            self.capacity_dead_reservoir[t+1] = 0

    def E_Sed_Curve_Adjustment(self, frac_act_stor_cap):
        # Method Description:
        # Re-adjust E_Sed matrix to account for uncertainty in sediment deposition curve if user has specified such. Requires that E_Sed has already been poulated in Data_Import section for this reservoir.

        # Method Inputs:
        #frac_act_stor_cap is the fraction of sediment that will be stored in the reservoir's active storage capacity. Typically will be set equal to Sampled_Parameter_Dict['10'][rz] in SedSim_v2.py.

        # Method outputs/results:
        # Resets self.E_Sed curve (for reservoir instance).

        # Step 1. Create list indicating what fraction of the total fraction each elevation stores.
        fraction = []  # Initialize empty list
        initial_dead_fract = self.E_Sed[1][self.E_Sed[0].index(self.Dead_Stor_Elev)]
        for elev in range(len(self.E_Sed[0])):
            if self.E_Sed[0][elev] < self.Dead_Stor_Elev:
                fraction.append(self.E_Sed[1][elev]/initial_dead_fract)

        # Step 2. Reset E_Sed[1] fraction values for elevations within the active storage range.
        for elev in range(len(self.E_Sed[0])):
            # Loop through elevations and store
            if self.E_Sed[0][elev] < self.Dead_Stor_Elev:
                pass  # Only looking to re-adjust E_Sed values in the active storage elevation range in this section of code. Will adjust lower values in next section (Step 3).
            elif self.E_Sed[0][elev] == self.Dead_Stor_Elev:
                self.E_Sed[1][elev] = (1-frac_act_stor_cap)
            else:
                if self.E_Sed[0][elev] < self.Act_Stor_Elev:
                    self.E_Sed[1][elev] = (1-frac_act_stor_cap) + (self.E_Sed[1][elev]-self.Dead_Stor_Elev)*((frac_act_stor_cap)/(self.Act_Stor_Elev-self.Dead_Stor_Elev))
                else:
                    # All values above active storage elevation set = 1
                    self.E_Sed[1][elev] = 1

        # Step 3. Down-adjust all values below the dead storage elevation to reflect the new fraction assumed to  be stored in active storage zone.
        for elev in range(len(self.E_Sed[0])):
            if self.E_Sed[0][elev] < self.Dead_Stor_Elev:
                self.E_Sed[1][elev] = fraction[elev]*((1-frac_act_stor_cap))

    def Reservoir_Volume_Reduction(self, t, E_Sed):
        # This method adjusts E-V-A curve to reflect sediment deposition. Can be done after Trapped_load() method completes.
        self.Settled_Volume_to_Distribute += self.Settled_volume[t]
        # Update E-V-A curve given new sediment deposition; perform this process once per month only, because this is very calculation intensive and becomes time consuming
        if self.current_date.day == 1:
            Extra_Sediment_Storage = 0
            for i in range(1, self.elev_counter + 1):
                self.E_V_A[1][i] -= self.Settled_Volume_to_Distribute * (E_Sed[i] - E_Sed[i-1])
                if (self.E_V_A[1][i] - self.E_V_A[1][i-1]) > 0:
                    # No additional action required, as more storage space for water is still available at this elevation
                    Incremental_Settled_Volume = self.Settled_Volume_to_Distribute * (E_Sed[i] - E_Sed[i-1])
                else:
                    # Additional sediment cannot be stored within the incremental storage space available at this elevation
                    Incremental_Extra_Sediment_Volume = abs(self.E_V_A[1][i] - self.E_V_A[1][i-1]) + i  # Add i to ensure some volume is left over for calculation/model execution purposes.
                    self.E_V_A[1][i] = self.E_V_A[1][i-1] + i  # To ensure that some storage is available at all elevations, to be realistic.
                    Extra_Sediment_Storage += Incremental_Extra_Sediment_Volume
                    Incremental_Settled_Volume = self.Settled_Volume_to_Distribute * (E_Sed[i] - E_Sed[i-1]) - Incremental_Extra_Sediment_Volume
                # Loop through all elements in list at higher elevations to be sure their cumulative volume is reduced. Must loop through entire list, not just up to the elev_counter (usually = active storage elevation).
                for x in range(i + 1, len(self.E_V_A[1])):
                    self.E_V_A[1][x] -= Incremental_Settled_Volume
            Incremental_Extra_Sediment_Volume = 0
            #  First, try to store the extra sediment at a higher elevation
            if Extra_Sediment_Storage > 0:
                for i in range(1, self.elev_counter + 1):
                    if self.E_V_A[1][i] - self.E_V_A[1][i-1] > i:
                        Incremental_Extra_Sediment_Volume = Extra_Sediment_Storage
                        self.E_V_A[1][i] = self.E_V_A[1][i] - Extra_Sediment_Storage
                        Extra_Sediment_Storage = 0
                        if self.E_V_A[1][i] - self.E_V_A[1][i-1] < 0:
                            Extra_Sediment_Storage = abs(self.E_V_A[1][i] - self.E_V_A[1][i-1]) + i
                            Incremental_Extra_Sediment_Volume = Incremental_Extra_Sediment_Volume - Extra_Sediment_Storage
                            self.E_V_A[1][i] = self.E_V_A[1][i-1] + i
                        else:
                            pass
                    else:
                        pass  # Elevation level of full of sediment, so move to next elevation
                    # Loop through all elements in list at higher elevations to be sure their cumulative volume is reduced. Must loop through entire list, not just up to the elev_counter (usually = active storage elevation).
                    for x in range(i + 1, len(self.E_V_A[1])):
                        self.E_V_A[1][x] -= Incremental_Extra_Sediment_Volume
                    if Extra_Sediment_Storage == 0:
                        break  # All extra sediment has been allocated to an elevation
                    else:
                        pass
            else:
                pass
            self.Settled_Volume_to_Distribute = 0
        else:
            pass  # Don't adjust E-V-A curve today; only do once per month.
        #self.E_V_A_TRACK[t+1] = deepcopy(self.E_V_A)  # In every time step, add a new E_V_A curve to the track list.

    def Reservoir_Volume_Addition(self, Load_Removed, E_Sed):
        '''
        Purpose: Method adjusts E-V-A curve (adds back volume) to reflect sediment removal via flushing or general
        sediment removal (e.g., dredging). Can be done after Trapped_load() method completes, or after various sediment management events have been completed.

        Args:
            Load_Removed: amount of sediment (kg) removed from sediment management. Will be equal to either
            Flushing_Load_Removed_Daily[t] or Sediment_Load_Removed_Daily[t]

            E_Sed: should be self.E_Sed for the reservoir, which allows us to drop the "self." notation for E_Sed
            within this method.

        Returns:

        '''

        for i in range(1, self.elev_counter + 1):
            self.E_V_A[1][i] += Load_Removed * (1/self.density_SS) * (E_Sed[1][i] - E_Sed[1][i-1])
            # Loop through all elements in list at higher elevations to be sure their cumulative volume is reduced.
            # Must loop through entire list, not just up to the elev_counter (usually = active storage elevation).
            for x in range(i + 1, len(self.E_V_A[1])):
                # Add load that was just added to elevation i to all elevations x above i, since this is a cumulative
                #  storage matrix.
                self.E_V_A[1][x] += Load_Removed * (1/self.density_SS) * (E_Sed[1][i] - E_Sed[1][i-1])

    def DPS_Set_WSE_Target(self, t):

        '''

        Method to determine a reservoir's daily release/target water level according to a specified operating policy.

        Before reservoir's mass balance is simulated for this time step, take externally specified reservoir
        operating policy (if user indicates such a RBF-based policy has been externally specified using DPS or other
        methods) and determine the required release/water level given system state variable values of interest. The
        method Set_DPS_Operations_Goal() is what actually determines the daily operational decision.

        Adapted from code originally written by Julianne Quinn, Cornell University (jdq2101@gmail.com)

        '''

        # Load inputs to DPS (RBF) policy function into an array (inputs), then normalize that array in range [0,
        # 1] as normInputs. Default assumption: Inputs are all at time t (beginning of time step).
        for i in range(self.op_policy_params['num_inputs']):
            try:
                if self.op_policy_params['ordered_input_names'][i] == 'current_date':
                    self.DPS_inputs[i] = getattr(self, self.op_policy_params['ordered_input_names'][i]).timetuple().tm_yday
                else:
                    delta_t = self.op_policy_params['ordered_input_time_index'][i]  # Will be 0 if user did not specify
                    if delta_t == 0:
                        self.DPS_inputs[i] = getattr(self, self.op_policy_params['ordered_input_names'][i])[t]
                    else:
                        if t > 0:
                            self.DPS_inputs[i] = getattr(self, self.op_policy_params['ordered_input_names'][i])[t-delta_t]
                        else:
                            # There is no (t-1) period on day 1.
                            self.DPS_inputs[i] = getattr(self, self.op_policy_params['ordered_input_names'][i])[t]

            except AttributeError:
                print("Error: User has specified DPS input variable that does not exist for %s" % self.name)

        # Temporarily set either (1) flow release, (2) target elevation, or (3) target water level. Convert to
        # elevation or storage target depending on decision type. Constraints may require reset.
        u = direct_policy_search.DPS_Policy_Decision(self.DPS_inputs, self.op_policy_params)  # Release decision
        self.DPS_Decision[t] = u[0]
        try:
            # Current valid entries for output variables: Storage_Target_Release_Goal, ELEVATION_TARGET, STORAGE_TARGET
            # Only set these targets if dealing with elevation or storage target. For flow, skip this.
            if self.op_policy_params['ordered_output_names'][0] == 'ELEVATION_TARGET':
                getattr(self, self.op_policy_params['ordered_output_names'][0])[t+1] = self.DPS_Decision[t]
                self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                self.ELEVATION_TARGET[t+1])
            elif self.op_policy_params['ordered_output_names'][0] == 'STORAGE_TARGET':
                getattr(self, self.op_policy_params['ordered_output_names'][0])[t+1] = self.DPS_Decision[t]
                self.ELEVATION_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation",
                                                                self.STORAGE_TARGET[t+1])
            elif self.op_policy_params['ordered_output_names'][0] == 'Storage_Target_Release_Goal':
                getattr(self, self.op_policy_params['ordered_output_names'][0])[t] = self.DPS_Decision[t]  # t, not t+1
                # Set a storage target associated with this release goal
                self.STORAGE_TARGET[t+1] = min(max(0, self.S[t] + (self.Q_in[t] - self.Evap[t] - self.Storage_Target_Release_Goal[t]) * self.dt * 86400 - self.Settled_volume[t]), Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.Act_Stor_Elev))
                self.ELEVATION_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation",
                                                                self.STORAGE_TARGET[t+1])
                self.Storage_Target_Release_Goal[t] = 0  # Reset to zero, as it will be reset later to same value.
        except AttributeError:
            print("Error: User has specified DPS output variable other than storage, WSE, or release flor rate for "
                  "%s" % self.name)

    def Sediment_Management_Daily_Import(self, t):
        '''
        Purpose: method loops through all sediment management types that exist for this reservoir and run the daily
        initiation check.

        Args:
            t: Current time period (integer, value range: [0,T])

        Returns:

        '''

        try:
            self.flushing_today = self.Operating_Policy["Flushing"].Current
            self.flushing_yesterday = self.Operating_Policy["Flushing"].previous_flushing
        except KeyError:
            self.flushing_today = 0
            self.flushing_yesterday = 0
        try:
            self.flushing_group_drawn_down_today = self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down
        except KeyError:
            self.flushing_group_drawn_down_today = 0
        try:
            self.sluicing_today = self.Operating_Policy["Sluicing"].Current
            self.sluicing_yesterday = self.Operating_Policy["Sluicing"].previous_sluicing
        except KeyError:
            self.sluicing_today = 0
            self.sluicing_yesterday = 0
        try:
            self.density_current_venting_today = self.Operating_Policy["Density Current Venting"].Current
            self.density_current_venting_yesterday = self.Operating_Policy["Density Current Venting"].previous_density_current_venting
        except KeyError:
            self.density_current_venting_today = 0
            self.density_current_venting_yesterday = 0

        # Loop through operating policy dictionary elements to call each sediment management type's daily initiation
        # check, which checks reservoir conditions (e.g., inflow) during time step t, before other reservoirs calcs
        # are done, to see if the particular sediment management type will be performed during time step t.
        if t != 0:
            if self.Q_jct[t-1] != 0:
                prev_day_inflow = self.Q_jct[t-1]  # Use upstream junction flow as trigger
            else:
                prev_day_inflow = self.Q_in[t-1]  # Use upstream junction flow as trigger
        else:
            if self.Q_jct[t] != 0:
                prev_day_inflow = self.Q_jct[t]  # Use upstream junction flow as trigger
            else:
                prev_day_inflow = self.Q_in[t]  # Use upstream junction flow as trigger

        for policy in self.Operating_Policy:
            if policy == "Flushing":
                self.Operating_Policy[policy].daily_initiation_check(self.current_date, prev_day_inflow,
                                                                     self.Flushing_Dictionary, self.sluicing_today,
                                                                     self.density_current_venting_today)
                self.flushing_today = self.Operating_Policy[policy].Current  # Keep track of activation of any sed mgmt operations.
            elif policy == "Dredging":
                self.BS_W[t+1] = self.Operating_Policy[policy].daily_initiation_check(self.current_date, t, self.BS_W[t+1])  # Update BS_W of the reservoir. Existing value is returned if no dredging occurs.
                Reservoir.Reservoir_Volume_Addition(self, self.Operating_Policy["Dredging"].Sediment_Load_Removed_Daily, self.E_Sed)  # Adjust E-V-A curve to reflect water storage volume that's been added due to sediment removal.
                Reservoir.Active_Dead_Storage_Capacity_Calculation(self, t)  # Compute new active storage capacity since dredged mass has been computed.
            elif policy == "Sluicing":
                self.Operating_Policy[policy].daily_initiation_check(self.current_date, self.Q_in[t], self.flushing_today, self.density_current_venting_today)
                self.sluicing_today = self.Operating_Policy[policy].Current  # Keep track of activation of any sed mgmt operations.
            elif policy == "Density Current Venting":
                [self.DCV_Flow] = self.Operating_Policy[policy].daily_initiation_check(prev_day_inflow, self.V_in[t], self.SS_W_in[t], self.flushing_today, self.sluicing_today)
                self.density_current_venting_today = self.Operating_Policy[policy].Current  # Keep track of activation of any sed mgmt operations.
        # Now that all daily initiation checks have been run establish the dominant technique, establish which sediment management techniques will be given priority in the event that two techniques need to be simulated
        # Send values back to the operating policy, in case any have been re-set.
        if self.flushing_yesterday == 1 and self.flushing_today == 1:
            if self.sluicing_today == 1:
                # Change sluicing start and stop days in table to next day, as sluicing must be delayed due to flushing
                self.Operating_Policy["Sluicing"].Specs[self.current_date]["Sluicing End Date"] += timedelta(days = 1)  # Delay end date by one day.
                self.Operating_Policy["Sluicing"].Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Update dictionary key so it reflects delay of sluicing by one day, as inflow conditions have not been met.
                self.Operating_Policy["Sluicing"].Current = 0
            try:
                self.Operating_Policy["Density Current Venting"].Current = 0
            except KeyError:
                pass
        elif self.sluicing_yesterday == 1 and self.sluicing_today == 1:
            if self.flushing_today == 1:
                self.Operating_Policy["Flushing"].Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Skip to next day; see explanation above for self.Non_Flushing_Group_Drawn_Down = 1            try:
                self.Operating_Policy["Flushing"].Current = 0
            try:
                self.Operating_Policy["Density Current Venting"].Current = 0
            except KeyError:
                pass
        elif self.density_current_venting_yesterday == 1 and self.density_current_venting_today == 1:
            if self.flushing_today == 1:
                self.Operating_Policy["Flushing"].Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Skip to next day; see explanation above for self.Non_Flushing_Group_Drawn_Down = 1            try:
                self.Operating_Policy["Flushing"].Current = 0
            if self.sluicing_today == 1:
                 # Change sluicing start and stop days in table to next day, as sluicing must be delayed due to flushing
                 self.Operating_Policy["Sluicing"].Specs[self.current_date]["Sluicing End Date"] += timedelta(days = 1)  # Delay end date by one day.
                 self.Operating_Policy["Sluicing"].Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Update dictionary key so it reflects delay of sluicing by one day, as inflow conditions have not been met.
                 self.Operating_Policy["Sluicing"].Current = 0

    def Elev_Stor_Target(self, t):
        '''
        Purpose: Method sets elevation and storage targets. Highly customizable, as it’s modified depending on what,
        if any, sediment management is taking place.

        Args:
            t: Current time period (integer, value range: [0,T])

        Returns:

        '''

        self.Environmental_Release_Goal[t] = self.Environmental_Flow_Data[self.current_date.month - 1]  # Establish the environmental release goal.
        self.hydro_multiplier = 1  # Power will be produced. Will be reset = 0 only if reservoir is sluicing and power is bring produced during sluicing.
        # If reservoir isn't full of sediment, and a sediment management activity is to occur, then set targets based on that activity.
        if self.reservoir_is_full == 0:
            try:
                if self.Operating_Policy["Flushing"].Current == 1 or self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down == 1:
                    # Flushing is happening (or reservoir is drawn down due to its inclusion in a group of reservoir that are to all be empty at the same time for flushing).
                    self.ELEVATION_TARGET[t+1] = self.Operating_Policy["Flushing"].Elev_Stor_Target(self.water_surface_elevation[t])
                    # If Flushing is to occur, STORAGE_TARGET must be set = 0, unless the storage target reflects an adjustment to control drawdown rate.
                    if self.ELEVATION_TARGET[t+1] == self.Orifices.Low_Level_Outlet_Invert_Elevation:
                        self.STORAGE_TARGET[t+1] = 0
                    else:
                        self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.ELEVATION_TARGET[t+1])  # User provided elevation targets for this reservoir. Convert elevation targets into storage targets.
                    return  # Do not execute other Storage/elevation target methods if this  method is accessed.
            except KeyError:
                pass
            try:
                if self.Operating_Policy["Sluicing"].Current == 1:
                    # Sluicing is happening, or reservoir is still partially drawn down from sluicing, so post_sluicing = 1
                    # Note: Post-sluicing is checked later.
                    if self.Operating_Policy["Sluicing"].Specs["Current Specs"]["Target drawdown elevation"] != "Inflow-based Sluicing":
                        self.ELEVATION_TARGET[t+1] = self.Operating_Policy["Sluicing"].Elev_Stor_Target(self.water_surface_elevation[t], ELEVATION_TARGET_PRE_SET = None)
                        #Proceed with regular water level determination.
                        self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.ELEVATION_TARGET[t+1])  # User provided elevation targets for this reservoir. Convert elevation targets into storage targets.
                        self.hydro_multiplier = self.Operating_Policy["Sluicing"].hydro_multiplier  # Will be reset = 0 only if reservoir is sluicing and power is being produced during sluicing.
                        return
                    else:
                        self.hydro_multiplier = self.Operating_Policy["Sluicing"].hydro_multiplier  # Will be reset = 0 only if reservoir is sluicing and power is being produced during sluicing.
            except KeyError:
                pass
            try:
                try:
                    if self.Operating_Policy["Sluicing"].Current == 1 and self.Operating_Policy["Sluicing"].Specs["Current Specs"]["Target drawdown elevation"] is "Inflow-based Sluicing":
                        return
                except KeyError:
                    pass
                if self.Operating_Policy["Density Current Venting"].Current == 1:
                    # Density Current venting is happening
                    [self.ELEVATION_TARGET[t+1], self.STORAGE_TARGET[t+1]] = self.Operating_Policy["Density Current Venting"].Elev_Stor_Target(self.water_surface_elevation[t], self.STORAGE_TARGET[t+1], self.MW_Capacity, self.Turbine_Efficiency, self.Q_in[t], self.S[t], self.V_in[t], self.SS_W_in[t], self.SS_C[t], self.Orifices.Capacity_low_level_Outlet[t], self.Orifices.Capacity_Controlled_Outlet[t], self.Orifices.Capacity_mid_level_Outlet[t], self.Orifices.Capacity_Spillway_Outlet[t], self.Orifices.Capacity_Hydropower_Outlet[t], self.Evap[t], self.Settled_volume[t], self.Dam_Tailwater_Elevation[t], self.Reservoir_Type, self.E_V_A)
                    return  # Do not execute other Storage/elevation target methods if this  method is accessed.
            except KeyError:
                pass
            # Still have not exited this method with "return", so either no sediment management exists, or none is being practiced right now
            if self.Reservoir_Operations_Goal_ID == 1:
                self.ELEVATION_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation", self.STORAGE_TARGET[t+1])  # User provided storage targets for this reservoir. Convert storage targets into elevation targets.
            elif self.Reservoir_Operations_Goal_ID == 2:
                self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.ELEVATION_TARGET[t+1])  # User provided elevation targets for this reservoir. Convert elevation targets into storage targets.
            elif self.Reservoir_Operations_Goal_ID == 3 or self.Reservoir_Operations_Goal_ID == 4:
                if self.Reservoir_Operations_Goal_ID == 3:
                    # Reservoir will be operated to achieve energy targets
                    self.Energy_Release_Goal[t] = self.Daily_Energy_Targets(self.current_date.month - 1) / (24 * 0.00981 * self.Turbine_Efficiency * max(0, self.water_surface_elevation[t] - self.Dam_Tailwater_Elevation[t]))  # Set a storage target, which will be minimum of the storage that would be achieved if energy flow were released and storage level when elevation is at the top of the active storage. Purpose is to prevent filling into flood zone when reservoir is operated for energy targets.
                    self.Storage_Target_Release_Goal[t] = self.Energy_Release_Goal[t]
                elif self.Reservoir_Operations_Goal_ID == 4:
                    # Reservoir will be operated to achieve daily release targets
                    self.Storage_Target_Release_Goal[t] = self.Daily_Release_Target(self.current_date.month - 1)
                    # Make sure that Environmental Release Goal is set if one exists.
                    if self.Reservoir_Type == "Diversion" or self.Reservoir_Type == "Diversion and Power":
                        # Reservoir diverts water away from downstream channel, so need to add environmental flow release to energy-based release
                        self.Storage_Target_Release_Goal[t] += self.Environmental_Release_Goal[t]
                    else:
                        self.Storage_Target_Release_Goal[t] = max(self.Storage_Target_Release_Goal[t], self.Environmental_Release_Goal[t])
                self.STORAGE_TARGET[t+1] = min(max(0, self.S[t] + (self.Q_in[t] - self.Evap[t] - self.Storage_Target_Release_Goal[t]) * self.dt * 86400 - self.Settled_volume[t]), Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.Act_Stor_Elev))
                self.ELEVATION_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation", self.STORAGE_TARGET[t+1])  # Store a time 0 elevation target corresponding to specified time zero storage target.
            elif self.Reservoir_Operations_Goal_ID == 6:
                # User wishes to define a target reservoir water surface elevation as a function of reservoir inflow.
                self.ELEVATION_TARGET[t+1] = Matrix_Interpolation(self.name, self.wse_inflow_target, "flow",
                                                         "elevation", self.Q_in[t])
                # Determine storage target corresponding to elevation taraget
                self.STORAGE_TARGET[t + 1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                                  self.ELEVATION_TARGET[t + 1])

            # If post sluicing is going on, then need to call that method. Being done here because post-sluicing requires that an elevation target has been pre-set.
            try:
                if self.Operating_Policy["Sluicing"].Post_Sluicing == 1:
                    if self.Operating_Policy["Sluicing"].Specs["Current Specs"]["Target drawdown elevation"] != "Inflow-based Sluicing":
                        self.ELEVATION_TARGET[t+1] = self.Operating_Policy["Sluicing"].Elev_Stor_Target(self.water_surface_elevation[t], self.ELEVATION_TARGET[t+1])
                        self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.ELEVATION_TARGET[t+1])  # User provided elevation targets for this reservoir. Convert elevation targets into storage targets.
            except KeyError:
                pass  # Sluicing doesn't exist for reservoir.

    def Set_total_release_goal(self, t):
        '''
        Purpose: Based on established storage targets, set storage target release goal, and then determine a final
        Q_out based on release goal and outlet capacity constraints

        Compute an initial estimate of what release must be made to meet storage target.

        Args:
            t: Current time period (integer, value range: [0,T])

        Returns:

        '''

        # Should include options for this module in future (e.g., release schedule)

        if self.reservoir_is_full == 0:
            if self.Reservoir_Operations_Goal_ID in [1, 2, 3, 4, 5, 6]:
                if self.Storage_Target_Release_Goal[t] == 0:
                    # If DPS being used, prevent ramping up and down of policies by resetting the Storage target
                    # release goal, and associated elevation/storage targets, if the elevation differential from one
                    # day to the next exceeds 2 m (in either direction).
                    if self.Reservoir_Operations_Goal_ID == 5:
                        reset_elev = 0
                        if (self.ELEVATION_TARGET[t+1] - self.water_surface_elevation[t]) > 2.0:
                            reset_elev = 1
                            self.ELEVATION_TARGET[t+1] = self.water_surface_elevation[t] + 2
                        elif (self.water_surface_elevation[t] - self.ELEVATION_TARGET[t+1]) > 2.0:
                            reset_elev = 1
                            self.ELEVATION_TARGET[t+1] = self.water_surface_elevation[t] - 2
                        if reset_elev == 1:
                            # Reset Storage target
                            self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation",
                                                                            "storage", self.ELEVATION_TARGET[t+1])
                    # Only set release goal if release hasn't already been specified.
                    self.S_delta[t] = self.S[t] + (self.Q_in[t] - self.Evap[t]) * self.dt * 86400 - self.Settled_volume[t] - self.STORAGE_TARGET[t+1]
                    if self.S_delta[t] < 0:
                        self.Storage_Target_Release_Goal[t] = 0
                    else:
                        self.Storage_Target_Release_Goal[t] = self.S_delta[t] / (self.dt * 86400)

            # Using the output from method set total release goal(), determine the total Q_out (release flow rate, m3/s) taking into consideration  all of the outlet’s capacities.
            if self.Reservoir_Type == "Power":
                # Power generation only. CASE in which MW capacity of plant cannot be exceeded to avoid wasting energy (you would need to move the Max_MW_Turbine_Flow to above this If/ElseIf loop for this to work: self.Q_out[t] = Min(Max(self.Storage_Target_Release_Goal[t], Environmental_Release_Goal[t]), self.Max_MW_Turbine_Flow[t], self.Orifices.Capacity_Hydropower_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t])
                v1 = min(max(self.Storage_Target_Release_Goal[t], self.Environmental_Release_Goal[t]),
                         self.hydro_multiplier * self.Orifices.Capacity_Hydropower_Outlet[t] +
                         self.Orifices.Capacity_Spillway_Outlet[t] + self.Orifices.Capacity_low_level_Outlet[t] +
                         self.Orifices.Capacity_mid_level_Outlet[t],
                         max(0, (self.S[t] + self.V_in[t] - self.Evap[t] * self.dt * 86400) / 86400))
                current_max_storage =Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage", self.Act_Stor_Elev)
                min_required_release = max(0, (self.S[t] + self.V_in[t] - self.Evap[t] * self.dt * 86400 - current_max_storage) / 86400)
                self.Q_out[t] = max(v1, min_required_release)
            elif self.Reservoir_Type == "Diversion":
                # Diversion only
                self.Q_out[t] = min(max(self.Storage_Target_Release_Goal[t], min(self.Orifices.Capacity_Controlled_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t] + self.Orifices.Capacity_low_level_Outlet[t] + self.Orifices.Capacity_mid_level_Outlet[t], self.Environmental_Release_Goal[t])), self.Orifices.Capacity_Controlled_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t] + self.Orifices.Capacity_Diverted_Outlet[t] + self.Orifices.Capacity_low_level_Outlet[t] + self.Orifices.Capacity_mid_level_Outlet[t], max(0, (self.S[t] + self.V_in[t] - self.Evap[t] * self.dt * 86400) / 86400))
            elif self.Reservoir_Type == "Diversion and Power":
                # Power generation and diversion. CASE in which MW capacity of plant cannot be exceeded to avoid wasting energy: self.Q_out[t] = Min(Max(self.Storage_Target_Release_Goal[t], self.Environmental_Release_Goal[t]), self.Max_MW_Turbine_Flow[t], Capacity_Controlled_Outlet[t](index_res(n)) + self.Orifices.Capacity_Hydropower_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t])
                self.Q_out[t] = min(max(self.Storage_Target_Release_Goal[t], min(self.Orifices.Capacity_Controlled_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t] + self.Orifices.Capacity_low_level_Outlet[t] + self.Orifices.Capacity_mid_level_Outlet[t], self.Environmental_Release_Goal[t])), self.Orifices.Capacity_Controlled_Outlet[t] + self.hydro_multiplier * self.Orifices.Capacity_Hydropower_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t] + self.Orifices.Capacity_low_level_Outlet[t] + self.Orifices.Capacity_mid_level_Outlet[t], max(0, (self.S[t] + self.V_in[t] - self.Evap[t] * self.dt * 86400) / 86400))
            elif self.Reservoir_Type == "Storage":
                # Storage/Release capabilities only
                self.Q_out[t] = min(max(self.Storage_Target_Release_Goal[t], self.Environmental_Release_Goal[t]), self.Orifices.Capacity_Controlled_Outlet[t] + self.Orifices.Capacity_Spillway_Outlet[t] + self.Orifices.Capacity_low_level_Outlet[t] + self.Orifices.Capacity_mid_level_Outlet[t], max(0, (self.S[t] + self.V_in[t] - self.Evap[t] * self.dt * 86400) / 86400))

    def Distribute_Outflow_Among_Outlets(self, t):
        '''
        This method takes a computed release self.Q_out[t] from method Set_total_release_goal and distributes that flow among the outlet works present at the dam.

        Args:
            t: Current time period (integer, value range: [0,T])

        Returns:

        '''
        if self.reservoir_is_full == 0:
            # If sluicing exists and is happening, call sluicing methods to allocate flow to outlets:
            if self.sluicing_today == 1:
                # Sluicing is happening.
                Input_Array = [self.Orifices.Capacity_Controlled_Outlet[t], self.Orifices.Capacity_mid_level_Outlet[t], self.Orifices.Capacity_low_level_Outlet[t], self.Orifices.Capacity_Diverted_Outlet[t], self.Orifices.Capacity_Hydropower_Outlet[t], self.Orifices.Capacity_Spillway_Outlet[t]]
                [self.Orifices.Q_controlled[t], self.Orifices.Q_mid_level_outlet[t], self.Orifices.Q_diversion[t], self.Orifices.Q_turbines[t], self.Orifices.Q_low_level_outlet[t], self.Orifices.Q_overflow[t], self.Orifices.Q_downstream[t]] = self.Operating_Policy["Sluicing"].Distribute_Outflow_Among_Outlets(Input_Array, t, self.Reservoir_Type, self.Environmental_Release_Goal[t], self.Q_out[t], self.Orifices.mid_level_outlet)

            # If DCV exists and is happening, call DCV methods to allocate flow to outlets:
            if self.density_current_venting_today == 1:
                # Density Current Venting is happening.
                Input_Array = [self.Orifices.Capacity_Controlled_Outlet[t], self.Orifices.Capacity_mid_level_Outlet[t], self.Orifices.Capacity_low_level_Outlet[t], self.Orifices.Capacity_Diverted_Outlet[t], self.Orifices.Capacity_Hydropower_Outlet[t], self.Orifices.Capacity_Spillway_Outlet[t]]
                [self.Orifices.Q_controlled[t], self.Orifices.Q_mid_level_outlet[t], self.Orifices.Q_diversion[t], self.Orifices.Q_turbines[t], self.Orifices.Q_low_level_outlet[t], self.Orifices.Q_overflow[t], self.Orifices.Q_downstream[t]] = self.Operating_Policy["Density Current Venting"].Distribute_Outflow_Among_Outlets(Input_Array, self.Reservoir_Type, self.Environmental_Release_Goal[t], self.Q_out[t], self.Q_in[t])

            # If sluicing or DCV is not proceeding, the sections of code below will allocate flow to outlets.
            # Even if Sluicing or DCV is occurring, flow allocation for some outlets is done the same regardless of whether these methods are proceeding, and is dependent on reservoir type. That code must be executed regardless of sediment management type currently ongoing.
            if self.Reservoir_Type in ["Diversion", "Diversion and Power"]:
                if self.sluicing_today == 0 and self.density_current_venting_today == 0:
                    # No targets set yet, so set them now.
                    self.Orifices.Q_controlled[t] = min(self.Q_out[t], self.Environmental_Release_Goal[t], self.Orifices.Capacity_Controlled_Outlet[t])
                    # Sluicing/venting are not being performed, so send water to diversion/turbine outlets, depending on which ones exist.
                    if self.Reservoir_Type == "Diversion":
                        self.Orifices.Q_diversion[t] = min(max(self.Q_out[t] - self.Orifices.Q_controlled[t], 0), self.Orifices.Capacity_Diverted_Outlet[t])
                        # self.Orifices.Q_turbines[t] = 0
                    else:
                        self.Orifices.Q_diversion[t] = min(max(self.Q_out[t] - self.Orifices.Q_controlled[t], 0), self.Orifices.Capacity_Hydropower_Outlet[t])
                        self.Orifices.Q_turbines[t] = self.Orifices.Q_diversion[t]

                # Code that must be executed regardless of sediment management scenario:
                if (self.Q_out[t] - (self.Orifices.Q_controlled[t] + self.Orifices.Q_diversion[t] + self.Orifices.Q_mid_level_outlet[t] + self.Orifices.Q_low_level_outlet[t])) > 0:
                    # Some flow remains and could possibly be distributed to the controlled outlet (assuming it was not fully allocated due to lower environmental release goal)
                    self.Orifices.Q_controlled[t] = min(self.Orifices.Q_controlled[t] + (self.Q_out[t] - (self.Orifices.Q_controlled[t] + self.Orifices.Q_diversion[t] + self.Orifices.Q_mid_level_outlet[t])), self.Orifices.Capacity_Controlled_Outlet[t])
                self.Orifices.Q_downstream[t] = max(self.Q_out[t] - self.Orifices.Q_diversion[t], 0)
                try:
                    if self.Operating_Policy["Flushing"].Current == 1 or self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down == 1:
                        self.Orifices.Q_low_level_outlet[t] = min(max(self.Q_out[t] - self.Orifices.Q_controlled[t] - self.Orifices.Q_diversion[t] - self.Orifices.Q_overflow[t] - self.Orifices.Q_mid_level_outlet[t], 0), self.Orifices.Capacity_low_level_Outlet[t])
                except KeyError:
                    pass  # Flushing doesn't exist for this reservoir
                self.Orifices.Q_overflow[t] = self.Orifices.Q_overflow[t] + min(max(self.Q_out[t] - self.Orifices.Q_diversion[t] - self.Orifices.Q_controlled[t] - self.Orifices.Q_mid_level_outlet[t] - self.Orifices.Q_low_level_outlet[t] - self.Orifices.Q_overflow[t], 0), self.Orifices.Capacity_Spillway_Outlet[t])
            elif self.Reservoir_Type == "Power":
                if self.sluicing_today == 0 and self.density_current_venting_today == 0:
                    # No targets set yet, so set them now.
                    # Sluicing and DCV are not being performed, so send water to turbine outlets (normal operations)
                    if self.Q_out[t] <= self.Orifices.Capacity_Hydropower_Outlet[t]:
                        # No spilling over the spillway needs to occur
                        self.Orifices.Q_turbines[t] = min(self.Q_out[t], self.Orifices.Capacity_Hydropower_Outlet[t])
                        self.Orifices.Q_overflow[t] = 0
                    else:
                        self.Orifices.Q_turbines[t] = self.Orifices.Capacity_Hydropower_Outlet[t]
                        # In a reservoir that only generates hydropower, if the release flow rate is higher than the hydropower plant capacity, this flow is probably being spilled through spillway, but need to see if water elevation is high enough for spilling
                else:
                    pass  # Targets already set

                # Code that must be executed either way:
                self.Orifices.Q_downstream[t] = self.Q_out[t]  # All water exiting the reservoir will go downstream for this type of reservoir
                self.Orifices.Q_diversion[t] = 0  # This reservoir does not have a diversion
                try:
                    if self.Operating_Policy["Flushing"].Current == 1 or self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down == 1:
                        self.Orifices.Q_low_level_outlet[t] = min(max(self.Q_out[t] - self.Orifices.Q_controlled[t] - self.Orifices.Q_turbines[t] - self.Orifices.Q_overflow[t], 0), self.Orifices.Capacity_low_level_Outlet[t])
                except KeyError:
                    pass  # Flushing doesn't exist for this reservoir
                self.Orifices.Q_overflow[t] = self.Orifices.Q_overflow[t] + min(max(self.Q_out[t] - self.Orifices.Q_turbines[t] - self.Orifices.Q_controlled[t] - self.Orifices.Q_mid_level_outlet[t] - self.Orifices.Q_low_level_outlet[t] - self.Orifices.Q_overflow[t], 0), self.Orifices.Capacity_Spillway_Outlet[t])
            elif self.Reservoir_Type == "Storage":
                if self.sluicing_today == 0 or self.density_current_venting_today == 0:
                    # No targets set yet, so set them now.
                    # Neither sluicing nor density current venting are being performed, so send water to regular outlets (normal operations)
                    if self.Q_out[t] <= self.Orifices.Capacity_Controlled_Outlet[t]:
                        self.Orifices.Q_controlled[t] = min(self.Q_out[t], self.Orifices.Capacity_Controlled_Outlet[t])
                        self.Orifices.Q_overflow[t] = 0
                    else:
                        self.Orifices.Q_controlled[t] = self.Orifices.Capacity_Controlled_Outlet[t]
                else:
                    pass  # Targets already set

                # Code that must be executed either way:
                try:
                    if self.Operating_Policy["Flushing"].Current == 1 or self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down == 1:
                        self.Orifices.Q_low_level_outlet[t] = min(max(self.Q_out[t] - self.Orifices.Q_controlled[t] - self.Orifices.Q_turbines[t] - self.Orifices.Q_overflow[t] - self.Orifices.Q_mid_level_outlet[t], 0), self.Orifices.Capacity_low_level_Outlet[t])
                except KeyError:
                    pass  # Flushing doesn't exist for this reservoir.
                self.Orifices.Q_overflow[t] = self.Orifices.Q_overflow[t] + min(max(self.Q_out[t] - self.Orifices.Q_controlled[t] - self.Orifices.Q_mid_level_outlet[t] - self.Orifices.Q_low_level_outlet[t] - self.Orifices.Q_overflow[t], 0), self.Orifices.Capacity_Spillway_Outlet[t])
            else:
                pass  # Currently there are no other reservoir types than the four covered above.

    def Terminal_WL_Calc(self, t):
        '''
        Determines the highest possible water level for a given predicted daily inflow rate, assuming there is an
        uncontrolled entry to bypass channel at the upstream end of the reservoir. This is a highly specific
        scenario; use with caution, as it assumes reservoir is not fully controllable by the dam's physical outlets
        and gates. Instead, water can spill from the upstream end of the reservoir.

        Args:
            t: Current time period (integer, value range: [0,T])

        Returns:

        '''

        self.Terminal_WSE[t+1] = Matrix_Interpolation(self.name, self.TERMINAL_WL_Rating_Curve, "flow", "elevation",
                                                      self.Q_in[t])

        # Regardless of the type of operations at the reservoir, entry to this method requires that the elevation
        # target be set to Terminal_WSE, and a new storage target be calculated.

        if self.ELEVATION_TARGET[t+1] > self.Terminal_WSE[t+1]:
            self.ELEVATION_TARGET[t+1] = self.Terminal_WSE[t+1]
            self.STORAGE_TARGET[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "elevation", "storage",
                                                            self.ELEVATION_TARGET[t+1])

    def Outflow_Storage_Check(self, t):
        '''
        Purpose: This method sets outflow volume, but then adjusts outflow volume (as well as flow/storage) if any
        errors have occurred in setting this value. This method also sets the final water surface elevation for the
        time step, as this is the last place where self.S[t+1] is updated in code.

        :param t: Current time period (integer, value range: [0,T])
        :return:
        '''

        self.V_out[t] = self.Q_out[t] * self.dt * 86400
        if self.S[t+1] <= 0:
            self.V_out[t] = max(0, self.S[t] + self.V_in[t] - self.Evap[t] * 86400)
            self.Q_out[t] = self.V_out[t] / (86400 * self.dt)
            self.S[t+1] = 0
        # Determine final water level elevation associated with final water storage, considering the impact of
        # sediment volume/elevation
        self.water_surface_elevation[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "elevation",
                                                                 self.S[t+1])
        # Determine final water level elevation associated with final water storage, considering the impact of
        # sediment volume/elevation
        self.water_surface_area[t+1] = Matrix_Interpolation(self.name, self.E_V_A, "storage", "area", self.S[t+1])

    def SedMgmt_Main(self, t):
        '''
        Purpose: Calls main mass balance routines for sediment management types if any exist and are being simulated.

        Args:
            t: Current time period (integer, value range: [0,T])

        Returns:

        '''

        try:
            if self.Operating_Policy["Flushing"].Current == 1:
                [self.BS_W[t+1], self.SS_W_out[t]] = self.Operating_Policy["Flushing"].SedMgmt_Main(t, self.S[t+1], self.E_V_A, self.water_surface_elevation[t], self.Q_out[t], self.BS_W[t+1], self.density_SS)
                if self.Operating_Policy["Flushing"].flushed_load[t] > 0:
                    self.res_flushed_load[t] = self.Operating_Policy["Flushing"].flushed_load[t]
                if self.Operating_Policy["Flushing"].Flushing_occurred == 1:
                    # Adjust E-V-A curve to reflect water storage volume that's been added due to sediment removal.
                    Reservoir.Reservoir_Volume_Addition(self,
                                                        self.Operating_Policy["Flushing"].Flushing_Load_Removed_Daily,
                                                        self.E_Sed)
                    # Compute new active storage capacity since flushed  mass has been computed.
                    Reservoir.Active_Dead_Storage_Capacity_Calculation(self, t)
        except KeyError:
            pass  # Flushing not defined for this reservoir
        try:
            if self.Operating_Policy["Density Current Venting"].Current == 1:
                [self.TE_avg[t], self.Settled_mass[t], self.Settled_volume[t], self.BS_W[t+1], self.BS_V[t+1], self.SS_W_out[t], self.DCV_Flow] = self.Operating_Policy["Density Current Venting"].SedMgmt_Main(self.BS_W[t+1], self.SS_W_in[t], self.V_in[t], self.Q_low_level_outlet[t], self.density_SS, self.capacity_active_reservoir[t+1], self.capacity_dead_reservoir[t+1])
        except KeyError:
            pass  # DCV not defined for this reservoir

    def Sediment_Mass_Balance(self, t):
        '''
        Purpose: Determines final sediment mass discharge from reservoir

        Accounts for all releases, including sediment management).

        :param t: Current time period (integer, value range: [0,T])
        :return:
        '''

        # Add sediment component released from outlet outflows. Components released via sediment management have
        # already been added. If DCV is happening, we don't want to add reservoir suspended sediment to the total
        # released sediment load.

        self.SS_W_out[t] += self.SS_C[t] * max(0, self.V_out[t] - self.DCV_Flow*self.dt*86400)
        # Now that we know the final storage self.S[t+1] (which includes evaporative losses), we can compute SS_C_adjusted
        if self.S[t+1] > 0:
            try:
                if self.flushing_today == 1:
                    # Flushing is occurring. Don't include flushed load in final reservoir storage concentration balance
                    outflow_replacement = self.SS_C[t] * max(0, self.V_out[t] - self.DCV_Flow*self.dt*86400)
                    self.SS_C_adjusted_for_evap[t] = (self.SS_W[t] + self.SS_W_in[t] - outflow_replacement -
                                                      self.Settled_mass[t]) / self.S[t+1]
                else:
                    # Flushing not currently taking place. Don't need to account for flushed load in concentration.
                    self.SS_C_adjusted_for_evap[t] = (self.SS_W[t] + self.SS_W_in[t] - self.SS_W_out[t] -
                                                      self.Settled_mass[t]) / self.S[t+1]
            except KeyError:
                # Flushing does not exist. Don't need to account for flushed load in storage concentration.
                self.SS_C_adjusted_for_evap[t] = (self.SS_W[t] + self.SS_W_in[t] - self.SS_W_out[t] -
                                                  self.Settled_mass[t]) / self.S[t+1]
        else:
            self.SS_C_adjusted_for_evap[t] = 0
        self.SS_W[t+1] = self.S[t + 1] * self.SS_C_adjusted_for_evap[t]

    def Hydropower_Calculations(self, t):

        # Determine dam tailwater elevation. Compute this regardless of whether it is a hydropower plant, as this value is used by other modules (e.g., DCV release component).
        if self.Dam_Tailwater_Elevation_1 is not None:
            self.Dam_Tailwater_Elevation[t+1] = self.Dam_Tailwater_Elevation_1  # User has specified constant value.
        else:
            # Determine elevation using tailwater rating curve.
            if self.re_calc_energy == 0:
                self.Dam_Tailwater_Elevation[t+1] = Matrix_Interpolation(self.name, self.Tailwater_Rating_Curve, "flow",
                                                                     "elevation", self.Q_out[t])
            else:
                # There is a junction downstream with multiple inflow elements. Use the flow at the junction to
                # determine tailwater elevation.
                self.Dam_Tailwater_Elevation[t+1] = Matrix_Interpolation(self.name, self.Tailwater_Rating_Curve, "flow",
                                                                     "elevation", self.Q_Downstream_Junction[t])
        # Compute hydropower output for period of interest
        if self.Reservoir_Type in ["Power", "Diversion and Power"]:
            # Reservoir has hydropower production capabilities
            self.hydropower_head_average[t] = max(0, 0.5 * ((self.water_surface_elevation[t+1] - self.Dam_Tailwater_Elevation[t]) + (self.water_surface_elevation[t+1] - self.Dam_Tailwater_Elevation[t+1])))
            # Determine the turbine flow that results in power production in excess of MW capacity of plant, given the expected head
            if self.MW_Capacity > 0:
                if self.hydropower_head_average[t] > 0:
                    self.Max_MW_Turbine_Flow[t] = (1 / 0.00981) * (self.MW_Capacity) * (1 / self.hydropower_head_average[t]) * (1 / self.Turbine_Efficiency)
                else:
                    self.Max_MW_Turbine_Flow[t] = 0
            else:
                self.Max_MW_Turbine_Flow[t] = 0
            # Compute hydropower/energy production
            if self.hydropower_head_average[t] >= self.Min_Net_Head:
                self.Hydropower_avg_MW[t] = min(self.MW_Capacity, 0.00981 * self.Orifices.Q_turbines[t] * self.hydropower_head_average[t] * self.Turbine_Efficiency)
            else:
                # Average hydraulic head for time period dropped below user-specified minimum. Hydropower/energy is
                # zeron on this day, though flow still released through turbines. Generators assumed to be off.
                self.Hydropower_avg_MW[t] = 0
            self.Hydropower_avg_MWH[t] = 24 * self.Hydropower_avg_MW[t]
            # computes flow rate of water released from the reservoir without generating hydropower (Q_spill).
            self.Orifices.Q_spill[t] = self.Orifices.Q_overflow[t] + self.Orifices.Q_controlled[t] + \
                                       self.Orifices.Q_low_level_outlet[t] + max(0, self.Orifices.Q_turbines[t] -
                                                                           self.Max_MW_Turbine_Flow[t])
            self.Curtailment_limited[t] = min(self.MW_Capacity, 0.00981 * self.Orifices.Q_spill[t] *
                                       self.hydropower_head_average[t] * self.Turbine_Efficiency)

            self.Curtailment_unlimited[t] = 0.00981 * self.Orifices.Q_spill[t] * self.hydropower_head_average[t] * self.Turbine_Efficiency
            # Compute the reservoir's number of hours that the powerhouse could theoretically operate at powerhouse capacity
            if self.Q_in[t] < self.Orifices.Capacity_Hydropower_Outlet[t]:
                #self.theoretical_peaking_capacity[t] = min(24,-(self.capacity_active_reservoir[
                # t+1]/self.Orifices.Capacity_Hydropower_Outlet[t])*(24/86400)/((self.Q_in[t]*3600*24/(
                # self.Orifices.Capacity_Hydropower_Outlet[t]*86400))-1))
                max_q_empty_act_stor = max(0, (
                self.S[t] + (self.Q_in[t] - self.Evap[t]) * self.dt * 86400 - self.Settled_volume[t] -
                self.capacity_dead_reservoir[t + 1]) / 86400)
                #self.theoretical_peaking_capacity[t] = min(24,
                # 24*(max_q_empty_act_stor/self.Orifices.Capacity_Hydropower_Outlet[t]))
                self.theoretical_peaking_capacity[t] = min(24,24*(self.Orifices.Q_turbines[t]/self.Orifices.Capacity_Hydropower_Outlet[t]))
            else:
                # Plant can peak for the entire day if inflow exceeds powerhouse capacity
                self.theoretical_peaking_capacity[t] = 24

        # Hydropower over/under flow threshold included here and in performance_measure_pre_calc section. Included
        # here in case hydropower calculations are being re-done in the event of tailwater depending on reservoir
        # operational decision.
        if self.Evaluate_Egg_Passage == 1:
            if self.Q_jct[t] >= self.larv_pass_flow_thold:
                self.Hydropower_avg_MWH_over_thold[t] = self.Hydropower_avg_MWH[t]
            else:
                self.Hydropower_avg_MWH_under_thold[t] = self.Hydropower_avg_MWH[t]

    def Performance_Measure_Pre_Calc(self, t):
        '''
        Purpose: This method does calculations helpful for evaluating performance of various measures with respect to
        reliability/resilience/vulnerability.

        More post-processing is done after simulation in separate file, but this prepares data for better use in
        post-processing.

        :param t: Current time period (integer, value range: [0,T])
        :return:
        '''
        if self.Evaluate_Egg_Passage == 1:
            self.Egg_Passage_Elevation_Target[t + 1] = Matrix_Interpolation(self.name, self.Egg_Passage_Rating_Curve,
                                                                            "flow", "elevation", self.Q_out[t])
            # Track the fraction of larvae actually released from a reservoir on a given day.
            # Accounts for mass balance, where in 1 unit of larvae is assumed to enter the reservoir every day,
            # or a fraction of 1 unit in the event where an upstream bypass exists that spills water before entering.
            # Hence, it is not a true mass balance (that includes concentrations and inflow/outflow rates). Instead,
            # it just attempts to account for larvae residence times.
            temp_larv_mass = self.larv_stor[t] + 1*self.fraction_Q_jct[t]  # Reservoir daily inflow of 1 unit of generic larvae mass.
            avail_stor = self.S[t] + self.V_in[t] - self.Evap[t] * self.dt * 86400 - self.Settled_volume[t]
            self.larv_mass_out[t] = temp_larv_mass*(self.V_out[t]/avail_stor)

            # In evaluating net passage, account for any bypassed larvae in calculating the daily fraction passing
            if (self.water_surface_elevation[t + 1] <= self.Egg_Passage_Elevation_Target[t + 1]) and (
                self.water_surface_elevation[t] <= self.Egg_Passage_Elevation_Target[t + 1]):
                # Successfully met reservoir's larvae passage velocity target
                self.larv_pass[t] = 1*self.fraction_Q_jct[t] + (1-self.fraction_Q_jct[t])
                self.larv_mass_out_surv[t] = self.larv_mass_out[t]
            else:
                # Failed to meet reservoir's larvae passage velocity target
                self.larv_pass[t] = 0*self.fraction_Q_jct[t] + (1-self.fraction_Q_jct[t])
                if t > 0 and self.larv_pass[t-1] == 0:
                    self.egg_pass_resilience[t] += 1
                self.larv_mass_out_surv[t] = 0
            self.larv_stor[t+1] = temp_larv_mass - self.larv_mass_out[t]
            self.larv_mass_out_surv_total[t] = self.larv_mass_out_surv[t] + (1-self.fraction_Q_jct[t])
            # Track passage success and residence time/reliability relative to inflow
            if self.Q_jct[t] >= self.larv_pass_flow_thold:
                self.larv_pass_over_thold[t] = self.larv_pass[t]
                self.larv_surv_over_thold[t] = self.larv_mass_out_surv_total[t]
                self.res_time_over_thold[t] = self.Residence_Time[t]
                self.Hydropower_avg_MWH_over_thold[t] = self.Hydropower_avg_MWH[t]
                if self.Residence_Time[t] > self.res_time_thold:
                    self.res_time_reliab_over_thold[t] = 0
                else:
                    self.res_time_reliab_over_thold[t] = 1
            else:
                self.larv_pass_under_thold[t] = self.larv_pass[t]
                self.larv_surv_under_thold[t] = self.larv_mass_out_surv_total[t]
                self.res_time_under_thold[t] = self.Residence_Time[t]
                self.Hydropower_avg_MWH_under_thold[t] = self.Hydropower_avg_MWH[t]
                if self.Residence_Time[t] > self.res_time_thold:
                    self.res_time_reliab_under_thold[t] = 0
                else:
                    self.res_time_reliab_under_thold[t] = 1

    def Element_Data_Transfer(self, t):
        '''
        Transfers data from reservoir element to other elements in the system whose methods may depend on the values
        of reservoir state variables.

        Method is performed at the end of every time step, storing information in the element_transfer dictionary
        that can then be sent to and used by other system objects. elem_xfer_output_dict is the primary vessel
        through which system elements (reservoirs, channels, etc.) communicate.
        '''

        # Passes current water surface elevation of reservoir to upstream junction if flow is split there and hence
        # depends on reservoir water level. Inflow element list for a reservoir can only contain a junction.
        for elems in self.Element_Sub_Dict["Inflow Elements"]:
            try:
                self.Element_Sub_Dict["Natural Bypass Water Level"] = {elems: self.water_surface_elevation[t+1]}
                self.elem_xfer_output_dict["Natural Bypass Water Level"] = {elems: self.water_surface_elevation[t+1]}
            except KeyError:
                pass

        # Passes diverted outlet flow to the downstream junction, with a flag that indicates this be distributed into the diversion
        # channel. Outflow element list for a reservoir can only contain a junction.
        if self.Reservoir_Type in ["Diversion", "Diversion and Power"]:
            for item_num in range(len(self.Element_Sub_Dict["Outflow Elements"])):
                self.elem_xfer_output_dict[self.Element_Sub_Dict["Outflow Elements"][item_num]] = {
                "Reservoir Diversion Flow": self.Orifices.Q_diversion[t]}
                #self.Element_Sub_Dict["Reservoir Diversion Flow"] = {elems: self.Orifices.Q_diversion[t]}
                #self.elem_xfer_output_dict["Reservoir Diversion Flow"] = {elems: self.Orifices.Q_diversion[t]}

    def Flow_Routing(self):
        return

    def Elevation_Calculation(self):
        return