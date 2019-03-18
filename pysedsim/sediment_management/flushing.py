# -*- coding: utf-8 -*-

'''

Module to define the Flushing class and asssociated methods.

'''

# Import relevant modules
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0)
from pysedsim.sediment_management.sediment_res_ops import Sediment_Res_Ops
import numpy as np
from pysedsim.data_processing.data_processing import *
from datetime import timedelta  # Used to add days/months/years to a datetime object
from datetime import datetime
from pysedsim.data_processing.matrix_interpolation import Matrix_Interpolation
import logging

class Flushing(Sediment_Res_Ops):
    '''
    This class represents the flushing technique for removing sediment from a reservoir. Class instances (objects)
    become stored in reservoir instances as part of its operating policy dictionary. The primary function of a
    dredging object is to perform a daily check to determine whether a reservoir will have flushing performed on a
    given day, and perform calculations for adjusting the reservoir's sediment mass balance if so.
    '''

    def __init__(self, element_name, T, Input_Data_File, density_SS, Flushing_Group_ID, Non_Flushing_Group_ID,
                 Act_Stor_Elev, Low_Level_Outlet_Invert_Elevation, Res_Avg_Area, Original_River_Bed_Elevation,
                 start_date, opt_params = None):
        self.T = T
        self.Sed_Mgmt_Type = 'Flushing'
        self.name = element_name
        self.density_SS = density_SS
        self.Flushing_Group_ID = Flushing_Group_ID
        self.Non_Flushing_Group_ID = Non_Flushing_Group_ID
        self.Res_Avg_Area = Res_Avg_Area  # Made an attribute because it's used in methods that get called frequently/externally from reservoir class, so don't want to constantly pass this argument.
        if hasattr(Sediment_Res_Ops, '__init__'):
            Sediment_Res_Ops.__init__(self, element_name, T, Input_Data_File, self.Sed_Mgmt_Type)   # Call parent constructor.
        Flushing.Array_Initialization(self, T)
        Flushing.Import_Data(self, Input_Data_File, element_name, Act_Stor_Elev, Low_Level_Outlet_Invert_Elevation,
                             Original_River_Bed_Elevation, start_date, opt_params)

    def Array_Initialization(self, T):
        # Initialize various variables, arrays, dictionaries, etc.
        self.water_surface_elevation_pre_flush = np.zeros(T + 1)
        self.Flushing_Group_Drawn_Down = 0  # Indicates whether any members of a flushing group are currently drawn down (1 = Yes; 0 = No, default)
        self.Flushing_Channel_Elevation = np.zeros(T + 1)  # Temporarily Stores elevation of flushing channel, which will be reset in reservoir.py using E-V-A data.
        self.Flushing_Event_Number = 0  # Stores each flushing event as a different number.
        self.Flushing_End_Date = 0  # Stores the end date of flushing. Gets updated.
        self.Daily_Sediment_layer_Depth = np.zeros(T)
        self.Sediment_Fraction_Available_Flushing = np.zeros(T+1)
        self.Flushing_Removal_Volume = 0
        self.Drawdown_End_Date = 0  # Variable used to see if drawdown has ended
        self.Time_Elapsed_Since_Drawdown = 0  # Used to track flushing duration (during actual full drawdown/flushing phase, not during drawdown phase).
        self.Time_Elapsed_During_Drawdown = 0  # Used to track duration of reservoir drawdown phase of flushing.
        self.Flushing_Output = []  # Stores summary of every flushing event (mass removed, average flow rate, etc.). A nested list. Each list within the main list will store 8 data points for a particular flushing event.
        self.Flushing_Total_Flows = []  # Used in preparing Flushing_Output.
        self.Flushing_Duration_TEMP = []  # Used in preparing Flushing_Output.
        self.previous_flushing = 0
        self.drawdown_terminated = 0  # Indicates when drawdown phase of flushing is completed.
        self.Flushing_Load_Removed_Daily = 0  # Zero out variable at time zero.
        self.flushed_load = np.zeros(T)  # Daily load released during flushing (series of: Flushing_Load_Removed_Daily)
        self.W_fb_coeff = 12.8  # Coefficient used to define function that describes flushing channel bottom width. This value will be used in deterministic simulations. Stochastic simulations will externally specify this value.

    def Import_Data(self, Input_Data_File, element_name, Act_Stor_Elev, Low_Level_Outlet_Invert_Elevation,
                    Original_River_Bed_Elevation, start_date, opt_params):
        # This method imports flushing relevant data from the Input Data File
        # Method inputs description:
        # Input_Data_File = file where flushing data are stored
        # element_name = name of reservoir at which flushing is to occur
        # Act_Stor_Elev from reservoir

        self.Flushing_Willingness_To_Wait = 7  # reservoir will try for x days (if inflows still aren't sufficient to initiate drawdown). After this long, flushing is cancelled for the year. change so this is user-specified.
        self.Flushing_Drawdown_Willingness_To_Wait = 14  # reservoir will wait x days before refilling while drawn down for another reservoir in a flushing group to finish flushing; change so this is user-specified.

        if 'Flushing' in Input_Data_File.sheetnames:
            self.Flushing_Data_Import = Excel_Data_Import(element_name, Input_Data_File, 'Flushing', 2,
														  16, max_distinct_data_types=None,
														  data_name_offset=2)
        else:
            self.error = 1
            logging.critical("Sediment management of type {0} does not have a corresponding and correctly named "
                             "worksheet in the input file.".format(self.Sed_Mgmt_Type))

        # Unpack flushing optimization preferences. User is not required to use all possible flushing optimization
        # variables (frequency, time of year, and trigger inflow), so here we unpack only those variables the user
        # included in the optimization.
        freq = 0
        day_of_year = 0
        inflow_trigger = 0
        end_date = start_date + timedelta(self.T-1)
        num_years = end_date.year - start_date.year + 1
        if opt_params is not None:
            opt_params = opt_params['Flushing Optimization']
            try:
                ct = 0
                for item in opt_params['Ordered input variable list']:
                    if item == 'Frequency':
                        freq = min(num_years-1, int(opt_params['Flushing Parameters'][ct]))
                    elif item == 'Day of year':
                        day_of_year = int(opt_params['Flushing Parameters'][ct])
                    elif item == 'Inflow trigger':
                        inflow_trigger = opt_params['Flushing Parameters'][ct]
                    ct += 1
            except KeyError:
                pass

        if freq == 0:
            # Use date list provided by user in Flushing sheet.
            date_list = self.Flushing_Data_Import[0]
        else:
            # A flushing frequency has been provided by the optimization routine. However, a day of year may or may
            # not have been provided by the optimization routine. If day_of_year isn't provided, this must be
            # provided in the input file by the user, and read in first here.
            date_list = self.Flushing_Data_Import[0]
            if day_of_year == 0:
                day_of_year = date_list[0].day
            n_flush_events = int(num_years/freq)
            date_list = [0 for i in range(n_flush_events)]
            start_year = start_date.year
            for i in range(n_flush_events):
                date_list[i] = datetime(start_year + freq + i*freq, 1, 1) + timedelta(day_of_year)
        # User-specified flushing data will be stored in a dictionary called Specs
        self.Specs = {}  # Initialize Dictionary
        i = 0  # Initialize counter
        for date_item in date_list:
            self.Specs[date_item] = {}  # Sub-dictionary inside this dictionary will store all relevant specs.
            try:
                self.Specs[date_item]["Minimum Discharge"] = self.Flushing_Data_Import[3][i]  # Minimum discharge target during the actual flushing process.
            except IndexError:
                self.Specs[date_item]["Minimum Discharge"] = self.Flushing_Data_Import[3][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Flushing Channel Bottom Width"] = self.Flushing_Data_Import[4][0]
            except IndexError:
                self.Specs[date_item]["Flushing Channel Bottom Width"] = self.W_fb_coeff * (self.Specs[date_item]["Minimum Discharge"]) ** 0.5  # User elects to have model determine flushing channel bottom width
            try:
                self.Specs[date_item]["Flushing Channel Side Slope"] = self.Flushing_Data_Import[5][0]
            except IndexError:
                self.Specs[date_item]["Flushing Channel Side Slope"] = (31.5 / 5) * (1 / 10) * ((self.density_SS / 1000) ** 4.7)  # User elects to have model determine flushing channel bottom side slope
            try:
                self.Specs[date_item]["Reservoir Bottom Width for Flushing"] = self.Flushing_Data_Import[8][0]
            except IndexError:
                 self.Specs[date_item]["Reservoir Bottom Width for Flushing"] = 0  # User elects to have model determine flushing channel bottom width
            try:
                self.Specs[date_item]["Reservoir Side Slope"] = self.Flushing_Data_Import[9][0]
            except IndexError:
                self.Specs[date_item]["Reservoir Side Slope"] = 1  # User elects to have model determine flushing channel bottom width. Set 45-degree side slope (1:1 ratio).
            if inflow_trigger > 0:
                self.Specs[date_item]["Minimum Drawdown Flow"] = inflow_trigger
            else:
                # Inflow trigger to initiate flushing not specified by optimization. Must import from user
                # preferences in flushing sheet instead.
                try:
                    self.Specs[date_item]["Minimum Drawdown Flow"] = self.Flushing_Data_Import[7][i] # Minimum discharge required to initiate drawdown for flushing.
                except IndexError:
                    try:
                        self.Specs[date_item]["Minimum Drawdown Flow"] = self.Flushing_Data_Import[7][0]  # Load default minimum inflow requirement from first flushing date (or zero if none specified)
                    except IndexError:
                        self.Specs[date_item]["Minimum Drawdown Flow"] = 0  # None spec'd. Inflow>0 triggers drawdown
            try:
                self.Specs[date_item]["Duration"] = self.Flushing_Data_Import[1][i]
            except IndexError:
                self.Specs[date_item]["Duration"] = self.Flushing_Data_Import[1][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            self.Specs[date_item]["Initial Duration"] = self.Specs[date_item]["Duration"]  # When self.specs for date_item gets added to Current Specs, you will decrement duration to track how long flushing has occurred, so it is necessary to store the original value before it's changed.
            try:
                self.Specs[date_item]["Max Drawdown WSE"] = self.Flushing_Data_Import[2][i]
            except IndexError:
                self.Specs[date_item]["Max Drawdown WSE"] = self.Flushing_Data_Import[2][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Max Drawdown Rate"] = self.Flushing_Data_Import[6][i]
            except IndexError:
                self.Specs[date_item]["Max Drawdown Rate"] = self.Flushing_Data_Import[6][0]  # Optional user input - flushing equation coefficient
            if self.Specs[date_item]["Max Drawdown Rate"] is None:
                self.Specs[date_item]["Max Drawdown Rate"] = 0  # Default: This value is used in calculations, so set it to zero if no value specified.
            try:
                self.Specs[date_item]["Removal Coefficient"] = self.Flushing_Data_Import[10][i]
            except IndexError:
                try:
                    self.Specs[date_item]["Removal Coefficient"] = self.Flushing_Data_Import[10][0]
                except IndexError:
                    self.Specs[date_item]["Removal Coefficient"] = 0  # User elects to have model determine flushing coefficient
            #   Optional user input - flushing equation exponent
            try:
                self.Specs[date_item]["Removal Exponent"] = self.Flushing_Data_Import[11][i]
            except IndexError:
                #  User elects to have model determine flushing Exponent
                try:
                    self.Specs[date_item]["Removal Exponent"] = self.Flushing_Data_Import[11][0]
                except IndexError:
                    self.Specs[date_item]["Removal Exponent"] = 0
            i += 1  # Increment counter
        self.Reservoir_Side_Slope = self.Specs[self.Specs.keys()[0]]["Reservoir Side Slope"]  # This variable will remain constant regardless of date. Set variable value with shorter name to use to make procedures below more concise.
        self.Reservoir_Bottom_Width_for_Flushing = self.Specs[self.Specs.keys()[0]]["Reservoir Bottom Width for Flushing"]  # This variable will remain constant regardless of date. Set variable value with shorter name to use to make procedures below more concise.
        self.Flushing_Channel_Bottom_Width = self.Specs[self.Specs.keys()[0]]["Flushing Channel Bottom Width"]  # This variable will remain constant regardless of date. Set variable value with shorter name to use to make procedures below more concise.
        self.Flushing_Channel_Side_Slope = self.Specs[self.Specs.keys()[0]]["Flushing Channel Side Slope"]  # This variable will remain constant regardless of date. Set variable value with shorter name to use to make procedures below more concise.

        # Initialize flushing channel elevation to original river bed elevation at each reservoir
        self.Original_Flushing_Channel_Elevation = Low_Level_Outlet_Invert_Elevation

        # Determine elevation at which flushing channel bank will hit the reservoir bank
        if self.Reservoir_Side_Slope > 0:
            self.W_Res = self.Reservoir_Bottom_Width_for_Flushing + (2 / (self.Reservoir_Side_Slope)) * (Low_Level_Outlet_Invert_Elevation - Original_River_Bed_Elevation)
        self.Flushing_Channel_Bottom_Width = min(self.Flushing_Channel_Bottom_Width, self.Reservoir_Bottom_Width_for_Flushing)
        self.W_tf = min(self.Flushing_Channel_Bottom_Width, self.W_Res) + (2 / self.Flushing_Channel_Side_Slope) * (Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation)
        self.W_t = self.Reservoir_Bottom_Width_for_Flushing + (2 / self.Reservoir_Side_Slope) * (Act_Stor_Elev - Original_River_Bed_Elevation)
        self.W = min(self.W_Res, self.Flushing_Channel_Bottom_Width)
        self.A_r = 0.5 * (self.W_t + self.Reservoir_Bottom_Width_for_Flushing) * (Act_Stor_Elev - Original_River_Bed_Elevation)

        if self.W_tf <= self.W_t and self.Flushing_Channel_Bottom_Width <= self.W_Res:
            #  Case 1: Flushing channel cross-section completely fits inside of reservoir cross-section.
            self.Flushing_Case = 1
            self.A_f = 0.5 * (self.W_tf + self.W) * (Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation)
            self.h_m = 0  # Value does not apply
        elif self.W_tf >= self.W_t and self.Flushing_Channel_Bottom_Width >= self.W_Res:
            #  Case 2: Flushing channel cross-section is completely outside of reservoir cross-section. Note that if flushing elevation not same as original channel elevation, LTCR will not = 1, otherwise it will.
            self.Flushing_Case = 2
            self.A_f = 0.5 * (self.W_Res + self.W_tf) * (Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation) # A_f = A_r ONLY if flushing outlet and original river bed elevation are identical.
        elif self.W_tf >= self.W_t and self.Flushing_Channel_Bottom_Width <= self.W_Res:
            # Case 3: Flushing channel cross-section is restricted by reservoir# s dimensions. Flushing channel is larger than reservoir at top (theoretically), but at bottom it's not.
            self.Flushing_Case = 3
            if self.Reservoir_Side_Slope > 0:
                self.h_m = 2 * (self.W_Res - self.Flushing_Channel_Bottom_Width) * (self.Reservoir_Side_Slope - self.Flushing_Channel_Side_Slope)
                New_Intersection_Flushing_Channel_Width = self.Flushing_Channel_Bottom_Width + (2 / self.Flushing_Channel_Side_Slope) * self.h_m
                self.h_l = Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation - self.h_m
                self.h_f = Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation
                self.A_f = self.W * self.h_f + (self.h_f + self.h_l) * self.h_m * self.Flushing_Channel_Side_Slope + (self.h_l ** 2) * self.Reservoir_Side_Slope
            else:
                pass  # Fix this section. Assume no slope for reservoir means it has a vertical side wall (rectangular cross-section).                            self.h_m = 0 #  User did not specify reservoir's side slope, so this intersection elevation cannot be determined.
        elif self.W_tf <= self.W_t and self.Flushing_Channel_Bottom_Width >= self.W_Res:
            # Case 4: Flushing channel is larger at the bottom of the reservoir, but smaller at the top. Must make adjustments to LTCR.
            self.Flushing_Case = 4
            if self.Reservoir_Side_Slope > 0:
                self.h_m = -self.Reservoir_Side_Slope * self.Flushing_Channel_Side_Slope * 0.5 * (self.Flushing_Channel_Bottom_Width - self.W_Res) / (self.Reservoir_Side_Slope - self.Flushing_Channel_Side_Slope)
                # self.h_m = 2 * (self.W_Res - self.Flushing_Channel_Bottom_Width) * (self.Reservoir_Side_Slope - self.Flushing_Channel_Side_Slope)
                New_Intersection_Flushing_Channel_Width = self.W_Res - self.Flushing_Channel_Side_Slope * (self.Flushing_Channel_Bottom_Width - self.W_Res) / (self.Reservoir_Side_Slope - self.Flushing_Channel_Side_Slope)
                self.h_l = Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation - self.h_m
                self.h_f = Act_Stor_Elev - Low_Level_Outlet_Invert_Elevation
                self.A_f = 0.5 * self.h_m * (self.W_Res + New_Intersection_Flushing_Channel_Width) + 0.5 * self.h_l * (New_Intersection_Flushing_Channel_Width + self.W_tf)
            else:
                pass  # Fix this section. Assume no slope for reservoir means it has a vertical side wall (rectangular cross-section).                            self.h_m = 0 ' User did not specify reservoir's side slope, so this intersection elevation cannot be determined.
        self.LTCR = self.A_f / self.A_r
        if self.LTCR > 1:
            self.LTCR = 1

    def Trapped_Load(self):
        # Computes TE and then determines how much sed. trapped and adjusts mass BS_W, BS_V, and SS_C
        return

    def daily_initiation_check(self, current_date, Daily_Inflow, flushing_dict, sluicing_now, density_current_venting_now):
        # This module checks reservoir conditions (e.g., inflow) during every time step t, before other reservoir's main calcs (e.g., setting/meeting storage or elevation targets by releasing water) are done, to see if the particular sediment management type will be performed during time step t.
        self.current_date = current_date
        self.sluicing_now = sluicing_now
        self.density_current_venting_now = density_current_venting_now
        self.Flushing_Dictionary = flushing_dict  # This dictionary is how the reservoirs communicate their current flushing/drawn down status to one another. Stuff is stored in this dictionary and passed around/updated.
        self.previous_flushing = self.Current  # self.Current is binary check on whether sediment management type is being done. Since it hasn't been set yet for t, here it represents whether this management type happened in day t-1.
        self.Daily_Inflow = Daily_Inflow
        self.Flushing_occurred = 0  # Variable that indicates whether flushing has occurred in a given day. Is zeroed out every day.
        try:
            # This section should only be entered on each first day of flushing. Imports all flushing data and makes flushing current, if there is a dictionary key with a name equal to today's date.
            self.Specs["Current Specs"] = self.Specs[self.current_date].copy()  # Create a new key and data in specs dictionary that sets data for current date to most current data.
            # Before loading all flushing specifications for this date, first need to see if minimum inflow requirement is met for flushing.
            if self.Specs["Current Specs"]["Minimum Drawdown Flow"] < self.Daily_Inflow and self.sluicing_now == 0 and self.density_current_venting_now == 0:
                # Only attempt drawdown in this time period if known reservoir inflow is high enough, AND if sluicing or density current venting is not continuing from previous day (would have been turned off in previous day if it was done).
                self.Non_Flushing_Group_Drawn_Down = 0  # reset maintain drawdown = 0, as we can reset = 1 below if we locate a reason to do so.
                if self.Non_Flushing_Group_ID is not None:
                    # If user specified Non_Flushing_Group_ID value, then there exist other reservoirs in the system that if they are being flushed at the same time, this reservoir (.self) can't be flushed until those ones finish.
                    for item in self.Flushing_Dictionary["Non Drawdown Flushing Group ID"][self.Non_Flushing_Group_ID]:
                        # Loop through reservoir names associated with the key = ID of the sub dictionary "Non Drawdown Flushing Group ID" within the flushing dictionary
                        if item is not self.name:
                            # if the located reservoir name is NOT the same as self (which obviously shares the same ID), then see if that reservoir with the same ID is currently being flushed.
                            if self.Flushing_Dictionary["Currently Flushing"][item] == 1 or self.Flushing_Dictionary["Currently Flushing Group Drawn Down"][item] == 1:
                                self.Non_Flushing_Group_Drawn_Down = 1  # other reservoir in group is being flushed (OR is being kept empty as its part of another flushing group), so we will prevent reservoir from drawing down and flushing today. wait until tomorrow.
                                break
                if self.Non_Flushing_Group_Drawn_Down == 1:
                    self.Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Skip to next day; see explanation above for self.Non_Flushing_Group_Drawn_Down = 1
                else:
                    self.Current = 1  # reservoir should have flushing performed for the specified number of days, starting now
                    self.Flushing_Dictionary["Currently Flushing"][self.name] = self.Current
                    self.Flushing_Event_Number += 1
                    self.Drawdown_Start_Date = self.current_date
                    self.Flushing_Output.append([0 for x in range(10)])  # Add 8 data points to the flushing data storage list. In this case zero-th element is left empty, though not for the others (e.g., Flushing_Duration_TEMP)
                    self.Flushing_Total_Flows.append([0 for x in range(1)])  # Add 1 data point to the flushing flow data storage sub-list.
                    self.Flushing_Duration_TEMP.append([0 for x in range(1)])  # Add 1 data point to the flushing duration temp data storage sub-list.
                    self.Flushing_Output[self.Flushing_Event_Number-1][1] = self.Drawdown_Start_Date  # Store data for this flushing event in its list.
                    self.Flushing_Output[self.Flushing_Event_Number-1][9] = self.LTCR  # Store data for this flushing event in its list.
                    self.Flushing_Duration_TEMP[self.Flushing_Event_Number-1][0] = self.Specs["Current Specs"]["Duration"]
            else:
                self.Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Change Specs dictionary key (date) to next day, as the inflow requirement for drawdown was not satisfied today. Check again tomorrow.
        except KeyError:
            pass  # current date did not exist in flushing dictionary
        if self.Flushing_Group_Drawn_Down == 1:  # If this reservoir is already being kept drawn down, then and only then do we check to see if we should keep it drawn down due to other reservoirs in the group being flushed.
            Flushing.Determine_if_Flushing_Group_Maintain_Drawdown(self)  # Call the method Flushing.Determine_if_Flushing_Group_Maintain_Drawdown(self). Need to set self.Flushing_Group_Drawn_Down = 1 if so, so this will be reflected in reservoir storage target.

    def Determine_if_Flushing_Group_Maintain_Drawdown(self):
        # Determine if we need to KEEP reservoir drawn down because it's in a group in which another reservoir is currently drawn down.
        # This code will be called several times, including twice in main flushing loop (SedMgmt_Main)
        self.Flushing_Group_Drawn_Down = 0  # reset maintain drawdown =0, as we can reset = 1 below if we locate a reason to do so.
        if self.Flushing_Group_ID is not None:
            # If user specified Non_Flushing_Group_ID value, then there exist other reservoirs in the system that if they are being flushed at the same time, this reservoir (.self) can't be flushed until those ones finish.
            for item in self.Flushing_Dictionary["Flushing Group ID"][self.Flushing_Group_ID]:
                # Loop through reservoir names associated with the key = ID of the sub dictionary "Flushing Group ID" within the flushing dictionary
                if item is not self.name:
                    # if the located reservoir name is NOT the same as self (which obviously shares the same ID), then see if that reservoir with the same ID is currently being flushed.
                    if self.Flushing_Dictionary["Currently Flushing"][item] == 1 or self.Flushing_Dictionary["Currently Flushing Group Drawn Down"][item] == 1:
                        self.Flushing_Group_Drawn_Down = 1  # other reservoir in group is being flushed (OR is being kept empty as its part of another flushing group), so we will prevent reservoir from drawing down and flushing today. wait until tomorrow.
                        break

    def Available_Flushing_Load(self, t, Settled_volume, capacity_active_reservoir, capacity_dead_reservoir, capacity_active_reservoir_init, capacity_dead_reservoir_init, Original_River_Bed_Elevation):
        '''

        Args:
            t:
            Settled_volume:
            capacity_active_reservoir:
            capacity_dead_reservoir:
            capacity_active_reservoir_init:
            capacity_dead_reservoir_init:
            Original_River_Bed_Elevation:

        Returns:

        '''
        # Having done TE, determine how much sediment is available to be flushed (flushing_removal_volume and sediment_fraction_available_flushing) and max. flushing load removed, then sets the amounts to actually be removed.
        # This method must be called in every time step, because we must track the evolution of the flushing channel geometry during times of non-flushing.

        # Method input description:
        # Settled_volume = Settled_volume[t] for reservoir
        # capacity_active_reservoir = capacity_active_reservoir[t+1] for reservoir
        # capacity_dead_reservoir = capacity_dead_reservoir[t+1] for reservoir
        # capacity_active_reservoir_init = capacity_active_reservoir[0] for reservoir
        # capacity_dead_reservoir_init = capacity_dead_reservoir[0] for reservoir

        # Before updating channel elevation, compute bottom width of trapezoidal section to be removed via flushing.
        if self.Flushing_Channel_Elevation[t] < self.Original_Flushing_Channel_Elevation:
            # Do not remove sediment until it has accumulated up to elevation of low-level outlet
            self.Flushing_Removal_Volume = 0
        else:
            # sediment layer is at or above elevation of low-level outlet; flushing can proceed now
            if self.Flushing_Case == 2:
                # All sediment can be removed
                self.Flushing_Removal_Volume += Settled_volume
            else:
                self.Daily_Sediment_layer_Depth[t] = Settled_volume / self.Res_Avg_Area
                self.Flushing_Channel_Elevation[t+1] = self.Flushing_Channel_Elevation[t] + self.Daily_Sediment_layer_Depth[t]
                self.Current_Flushing_Channel_Top_Width = min(self.W_tf, self.Flushing_Channel_Bottom_Width + 2 * (self.Flushing_Channel_Elevation[t+1] - self.Original_Flushing_Channel_Elevation) / self.Flushing_Channel_Side_Slope)
                self.Current_Reservoir_Sediment_Top_Width = min(self.W_t, self.Reservoir_Bottom_Width_for_Flushing + 2 * (self.Flushing_Channel_Elevation[t+1] - Original_River_Bed_Elevation) / self.Reservoir_Side_Slope)
                if ((capacity_active_reservoir + capacity_dead_reservoir) / (capacity_active_reservoir_init + capacity_dead_reservoir_init)) >= self.LTCR:
                    self.Flushing_Removal_Volume += min(Settled_volume, Settled_volume * (self.Current_Flushing_Channel_Top_Width / self.Current_Reservoir_Sediment_Top_Width))
                else:
                    self.Flushing_Removal_Volume += Settled_volume
                self.Sediment_Fraction_Available_Flushing[t+1] = self.Current_Flushing_Channel_Top_Width / self.Current_Reservoir_Sediment_Top_Width

        # If Flushing = 1 now, compute actual volumes/loads to be removed.
        if self.Current == 1:
            if self.Specs["Current Specs"]["Duration"] == self.Specs["Current Specs"]["Initial Duration"]:
                # Enter this portion of code only if we are on the first day of flushing (don't want to do this multiple times, as we would then incorrectly use the BS_V(t+1) that has been updated after flushing has taken place
                # Check to see whether reservoir has reached the point where it is full up to the LTCR. Location of this code is important, as it follows removal of sediment from management techniques, and sediment trapping with Brune's curve.
                    # Determine sediment load discharge based on what information the user has specified.
                    if self.Specs["Current Specs"]["Removal Coefficient"] == 0:
                        # self.Max_Flushing_Load_Removed = LTCR(index_res(n)) * (BS_V(t + 1, n) - Sed_Vol_at_end_of_last_flushing(index_res(n))) * density_SS / self.Specs["Current Specs"]["Duration"]
                        self.Max_Flushing_Load_Removed = self.density_SS * self.Flushing_Removal_Volume / self.Specs["Current Specs"]["Initial Duration"]
                        self.counter_23 = 1
                    else:
                        self.Max_Flushing_Load_Removed = self.Specs["Current Specs"]["Removal Coefficient"] * (self.Daily_Inflow ** self.Specs["Current Specs"]["Removal Exponent"]) * self.Daily_Inflow * self.dt * 86400
                        self.counter_23 = 0
            else:
                # Not the first day of flushing, but the reservoir is still scheduled to be flushed.
                if self.counter_23 == 0:
                    self.Max_Flushing_Load_Removed = self.Specs["Current Specs"]["Removal Coefficient"] * (self.Daily_Inflow ** self.Specs["Current Specs"]["Removal Exponent"]) * self.Daily_Inflow * self.dt * 86400

    def Elev_Stor_Target(self, water_surface_elevation):
        # Determine the TARGET water level elevation at the end of the time period, taking into account how much sediment has accumulated in the reservoir and what storage(t+1) the user specified.
        ELEVATION_TARGET = self.Original_Flushing_Channel_Elevation
        # Make sure that achieving elevation target does not result in exceeding the user-specified maximum drawdown rate
        if self.Specs["Current Specs"]["Max Drawdown Rate"] > 0:
            if (water_surface_elevation - ELEVATION_TARGET) > self.Specs["Current Specs"]["Max Drawdown Rate"]:
                ELEVATION_TARGET = water_surface_elevation - self.Specs["Current Specs"]["Max Drawdown Rate"]
        # Sets elevation and storage targets. Highly customizable, as it’s modified depending on what sediment management is taking place.
        return ELEVATION_TARGET

    def SedMgmt_Main(self, t, S, EVA, water_surface_elevation, Q_out, BS_W, density_SS):
        '''
        Purpose: Main method for flushing mass balance. Removes sediment then adjusts mass balance.

        Args:
            t: Current time period (integer, value range: [0,T])
            S: self.S[t+1] for reservoir
            EVA: self.E_V_A for reservoir
            water_surface_elevation: self.water_surface_elevation[t] for reservoir
            Q_out: self.Q_out[t] for reservoir
            BS_W: self.BS_W[t+1] for reservoir
            density_SS: self.density_SS for reservoir

        Returns:

        '''

        # Initialize local variables. Their values should only be returned if > 0
        SS_W_out = 0
        BS_W = BS_W

        # Determine final water level elevation associated with final water storage, before flushing is successful in removing any sediment.
        if self.Current == 1 or self.Flushing_Group_Drawn_Down == 1:
            self.water_surface_elevation_pre_flush[t+1] = Matrix_Interpolation(self.name, EVA, "storage", "elevation", S)

        # Main flushing sediment mass balance routine.
        if self.Drawdown_End_Date != 0:
            # Drawdown was previously completed, we are in the flushing phase now, so keep track of how long we are flushing
            self.Time_Elapsed_Since_Drawdown = (self.current_date - self.Drawdown_End_Date).days + 1  # track time Once reservoir is drawn down, track how long it stays drawn down
        else:
            # Drawdown is still proceeding; need to track its duration and terminate flushing if willingness to wait for drawdown is exceeded
            self.Time_Elapsed_During_Drawdown = (self.current_date - self.Drawdown_Start_Date).days + 1  # Track time elapsed since drawdown was initiated.
        if ((water_surface_elevation <= self.Specs["Current Specs"]["Max Drawdown WSE"]) and (self.water_surface_elevation_pre_flush[t+1] <= self.Specs["Current Specs"]["Max Drawdown WSE"])):
            if self.Drawdown_End_Date == 0:
                # Record that drawdown is completed so we know how long we have been drawn down, but only record once, so we don# t write over value each time
                self.Drawdown_End_Date = self.current_date
                self.Time_Elapsed_During_Drawdown = 0 # Variable is just used to make sure willingness to wait for drawdown has not been exceeded, so here we eliminate that constraint since drawdown has been achieved.
            if (Q_out - self.Specs["Current Specs"]["Minimum Discharge"] >= 0) or (self.Specs["Current Specs"]["Minimum Discharge"] - Q_out <= 0.2 * self.Specs["Current Specs"]["Minimum Discharge"]):
                # This time period resulted in successful flushing; reduce number of remaining required flushing days
                self.Flushing_occurred = 1
                self.Specs["Current Specs"]["Duration"] = self.Specs["Current Specs"]["Duration"] - 1
                self.Flushing_Output[self.Flushing_Event_Number - 1][4] = 1 + self.Flushing_Output[self.Flushing_Event_Number - 1][4] # Increment Number of Flushing days by 1
                self.Flushing_Load_Removed_Daily = min(BS_W, self.Max_Flushing_Load_Removed)
                self.flushed_load[t] = self.Flushing_Load_Removed_Daily
                BS_W -= self.Flushing_Load_Removed_Daily  # This is simply updating the BS_W(t+1,n) that has already been updated in the sediment trapping section above.
                BS_V = BS_W / density_SS
                if self.Flushing_Output[self.Flushing_Event_Number - 1][2] == 0:
                    self.Flushing_Start_Date = self.current_date
                    self.Flushing_Output[self.Flushing_Event_Number - 1][2] = self.Flushing_Start_Date
                self.Flushing_Output[self.Flushing_Event_Number - 1][8] = self.Flushing_Output[self.Flushing_Event_Number - 1][8] + self.Flushing_Load_Removed_Daily
                self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] += Q_out
                if self.Specs["Current Specs"]["Duration"] == 0:
                    # Flushing is now complete, though the reservoir may need to remain drawn down if upstream grouped reservoirs are still being flushed.
                    self.Current = 0
                    self.Flushing_Dictionary["Currently Flushing"][self.name] = self.Current
                    Flushing.Determine_if_Flushing_Group_Maintain_Drawdown(self)  # Call the method Flushing.Determine_if_Flushing_Group_Maintain_Drawdown(self). Need to set self.Flushing_Group_Drawn_Down = 1 if so, so this will be reflected in reservoir storage target.
                    self.Flushing_End_Date = self.current_date
                    self.Flushing_Output[self.Flushing_Event_Number - 1][3] = self.Flushing_End_Date
                    counter_1 = (self.Flushing_Start_Date - self.Drawdown_Start_Date).days + 1
                    self.Flushing_Output[self.Flushing_Event_Number - 1][4] = counter_1
                    counter_1 = (self.current_date - self.Flushing_Start_Date).days + 1
                    self.Flushing_Output[self.Flushing_Event_Number - 1][5] = counter_1
                    self.Flushing_Output[self.Flushing_Event_Number - 1][7] = self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] / self.Flushing_Duration_TEMP[self.Flushing_Event_Number - 1][0]
                    self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] = 0  # Reset so array can be used next time for Flushing
                    self.counter_23 = 0
                    self.Flushing_Removal_Volume = 0
                    self.Flushing_Start_Date = 0
                    self.Drawdown_Start_Date = 0
                    self.Drawdown_End_Date = 0
                    self.Time_Elapsed_Since_Drawdown = 0
                SS_W_out = self.Flushing_Load_Removed_Daily  # This is only a sub-component of SS_W_out. Sediment suspended in released water from other outlets will be added separately later.
            else:
                # The minimum flow rate required to achieve flushing was not satisfied.
                #  if time elapsed during drawdown has exceeded willingness to wait, terminate flushing
                if self.Time_Elapsed_During_Drawdown >= self.Flushing_Drawdown_Willingness_To_Wait and self.Flushing_Drawdown_Willingness_To_Wait > 0 and self.Time_Elapsed_During_Drawdown > 0:
                    self.drawdown_terminated = 1
                    self.Current = 0
                    self.Flushing_Dictionary["Currently Flushing"][self.name] = self.Current
        else:
            # Storage was not sufficiently drained for flushing to occur, so try in next time period to finish draining.
            #  if time elapsed during drawdown has exceeded willingness to wait, terminate flushing
            if self.Time_Elapsed_During_Drawdown >= self.Flushing_Drawdown_Willingness_To_Wait and self.Flushing_Drawdown_Willingness_To_Wait > 0 and self.Time_Elapsed_During_Drawdown > 0:
                self.drawdown_terminated = 1
                self.Current = 0
                self.Flushing_Dictionary["Currently Flushing"][self.name] = self.Current

        # Determine whether maximum willingness to wait while drawn down has been exceeded.
        if self.Current == 1:
            if self.Time_Elapsed_Since_Drawdown >= self.Flushing_Willingness_To_Wait and self.Flushing_occurred == 1:
                pass # if we have waited the maximum duration for flushing to be successful, but it was successful in this time period, then we can wait one extra day. Extend flushing by one more day.
            else:
                if self.Time_Elapsed_Since_Drawdown > 0 and self.Time_Elapsed_Since_Drawdown >= self.Flushing_Willingness_To_Wait:
                    # Before we end flushing to begin refill because we have exceeded willingness to wait, we must make sure no other reservoirs in our flushing group are still being flushed.
                    # In this scenario, we have exceeded flushing willingness to wait, and the reservoir was NOT flushed successfully today, so we should turn flushing off at ROI if flushing at upstream reservoir is ending.
                    Flushing.Determine_if_Flushing_Group_Maintain_Drawdown(self)
                    # Willingness to stay drawn down has been exceeded, so it is now time to end flushing; the question is whether to keep the group drawn down.
                    if self.Flushing_Group_Drawn_Down == 1:
                        pass  # Do not end flushing
                    else:
                        # Flushing must now be stopped, as (1) time elapsed since drawdown has exceeded the dam operator's willingness to wait for sufficient flushing inflow, and (2) no reservoirs upstream within this flushing group are currently being flushed.
                        self.Current = 0
                        self.Flushing_Dictionary["Currently Flushing"][self.name] = self.Current
                        self.Time_Elapsed_Since_Drawdown = 0
                        self.Drawdown_End_Date = 0
                        self.Flushing_End_Date = self.current_date
                        self.Flushing_Output[self.Flushing_Event_Number - 1][3] = self.Flushing_End_Date

                        if (self.Flushing_Duration_TEMP[self.Flushing_Event_Number - 1][0] - self.Specs["Current Specs"]["Duration"]) > 0:
                            # Some flushing has occurred, so there is technically a flushing end date now
                            counter_1 = (self.Flushing_Start_Date - self.Drawdown_Start_Date).days + 1
                            self.Flushing_Output[self.Flushing_Event_Number - 1][4] = counter_1
                            counter_1 = (self.current_date - self.Flushing_Start_Date).days + 1
                            self.Flushing_Output[self.Flushing_Event_Number - 1][5] = counter_1
                        else:
                            counter_1 = (self.current_date - self.Drawdown_Start_Date).days + 1
                            self.Flushing_Output[self.Flushing_Event_Number - 1][4] = counter_1
                        if self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] > 0:
                            self.Flushing_Output[self.Flushing_Event_Number - 1][7] = self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] / (self.Flushing_Duration_TEMP[self.Flushing_Event_Number - 1][0] - self.Specs["Current Specs"]["Duration"])
                        self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] = 0  # Reset so array can be used next time for Flushing
                        self.counter_23 = 0
                        self.Flushing_Removal_Volume = self.Flushing_Removal_Volume - (1 / density_SS) * (self.Flushing_Duration_TEMP[self.Flushing_Event_Number - 1][0] - self.Specs["Current Specs"]["Duration"]) * self.Flushing_Load_Removed_Daily  # Not all of the sediment in the flushing channel was removed, since flushing was cut short by operator.
                        #  Reset remaining variables to zero now that they have been used.
                        self.Flushing_Start_Date = 0
                        self.Drawdown_Start_Date = 0
        if self.drawdown_terminated == 1 and self.Current == 0:
            # Check to see if flushing drawdown willinginess to wait has been exceeded and flushing terminated
            # Drawdown (and therefore flushing, which never started by definition) was terminated due to its duration exceeding the specified willingness to wait for drawdown.
            #  Reset all variables to zero
            self.Time_Elapsed_Since_Drawdown = 0
            self.Time_Elapsed_During_Drawdown = 0
            self.Drawdown_End_Date = 0
            self.Flushing_End_Date = 0
            self.Flushing_Total_Flows[self.Flushing_Event_Number - 1][0] = 0  # Reset so array can be used next time for Flushing
            self.counter_23 = 0
            self.Flushing_Start_Date = 0
            self.Drawdown_Start_Date = 0
        self.drawdown_terminated = 0  # zero out variable for use next time
        return BS_W, SS_W_out

    def Distribute_Outflow_Among_Outlets(self):
        # Takes a computed release from function Total_Outflow_Calc() and distributes that flow among the outlet works present at the dam.
        return

    def mass_balance(self, t=None):
        return