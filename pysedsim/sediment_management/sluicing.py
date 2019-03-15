# -*- coding: utf-8 -*-

'''
This module contains the sluicing class definition.

Sluicing is a type of sediment management that releases sediment from the reservoir through mid- and/or low-level
outlets in the reservoir before it can settle. This is frequently done through a partial drawdown  of the reservoir's
water level (though no drawdown is required) for some period of time. This will ultimately be part of a reservoir's
operating policy if the user indicates sluicing is to occur at the reservoir.

'''

# Import relevant modules
from __future__ import division  # This ensures result of quotient of two integers will be a float
from pysedsim.sediment_management.sediment_res_ops import Sediment_Res_Ops
from pysedsim.data_processing.data_processing import *
from datetime import timedelta  # Used to add days/months/years to a datetime object

class Sluicing(Sediment_Res_Ops):
    def __init__(self, element_name, T, Input_Data_File, start_date):
        self.Sed_Mgmt_Type = 'Sluicing'
        self.start_date = start_date
        if hasattr(Sediment_Res_Ops, '__init__'):
            Sediment_Res_Ops.__init__(self, element_name, T, Input_Data_File, self.Sed_Mgmt_Type)  # Parent constructor.
        Sluicing.Array_Initialization(self, T)
        Sluicing.Import_Data(self, Input_Data_File, element_name)

    def Array_Initialization(self, T):
        self.Post_Sluicing = 0
        self.Post_Sluicing_1 = 0
        self.previous_sluicing = 0

    def Import_Data(self, Input_Data_File, element_name):
        # See if sluicing is scheduled to occur on this date. If so, import data and set Sluicing = 1.
        if 'Sluicing' in Input_Data_File.sheetnames:
            self.Sluicing_Data_Import = data_processing.Excel_Data_Import(element_name, Input_Data_File, 'Sluicing', 2,
                                                                          8, max_distinct_data_types=None,
                                                                          data_name_offset=2)
        else:
            self.error = 1
            print "Error: Sediment management of type %s does not have a corresponding and correctly named worksheet " \
                  "in the input file." % self.Sed_Mgmt_Type

        # User-specified sluicing data will be stored in a dictionary called Specs
        self.Specs = {}  # Initialize Dictionary
        try:
            if self.Sluicing_Data_Import[0][0] is not None:
                self.Non_date_based_sluicing = 0  # Sluicing data should be specified for different dates.
        except IndexError:
            # Only one row of data are specified, as sluicing will be initiated based on inflow conditions,
            # not according to particular dates. No dates or durations are specified.
            self.Non_date_based_sluicing = 1
            try:
                self.Sluicing_Data_Import[0][0] = self.start_date  # Only one date will exist; the simulation start date. hence there is only one sub-dictionary. Will be updated throughout simulation.
            except IndexError:
                self.Sluicing_Data_Import[0].append(self.start_date)
        i = 0  # Initialize counter
        for date_item in self.Sluicing_Data_Import[0]:
            self.Specs[self.start_date] = {}  # Sub-dictionary inside this dictionary will store all relevant specs.
            try:
                self.Specs[date_item]["Drawdown trigger inflow"] = self.Sluicing_Data_Import[1][i]
            except IndexError:
                self.Specs[date_item]["Drawdown trigger inflow"] = self.Sluicing_Data_Import[1][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Refill trigger inflow"] = self.Sluicing_Data_Import[3][i]  # if Q_in exceeds this, sluicing will continue, otherwise not.
            except IndexError:
                self.Specs[date_item]["Refill trigger inflow"] = self.Sluicing_Data_Import[3][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Target drawdown elevation"] = self.Sluicing_Data_Import[4][i]
            except IndexError:
                try:
                    self.Specs[date_item]["Target drawdown elevation"] = self.Sluicing_Data_Import[4][0]  # User did not specify an appropriate value; Load the default value stored in the first row
                except IndexError:
                    # No first row value specified. Sluicing will proceed at any elevation when Q>Q_trigger
                    self.Specs[date_item]["Target drawdown elevation"] = 'Inflow-based Sluicing'
            try:
                self.Specs[date_item]["Max drawdown rate"] = self.Sluicing_Data_Import[5][i]
            except IndexError:
                self.Specs[date_item]["Max drawdown rate"] = self.Sluicing_Data_Import[5][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Max refill rate"] = self.Sluicing_Data_Import[6][i]
            except IndexError:
                self.Specs[date_item]["Max refill rate"] = self.Sluicing_Data_Import[6][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Hydropower during sluicing"] = self.Sluicing_Data_Import[7][i]
            except IndexError:
                self.Specs[date_item]["Hydropower during sluicing"] = self.Sluicing_Data_Import[7][0]  # User did not specify an appropriate value; Load the default value stored in the first row
            if self.Specs[date_item]["Hydropower during sluicing"] == "Yes":
                self.Sluicing_power_yes_no = 1
            else:
                self.Sluicing_power_yes_no = 0
            if self.Non_date_based_sluicing == 0:
                try:
                    self.Specs[date_item]["Duration"] = self.Sluicing_Data_Import[2][i]
                except IndexError:
                    self.Specs[date_item]["Duration"] = self.Sluicing_Data_Import[2][0] #  User did not specify an appropriate value; Load the default value stored in the first row
            else:
                self.Specs[date_item]["Duration"] = 0
            self.Specs[date_item]["Sluicing End Date"] = date_item + timedelta(days = self.Specs[date_item]["Duration"])
            i += 1  # Increment counter

    def Trapped_Load(self):
        # Computes TE and then determines how much sed. trapped and adjusts mass BS_W, BS_V, and SS_C
        return

    def daily_initiation_check(self, current_date, Daily_Inflow, flushing_now, density_current_venting_now):
        # This module checks reservoir conditions (e.g., inflow) during time step t, before other reservoirs calcs are done, to see if the particular sediment management type will be performed during time step t.
        # Store method inputs.
        self.current_date = current_date
        self.flushing_now = flushing_now
        self.density_current_venting_now = density_current_venting_now

        # Before sluicing loop, see if sluicing completed in last loop, to see whether or not to begin post-sluicing refill routine.
        if self.Post_Sluicing_1 == 1:
            self.Post_Sluicing = 1
        self.previous_sluicing = self.Current

        # Attempt to locate scheduled sluicing event. If it's today, then our "Current" sluicing data gets switched to data for sluicing event scheduled for this date.
        try:
            if Daily_Inflow > self.Specs[self.current_date]["Drawdown trigger inflow"] and self.flushing_now == 0 and self.density_current_venting_now == 0:
                # This is the first day of sluicing. This section will not be entered when sluicing is ongoing after the first day, until a new sluicing event is initiated.
                # Imports all sluicing data and makes sluicing current, if there is a dictionary key with a name equal to today's date.
                self.Current = 1
                self.Specs["Current Specs"] = self.Specs[self.current_date].copy()  # Create a new key and data in specs dictionary that sets data for current date to most current data.
                self.Specs["Current Specs"]["Start date"] = self.current_date
            else:
                self.Specs[self.current_date]["Sluicing End Date"] += timedelta(days = 1)  # Delay end date by one day.
                self.Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # Update dictionary key so it reflects delay of sluicing by one day, as inflow conditions have not been met.
                self.Current = 0  # Sluicing is 0 by definition in this case
        except KeyError:
            pass  # current date did not exist in sluicing dictionary, so this is not the first day. Reservoir could be currently sluicing, though.
        try:
            if self.current_date >= self.Specs["Current Specs"]["Start date"] and self.current_date <= self.Specs["Current Specs"]["Sluicing End Date"]:
                if self.current_date == self.Specs["Current Specs"]["Sluicing End Date"]:
                    # reservoir should be sluiced during this time period. Last date of sluicing, so check to make sure inflow conditions are met, otherwise extend sluicing by another day.
                    if self.Specs["Current Specs"]["Refill trigger inflow"] is not None:
                        if Daily_Inflow > self.Specs["Current Specs"]["Refill trigger inflow"]:
                            self.Specs["Current Specs"]["Sluicing End Date"] += timedelta(days = 1)  # Extend sluicing by an additional day in CURRENT dictionary, as inflow is still higher than threshold.
                            if self.Non_date_based_sluicing == 1:
                                self.Specs[self.current_date]["Sluicing End Date"] += timedelta(days = 1)  # Delay end date by one day in ORIGINAL user input.
                                self.Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # REPLACE ORIGINAL dictionary key so it reflects delay of sluicing by one day, as inflow conditions have not been met.
                        else:
                            self.Post_Sluicing_1 = 1
                            if self.Non_date_based_sluicing == 1:
                                # If sluicing is only triggered by inflow, not dates, then reset the next sluicing start date to tomorrow.
                                self.Specs[self.current_date]["Sluicing End Date"] += timedelta(days = 1)  # Delay end date by one day in ORIGINAL user input.
                                self.Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # REPLACE ORIGINAL dictionary key so it reflects delay of sluicing by one day, as inflow conditions have not been met.
                    else:
                        # no minimum sluicing inflow value for refill was set by user
                        self.Post_Sluicing_1 = 1
                        if self.Non_date_based_sluicing == 1:
                            # If sluicing is only triggered by inflow, not dates, then reset the next sluicing start date to tomorrow.
                            self.Specs[self.current_date]["Sluicing End Date"] += timedelta(days = 1)  # Delay end date by one day in ORIGINAL user input.
                            self.Specs[self.current_date + timedelta(days = 1)] = self.Specs.pop(self.current_date)  # REPLACE ORIGINAL dictionary key so it reflects delay of sluicing by one day, as inflow conditions have not been met.
                self.Current = 1  # Not the first or last sluicing day, but still a sluicing day. reservoir should be sluiced during this time period.
            else:
                self.Current = 0  # The current specs dictionary doesn't have a sluicing time frame within which the current date falls.
        except KeyError:
            pass
        # Set a variable value that will allow capacity of hydropower (and/or hydropower/diversion outlets) to be set to zero when sluicing is occurring if no power production is allowed during sluicing.
        if self.Current == 1:
            if self.Sluicing_power_yes_no == 1:
                self.hydro_multiplier = 1
            else:
                self.hydro_multiplier = 0
        else:
            self.hydro_multiplier = 1

    def Elev_Stor_Target(self, water_surface_elevation, ELEVATION_TARGET_PRE_SET):
        # Sets elevation and storage targets. Highly customizable, as it’s modified depending on what sediment management is taking place.
        # Convert the sluicing target elevation into a storage target that can be used for reservoir operations.
        if self.Current == 1:
            ELEVATION_TARGET = self.Specs["Current Specs"]["Target drawdown elevation"]
            # Make sure that achieving elevation target does not result in exceeding the user-specified maximum drawdown rate
            if (water_surface_elevation - ELEVATION_TARGET) > self.Specs["Current Specs"]["Max drawdown rate"]:
                ELEVATION_TARGET = water_surface_elevation - self.Specs["Current Specs"]["Max drawdown rate"]
            else:
                pass  # No action required; already set the sluicing elevation target above
        elif self.Post_Sluicing == 1:
            ELEVATION_TARGET = ELEVATION_TARGET_PRE_SET
            if ELEVATION_TARGET_PRE_SET - water_surface_elevation > self.Specs["Current Specs"]["Max refill rate"]:
                ELEVATION_TARGET = water_surface_elevation + self.Specs["Current Specs"]["Max refill rate"]
            else:
                pass  # No action required; already set the sluicing elevation target above
        return ELEVATION_TARGET

    def Post_Sluicing_Check(self, STORAGE_TARGET, S):
        # This method checks to see if sluicing is completed (elevation returned to within 5% of target WSE for that day). Means operating policy is back on its normal path.

        # Method inputs:
        # STORAGE_TARGET = self.STORAGE_TARGET[t+1] for reservoir
        # S = self.S[t + 1] for reservoir

        # Main method:
        if self.Post_Sluicing == 1:
            if STORAGE_TARGET > 0:
                if 100 * ((S - STORAGE_TARGET) / STORAGE_TARGET) >= -5:
                    self.Post_Sluicing = 0
                    self.Post_Sluicing_1 = 0

    def Distribute_Outflow_Among_Outlets(self, Input_Array, t, Reservoir_Type, Environmental_Release_Goal, Q_out, mid_level_outlet):
        # Takes a computed release from function Total_Outflow_Calc() and distributes that flow among the outlet works present at the dam.

        # Method inputs:
        # t = current time step in main simulation
        # Orifices = self.Orifices (from reservoir) object containing all daily outlet flows and capacities for all outlets, as well as the outlet data dictionary.
        # Environmental_Release_Goal = self.Environmental_Release_Goal[t] for reservoir
        # Q_out = self.Q_out[t] for reservoir
        # Reservoir_Type = self.Reservoir_Type for reservoir
        # mid_level_outlet = self.Orifices.mid_level_outlet for reservoir (indicates whether a mid level outlet exists)
        # Input_Array contains all outlet capacities (self.Orifices.Outlet_Capacity...) during time step t for reservoir.
        Capacity_Controlled_Outlet = Input_Array[0]
        Capacity_mid_level_Outlet = Input_Array[1]
        Capacity_low_level_Outlet = Input_Array[2]
        Capacity_Diverted_Outlet = Input_Array[3]
        Capacity_Hydropower_Outlet = Input_Array[4]
        Capacity_Spillway_Outlet = Input_Array[5]

        # Zero out all important local variables in every time step
        Q_controlled = 0
        Q_mid_level_outlet = 0
        Q_diversion = 0
        Q_turbines = 0
        Q_low_level_outlet = 0
        Q_overflow = 0
        Q_downstream = 0

        if Reservoir_Type in ["Diversion", "Diversion and Power"]:
            Q_controlled = min(Q_out, Environmental_Release_Goal, Capacity_Controlled_Outlet)
            if self.Sluicing_power_yes_no == 0:
                # Sluicing will be performed, but no power will be produced during this process, so water should be released primarily from mid-level outlets
                Q_mid_level_outlet = min(max(Q_out - Q_controlled, 0), Capacity_mid_level_Outlet)
                if mid_level_outlet == 0:
                    # No mid-level outlets, so can use LL outlet for sluicing
                    Q_low_level_outlet = min(max(Q_out - Q_mid_level_outlet, 0), Capacity_low_level_Outlet)
            else:
                # Sluicing will be performed, but power will be produced during this process, so water should be first released through power/diversion outlets, and then through mid-levels
                if Reservoir_Type == "Diversion":
                    Q_diversion = min(max(Q_out - Q_controlled, 0), Capacity_Diverted_Outlet)
                    Q_turbines = 0
                else:
                    Q_diversion = min(max(Q_out - Q_controlled, 0), Capacity_Hydropower_Outlet)
                    Q_turbines = Q_diversion
                # Having allocated water to diversion and/or power outlets, now allocate water to mid-level outlets
                Q_mid_level_outlet = min(max(Q_out - (Q_controlled + Q_diversion), 0), Capacity_mid_level_Outlet)
                if mid_level_outlet == 0:
                    # No mid-level outlets, so can use LL outlet for sluicing
                    Q_low_level_outlet = min(max(Q_out - (Q_controlled + Q_diversion + Q_mid_level_outlet), 0), Capacity_low_level_Outlet)
        elif Reservoir_Type == "Power":
            if self.Sluicing_power_yes_no == 0:
                # Sluicing will be performed, but no power will be produced during this process, so water should be released primarily from mid-level outlets
                Q_mid_level_outlet = min(Q_out, Capacity_mid_level_Outlet)
                if mid_level_outlet == 0:
                    # No mid-level outlets, so can use LL outlet for sluicing
                    Q_low_level_outlet = min(max(Q_out - Q_mid_level_outlet, 0), Capacity_low_level_Outlet)
            else:
                # Sluicing will be performed, and power will be produced during this process, so water should be first released through power outlets, and then through mid-levels
                if Q_out <= Capacity_Hydropower_Outlet:
                    # No spilling over the spillway needs to occur
                    Q_turbines = min(Q_out, Capacity_Hydropower_Outlet)
                else:
                    Q_turbines = Capacity_Hydropower_Outlet
                    # In a reservoir that only generates hydropower, if the release flow rate is higher than the hydropower plant capacity, this flow is probably being spilled through spillway, but need to see if water elevation is high enough for spilling
                # Having allocated water to power outlets, now allocate water to mid-level outlets
                Q_mid_level_outlet = min(max(Q_out - Q_turbines, 0), Capacity_mid_level_Outlet)
                if mid_level_outlet == 0:
                    # No mid-level outlets, so can use LL outlet for sluicing
                    Q_low_level_outlet = min(max(Q_out - Q_turbines - Q_mid_level_outlet, 0), Capacity_low_level_Outlet)
        elif Reservoir_Type == "Storage":
            Q_turbines = 0
            Q_diversion = 0
            Q_downstream = Q_out
            Q_mid_level_outlet = min(Q_out, Capacity_mid_level_Outlet)
            if mid_level_outlet == 0:
                # No mid-level outlets, so can use LL outlet for sluicing
                Q_low_level_outlet = min(max(Q_out - Q_mid_level_outlet, 0), Capacity_low_level_Outlet)
            if Q_out - Q_mid_level_outlet - Q_low_level_outlet > 0:
                # Some flow remains to be allocated. Attempt allocation to controlled outlets.
                if Q_out - Q_mid_level_outlet <= Capacity_Controlled_Outlet:
                    Q_controlled = min(Q_out - Q_mid_level_outlet, Capacity_Controlled_Outlet)
                    Q_overflow = 0
                else:
                    Q_controlled = Capacity_Controlled_Outlet
                    Q_overflow = min(max(Q_out - Q_controlled - Q_mid_level_outlet, 0), Capacity_Spillway_Outlet)
        else:
            pass  # Currently there are no other reservoir types than the four covered above.
        return Q_controlled, Q_mid_level_outlet, Q_diversion, Q_turbines, Q_low_level_outlet, Q_overflow, Q_downstream

    def mass_balance(self, constant, t = None):
        self.att = constant
        self.removed_load = 2
        return self.removed_load
