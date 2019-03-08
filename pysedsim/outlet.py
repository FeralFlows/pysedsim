# -*- coding: utf-8 -*-
'''

Module to define the Dam Orifice (outlet) class and methods.

A reservoir can contain these outlet objects, which compute the reservoir's maximum capacity to release flow given
the reservoir's water level.

'''
# Description: This file creates class definitions for several dam outlet types

# import relevant libraries
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0)
import numpy as np
import data_processing
import calendar  # Used to determine number of days in a given month
from matrix_interpolation import Matrix_Interpolation

# Ideas:
# Since the imports here (openpyxl modules) get executed every time this file gets imported in reservoir.py, you should just feed in the entire list of reservoirs here, maybe.

class Outlet:
    def __init__(self, storage_element_name, T, Input_Data_File):
        # storage_element_name = reservoir name, dam name, channel name, etc.
        # Input_Data_File = .xlsx input file (workbook) name
        # Load Outlet Capacity Table
        [self.Outlet_Capacity_Dict, self.num_outlets] = data_processing.Excel_Data_Import(storage_element_name,
                                                                                          Input_Data_File,
                                                                                          "Outlet Capacity Data", 2, 2,
                                                                                          max_distinct_data_types=7,
                                                                                          data_name_offset=4)
        # Initialize Arrays
        Outlet.Array_Initialization(self, storage_element_name, T)

    def Array_Initialization(self, storage_element_name, T):
        self.dt = 1  # 1-day time step
        self.Q_turbines = np.zeros(T)
        self.Q_spill = np.zeros(T)
        self.Q_overflow = np.zeros(T)
        self.Q_controlled = np.zeros(T)
        self.Q_diversion = np.zeros(T)
        self.Q_low_level_outlet = np.zeros(T)
        self.Q_mid_level_outlet = np.zeros(T)
        self.Q_downstream = np.zeros(T)
        self.mid_level_outlet = 0  # Indicates whether a mid level outlet exists. Will be reset later if one exists.
        self.Capacity_Controlled_Outlet = np.zeros(T)
        self.Capacity_Diverted_Outlet = np.zeros(T)
        self.Capacity_Hydropower_Outlet = np.zeros(T)
        self.Capacity_low_level_Outlet = np.zeros(T)
        self.Capacity_Spillway_Outlet = np.zeros(T)
        self.Capacity_mid_level_Outlet = np.zeros(T)

        # Loop through outlet dictionary keys, and initialize all appropriate arrays.
        for key in self.Outlet_Capacity_Dict:
            if key == 'Hydropower Outlet':
                # Determine elevation of outlet (at zero capacity)
                self.hydropower_outlet_elevation = self.Outlet_Capacity_Dict['Hydropower Outlet'][0][0]
            elif key == 'Spillway Outlet':
                # Determine elevation of outlet (at zero capacity)
                self.spillway_outlet_elevation = self.Outlet_Capacity_Dict['Spillway Outlet'][0][0]
            elif key == 'Controlled Outlet':
                # Determine elevation of outlet (at zero capacity)
                self.controlled_outlet_elevation = self.Outlet_Capacity_Dict['Controlled Outlet'][0][0]
            elif key == 'Diverted Outlet':
                # Determine elevation of outlet (at zero capacity)
                self.diverted_outlet_elevation = self.Outlet_Capacity_Dict['Diverted Outlet'][0][0]
            elif key == 'Hydropower/Diversion Outlet':
                # Determine elevation of outlet (at zero capacity)
                self.hydropower_outlet_elevation = self.Outlet_Capacity_Dict['Hydropower/Diversion Outlet'][0][0]
            elif key == 'Low Level Outlet':
                # Low level outlet elevation
                self.Low_Level_Outlet_Invert_Elevation = self.Outlet_Capacity_Dict['Low Level Outlet'][0][0]
                self.low_level_outlet_elevation = self.Low_Level_Outlet_Invert_Elevation
            elif key == 'Mid Level Outlet':
                # Determine existence and elevation of mid level outlet (at zero capacity)
                self.mid_level_outlet = 1  # Indicates whether a mid level outlet exists.
                self.mid_level_outlet_elevation = self.Outlet_Capacity_Dict['Mid Level Outlet'][0][0]

    def Outlet_RK_Capacity_Estimation(self, t):

        # Runge Kutta Loop
        # For now just set to high value until RK loop is coded.
        self.Capacity_Controlled_Outlet = 0
        self.Capacity_Diverted_Outlet = 0
        self.Capacity_Hydropower_Outlet = 17668
        self.Capacity_low_level_Outlet = 100000
        self.Capacity_Spillway_Outlet = 60000
        self.Capacity_mid_level_Outlet = 0

        return

    def Outlet_Euler_Capacity_Estimation(self, res_name, t, current_date, Monthly_Evap_Data, S, water_surface_elevation,
                                         water_surface_area, Environmental_Release_Goal, Flushing,
                                         Flushing_Group_Drawn_Down, Sluicing, Density_Current_Venting, Reservoir_Type,
                                         Settled_volume, Q_in, E_V_A):

        # Inputs:
        # S = self.S[t] for reservoir
        # current_date = self.current_date for reservoir
        # Monthly_Evap_Data = self.Monthly_Evap_Data for reservoir
        # water_surface_elevation = self.water_surface_elevation[t] for reservoir
        # water_surface_area = self.water_surface_area[t] for reservoir
        # Environmental_Release_Goal = self.Environmental_Release_Goal[t] for reservoir
        # res_name = self.name for reservoir. Used only to produce error messages from interpolation process in outlet capacity elevation table.
        # E_V_A = self.E_V_A for reservoir
        # Reservoir_Type = self.Reservoir_Type for reservoir
        # Flushing = self.Operating_Policy["Flushing"].Current for reservoir
        # Flushing_Group_Drawn_Down = self.Operating_Policy["Flushing"].Flushing_Group_Drawn_Down for reservoir
        # Sluicing = self.Operating_Policy["Sluicing"].Current
        # Density_Current_Venting = self.Operating_Policy["Density Current Venting"].Current
        # Q_in = self.Q_in[t] for reservoir
        # Settled_volume = self.Settled_volume[t] for reservoir
        # E_V_A = self.E_V_A for reservoir
        # Determine release capacities of all outlets based on current storage estimate and one projected future estimate.

        # Estimate evaporation rate (m3/s) that is expected during time period, given current storage S(t,n) and expected storage STORAGE_TARGET(t+1,index_res(n))
        RK_Evap = (water_surface_area * 10000 * (Monthly_Evap_Data[current_date.month - 1] / 1000) * (
        1 / (calendar.monthrange(current_date.year, current_date.month)[1]))) * (1 / 86400)
        Max_Release_Storage_Level = max(0, S + (Q_in - RK_Evap) * self.dt * 86400 - Settled_volume)

        # For each outlet perform search and interpolation to produce a rough estimate of the outlet capacity for the average storage during each time period.
        # Then, Select the minimum of the outlet capacity and the volume above each outlet divided by seconds per day

        for key in self.Outlet_Capacity_Dict:
            if key == "Controlled Outlet":
                self.Capacity_Controlled_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Controlled Outlet'], "elevation", "flow", water_surface_elevation)
                # Make sure controlled outlet doesn't release more than environmental flows, since that is all they are used to accomplish in SedSim.
                if Reservoir_Type == "Power" or Reservoir_Type == "Diversion" or Reservoir_Type == "Diversion and Power":
                    self.Capacity_Controlled_Outlet[t] = min(self.Capacity_Controlled_Outlet[t], Environmental_Release_Goal)
                self.Capacity_Controlled_Outlet[t] = min(self.Capacity_Controlled_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.controlled_outlet_elevation)))
            elif key == "Diverted Outlet":
                self.Capacity_Diverted_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Diverted Outlet'], "elevation", "flow", water_surface_elevation)
                self.Capacity_Diverted_Outlet[t] = min(self.Capacity_Diverted_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.diverted_outlet_elevation)))
            elif key == "Spillway Outlet":
                self.Capacity_Spillway_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Spillway Outlet'], "elevation", "flow", water_surface_elevation)
                self.Capacity_Spillway_Outlet[t] = min(self.Capacity_Spillway_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.spillway_outlet_elevation)))
            elif key == "Hydropower Outlet":
                self.Capacity_Hydropower_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Hydropower Outlet'], "elevation", "flow", water_surface_elevation)
                self.Capacity_Hydropower_Outlet[t] = min(self.Capacity_Hydropower_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.hydropower_outlet_elevation)))
            elif key == "Hydropower/Diversion Outlet":
                self.Capacity_Hydropower_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Hydropower/Diversion Outlet'], "elevation", "flow", water_surface_elevation)
                self.Capacity_Hydropower_Outlet[t] = min(self.Capacity_Hydropower_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.hydropower_outlet_elevation)))
            elif key == "Low Level Outlet":
                if Flushing == 0 and Flushing_Group_Drawn_Down == 0 and Sluicing == 0 and Density_Current_Venting == 0:
                    #No flushing/sluicing/venting during this time period; do not use low level gates just to achieve new storage/elevation target
                    self.Capacity_low_level_Outlet[t] = 0
                else:
                    self.Capacity_low_level_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Low Level Outlet'], "elevation", "flow", water_surface_elevation)
                self.Capacity_low_level_Outlet[t] = min(self.Capacity_low_level_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.low_level_outlet_elevation)))
            elif key == "Mid Level Outlet":
                if Sluicing == 0 and Density_Current_Venting == 0: #No sluicing during this time period; do not use low level gates just to achieve new storage/elevation target
                    #Do nothing; no capacity
                    self.Capacity_mid_level_Outlet[t] = 0
                else:
                    self.Capacity_mid_level_Outlet[t] = Matrix_Interpolation(res_name, self.Outlet_Capacity_Dict['Mid Level Outlet'], "elevation", "flow", water_surface_elevation)
                self.Capacity_mid_level_Outlet[t] = min(self.Capacity_mid_level_Outlet[t], (1 / 86400) * max(0, Max_Release_Storage_Level - Matrix_Interpolation(res_name, E_V_A, "elevation", "flow", self.mid_level_outlet_elevation)))

        #Estimate evaporation that is expected during time period, given current storage S(t,n) and expected storage STORAGE_TARGET(t+1,index_res(n))
        Evap = RK_Evap
        return Evap

    def outlet_capacity(self, cap):
        self.outlet_cap = cap
        return