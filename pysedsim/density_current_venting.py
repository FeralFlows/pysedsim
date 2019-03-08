# -*- coding: utf-8 -*-

'''

Module to define the Density Current Venting class and asssociated methods.

This class represents the density current venting technique for routing sediment through a reservoir. Class instances (objects) become
stored in reservoir instances as part of its operating policy. Hence, this class contains methods for determining whether density current
venting will happen at a reservoir on a given day, definition of required reservoir water level targets, distribution of flow among
outlets to achieve those elevation targets, and mass balance calculations for the reservoir.

'''

# This module contains the density current class definition

# Import relevant modules
from __future__ import division
from sediment_res_ops import Sediment_Res_Ops
import data_processing
from matrix_interpolation import Matrix_Interpolation
from math import log

class Density_Current_Venting(Sediment_Res_Ops):
    def __init__(self, element_name, T, Input_Data_File, res_length):
        self.Sed_Mgmt_Type = 'Density Current Venting'
        self.reservoir_name = element_name
        if hasattr(Sediment_Res_Ops, '__init__'):
            Sediment_Res_Ops.__init__(self, element_name, T, Input_Data_File, self.Sed_Mgmt_Type)  # Call parent constructor.
        Density_Current_Venting.Array_Initialization(self, T)
        Density_Current_Venting.Import_Data(self, Input_Data_File, element_name)
        self.Reservoir_Length = res_length

    def Array_Initialization(self, T):
        return

    def Import_Data(self, Input_Data_File, element_name):
        # See if DCV is scheduled to occur on this date. If so, import data and set Density Current Venting = 1.
        if 'Density Current Venting' in Input_Data_File.sheetnames:
            [self.Min_venting_efficiency, self.Min_venting_stor_elev_target, self.Venting_downstream_max_conc,
             self.continue_venting_despite_max_exceedance, self.reservoir_bed_slope, self.reservoir_bottom_width,
             self.minimum_power_during_venting] = data_processing.Excel_Data_Import(element_name, Input_Data_File,
                                                                                    'Density Current Venting', 1, 7,
                                                                                    max_distinct_data_types=None,
                                                                                    data_name_offset=None)  # Optional
        else:
            self.error = 1
            print "Error: Sediment management of type %s does not have a corresponding and correctly named worksheet in the input file." \
                  % self.Sed_Mgmt_Type
        # Set defaults if none specified
        if self.Min_venting_efficiency is not None:
            pass  # do nothing, value set already by user
        else:
            if self.Reservoir_Length is not None:
                self.Min_venting_efficiency = -0.08 * log(self.Reservoir_Length) + 0.5384
            else:
                self.Min_venting_efficiency = 101
        if self.continue_venting_despite_max_exceedance == "Yes":
            self.continue_venting_despite_max_exceedance = 1
        elif self.continue_venting_despite_max_exceedance == "No":
            self.continue_venting_despite_max_exceedance = 0
        else:
            self.continue_venting_despite_max_exceedance = 1  # Set default value

    def daily_initiation_check(self, Daily_Inflow, Daily_Volume_in, Daily_Sediment_in, flushing_now, sluicing_now):
        # This module checks reservoir conditions (e.g., inflow) during time step t, before other reservoirs calcs are done,
        # to see if the particular sediment management type will be performed during time step t.
        # Import method inputs
        self.flushing_now = flushing_now
        self.sluicing_now = sluicing_now
        self.previous_density_current_venting = self.Current
        self.percent_passing_d_90 = 100
        previous_percent_passing_d_90 = 0.1
        # Set/initialize values (ones that must be initialized every day)
        self.DCV_Flow = 0  # Set = 0 every day. Gets set = Q_low_level_outlet later if DCV happens, for use in preventing overestimate of sediment released during DCV.
        # Main loop
        while abs(self.percent_passing_d_90 - previous_percent_passing_d_90) / previous_percent_passing_d_90 > 0.005:
            # 1. Determine density
            water_density = 0.0006 * (self.percent_passing_d_90 / 100) * Daily_Sediment_in / Daily_Volume_in + 0.998232
            # 2. Determine water depth at plunge point of DC
            # water_depth_plunge_point = ((Q_in / (0.78 * self.reservoir_bottom_width)) ^ (2 / 3)) * ((water_density - 0.998232) * 9.81 / 0.998232) ^ (-1 / 3)
            # 3. Calculate the velocity of the current based on the solids concentration in the inflow
            velocity_density_current = ((8 / 0.025) * ((water_density - 0.998232) * 9.81 / 0.998232) * (Daily_Inflow / self.reservoir_bottom_width) * self.reservoir_bed_slope) ^ (1 / 3)
            # 4. Based on the velocity, determine the maximum grain size that can be transported
            d_90_density_current = -0.0074 * (velocity_density_current ^ 2) + 0.0369 * velocity_density_current + 0.0007
            # 5. Remove all larger grains from transport (they settle), and return to step 3 to recompute the carrying capacity of the solids that remain in suspension.
            previous_percent_passing_d_90 = self.percent_passing_d_90
            self.percent_passing_d_90 = 15.226 * log(d_90_density_current) + 95.839
            if self.percent_passing_d_90 <= 0:
                self.percent_passing_d_90 = 0
                break

        if self.percent_passing_d_90 > self.Min_venting_efficiency:
            self.Current = 1  # reservoir should be vented during this time period
        else:
            self.Current = 0

    def Elev_Stor_Target(self, water_surface_elevation, STORAGE_TARGET, MW_Capacity, Turbine_Efficiency, Q_in, S, V_in, SS_W_in, SS_C, Capacity_low_level_Outlet, Capacity_Controlled_Outlet, Capacity_mid_level_Outlet, Capacity_Spillway_Outlet, Capacity_Hydropower_Outlet, Evap, Settled_volume, Dam_Tailwater_Elevation, Reservoir_Type, E_V_A):
        # Sets elevation and storage targets. Highly customizable, as it’s modified depending on what sediment management is taking place.
        # Inputs:
        # water_surface_elevation = self.water_surface_elevation[t] for reservoir (set in previous day)
        # MW_Capacity = self.MW_Capacity for reservoir
        # Turbine_Efficiency = self.Turbine_Efficiency for reservoir
        # Q_in = Q_in[t] for reservoir
        # STORAGE_TARGET = self.STORAGE_TARGET[t + 1] for reservoir
        # S = self.S[t] for reservoir
        # V_in = self.V_in[t] for reservoir
        # SS_W_in = self.SS_W_in[t] for reservoir
        # SS_C = self.SS_C[t] for reservoir
        # Capacity_low_level_Outlet = self.Capacity_low_level_Outlet[t-1] for reservoir
        # Capacity_Controlled_Outlet = self.Capacity_Controlled_Outlet[t-1] for reservoir
        # Capacity_mid_level_Outlet = self.Capacity_mid_level_Outlet[t-1] for reservoir
        # Capacity_Spillway_Outlet = self.Capacity_Spillway_Outlet[t - 1] for reservoir
        # Capacity_Hydropower_Outlet = Capacity_Hydropower_Outlet[t - 1] for reservoir
        # Evap = self.Evap[t] for reservoir
        # Settled_volume = self.Settled_volume[t]
        # Dam_Tailwater_Elevation = self.Dam_Tailwater_Elevation[t] for reservoir
        # Reservoir_Type = self.Reservoir_Type for reservoir.
        # E_V_A = self.E_V_A for reservoir

        # Determine the maximum turbine flow based on capacity. Only used internally to this method to set release targets during DCV.

        hydropower_head_average = max(0, water_surface_elevation - Dam_Tailwater_Elevation)
        if MW_Capacity > 0:
            if hydropower_head_average > 0:
                Max_MW_Turbine_Flow = (1 / 0.00981) * (MW_Capacity) * (1 / hydropower_head_average) * (1 / Turbine_Efficiency)
            else:
                Max_MW_Turbine_Flow = 0
        else:
            Max_MW_Turbine_Flow = 0

        # Capacity_low_level_Outlet is at time t-1 (since this has not yet been computed for time t). This results in an estimation error.
        original_storage_target_before_venting = STORAGE_TARGET
        Density_Current_low_level_flow = min(Q_in, Capacity_low_level_Outlet)
        self.Density_Current_Turbine_Flow = min(Max_MW_Turbine_Flow, self.minimum_power_during_venting / (0.00981 * (water_surface_elevation - Dam_Tailwater_Elevation) * Turbine_Efficiency))
        if Reservoir_Type == "Diversion" or Reservoir_Type == "Diversion and Power":
            # Reservoir is a diversion, the hydropower flow from self.Density_Current_Turbine_Flow cannot be used to satisfy the downstream concentration requirement. Situations for the other 2 reservoir types are satisfied below.
            density_current_concentration = ((self.percent_passing_d_90 / 100) * (SS_W_in / V_in) * Density_Current_low_level_flow) / (Density_Current_low_level_flow)
        else:
            density_current_concentration = ((self.percent_passing_d_90 / 100) * (SS_W_in / V_in) * Density_Current_low_level_flow + SS_C * self.Density_Current_Turbine_Flow) / (self.Density_Current_Turbine_Flow + Density_Current_low_level_flow)
        if density_current_concentration < (self.Venting_downstream_max_conc / 1000):
            # Additional releases not required, so storage target can now be established based on known releases.
            STORAGE_TARGET = S + V_in - (Density_Current_low_level_flow + self.Density_Current_Turbine_Flow + Evap) * self.dt * 86400 - Settled_volume
        else:
            # Additional releases are required to satisfy downstream concentration target
            # First, determine whether the extra releases needed to achieve the desired concentration target are even possible.
            required_extra_release_for_max_conc = Density_Current_low_level_flow * ((self.Venting_downstream_max_conc / 1000) - (self.percent_passing_d_90 / 100) * (SS_W_in / V_in)) / (SS_C - self.Venting_downstream_max_conc / 1000)  # This represents whatever amount in excess of the low level outlets must be released; this includes water from hydro (if not diverted)
            if self.continue_venting_despite_max_exceedance == 1:
                pass  # No action required (don't need to check to see if we can really meet target given capacities of available outlets that release flow downstream)
            else:
                # Check to see if available outlets that can release flow downstream have enough capacity to meet conc. target. if not, cancel DCV.
                # First, determine actual capacity of mid-level outlets in previous time step. Frequently this is set to zero (in code below) if DCV=0 or Sluicing=0 to avoid errors, so it must be recomputed here to make sure we are using the right value.
               if Reservoir_Type == "Diversion" or Reservoir_Type == "Diversion and Power":
                   if required_extra_release_for_max_conc < (Capacity_Controlled_Outlet + Capacity_mid_level_Outlet + Capacity_Spillway_Outlet):
                       pass  # No action required, as it appears concentration target will be met
                   else:
                       self.Current = 0  # Downstream concentration target cannot likely be met. Cancel DCV.
               else:
                   if required_extra_release_for_max_conc < (Capacity_Controlled_Outlet + Capacity_Hydropower_Outlet + Capacity_mid_level_Outlet + Capacity_Spillway_Outlet):
                       pass  # No action required, as it appears concentration target will be met
                   else:
                       self.Current = 0  # Downstream concentration target cannot likely be met. Cancel DCV.

            if Reservoir_Type == "Diversion" or Reservoir_Type == "Diversion and Power":
                # water is diverted, so must explicitly include that hydro-produced water in the storage target
                STORAGE_TARGET = S + V_in - (required_extra_release_for_max_conc + self.Density_Current_Turbine_Flow + Density_Current_low_level_flow + Evap) * self.dt * 86400 - Settled_volume
            else:
                # water is not diverted, so don't include the hydro water in the storage target because it's already included in the required_extra_release variable
                STORAGE_TARGET = S + V_in - (required_extra_release_for_max_conc + self.Density_Current_Turbine_Flow + Density_Current_low_level_flow + Evap) * self.dt * 86400 - Settled_volume

        # Pick the least of the original elevation/storage target and the new target; no problem to go lower than DCV target to get closer to original target, but going higher to hit original target would reduce outflow and thus prevent satisfying DCV objectives (power, concentration, etc.)
        STORAGE_TARGET = min(STORAGE_TARGET, original_storage_target_before_venting)
        # Determine water surface elevation associated with new storage target
        ELEVATION_TARGET = Matrix_Interpolation(self.reservoir_name, E_V_A, "storage", "elevation", STORAGE_TARGET)

        self.power_priority = 0
        if ELEVATION_TARGET > self.Min_venting_stor_elev_target:
            pass  # Did not go below minimum water surface elevation, so no action required
        else:
            # Set elevation target to minimum value, and then reset storage target to corresponding value.
            self.power_priority = 1  # Give power production the priority in this situation
            ELEVATION_TARGET = self.Min_venting_stor_elev_target
            STORAGE_TARGET = Matrix_Interpolation(self.reservoir_name, E_V_A, "elevation", "storage", ELEVATION_TARGET)
        return ELEVATION_TARGET, STORAGE_TARGET

    def Distribute_Outflow_Among_Outlets(self, Input_Array, Reservoir_Type, Environmental_Release_Goal, Q_out, Q_in):
        # Takes a computed release from function Total_Outflow_Calc() and distributes that flow among the outlet works present at the dam.

        # Method inputs:
        # t = current time step in main simulation
        # Orifices = self.Orifices (from reservoir) object containing all daily outlet flows and capacities for all outlets, as well as the outlet data dictionary.
        # Environmental_Release_Goal = self.Environmental_Release_Goal[t] for reservoir
        # Q_out = self.Q_out[t] for reservoir
        # Reservoir_Type = self.Reservoir_Type for reservoir
        # Input_Array contains all outlet capacities (self.Orifices.Outlet_Capacity...) during time step t for reservoir.

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

        # Takes a computed release from function Total_Outflow_Calc() and distributes that flow among the outlet works present at the dam.
        if Reservoir_Type in ["Diversion", "Diversion and Power"]:
            if self.power_priority == 0:
                # Reservoir level is drawn down close to minimum user-specified level; power production must receive priority (DCV no longer a priority due to low res. levels)
                Q_low_level_outlet = min(Q_in, Q_out, Capacity_low_level_Outlet)
                if Reservoir_Type == "Diversion":
                    Q_diversion = min(max(Q_out - Q_low_level_outlet, 0), Capacity_Diverted_Outlet)
                    Q_turbines = 0
                else:
                    Q_diversion = min(max(Q_out - Q_low_level_outlet, 0), Capacity_Hydropower_Outlet)
                    Q_turbines = Q_diversion
            else:
                if Reservoir_Type == "Diversion":
                    Q_low_level_outlet = min(Q_in, Q_out, Capacity_low_level_Outlet)
                    Q_diversion = min(Q_out, self.Density_Current_Turbine_Flow, Capacity_Diverted_Outlet)
                    Q_turbines = 0
                else:
                    Q_diversion = min(Q_out, self.Density_Current_Turbine_Flow, Capacity_Hydropower_Outlet)
                    Q_turbines = Q_diversion
                    Q_low_level_outlet = min(Q_in, max(0, Q_out - Q_diversion), Capacity_low_level_Outlet)
            # Now, re-assign additional flow to the turbines (and/or diversion outlets) and low level outlet; then proceed to mid-level and overflow outlets
            if Reservoir_Type == "Diversion":
                Q_diversion = Q_diversion + min(max(Q_out - Q_low_level_outlet - Q_diversion, 0), Capacity_Diverted_Outlet - Q_diversion)
            else:
                Q_diversion = Q_diversion + min(max(Q_out - Q_low_level_outlet - Q_diversion, 0), Capacity_Hydropower_Outlet - Q_diversion)
                Q_turbines = Q_diversion
            Q_low_level_outlet = Q_low_level_outlet + min(max(Q_out - Q_low_level_outlet - Q_diversion, 0), Capacity_low_level_Outlet - Q_low_level_outlet)
            Q_mid_level_outlet = min(max(Q_out - Q_diversion - Q_low_level_outlet, 0), Capacity_mid_level_Outlet)
            Q_controlled = min(max(Q_out - Q_low_level_outlet - Q_diversion - Q_mid_level_outlet, 0), Capacity_Controlled_Outlet)
            Q_overflow = min(max(Q_out - Q_low_level_outlet - Q_diversion - Q_mid_level_outlet - Q_controlled, 0), Capacity_Spillway_Outlet)
        elif Reservoir_Type == "Power":
            Q_controlled = 0
            if self.power_priority == 0:
                Q_low_level_outlet = min(Q_in, Q_out, Capacity_low_level_Outlet)
                Q_turbines = min(Q_out - Q_low_level_outlet, Capacity_Hydropower_Outlet)
            else:
                Q_turbines = min(Q_out, self.Density_Current_Turbine_Flow, Capacity_Hydropower_Outlet)
                Q_low_level_outlet = min(Q_in, Q_out - Q_turbines, Capacity_low_level_Outlet)
            # Now, re-assign additional flow to the hydro, then low level outlet; then proceed to mid-level and overflow outlets
            Q_turbines = Q_turbines + min(max(Q_out - Q_low_level_outlet - Q_turbines, 0), Capacity_Hydropower_Outlet - Q_turbines)
            Q_low_level_outlet = Q_low_level_outlet + min(max(Q_out - Q_low_level_outlet - Q_turbines, 0), Capacity_low_level_Outlet - Q_low_level_outlet)
            Q_mid_level_outlet = min(max(Q_out - Q_low_level_outlet - Q_turbines, 0), Capacity_mid_level_Outlet)
            Q_overflow = min(max(Q_out - Q_low_level_outlet - Q_turbines - Q_mid_level_outlet, 0), Capacity_Spillway_Outlet)
        elif Reservoir_Type == "Storage":
            Q_turbines = 0
            Q_diversion = 0
            Q_downstream = Q_out
            Q_controlled = 0
            Q_low_level_outlet = min(Q_out, Capacity_low_level_Outlet)
            Q_mid_level_outlet = min(max(Q_out - Q_low_level_outlet, 0), Capacity_mid_level_Outlet)
            Q_overflow = min(max(Q_out - Q_low_level_outlet - Q_mid_level_outlet, 0), Capacity_Spillway_Outlet)
        else:
            pass  # Currently there are no other reservoir types than the four covered above.
        return Q_controlled, Q_mid_level_outlet, Q_diversion, Q_turbines, Q_low_level_outlet, Q_overflow, Q_downstream

    def SedMgmt_Main(self, BS_W, SS_W_in, V_in, Q_low_level_outlet, density_SS, capacity_active_reservoir, capacity_dead_reservoir):
        # Main method for DCV mass balance. Removes sediment then adjusts mass balance.
        # Computes TE and then determines how much sed. trapped and adjusts mass BS_W, BS_V, and SS_C

        # Method inputs:
        # V_in = self.V_in[t] for reservoir
        # Q_low_level_outlet = self.Q_low_level_outlet[t] for reservoir
        # SS_W_in = self.SS_W_in[t] for reservoir
        # TE_avg = self.TE_avg[t] for reservoir
        # density_SS = self.density_SS for reservoir
        # capacity_active_reservoir = self.capacity_active_reservoir[t+1] for reservoir.
        # capacity_dead_reservoir = self.capacity_dead_reservoir[t+1] for reservoir.

        # Determine new TE and settled mass,volume
        self.DCV_Flow = Q_low_level_outlet  # For use in preventing overestimate of sediment released during DCV.
        BS_W = BS_W
        TE_avg = 1 - (Q_low_level_outlet * self.dt * 86400 * (SS_W_in / V_in) * self.percent_passing_d_90 / 100) / SS_W_in  # Never computed this above, because percent_passing_d_90 doesn't reflect how much sediment mass was actually released (the latter is impacted by release capacity, among other things)
        Settled_mass = TE_avg * SS_W_in # Sediment in suspension from previous time step does not get trapped
        Settled_volume = Settled_mass / density_SS
        # Ensure settled volume isn't more than available capacity
        if Settled_volume > capacity_active_reservoir + capacity_dead_reservoir:
            Settled_volume = capacity_active_reservoir + capacity_dead_reservoir
            Settled_mass = Settled_volume * density_SS
        # Determine new sediment mass, volume balances in reservoir
        BS_W += Settled_mass
        BS_V = BS_W / density_SS
        SS_W_out = SS_W_in - Settled_mass  # This is only a sub-component of SS_W_out. Sediment suspended in released water from other outlets will be added separately later.
        return TE_avg, Settled_mass, Settled_volume, BS_W, BS_V, SS_W_out, self.DCV_Flow

    def SSC_Adjustment_Evap(self):
        # Adjusts SSC concentration to account for evaporation.
        return

    def Distribute_Outflow_Among_Outlets(self):
        # Takes a computed release from function Total_Outflow_Calc() and distributes that flow among the outlet works present at the dam.
        return