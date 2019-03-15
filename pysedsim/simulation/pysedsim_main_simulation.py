# -*- coding: utf-8 -*-

'''

This module houses the main PySedSim simulation loop.

This includes creation of all the simulation objects, which occurs in every simulation, and the balance
component as well, which consists primarily of calls to the Master_Method_Caller of the class the object being
simulated was instantiated from. The class methods actually contain the mass balance simulation code.

'''

# Main imports
from pysedsim.river_basin_elements.system_element_creation import System_Object_Creation
import pandas as pd
import numpy as np

def SedSim_Main_Simulation(Num_Realizations, T, Input_Data_File, element_stochastic_components, SystemObjects, Element_Dict,
                           Flushing_Group_Dict, Stochastic_Sim, Parameter_Input_Dictionary, simulation_dates_no_leap,
                           Time_Series_Output_Dictionary, Sim_Dur, simulation_dates, Col_Names, Simulation_mode,
                           Output_Object_Dict, var_sub_list, element_export_list, Sampled_Parameter_Dict,
                           Synthetic_Inflow_dataframe_name_LIST, Synthetic_Inflows_dictionary, op_policy_params=None):

    # Inputs:
    # Notes: All inputs are generated automatically in the main top-level pysedsim.py file through various function calls there.

    try:
        stochastic_flow_list = Parameter_Input_Dictionary['1']['Locations']
    except KeyError:
        stochastic_flow_list = []

    distribution_name = {}
    for i in SystemObjects["Ordered Simulation List"]:
        distribution_name[i] = None

    for rz in range(Num_Realizations):
        # Create all objects (instances of classes) to be simulated (reservoirs, reaches, etc.). SystemObjects is fed in and then updated.
        [SystemObjects, Flushing_Group_Dict, element_stochastic_components] = System_Object_Creation(T, Input_Data_File,
                                                                                                     element_stochastic_components,
                                                                                                     SystemObjects,
                                                                                                     Element_Dict,
                                                                                                     Flushing_Group_Dict,
                                                                                                     stochastic_flow_list,
                                                                                                     op_policy_params=op_policy_params)

        # If this is a monte carlo, set relevant values for all objects (e.g., daily inflows at junction for this simulation) that were set
        # in the monte carlo function.
        if Stochastic_Sim == 1:
            for i in SystemObjects["Ordered Simulation List"]:
                # Loop through System Objects. Every time you hit a junction, set the junction incremental flow variable equal to the next
                # Synthetic flow realization column from the Synthetic_Inflows_dictionary.
                try:
                    distribution_name[i] = Sampled_Parameter_Dict[i]['2']  # Need to feed this to junction calibration routine if it exists.
                except KeyError:
                    pass
                if SystemObjects[i].Element_Sub_Dict["Type"] == "Junction":
                    for keys in Synthetic_Inflows_dictionary:
                        if i == keys:
                            # We have a match: system object name is identical to the name of the dataframe for which synthetic flows
                            # exist in Synthetic_Inflows_dictionary. Take the rz-th (e.g., 1 through 100th) member of the list
                            # corresponding to junction name of interest. Use .values so existing variable is not converted to pandas DF.
                            SystemObjects[i].Q_incremental = Synthetic_Inflows_dictionary[keys][
                                Synthetic_Inflow_dataframe_name_LIST[keys][rz]].values
                            SystemObjects[i].Q_incremental = SystemObjects[i].Q_incremental.astype(float)  # Cast as float avoids errors.
                            SystemObjects[i].Q_out_unreg_for_param_calib = SystemObjects[i].Q_incremental  # Must follow above step.
                            break
                    if '5' in Sampled_Parameter_Dict[i].keys():
                        SystemObjects[i].sed_beta = Sampled_Parameter_Dict[i]['5'][rz]
                    if '3' in Sampled_Parameter_Dict[i].keys():
                        # Set both the annual sediment load and cumulative annual sediment load. The user provides the cumulative load,
                        # but it is necessary to temporarily set the annual load before an adjustment takes place below.
                        SystemObjects[i].Annual_SED_LOAD = Sampled_Parameter_Dict[i]['3'][rz]
                        SystemObjects[i].Cum_Annual_SED_LOAD = Sampled_Parameter_Dict[i]['3'][rz]
                elif SystemObjects[i].Element_Sub_Dict["Type"] == "Reservoir":
                    if '7' in Sampled_Parameter_Dict[i].keys():
                        SystemObjects[i].Sed_Trapping_Curve_Spec = Sampled_Parameter_Dict[i]['7'][rz]
                    if '8' in Sampled_Parameter_Dict[i].keys():
                        SystemObjects[i].density_SS = Sampled_Parameter_Dict[i]['8'][rz]
                    if '10' in Sampled_Parameter_Dict[i].keys():
                        SystemObjects[i].E_Sed_Curve_Adjustment(Sampled_Parameter_Dict[i]['10'][rz])
                    if '12' in Sampled_Parameter_Dict[i].keys():
                        if 'Sluicing' in SystemObjects[i].Operating_Policy.keys():
                            SystemObjects[i].stoch_trap_adjust = Sampled_Parameter_Dict[i]['12'][rz]
                    if '13' in Sampled_Parameter_Dict[i].keys():
                        if 'Flushing' in SystemObjects[i].Operating_Policy.keys():
                            SystemObjects[i].Operating_Policy["Flushing"].W_fb_coeff = Sampled_Parameter_Dict[i]['13'][rz]
                elif SystemObjects[i].Element_Sub_Dict["Type"] == "Bypass Structure":
                    if '11' in Sampled_Parameter_Dict[i].keys():
                        SystemObjects[i].Bypass_Fraction = Sampled_Parameter_Dict[i]['11'][rz]

        i = 0  # initialize counter

        # Only execute this portion during the first realization.
        # Check to see that incremental sediment load calibration is required anywhere.
        incremental_calibration = 0  # Default. Will be reset below if applicable.
        for i in SystemObjects["Ordered Simulation List"]:
            if SystemObjects[i].Element_Sub_Dict["Type"] == "Junction":
                if SystemObjects[i].calibration_preference == 1:
                    incremental_calibration = 1
                    break

        # Check to see that channel carrying capacity sediment load calibration is required anywhere.
        carrying_capacity_calibration = 0  # Set default. Will be reset below if applicable.
        for i in SystemObjects["Ordered Simulation List"]:
            if SystemObjects[i].Element_Sub_Dict["Type"] == "Reach":
                if SystemObjects[i].calibration_preference == 1:
                    carrying_capacity_calibration = 1
                    break
                else:
                    carrying_capacity_calibration = 0

        # Loop to store incremental annual sediment loads.
        if incremental_calibration == 1:
            for i in SystemObjects["Ordered Simulation List"]:
                if SystemObjects[i].Cum_Annual_SED_LOAD == 0:
                    # Not a junction, so by definition no incremental load exists. However, cumulative sediment load data need to be stored
                    # for use in channel carrying capacity calibration.
                    for item in SystemObjects[i].Element_Sub_Dict["Inflow Elements"]:
                        # Determine daily element inflows, use it in master caller.
                        SystemObjects[i].Cum_Annual_SED_LOAD += SystemObjects[item].Cum_Annual_SED_LOAD
                else:
                    # Element is a junction. Do nothing, as cumulative sediment data were provided by user for this junction. Still need to
                    # determine the incremental annual sediment load, though.
                    for item in SystemObjects[i].Element_Sub_Dict["Inflow Elements"]:
                        SystemObjects[i].Annual_SED_LOAD -= SystemObjects[item].Cum_Annual_SED_LOAD
                    SystemObjects[i].Calibration_Incremental_Sediment_Load(element_stochastic_components[i], distribution_name[i])
        # Carrying capacity calibration - determination of cumulative daily flow passing each point.
        if carrying_capacity_calibration == 1:
            for i in SystemObjects["Ordered Simulation List"]:
                # Do the reach carrying capacity calibration. First determine cumulative daily element
                for item in SystemObjects[i].Element_Sub_Dict["Inflow Elements"]:
                    SystemObjects[i].Q_out_unreg_for_param_calib += SystemObjects[item].Q_out_unreg_for_param_calib
                if SystemObjects[i].Element_Sub_Dict["Type"] == "Reach":
                    SystemObjects[i].Calibration_Reach_Sediment_Carrying_Capacity()

        # Add dredging inflow elements to reservoirs that are the "Dredging Outflow Element" of any reservoirs being dredged.
        for i in SystemObjects["Ordered Simulation List"]:
            try:
                if SystemObjects[i].Element_Sub_Dict["Dredging Outflow Element"] is not None:
                    # Dredging exists, so add it to destination elements list.
                    SystemObjects[SystemObjects[i].Element_Sub_Dict["Dredging Outflow Element"]].Element_Sub_Dict[
                        "Dredging Inflow Elements"].append(i)
            except KeyError:
                pass  # No dredging exists for reservoir

        # Main PySedSim simulation loop. Loops through objects in order at each time step. Simulation begins at time t = 0, though time t=0
        # values for storage-type variables are loaded during instantiation of classes (reservoirs, reaches, etc. have initial/time zero
        # values loaded there).
        for t in range(0, T):
            for i in SystemObjects["Ordered Simulation List"]:
                # Before simulating next object, locate data that is to be shared among elements, stored in each
                # element's self.elem_xfer_output_dict dictionary. Then store in destination element's input dictionary.
                for j in SystemObjects["Ordered Simulation List"]:
                    if j in SystemObjects["Ordered Simulation List"][0:SystemObjects["Ordered Simulation List"].index(i) + 1]:
                        # If element j has been simulated in time step t already.
                        try:
                            # See if object j has data to transfer to object i. Transfer if so. In case element that needs to transfer
                            # data is downstream and hence hasn't yet been simuluated, need to check this again later.
                            SystemObjects[i].elem_xfer_input_dict[j] = SystemObjects[j].elem_xfer_output_dict[i]
                        except KeyError:
                            pass
                Flow_in = 0
                Sed_in = 0
                for item in SystemObjects[i].Element_Sub_Dict["Inflow Elements"]:
                    # Determine daily element inflows, use it in master caller.  Locate the outflows from all the upstream elements,
                    # whose names are contained in Element_Dict.
                    Flow_in += SystemObjects[item].Element_Sub_Dict["Daily Water Outflows"][i]
                    Sed_in += SystemObjects[item].Element_Sub_Dict["Daily Sediment Outflows"][i]
                    if SystemObjects[item].Element_Sub_Dict["Type"] == "Bypass Structure":
                        # If "i" is a junction for which "item" is a bypass upstream, set water and sediment inflow rates for
                        # reservoir/bypass
                        SystemObjects[i].Qbypass = SystemObjects[item].Q_bypass[t]
                        SystemObjects[i].SSWbypass = SystemObjects[item].SS_W_Bypass[t]
                        SystemObjects[i].QReservoir = SystemObjects[item].Res_Q_in[t]
                        SystemObjects[i].SedReservoir = SystemObjects[item].Res_SS_W_in[t]
                        SystemObjects[i].VinReservoir = SystemObjects[item].Res_V_in[t]
                if SystemObjects[i].Element_Sub_Dict["Type"] == "Junction":
                    # If junction splits to 2+ downstream elements, if any are reservoirs, send back reservoir water level to junction
                    # for purposes of computing distribution of flow as a function of flow and downstream reservoir water level.
                    for item in SystemObjects[i].Element_Sub_Dict["Outflow Elements"]:
                        try:
                            SystemObjects[i].DS_res_WSE = SystemObjects[item].Element_Sub_Dict["Natural Bypass Water Level"][i]
                        except KeyError:
                            pass
                    Flow_in += SystemObjects[i].Q_incremental[t]  # Add incremental flows only if element is a junction
                    Sed_in += SystemObjects[i].Incremental_Sed_Load_Junction[t]
                if SystemObjects[i].Element_Sub_Dict["Type"] == "Reservoir":
                    # Execute all main routines for time t for reservoir from here.
                    SystemObjects[i].Master_Method_Caller(t, Flow_in, Sed_in, Flushing_Group_Dict)
                    Flushing_Group_Dict = SystemObjects[i].Flushing_Dictionary  # Update Flushing Group Dict [t] in case reservoirs are grouped.
                else:
                    # For all other elements, execute all primary routines of the element for time t from here.
                    SystemObjects[i].Master_Method_Caller(t, Flow_in, Sed_in)

                # Now that element's sediment mass balance has been simulated, if there's any dredged sediment, send it to destination
                # element's BS_W
                try:
                    SystemObjects[SystemObjects[i].Element_Sub_Dict["Dredging Outflow Element"]].BS_W[t + 1] += \
                    SystemObjects[i].Operating_Policy["Dredging"].Sediment_Load_Removed_Daily[t]
                except KeyError:
                    pass  # no dredging exists for reservoir i

            for i in SystemObjects["Ordered Simulation List"]:
                for j in SystemObjects["Ordered Simulation List"]:
                    try:
                        # See if object j has data to transfer to object i. Transfer if so. In case element that needs to transfer
                        # data is downstream and hence hasn't yet been simuluated, need to check this again later.
                        SystemObjects[i].elem_xfer_input_dict[j] = SystemObjects[j].elem_xfer_output_dict[i]
                        # Before proceeding to next time step, for all reservoirs, re-run energy calculations to account for flow
                        # at downstream junction impacting tailwater, if this reservoir is one that has a downstream junction
                        # with multiple inflow elements.
                        if SystemObjects[i].Element_Sub_Dict["Type"] == "Reservoir":
                            if SystemObjects[i].re_calc_energy == 1:
                                SystemObjects[i].Import_External_Element_State(t)
                                SystemObjects[i].Hydropower_Calculations(t)
                    except KeyError:
                        pass

        # If this is the first simulation, can now initialize an output dictionary of dataframes.
        if rz == 0:
            num_states = {}  # Stores num. state variables that are stored for each system element. Each element type has diff. number.
            state_list = {}  # Stores list of state var names for each object instance. Can only be done upon completion of simulation
            state_list_excel = {}  # Stores list of state var names for each object instance, truncated so name will fit into a worksheet name.
            if element_export_list is None:
                # User did not indicate for which elements to export output, so export all elements.
                element_export_list = SystemObjects["Ordered Simulation List"]
            else:
                # User did indicate for which elements to export output. To avoid misspellings, only export those
                # elements that were spelled correctly (according to names defined in Network Connectivity sheet).
                element_export_list = list(set(element_export_list) & set(SystemObjects["Ordered Simulation List"]))
            for i in element_export_list:
                outlets_list = []
                if var_sub_list is None:
                    # User did not indicate for which elements to export output, so export all time series state
                    # variables that apply to each element.
                    user_export_indicated = 0  # User did not indicate what to export.
                    var_list = SystemObjects[i].__dict__.keys()
                else:
                    # User did indicate for which variables to export output. To avoid misspellings, only store in
                    # Time_Series_Output_Dictionary those variables that were spelled correctly.
                    user_export_indicated = 1  # User did indicate what variables to export.
                    var_list = list(set(SystemObjects[i].__dict__.keys()) & set(var_sub_list))
                    # Account for additional time series variables user may want to export not stored in __dict__,
                    # and add back into var list.
                    try:
                        outlets_list = list(set(SystemObjects[i].Orifices.__dict__.keys()) & set(var_sub_list))
                        var_list += outlets_list
                    except AttributeError:
                        pass  # Element has no orifice attribute.
                num_states[i] = 0  # Initialize before counting how many time series variables exist.
                state_list[i] = []  # Initialize empty list to store state variables for each object
                state_list_excel[i] = []  # Initialize empty list to store state variables for each object
                for item in range(len(var_list)):
                    # Use correct object to export time series (orifices require opening up a reservoir object to
                    # grab underlying time series)
                    if var_list[item] not in outlets_list:
                        object_switcher = SystemObjects[i]
                    else:
                        object_switcher = SystemObjects[i].Orifices

                    if type(getattr(object_switcher, var_list[item])) == np.ndarray:
                        if user_export_indicated == 1:
                            # Assume user has correctly specified a time series style variable's name.
                            num_states[i] += 1
                            state_list[i].append(var_list[item])  # Add variable name to time series state list
                            state_list_excel[i].append(var_list[item][0:30])  # Only keep first 30 charac, as name will be exported to excel.
                        else:
                            # PySedSim is internally selecting time series variables for export. Need to make sure
                            # those variables selected are actually time series, and not just numpy arrays.
                            if (len(getattr(object_switcher, var_list[item])) == Sim_Dur or len(
                                    getattr(object_switcher, var_list[item])) == Sim_Dur + 1):
                                num_states[i] += 1
                                state_list[i].append(var_list[item])  # Add variable name to time series state list
                                state_list_excel[i].append(var_list[item][0:30])  # Only keep first 30 charac, as name will be exported to excel.
                Time_Series_Output_Dictionary[i] = pd.Panel(np.zeros((num_states[i], Sim_Dur, Num_Realizations)), items=state_list[i],
                                                            major_axis=simulation_dates, minor_axis=Col_Names)

        # Simulation of length Sim_Dur is now complete for Realization_i. Store data in output dict/dataframe.
        # Loop through SystemObjects and states/variables to store for this realization.
        if Simulation_mode == 'debug':
            Output_Object_Dict[rz] = {}  # Initialize sub-dict for each realization, if in debug mode.
        # Loop through each object's attributes, identify time series arrays, and store them in the Time Series Output Dictionary.
        outlets_list = ['Q_downstream', 'Q_overflow', 'Q_diversion', 'Q_low_level_outlet', 'Q_controlled',
                        'Q_turbines', 'Q_mid_level_outlet']
        for element in element_export_list:
            for i in state_list[element]:
                # Use correct object to export time series (orifices require opening up a reservoir object to
                # grab underlying time series)
                if i not in outlets_list:
                    object_switcher = SystemObjects[element]
                else:
                    object_switcher = SystemObjects[element].Orifices

                if len(getattr(object_switcher, i)) == Sim_Dur:
                    # Attribute is an array of length Sim_Dur, so handle it accordingly.
                    Time_Series_Output_Dictionary[element][i][Col_Names[rz]] = getattr(object_switcher, i)[0:Sim_Dur]
                elif len(getattr(object_switcher, i)) == (Sim_Dur + 1):
                    # Attribute is an array of length Sim_Dur + 1, so handle it accordingly.
                    Time_Series_Output_Dictionary[element][i][Col_Names[rz]] = getattr(object_switcher, i)[1:Sim_Dur+1]
                else:
                    pass  # User does not want to store this array in output dictionary.

            # Store SystemObjects output in a dictionary, if in debug mode
            if Simulation_mode == 'debug':
                Output_Object_Dict[rz][element] = SystemObjects[element]

        print("Simulation Realization %s Complete.") % (rz+1)

    return state_list_excel, Time_Series_Output_Dictionary, Output_Object_Dict