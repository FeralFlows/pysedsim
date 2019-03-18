# -*- coding: utf-8 -*-

'''
Module creates the network of system elements (reaches, reservoirs, and junctions) that will be simulated.
'''

# Import relevant libraries
from __future__ import division
from pysedsim.river_basin_elements.reservoir import Reservoir
from pysedsim.river_basin_elements.channel import Channel
from pysedsim.river_basin_elements.channel import Diversion_Channel
from pysedsim.river_basin_elements.junction import Junction
from pysedsim.river_basin_elements.bypassing import Bypass_Structure
import logging

def Simulation_Network_Creation(Input_Data_File):

    # Purpose: To import the user's simulation preferences.
    # Inputs:
    # (1) Input_Data_File is the input data file that was created in the import_simulation_preferences module.

    # Create a network connectivity structure from user-defined elements in "Network connectivity" worksheet of input file
    [Element_Dict,Ordered_Object_List] = Network_Connectivity(Input_Data_File)

    # Create system object dictionary in which time series output will be stored for each realization (e.g., a 100-year simulation)
    Time_Series_Output_Dictionary = {}  # All relevant time series output from each simulation (resulting from each realization) will be
    # stored in a dataframe column for that realization. Each key points to a dataframe.
    Output_Object_Dict = {}  # Stores all SystemObject information for each realization.

    # List of system objects to be simulated. Ordered list of elements for simulation is the only non- simulation object member of the
    # dictionary. Rest of dictionary entries are the system simulation objects.
    SystemObjects = {}
    SystemObjects["Ordered Simulation List"] = Ordered_Object_List

    element_stochastic_components = {}  # This is really used for MC simulations, but it needs to be defined regardless of whether
    # simulation is stochastic or deterministic.  Will store for each element a list of MC sheets to be imported from Paramaters_for_MC.csv
    for item in SystemObjects["Ordered Simulation List"]:
        # For each element, stores what parameters are going to be stochastically sampled. Keys are element names, known only if they are defined
        # in the Stochastic_Sim section below.
        element_stochastic_components[item] = []  # Will be populated later if user has actually specified stochastic parameters.

    Flushing_Group_Dict = {"Flushing Group ID": {}, "Non Drawdown Flushing Group ID": {}, "Currently Flushing": {}}

    Parameter_Input_Dictionary = {}  # Define dictionary

    return SystemObjects, Time_Series_Output_Dictionary, Output_Object_Dict, Element_Dict, Flushing_Group_Dict, \
           element_stochastic_components, Parameter_Input_Dictionary

def Network_Connectivity(Input_Data_File):
    NetConMat_Sheet = Input_Data_File['Network connectivity']  # sheet to search: Network connectivity matrix
    Element_Dict = {}  # Dictionary that stores all basic network connectivity data for all user-defined system elements
    # Permit a wide variety of user inputted names to be allowed:
    Valid_Element_Names = ["ReachElement", "Reach Element", "Reach", "Channel", "ReservoirElement", "Reservoir Element", "Reservoir", "JunctionElement", "Junction Element", "Junction", "DivertedOutletElement", "Diverted Outlet Element", "Diversion Element", "Diversion", "Bypass Structure", "Bypass", "BypassElement", "Bypass Element"]
    Valid_Element_Names_Dict = {'Reaches': ["ReachElement", "Reach Element", "Reach", "Channel"], 'Reservoirs': ["ReservoirElement", "Reservoir Element", "Reservoir"], 'Junctions': ["JunctionElement", "Junction Element", "Junction"], 'Diversions': ["DivertedOutletElement", "Diverted Outlet Element", "Diversion Element", "Diversion"], 'Bypass Structures': ["Bypass Structure", "Bypass", "BypassElement", "Bypass Element"]}
    start_date = Input_Data_File['Simulation Specifications']['B3'].value
    Num_Res = 0
    Num_Reach = 0
    Num_Jct = 0
    Num_Div = 0
    Num_Byp = 0
    i = 1
    assign_element_ID = 1 # assume model should assign element ID values unless spreadsheet check reveals user has already specified them.
    while Input_Data_File['Network connectivity'].cell(row = i, column = 1).value is not None:
        current_spreadsheet_value = NetConMat_Sheet.cell(row = i, column = 1).value
        if current_spreadsheet_value in Valid_Element_Names:
            Element_Name = NetConMat_Sheet.cell(row = i, column = 3).value
            Element_Dict[Element_Name] = {} # Create a nested dictionary for this element, to store its network connectivity data
            Element_Dict[Element_Name]["Simulation Start Date"] = start_date # Create a nested dictionary for this element, to store its network connectivity data
            if i == 1 and NetConMat_Sheet.cell(row = i, column = 2).value is not None:
                # User intends to assign ID values to all elements
                assign_element_ID = 0
            if assign_element_ID == 1:
                Element_Dict[Element_Name]["ID"] = i # Assign ID value to element.
            else:
                Element_Dict[Element_Name]["ID"] = NetConMat_Sheet.cell(row = i, column = 2).value # Use user assigned ID value for element
            Element_Dict[Element_Name]["Input File Row Location"] = i # Store location of element in worksheet
            # Store element type:
            if current_spreadsheet_value in Valid_Element_Names_Dict['Reservoirs']:
                Element_Dict[Element_Name]["Type"] = "Reservoir"
                Num_Res += 1
            elif current_spreadsheet_value in Valid_Element_Names_Dict['Reaches']:
                Element_Dict[Element_Name]["Type"] = "Reach"
                Num_Reach += 1
            elif current_spreadsheet_value in Valid_Element_Names_Dict['Junctions']:
                Element_Dict[Element_Name]["Type"] = "Junction"
                Num_Jct += 1
            elif current_spreadsheet_value in Valid_Element_Names_Dict['Diversions']:
                Element_Dict[Element_Name]["Type"] = "Diversion"
                Num_Reach += 1  # A diversion is a sub-type of reach
                Num_Div += 1
            elif current_spreadsheet_value in Valid_Element_Names_Dict['Bypass Structures']:
                Element_Dict[Element_Name]["Type"] = "Bypass Structure"
                Num_Byp += 1
            # Initialize element lists (need to initialize lists before you can store data in them)
            Element_Dict[Element_Name]["Inflow Elements"] = []  # List of all inflow elements of current element
            Element_Dict[Element_Name]["Inflow Elements Type"] = []  # List of all inflow element types
            Element_Dict[Element_Name]["Outflow Elements"] = []  # List of all outflow elements for current element
            Element_Dict[Element_Name]["Outflow Elements Type"] = []  # List of all outflow element types (same order as list of outflow elements)
            Element_Dict[Element_Name]["Exception Dictionary"] = {}  # For junctions; information will be stored regarding upstream diversion dams and/or bypass structures for use during outflow element allocation.
            Element_Dict[Element_Name]["Dredging Inflow Elements"] = []
            Element_Dict[Element_Name]["Simulation Spot"] = 0
        i += 1
    Num_elem = Num_Res + Num_Reach

    # Run a second loop to store inflow and outflow elements for each element. Did not perform this above because all junction elements needed to be defined first to run the loop below.
    # Store the inflow/outflow elements for every element type (including junction element types).
    for key in Element_Dict:
        row_start_location = Element_Dict[key]["Input File Row Location"]
        exit_loop = 0
        counter = 0
        if Element_Dict[key]["Type"] in ["Reach", "Reservoir", "Diversion", "Bypass Structure"]:
            # Not a junction so inflow and outflow elements are required to be specified by user in network connectivity matrix.
            while exit_loop == 0:
                counter = counter + 1
                if NetConMat_Sheet.cell(row = row_start_location + counter, column = 1).value in ["Inflow Node", "Inflow Junction"]:
                    if type(NetConMat_Sheet.cell(row = row_start_location + counter, column = 4).value) in [str, unicode]:
                        # If inflow junction is a called out by its name (string or unicode),
                        Element_Dict[key]["Inflow Elements"].append(NetConMat_Sheet.cell(row = row_start_location + counter, column = 4).value)  # junction is upstream of element of interest
                        Element_Dict[NetConMat_Sheet.cell(row = row_start_location + counter, column = 4).value]["Outflow Elements"].append(key)  # by definition, element of interest is thus downstream of junction
                        Element_Dict[NetConMat_Sheet.cell(row = row_start_location + counter, column = 4).value]["Outflow Elements Type"].append(Element_Dict[key]["Type"])  # by definition, element of interest is thus downstream of junction
                    else:
                        # Inflow Junction for element key is called out by its ID number (int) rather than its name.
                        for elem in Element_Dict:
                            if Element_Dict[elem]["ID"] == NetConMat_Sheet.cell(row = row_start_location + counter, column = 4).value:
                                Element_Dict[key]["Inflow Elements"].append(elem)
                                Element_Dict[elem]["Outflow Elements"].append(key)  # by definition, element of interest is thus downstream of junction
                                Element_Dict[elem]["Outflow Elements Type"].append(Element_Dict[key]["Type"])
                    Element_Dict[key]["Inflow Elements Type"].append("Junction")  # Only possible upstream type
                elif NetConMat_Sheet.cell(row = row_start_location + counter, column = 1).value in ["Outflow Node", "Outflow Junction"]:
                    if type(NetConMat_Sheet.cell(row = row_start_location + counter, column = 5).value) in [str, unicode]:  # If outflow element is a called out by its name (string),
                        Element_Dict[key]["Outflow Elements"].append(NetConMat_Sheet.cell(row = row_start_location + counter, column = 5).value)  # junction is downstream of element of interest
                        Element_Dict[key]["Outflow Elements Type"].append(Element_Dict[NetConMat_Sheet.cell(row = row_start_location + counter, column = 5).value]["Type"])  # junction is downstream of element of interest
                        Element_Dict[NetConMat_Sheet.cell(row = row_start_location + counter, column = 5).value]["Inflow Elements"].append(key)  # by definition, element of interest is thus upstream of junction
                        Element_Dict[NetConMat_Sheet.cell(row = row_start_location + counter, column = 5).value]["Inflow Elements Type"].append(Element_Dict[key]["Type"])
                    else:
                        # Outflow Junction for element key is called out by its ID number (int) rather than its name.
                        for elem in Element_Dict:
                            if Element_Dict[elem]["ID"] == NetConMat_Sheet.cell(row = row_start_location + counter, column = 5).value:
                                Element_Dict[key]["Outflow Elements"].append(elem)
                                Element_Dict[key]["Outflow Elements Type"].append(Element_Dict[elem]["Type"])
                                Element_Dict[elem]["Inflow Elements"].append(key)  # by definition, element of interest is thus downstream of junction
                                Element_Dict[elem]["Inflow Elements Type"].append(Element_Dict[key]["Type"])  # by definition, element of interest is thus downstream of junction
                else:
                    exit_loop = 1

    # For every junction, if there's a diversion at its outlet, then the reservoir at the inflow needs to be labeled a diversion. This is to avoid the user needing to specify in the "Reservoir Specifications" sheet what the reservoir's specific functionality is (i.e., Diversion). If it's a diversion reservoir, we can tell from netconmat, so don't need user to specify this.
    for key in Element_Dict:
        if Element_Dict[key]["Type"] == "Junction":
            for item in Element_Dict[key]["Outflow Elements"]:
                if len(Element_Dict[key]["Outflow Elements"])>1:
                    if Element_Dict[item]["Type"] == "Reservoir":
                        Element_Dict[item]["Natural Bypass Water Level"] = {key: 0}
                if Element_Dict[item]["Type"] == "Diversion":
                    for item2 in Element_Dict[key]["Inflow Elements"]:
                        if Element_Dict[item2]["Type"] == "Reservoir":
                            Element_Dict[item2]["Reservoir Type"] = "Diversion"

    # For every junction, if there's a bypass structure or a diversion dam at its inlet, store that dam's name in a dictionary so it can be looked up when the junction is later allocating outflows to its downstream elements.
    for key in Element_Dict:
        if Element_Dict[key]["Type"] == "Junction":
            for item in Element_Dict[key]["Inflow Elements"]:
                if Element_Dict[item]["Type"] == "Diversion":
                    Element_Dict[key]["Exception Dictionary"]["Diversion"] = item
                elif Element_Dict[item]["Type"] == "Bypass Structure":
                    Element_Dict[key]["Exception Dictionary"]["Bypass"] = item

    # Call Element_Simulation_Order method to give each element an order in the simulation, to be stored in the Element_Dict
    [Element_Dict, ordered_element_list] = Element_Simulation_Order(Element_Dict, Num_elem, Num_Jct, Num_Byp)
    return Element_Dict, ordered_element_list

def Element_Simulation_Order(Element_Dict, Num_elem, Num_Jct, Num_Byp):
    # Description: This function determines the order in which elements will be simulated, to be sure each element is not simulated before all of its upstream elements are simulated.
    # Function inputs:
    # Element_Dict: Dictionary (creation starts in Network_Connectivity function) that stores basic network connectivity information for each system element so they can all be simulated as a connected network of elements.
    # Num_elem: Number of reaches and reservoirs contained in Element_Dict
    # Num_Jct: Number of junctions (reaches and reservoirs) contained in Element_Dict

    ordered_element_list = [0 for x in range(Num_elem + Num_Jct + Num_Byp)]  # Initialize matrix that indicates whether calcs have been done for the reach/reservoir of interest in time period.
    for m in range(1, Num_elem + Num_Jct + Num_Byp + 1):
        # Run this loop with Num_elem + Num_Jct + Num_Byp iterations, so each element gets a spot in the simulation order.
        for key1 in Element_Dict:
            if Element_Dict[key1]["Simulation Spot"] == 0:  # we are on an element that has not been done yet.
                skip = 0
                for key2 in Element_Dict:  # Now see if all elements that run into the element key1 have been done. If not, have to move on.
                    # if outflow junctions of ANY key2 is in outflow junctions of any key1, see if the element key2 has a simulation spot assigned. if any dont have spots assigned, break the loop.
                    if key1 not in Element_Dict[key2]["Outflow Elements"]:
                        pass  # key1 not in key2 Element's dictionary, so skip it.
                    else:
                        if Element_Dict[key2]["Simulation Spot"] != 0:
                            pass  # Continue
                        else:  # This other element key2 that flows into key1 hasn't been done yet. Break loop.
                            skip = 1
                            break
                if skip == 0:  # if skip = 0, this will be the next element added to the simulation list
                    Element_Dict[key1]["Simulation Spot"] = m
                    ordered_element_list[m-1] = key1  # Use this list to create instances of each reach or reservoir class, so the objects are all alredy in the correct order.
                    break
    return Element_Dict, ordered_element_list

def System_Object_Creation(T, Input_Data_File, element_stochastic_components, SystemObjects, Element_Dict, Flushing_Group_Dict,
                           stochastic_flow_list, op_policy_params=None):

    # Purpose: To create all reservoir/reach/junction/diversion objects (instances of classes, e.g., Reservoir)

    # Inputs:
    # (1) All inputs created automatically within SedSim_v2.py using various functions called there.

    # Create all system objects (reservoirs, channels, diversions, junctions), already loaded in their correct simulation order.
    for i in SystemObjects["Ordered Simulation List"]:
        if Element_Dict[i]["Type"] == "Reservoir":
            SystemObjects[i] = Reservoir(i, T, Input_Data_File, Element_Dict[i],
                                         stochastic_components=element_stochastic_components[i],
                                         op_policy_params=op_policy_params)
            # If the simulated system will have reservoirs with flushing performed, need to create a dictionary storing flushing group ID
            #  and non_drawdown group ID, to be fed back to each reservoir and stored there.
            if SystemObjects[i].Flushing_Group_ID is not None:
                # If reservoir has a flushing group ID, create such a key in the Flushing ID dictionary.
                if SystemObjects[i].Flushing_Group_ID not in Flushing_Group_Dict["Flushing Group ID"]:
                    # Store reservoir names with this ID in a list.
                    Flushing_Group_Dict["Flushing Group ID"][SystemObjects[i].Flushing_Group_ID] = []
                    Flushing_Group_Dict["Flushing Group ID"][SystemObjects[i].Flushing_Group_ID].append(SystemObjects[i].name)
                else:
                    # Key ID already exists. Add Reservoir name to list.
                    Flushing_Group_Dict["Flushing Group ID"][SystemObjects[i].Flushing_Group_ID].append(SystemObjects[i].name)
            if SystemObjects[i].Non_Flushing_Group_ID is not None:
                # If reservoir has a flushing group NON-DRAWDOWN ID, create such a key in the Non-drawdown Flushing ID dictionary.
                if SystemObjects[i].Non_Flushing_Group_ID not in Flushing_Group_Dict["Non Drawdown Flushing Group ID"]:
                    # Store reservoir names with this ID in a list.
                    Flushing_Group_Dict["Non Drawdown Flushing Group ID"][SystemObjects[i].Non_Flushing_Group_ID] = []
                    Flushing_Group_Dict["Non Drawdown Flushing Group ID"][SystemObjects[i].Non_Flushing_Group_ID].append(
                        SystemObjects[i].name)
                else:
                    # Key ID already exists. Add Reservoir name to list.
                    Flushing_Group_Dict["Non Drawdown Flushing Group ID"][SystemObjects[i].Non_Flushing_Group_ID].append(
                        SystemObjects[i].name)
        elif Element_Dict[i]["Type"] == "Reach":
            SystemObjects[i] = Channel(i, T, Input_Data_File, Element_Dict[i], element_stochastic_components[i])
        elif Element_Dict[i]["Type"] == "Junction":
            SystemObjects[i] = Junction(i, T, Input_Data_File, Element_Dict[i], element_stochastic_components[i], stochastic_flow_list)
        elif Element_Dict[i]["Type"] == "Diversion":
            SystemObjects[i] = Diversion_Channel(i, T, Input_Data_File, Element_Dict[i], element_stochastic_components[i])
        elif Element_Dict[i]["Type"] == "Bypass Structure":
            SystemObjects[i] = Bypass_Structure(i, T, Input_Data_File, Element_Dict[i])

    # Store the Flushing Group Dictionary as part of each reservoir's attributes
    for i in SystemObjects["Ordered Simulation List"]:
        if Element_Dict[i]["Type"] == "Reservoir":
            SystemObjects[i].Flushing_Dictionary = Flushing_Group_Dict

    return SystemObjects, Flushing_Group_Dict, element_stochastic_components