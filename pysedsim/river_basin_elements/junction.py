# -*- coding: utf-8 -*-
'''

Module to define the Junction class and methods.

'''

# import relevant modules
from __future__ import division  # Ensures result of quotient of two integers will be a float, not an integer.
import numpy as np
import pysedsim.data_processing.data_processing
from pysedsim.river_basin_elements.storage_element import Storage_Element
from pysedsim.data_processing.matrix_interpolation import Matrix_Interpolation
from pysedsim.data_processing.matrix_interpolation import Return_Closest_Values
#from sklearn import linear_model

class Junction(Storage_Element):
    def __init__(self, name, T, Input_Data_File, Element_Sub_Dict, stochastic_components, stochastic_flow_list):
        '''
        Purpose: Constructor for Junction class

        :param name: user-defined name of element
        :param T: Simulation Duration (integer, days)
        :param Input_Data_File: openpyxl MS Excel (.xlsx) object input file for current simulation run
        :param Element_Sub_Dict: Dictionary for element that contains basic data (e.g., name, type, connectivity)
        :param stochastic_components: Dictionary of Monte Carlo parameter assumptions relevant to class
        :param stochastic_flow_list: List of junction names (strings) for which stochastic inflows will be provided
        for simulations in a Monte Carlo simulation setting.
        :return:
        '''

        # Note: The order in which method calls appear below is important. They should not be rearranged.
        if hasattr(Storage_Element, '__init__'):
            Storage_Element.__init__(self, name, T, Element_Sub_Dict)   # Call parent constructor.
        self.name = name
        self.Input_Data_File = Input_Data_File  # for use of file name outside of data import and array initialization methods.
        # Initialize as arrays critical reservoir attributes.
        Junction.Array_Initialization(self, T)
        self.ID = None
        Junction.Import_Data(self, T, Input_Data_File, stochastic_components, self.Element_Sub_Dict["Simulation Start Date"],
                             stochastic_flow_list)

    def Master_Method_Caller(self, t, Flow_in, Sed_in):
        '''
        Purpose: Called in every time step, this method calls all relevant methods in this class in the correct order.

        :param t: Current time period (integer, value range: [0,T])
        :param Flow_in: Inflow rate (m^3/s) during time period t
        :param Sed_in: Sediment inflow (kg) during time period t
        :return: Nothing (calls class methods that set attribute values)
        '''

        Storage_Element.Master_Method_Caller(self, t)  # Increments date
        Junction.Import_External_Element_State(self, t)
        Storage_Element.Element_Inflows(self, t, Flow_in, Sed_in)
        Junction.Mass_Balance(self, t)
        Junction.Outflow_Element_Allocation(self, t)
        Junction.Element_Data_Transfer(self, t)

    def Array_Initialization(self, T):
        '''
        Purpose: Initializes variable arrays, called once (as part of class constructor).

        :param T: Simulation Duration (days)
        :return:
        '''

        self.Q_in = np.zeros(T)  # Total daily junction inflow rate (m3/s), including all upstream elements and incremental flows
        self.Q_incremental = np.zeros(T)
        self.Incremental_Sed_Load_Junction = np.zeros(T)
        self.random_concentration_array = np.zeros(T)  # Store normally or lognormally dist. r.v. values for sediment load generation.
        self.concentration_deterministic = np.zeros(T)  # Store normally or lognormally dist. r.v. values for sediment load generation.
        self.Incremental_Sed_Load_Junction_Deterministic = np.zeros(T)  # Store daily sediment loads if no uncertainty in rating-curve
        # based sediment load prediction.
        self.Cum_Annual_SED_LOAD = 0
        self.Annual_SED_LOAD = 0
        self.sed_alpha = 0
        self.sed_beta = 0
        self.DS_res_WSE = 0  # Stores downstream reservoir's water level at time t (end of previous time period)
        self.US_res_Q_Diversion = 0  # Stores upstream reservoir's releases into any diverted outlets.
        self.fraction = np.zeros(T)
        self.cum_fraction = np.zeros(T)  # Stores cumulative fraction of junction inflow that has been assigned to downstream elements.

    def Import_Data(self, T, Input_Data_File, stochastic_components, start_date, stochastic_flow_list):
        '''
        Purpose: Imports all user-specified data relevant to the channel segment of interest.

        :param T: Simulation Duration (days)
        :param Input_Data_File: openpyxl MS Excel (.xlsx) object input file for current simulation run
        :param start_date: date [string] in MM/DD/YYYY format, used to specify start date of time series data import.
        :param stochastic_flow_list: List of junction names (strings) for which stochastic inflows will be provided
        for simulations in a Monte Carlo simulation setting.
        :param stochastic_components: Dictionary of Monte Carlo parameter assumptions relevant to junction class
        :return: Nothing (sets self attributes)
        '''

        if ('Incremental Flows' in Input_Data_File.sheetnames):
            if self.name not in stochastic_flow_list:
                # This junction does not have stochastic flows specified in a file. Therefore, take any provided deterministic time
                # series of inflows and use those if they exist.
                self.Q_incremental = data_processing.Excel_Data_Import(self.name, Input_Data_File, 'Incremental Flows',
                                                                       0, T, max_distinct_data_types=None,
                                                                       data_name_offset=None,
                                                                       start_date=start_date)  # Required input
                self.Q_out_unreg_for_param_calib = self.Q_incremental
        if 'Sediment Loads' in Input_Data_File.sheetnames:
            [self.Cum_Annual_SED_LOAD, self.sed_alpha, self.sed_beta,
             self.calibration_preference] = data_processing.Excel_Data_Import(self.name, Input_Data_File,
                                                                              'Sediment Loads', 1, 4,
                                                                              max_distinct_data_types=None,
                                                                              data_name_offset=None)
            # Call of Excel_Data_Import if junction name isn't listed will return "None". This instead needs to be a zero value.
            if self.Cum_Annual_SED_LOAD is None:
                self.Cum_Annual_SED_LOAD = 0
            self.Annual_SED_LOAD = self.Cum_Annual_SED_LOAD # Set initial value to this, but later it will be reduced by any and all upstream cumulative loads that exist.
        else:
            if 'Incremental Sediment Loads' not in Input_Data_File.sheetnames:
                print "Error: No 'Sediment Loads' worksheet is provided. This sheet is required since incremental " \
                      "sediment loads are provided iun the 'Incremental Sediment Loads' worksheet"

        # Incremental sediment load calibration preferences
        if self.calibration_preference == "Calibrate coefficient":
            self.calibration_preference = 1
        elif self.calibration_preference == "Specify coefficient":
            self.calibration_preference = 2  # Do nothing else; alpha and beta values have already been imported.
            if self.sed_alpha >= 0 and self.sed_beta >= 0:
                pass
            else:
                print "Error: Incremental sediment load parameters not properly specified for element %s in Sediment " \
                      "Loads worksheet" % self.name
        elif self.calibration_preference == "Specify daily incremental sediment loads (kg/day)":
            self.calibration_preference = 3  # sediment loads have already been specified
        else:
            # User did not specify calibration preferences
            if self.sed_alpha > 0 and self.sed_beta > 0:
                self.calibration_preference = 2  # User specified parameter values, so use those.
            else:
                if self.Annual_SED_LOAD > 0:
                    self.calibration_preference = 1
                else:
                    self.calibration_preference = 4  # No preferences or parameters specified. No sediment loads will be produced here.

        # If user wishes to specify sediment loads externally, import them now.
        if self.calibration_preference == 3 or self.calibration_preference == 4:
            self.sed_alpha = 0  # No cumulative annual sediment load specified. Set alpha = beta = 0.
            self.sed_beta = 0  # No cumulative annual sediment load specified. Set alpha = beta = 0.
            if self.calibration_preference == 3:
                if 'Incremental Sediment Loads' in self.Input_Data_File.sheetnames:
                    self.Incremental_Sed_Load_Junction = data_processing.Excel_Data_Import(self.name,
                                                                                           self.Input_Data_File,
                                                                                           'Incremental Sediment Loads',
                                                                                           0, T,
                                                                                           max_distinct_data_types=None,
                                                                                           data_name_offset=None,
                                                                                           start_date=start_date)
                else:
                    print "Error: No 'Incremental Sediment Loads' worksheet is provided. This sheet is required given " \
                          "the sediment load preferences specified in the 'Sediment Loads' worksheet."

        self.num_branches = len(self.Element_Sub_Dict['Outflow Elements'])  # number of outflow elements per junction
        # Set number of effective branches as number of branches minus 1 if there is a diversion channel that has been allocated flow
        # from the reservoir before it can reach the downstream channel.
        if 'Diversion' in self.Element_Sub_Dict['Outflow Elements Type']:
            self.num_effective_branches = self.num_branches - 1
        else:
            self.num_effective_branches = self.num_branches

        self.flow_fraction_specified = 0  # 0= No (Default value), 1=yes. Gets reset below if user specifies fractions.

        # If only 1 downstream element, do not need to distribute flow among downstream elements.
        # Initialize dictionary storing fraction of flow that makes it downstream.
        if self.num_branches == 1:
            self.fraction = {ds_elem_name: np.ones(T) for ds_elem_name in self.Element_Sub_Dict['Outflow Elements']}
        else:
            self.fraction = {ds_elem_name: np.zeros(T) for ds_elem_name in self.Element_Sub_Dict['Outflow Elements']}

        if self.num_effective_branches > 1:
            if 'Reservoir' not in self.Element_Sub_Dict['Outflow Elements Type']:
                Junction.Junction_Flow_Dist_Sheet(self, Input_Data_File)  # See if user has specified jct fraction.
            else:
                # A reservoir is downstream of junction, along with at least one other element (e.g., channel).
                if 'Reservoir Natural Bypass' in Input_Data_File.sheetnames:
                    Junction.Res_Nat_Byp_Sheet(self, Input_Data_File)
                    if all(item == [] for item in self.Jct_Flow_Distribution):
                        Junction.Junction_Flow_Dist_Sheet(self, Input_Data_File)  # See if user has specified jct fraction.
                else:
                    Junction.Junction_Flow_Dist_Sheet(self, Input_Data_File)  # See if user has specified jct fraction.

    def Junction_Flow_Dist_Sheet(self, Input_Data_File):
        if 'Junction Flow Distribution' in Input_Data_File.sheetnames:
            # Load flow distribution data for junction. Then parse the data so only the required portions remain.
            self.Jct_Flow_Distribution = data_processing.Excel_Data_Import(self.name, Input_Data_File,
                                                                            'Junction Flow Distribution', 2,
                                                                            self.num_effective_branches + 1,
                                                                            max_distinct_data_types=None,
                                                                            data_name_offset=1)
            if all(item != [] for item in self.Jct_Flow_Distribution):
                self.Jct_Flow_Dist_Dict = {key: [[] for i in range(2)] for key in self.Element_Sub_Dict['Outflow Elements']}
                for loc in range(self.num_effective_branches):
                    try:
                        self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 1][0]][0] = self.Jct_Flow_Distribution[0][2:]
                        self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 1][0]][1] = self.Jct_Flow_Distribution[loc + 1][2:]
                        self.flow_fraction_specified = 1  # User has specified junction fraction distribution.
                    except KeyError:
                        print "Error: provided outflow element for Junction %s in Junction Flow Distribution " \
                              "worksheet doesn't match a downstream element." % self.name
                        raise KeyError('Provided outflow element for Junction %s in Junction Flow. '
                                       'Distribution worksheet doesnt match a downstream element.' % self.name)

    def Res_Nat_Byp_Sheet(self, Input_Data_File):
        self.Jct_Flow_Distribution = data_processing.Excel_Data_Import(self.name, Input_Data_File,
                                                                       'Reservoir Natural Bypass', 2,
                                                                       2+self.num_effective_branches,
                                                                       max_distinct_data_types=None,
                                                                       data_name_offset=1) # 4 lists of data

        if all(item != [] for item in self.Jct_Flow_Distribution):
            self.Jct_Flow_Dist_Dict = {key: [[] for i in range(3)] for key in self.Element_Sub_Dict['Outflow Elements']}
            self.use_regr_for_frac = 0  # TEMPORARY
            if self.use_regr_for_frac == 1:
                self.regr = {}
            for loc in range(self.num_effective_branches):
                try:
                    self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 2][0]][0] = self.Jct_Flow_Distribution[0][2:]
                    self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 2][0]][1] = self.Jct_Flow_Distribution[1][2:]
                    self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 2][0]][2] = self.Jct_Flow_Distribution[loc + 2][2:]
                    self.flow_fraction_specified = 1  # User has specified junction fraction distribution.

                    if self.use_regr_for_frac == 0:
                        # Loop through downstream locations, and flow values for those locations,
                        # to store elevation-fraction values for that location as two lists.
                        self.Flow_Elev_Frac_Dict = {}
                        for keys in self.Jct_Flow_Dist_Dict:
                            counter = 0
                            self.Flow_Elev_Frac_Dict[keys] = {}
                            # Add flow values as keys to dictionary
                            for list_items in self.Jct_Flow_Dist_Dict[keys][0]:
                                if list_items not in self.Flow_Elev_Frac_Dict[keys].keys():
                                    # Create key for this flow that stores empty elevation-fraction lists
                                    self.Flow_Elev_Frac_Dict[keys][list_items] = [[], []]
                                # Store Elevation and fraction for this downstream element
                                self.Flow_Elev_Frac_Dict[keys][list_items][0].append(self.Jct_Flow_Dist_Dict[keys][1][counter])
                                self.Flow_Elev_Frac_Dict[keys][list_items][1].append(self.Jct_Flow_Dist_Dict[keys][2][counter])
                                counter += 1  # Step to next line of flow values
                    else:
                        x = np.vstack([np.array(self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 2][0]][0]),
                                       np.array(self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 2][0]][1])]).T
                        self.regr[self.Jct_Flow_Distribution[loc + 2][0]] = linear_model.LinearRegression()
                        self.regr[self.Jct_Flow_Distribution[loc + 2][0]].fit(x,self.Jct_Flow_Dist_Dict[self.Jct_Flow_Distribution[loc + 2][0]][2])
                except KeyError:
                    print "Error: provided outflow element for Junction %s in Junction Flow Distribution " \
                          "worksheet doesn't match a downstream element." % self.name
                    raise KeyError

    def Outflow_Element_Allocation(self, t):
        '''
        Purpose: Allocates the element outflow to each of the downstream elements

        :param t: Current time period (integer, value range: [0,T])
        :return:
        '''

        if len(self.Element_Sub_Dict["Daily Water Outflows"]) <= 1:
            # Allocate all outflow to next element's inflow. If only one downstream element, use Storage_Element()
            # class method.
            Storage_Element.Outflow_Element_Allocation(self, t)
        else:
            # There are multiple outflow elements. This is therefore a junction by definition (because it has multiple outflow
            # elements), so we must fractionally allocate element outflow to each downstream element. This code will be contained in a
            # junction.py class method.

            # If there is a bypass structure upstream of this junction, distribute flow into the bypass, and the rest
            # into the remaining channel/reservoir.
            try:
                if self.Element_Sub_Dict["Exception Dictionary"]["Bypass"] is not None:
                    # Loop through outflow elements; locate the channel; send bypass flows there; rest to reservoir. These values are set
                    #  for a junction in pysedsim_main_simulation.SedSim_Main_Simulation()
                    for i in range(len(self.Element_Sub_Dict["Outflow Elements"])):
                        if self.Element_Sub_Dict["Outflow Elements Type"][i] == "Reach":
                            # We assume here that a bypass structure is erected to bypass flow into a channel FROM a reservoir. We don't
                            # allow a bypass structure to divert flow from a main channel to another channel at this time.
                            self.Element_Sub_Dict["Daily Water Outflows"][self.Element_Sub_Dict["Outflow Elements"][i]] = self.Qbypass
                            self.Element_Sub_Dict["Daily Sediment Outflows"][self.Element_Sub_Dict["Outflow Elements"][i]] = self.SSWbypass
                        elif self.Element_Sub_Dict["Outflow Elements Type"][i] == "Reservoir":
                            # Distribute rest of bypass flow into reservoir
                            self.Element_Sub_Dict["Daily Water Outflows"][self.Element_Sub_Dict["Outflow Elements"][i]] = self.QReservoir
                            self.Element_Sub_Dict["Daily Sediment Outflows"][self.Element_Sub_Dict["Outflow Elements"][i]] = self.SedReservoir
                    return  # Exit method. Do not enter junction distribution code below.
            except KeyError:
                pass  # This junction does not have an upstream bypass

            # Loop through Junction Flow Distribution Dictionary keys, allocate flow to each element via self.fraction
            for keys in self.Element_Sub_Dict['Outflow Elements']:
                if self.flow_fraction_specified == 1:
                    if 'Reservoir' not in self.Element_Sub_Dict['Outflow Elements Type']:
                        self.fraction[keys][t] = Matrix_Interpolation(self.name, self.Jct_Flow_Dist_Dict[keys], "flow", "fraction",
                                                        self.Q_out[t])
                    else:
                        if self.use_regr_for_frac == 0:
                            # Find closest flow values in dictionary for each downstream junction.
                            # Loop through keys (in order) to interpolate.
                            flow_list = Return_Closest_Values(self.Flow_Elev_Frac_Dict[keys].keys(),
                                                              self.Q_out[t])
                            flow_list = sorted(flow_list)  # Sort flows in ascending order
                            # Interpolate within each of those flows.

                            interp_frac = []
                            for flow_val in flow_list:
                                interp_frac.append(Matrix_Interpolation(self.name, self.Flow_Elev_Frac_Dict[keys][
                                    flow_val], "elevation", "fraction", self.DS_res_WSE))
                            if len(flow_list) == 1:
                                # Flow identical to key value
                                self.fraction[keys][t] = interp_frac[0]
                            else:
                                # Need to interpolate across flow values.
                                self.fraction[keys][t] = interp_frac[0] + (self.Q_out[t] - flow_list[0]) * (
                                interp_frac[1] - interp_frac[0]) / (flow_list[1] - flow_list[0])
                        else:
                            self.fraction[keys][t] = self.regr[keys].predict([self.Q_out[t],self.DS_res_WSE])  # WILL BE AN ARRAY.
                else:
                    self.fraction[keys][t] = 1/self.num_effective_branches  # Divide flow equally among all branches

                if self.fraction[keys][t] > 1:
                    self.fraction[keys][t] = 1
                elif self.fraction[keys][t] < 0:
                    self.fraction[keys][t] = 0
                if self.cum_fraction[t] == 1:
                    self.fraction[keys][t] = 0
                else:
                    self.cum_fraction[t] += self.fraction[keys][t]
                    if self.cum_fraction[t] > 1:
                        self.fraction[keys][t] = 1-(self.cum_fraction[t]-self.fraction[keys][t])  # Junction flow has been overallocated to downstream elements. Adjust.

                # If there's a diversion, allocate all of the diversion water to the diversion. In this case, if there are downstream
                # channels, adjust available flow for those channels to reflect portion taken by diversion.
                if 'Diversion' in self.Element_Sub_Dict['Outflow Elements Type']:
                    # There is a diversion. See if this element is the diversion channel. Allocate water accordingly.
                    if self.Element_Sub_Dict['Outflow Elements Type'][
                        self.Element_Sub_Dict['Outflow Elements'].index(keys)] == 'Diversion':
                        self.Element_Sub_Dict["Daily Water Outflows"][keys] = self.US_res_Q_Diversion
                        self.Element_Sub_Dict["Daily Sediment Outflows"][keys] = self.SS_W_out[t] * (
                        self.US_res_Q_Diversion / self.Q_out[t])
                    else:
                        self.Element_Sub_Dict["Daily Water Outflows"][keys] = (self.Q_out[t]-self.US_res_Q_Diversion)*self.fraction[keys][t]
                        self.Element_Sub_Dict["Daily Sediment Outflows"][keys] = ((self.Q_out[t]-self.US_res_Q_Diversion)/self.Q_out[t])*self.SS_W_out[t]*self.fraction[keys][t]
                else:
                    self.Element_Sub_Dict["Daily Water Outflows"][keys] = self.Q_out[t]*self.fraction[keys][t]
                    self.Element_Sub_Dict["Daily Sediment Outflows"][keys] = self.SS_W_out[t]*self.fraction[keys][t]

    def Calibration_Incremental_Sediment_Load(self, stochastic_components, distribution_name=None):

        # This call cannot occur as part of __init, as ann sed loads arent available yet. So it is called from SedSim_v2.py.

        # Inputs:
        # stochastic_components: a list of what parameters/model components will be stochastic for this junction.
        # distribution_name = string name of probability distribution from which random sediment noise is to be generated (e.g., 'normal')
        #  If calibration_preference = 1, perform calibration of incremental sediment load parameters, and compute daily incremental
        # sediment loads.
        if self.calibration_preference == 1:
            if self.Annual_SED_LOAD > 0:
                self.sed_alpha = self.Annual_SED_LOAD / (sum(86400 * np.power(self.Q_incremental, (self.sed_beta + 1))) / (self.T / 365))
            else:
                print "Error: Sediment calibration requested for Junction %s, but no annual sediment load provided in Sediment Loads worksheet" % self.name
                self.sed_alpha = 0  # No cumulative annual sediment load specified. Set alpha = beta = 0.
                self.sed_beta = 0  # No cumulative annual sediment load specified. Set alpha = beta = 0.
        if self.calibration_preference in [1, 2]:
            if '2' in stochastic_components:
                self.random_concentration_array = np.maximum(0, getattr(np.random, distribution_name)(self.sed_alpha*np.power(self.Q_incremental, self.sed_beta), .020))
                self.concentration_deterministic = self.sed_alpha*np.power(self.Q_incremental, self.sed_beta)
                self.Incremental_Sed_Load_Junction_Deterministic = self.dt*self.SPD*self.Q_incremental*(self.sed_alpha*np.power(self.Q_incremental, self.sed_beta))  # Using provided mean and standard deviation, produces array of RVs of len=Sim_Dur.
                self.Incremental_Sed_Load_Junction = self.dt*self.SPD*self.Q_incremental*(self.random_concentration_array)  # Using provided mean and standard deviation, produces array of RVs of len=Sim_Dur.
            else:
                self.Incremental_Sed_Load_Junction = self.sed_alpha*self.dt*self.SPD*np.power(self.Q_incremental, self.sed_beta + 1)
        else:
            pass  # junction incremental sediment loads were already separtely specified and imported

    def Import_External_Element_State(self, t):
        '''
        Imports data from other elements being simulated (e.g., junctions, reaches and reservoirs) for use in the
        simulation of this element.

        Input data are stored in the elem_xfer_input_dict (which originally came from other elements'
        elem_xfer_output_dict). This module simply takes this information and converts it into a specific use for
        this element. For this reason, each element will have their own unique implementation of this method. This
        method is executed at the beginning of every time step, before any processes are simulated. The module must be run a second time
        for reservoir elements that have a junction downstream with tailwater impacted by multiple elements flowing into junction.

        :param t: Current time period (integer, value range: [0,T])
        :returns:
        '''

        # Loop through elements sending data, and unpack it to provide use for this element.
        for keys in self.elem_xfer_input_dict:
            for sub_keys in self.elem_xfer_input_dict[keys]:
                if sub_keys == "Reservoir Diversion Flow":
                    self.US_res_Q_Diversion = self.elem_xfer_input_dict[keys][sub_keys]

    def Element_Data_Transfer(self, t):
        '''
        Transfers data from reservoir element to other elements in the system whose methods may depend on the values
        of junction state variables.

        Method is performed at the end of every time step, storing information in the element_transfer dictionary
        that can then be sent to and used by other system objects. elem_xfer_output_dict is the primary vehicle
        through which system elements (reservoirs, channels, etc.) communicate.

        :param t: Current time period (integer, value range: [0,T])

        '''

        # Passes current junction flow to upstream element, if that upstream element is a reservoir.
        for item_num in range(len(self.Element_Sub_Dict["Inflow Elements"])):
            if self.Element_Sub_Dict["Inflow Elements Type"][item_num] == "Reservoir":
                self.elem_xfer_output_dict[self.Element_Sub_Dict["Inflow Elements"][item_num]] = {"Downstream Junction Flow": self.Q_out[t]}
                self.elem_xfer_output_dict[self.Element_Sub_Dict["Inflow Elements"][item_num]]["Num Inflow Elements"] = len(
                    self.Element_Sub_Dict["Inflow Elements"])
        for item_num in range(len(self.Element_Sub_Dict["Outflow Elements"])):
            if self.Element_Sub_Dict["Outflow Elements Type"][item_num] in ["Reservoir", "Reach"]:
                self.elem_xfer_output_dict[self.Element_Sub_Dict["Outflow Elements"][item_num]] = {"Upstream Junction Inflow Fraction": self.fraction[self.Element_Sub_Dict["Outflow Elements"][item_num]][t], "Upstream Junction Inflow": self.Q_in[t]}

    def Mass_Balance(self, t):
        '''
        Purpose: Performs a mass balance, to track net inflows/outflows of water/sediment.

        At a junction, no water or sediment mass is stored. Junctions are used only to distribute flows among and
        between reservoirs and river channel segments.

        :param t: Current time period (integer, value range: [0,T])
        :return:
        '''

        self.Q_out[t] = self.Q_in[t]
        self.SS_W_out[t] = self.SS_W_in[t]