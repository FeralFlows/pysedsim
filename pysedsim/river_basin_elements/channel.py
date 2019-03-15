# -*- coding: utf-8 -*-

'''
Module to define the (river) Channel class and sub-class diversion channel (sub-class of Channel), which are
sub-classes of the superclass Storage_Element.

The Channel class applies to any stretch of pressure conveyance or open channel flow that that conveys water
and/or sediment from one point to another in the modeled network. The diversion class is a type of channel class that
refers to diversions at the dam itself (where diversions can divert water elsewhere with or without first going
through a power house, and for bypass diversions at the upstream end of reservoir).
'''

# Import relevant modules
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer.
import numpy as np
from pysedsim.river_basin_elements.storage_element import Storage_Element
import pysedsim.data_processing.data_processing
import math

class Channel(Storage_Element):
    def __init__(self, name, T, Input_Data_File, Element_Sub_Dict, stochastic_components=None):
        '''
        Purpose: Constructor for Channel class

        :param name: User-assigned element name (string)
        :param T: Simulation Duration (integer, days)
        :param Input_Data_File: openpyxl MS Excel (.xlsx) object input file for current simulation run
        :param Element_Sub_Dict: Dictionary for element that contains basic data (e.g., name, type, connectivity)
        :param stochastic_components: Dictionary of Monte Carlo parameter assumptions relevant to channel class
        :return: Nothing (sets self attributes)
        '''
        # Note: The order in which method calls appear below is important. They should not be rearranged.
        if hasattr(Storage_Element, '__init__'):
            Storage_Element.__init__(self, name, T, Element_Sub_Dict)  # Call parent constructor.
        Channel.Array_Initialization(self, T)  # Initialize arrays channels will have.
        Channel.Import_Data(self, T, Input_Data_File)
        Storage_Element.Time_Zero_Initialization(self)

    def Master_Method_Caller(self, t, Flow_in, Sed_in):
        '''
        Purpose: Called in every time step, this method calls all relevant methods in this class in the correct order.

        :param t: Current time period (integer, value range: [0,T])
        :param Flow_in: Inflow rate (m^3/s) during time period t
        :param Sed_in: Sediment inflow (kg) during time period t
        :return: Nothing (calls class methods that set attribute values)
        '''

        Storage_Element.Master_Method_Caller(self, t)
        Channel.Import_External_Element_State(self, t)  # Imports relevant real-time data from other elements.
        Storage_Element.Element_Inflows(self, t, Flow_in, Sed_in)
        Channel.Flow_Routing(self, t)
        Channel.Sediment_Routing(self, t)
        Storage_Element.Outflow_Element_Allocation(self, t)  # Flow allocated to single downstream junction.
        Storage_Element.Post_Mass_Balance_Calcs(self, t)

    def Array_Initialization(self, T):
        '''
        Purpose: Initializes variable arrays, called once (as part of class constructor).

        :param T: Simulation duration (days)
        :return: Nothing (sets self attributes)
        '''
        self.TS_surplus_deficit = np.zeros(T+1)  # Surplus/Deficit of sediment mass in element.
        self.fraction_Q_jct = np.zeros(T)  # Fraction of flow present at upstream junction that enters reservoir.
        self.Q_jct = np.zeros(T)  # Upstream junction inflow, before any flow is distributed among outflow elements

    def Import_Data(self, T, Input_Data_File):
        '''
        Purpose: Imports all user-specified data relevant to the channel segment of interest.

        :param T: Simulation Duration (days)
        :param Input_Data_File: openpyxl MS Excel (.xlsx) object input file for current simulation run
        :return: Nothing (sets self attributes)
        '''
        if 'Reach Specifications' in Input_Data_File.sheetnames:
            [self.Routing_Coefficient, self.Routing_Exponent, self.Pool_Volume, self.Initial_Storage, self.sed_alpha,
             self.sed_beta, self.Initial_Sediment_Mass, channel_routing, self.density_SS,
             self.calibration_preference, self.L, self.W_bkfl, self.z_ch, self.depth_bkfl, self.ch_slope,
             self.mann_n, self.PRF, self.CSP, self.SP_EXP,
             self.flow_thold] = Excel_Data_Import(self.name, Input_Data_File,
                                                  'Reach Specifications', 1, 20,
                                                   max_distinct_data_types=None,
                                                   data_name_offset=None)

            # Read in a flow threshold, above/below which to track performance of particular metrics.
            if self.flow_thold is None:
                self.flow_thold = 0  # Track passage for all inflows

            # Energy production above/below threshold
            self.Q_out_over_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
            self.Q_out_over_thold[:] = 'NaN'  # Fill array with 'NaN'
            self.Q_out_under_thold = np.zeros(T)  # Binary (yes=1/no=0), Qin < provided flow threshold.
            self.Q_out_under_thold[:] = 'NaN'  # Fill array with 'NaN'

            # Channel routing
            if channel_routing == "Null Routing (Flow in = Flow Out)":
                self.Channel_Routing_Method = 0
            elif channel_routing == "Storage-Outflow Routing":
                self.Channel_Routing_Method = 1
            elif channel_routing == 'Manning Routing':
                self.Channel_Routing_Method = 2
                self.W_btm = max(0,self.W_bkfl - 2*self.z_ch*self.depth_bkfl)  # Channel bottom width
            else:
                self.Channel_Routing_Method = 0  # Set default = "Null Routing (Flow in = Flow Out)"

            if self.density_SS > 0:
                pass
            else:
                self.density_SS = 1000  # default (if no value provided by user): 1000 kg/m^3 (density of water).
        # Reach sediment carrying capacity calibration preferences. Reset string value to a number.
        if self.calibration_preference == "Calibrate coefficient":
            self.calibration_preference = 1  # Perform Calibration in separate method.
        elif self.calibration_preference == "Specify coefficient":
            self.calibration_preference = 2  # Do nothing else; alpha and beta values have already been imported.
            if self.sed_alpha >= 0 and self.sed_beta >= 0:
                pass
            else:
                print "Error: Sediment carrying capacity parameters not properly specified for element %s in Reach " \
                      "Specifications worksheet" % self.name
        elif self.calibration_preference == "Sediment mass out (kg) = Sediment mass in (kg)":
            self.calibration_preference = 3  # Sediment production parameters are not needed
            self.sed_alpha = 0  # Set alpha = beta = 0.
            self.sed_beta = 0  # Set alpha = beta = 0.
        elif self.calibration_preference == "Specify Bagnold Equation Parameters":
            self.calibration_preference = 4  # Sediment production parameters are not needed
            self.sed_alpha = 0  # Set alpha = beta = 0.
            self.sed_beta = 0  # Set alpha = beta = 0.
        else:
            # User did not specify calibration preferences
            if self.sed_alpha > 0 and self.sed_beta > 0:
                self.calibration_preference = 2  # User specified parameter values, so use those.
            else:
                self.calibration_preference = 3  # No values specified. Default: "Sediment mass out (kg) = Sediment
                #  mass in (kg)"
                self.sed_alpha = 0
                self.sed_beta = 0

    def Calibration_Reach_Sediment_Carrying_Capacity(self):
        '''
        Purpose: Performs calibration of rating curve parameters for function defining channel segment sediment
        carrying capacity.

        This call cannot occur as part of __init__, as annual sediment loads are not available yet. So this method is
        called from the main pysedsim simulation file (pysedsim_main_simulation.py)

        :return: Nothing (sets self attributes)
        '''

        if self.calibration_preference == 1:
            if self.Cum_Annual_SED_LOAD > 0 and self.sed_beta > 0:
                self.sed_alpha = self.Cum_Annual_SED_LOAD / (
                sum(86400 * np.power(self.Q_out_unreg_for_param_calib, (self.sed_beta + 1))) / (
                self.T / 365))  # Perform calibration
            else:
                print "Error: Two errors are possible (1) Sediment calibration requested for channel element %s, " \
                      "but no cumulative sediment load exists for this channel segment, per Sediment Loads worksheet. " \
                      "Or (2): No exponent value is provided for the carrying capacity calculation in the Reach " \
                      "specifications worksheet." % self.name
                self.sed_alpha = 0  # No annual sediment load specified. Set alpha = beta = 0.
                self.sed_beta = 0  # No annual sediment load specified. Set alpha = beta = 0.

    def Flow_Routing(self, t):
        '''
        Purpose: Contains channel flow routing routines. User preferences in input file dictate which method(s) are
        used.

        :param t: Current time period (integer, value range: [0,T])
        :return: Nothing (sets self attributes)
        '''
        # Use user-specified routing technique to determine channel outflow
        if self.Channel_Routing_Method == 0:
            # Implement Null Channel Routing (Q_in = Q_out)
            self.Q_out[t] = self.Q_in[t]
            self.V_out[t] = self.V_in[t]
        elif self.Channel_Routing_Method == 1:
            # Determine Channel outflow as a function of reach storage
            if max(self.S[t] + self.V_in[t] - self.Pool_Volume, 0) >= 1:
                routing_exponent = self.Routing_Exponent
            else:
                routing_exponent = 1 # To avoid raising a fraction to a fraction
            # Storage used to compute outflow cannot be less than zero, and outflow cannot reduce storage[t+1] to
            # below the pool volume.
            self.Q_out[t] = min((self.S[t] + self.V_in[t]) / 86400, self.Routing_Coefficient * (
            max(self.S[t] + self.V_in[t] - self.Pool_Volume, 0)) ** routing_exponent)
            self.V_out[t] = self.Q_out[t] * self.dt * 86400
        elif self.Channel_Routing_Method == 2:
            # Implement Manning Routing methods
            # Methods are borrowed from SWAT 2009 documentation

            self.A_cx = self.S[t]/self.L # Compute cross-sectional flow area at beginning of time step t
            # Water depth at time t:
            self.depth = math.sqrt((self.A_cx/self.z_ch) + np.power((self.W_btm/(2*self.z_ch)),2)) - self.W_btm/(
                2*self.z_ch)
            self.wetted_perim = self.W_btm + 2*self.depth*math.sqrt(1 + np.power(self.z_ch,2))
            self.hydraul_rad = self.A_cx/self.wetted_perim
            self.Q_out[t] = self.A_cx*(np.power(self.hydraul_rad, 2/3))*math.sqrt(self.ch_slope)/self.mann_n
        self.S[t+1] = self.S[t] + self.V_in[t] - self.V_out[t]

        if self.Q_jct[t] >= self.flow_thold:
            self.Q_out_over_thold[t] = self.Q_out[t]
        else:
            self.Q_out_under_thold[t] = self.Q_out[t]

    def Sediment_Routing(self, t):
        '''
        Purpose: contains channel sediment routing routines. User preferences in input file dictate which method(s)
        are used.

        :param t: Current time period (integer, value range: [0,T])
        :return: Nothing (sets self attributes)
        '''
        # Begin sediment routing calculations
        if self.calibration_preference == 3:
            SS_W_carycap = self.SS_W_in[t]
            self.SS_W_out[t] = self.SS_W_in[t]
        elif self.calibration_preference == 4:
            self.q_ch_pk = (self.PRF)*self.Q_out[t]
            self.v_ch_pk = self.q_ch_pk/self.A_cx
            SS_W_carycap = self.CSP * (self.v_ch_pk ** (self.SP_EXP)) * (1/1000) * self.dt * 86400
        else:
            SS_W_carycap = self.sed_alpha * (self.Q_out[t] ** (self.sed_beta + 1)) * self.dt * 86400
        self.SS_W_out[t] = min(SS_W_carycap, self.BS_W[t] + self.SS_W[t] + self.SS_W_in[t])
        if (self.SS_W[t] + self.SS_W_in[t]) >= SS_W_carycap:
            # If more sediment exists in suspension than the carrying capacity, then some settling/deposition will
            # occur.
            if self.SS_W[t] + self.SS_W_in[t] - self.SS_W_out[t] - SS_W_carycap >= 0:
                # In case dredging placed any sediment here, begin with t+1 value here
                BS_W_settle = (self.SS_W[t] + self.SS_W_in[t]) - self.SS_W_out[t] - SS_W_carycap
                self.BS_W[t + 1] = self.BS_W[t + 1] + self.BS_W[t] + BS_W_settle
                self.SS_W[t+1] = SS_W_carycap
            else:
                BS_W_resuspend = min(self.BS_W[t], abs((self.SS_W[t] + self.SS_W_in[t]) - self.SS_W_out[t] - SS_W_carycap))
                if self.BS_W[t] > abs((self.SS_W[t] + self.SS_W_in[t]) - self.SS_W_out[t] - SS_W_carycap):
                    # In case dredging placed any sediment here, begin with t+1 value here
                    self.BS_W[t+1] = self.BS_W[t+1] + self.BS_W[t] - BS_W_resuspend
                    self.SS_W[t+1] = SS_W_carycap
                else:
                    self.BS_W[t+1] = 0 + self.BS_W[t+1]  # In case dredging placed any sediment here, keep t+1 here
                    self.SS_W[t+1] = (self.SS_W[t] + self.SS_W_in[t]) - SS_W_carycap + BS_W_resuspend
        else:
            # If less sediment exists in suspension than the carrying capacity, then some scour from available BS_W
            # will occur.
            if (SS_W_carycap <= (self.BS_W[t] + (self.SS_W[t] + self.SS_W_in[t]))):
                BS_W_resuspend = min(self.BS_W[t], abs((self.SS_W[t] + self.SS_W_in[t]) - self.SS_W_out[t] - SS_W_carycap))
                if abs((self.SS_W[t] + self.SS_W_in[t]) - self.SS_W_out[t] - SS_W_carycap) <= self.BS_W[t]:
                    self.SS_W[t+1] = SS_W_carycap
                    # In case dredging placed any sediment here, begin with t+1 value here
                    self.BS_W[t+1] = self.BS_W[t+1] + self.BS_W[t] - BS_W_resuspend
                else:
                    self.SS_W[t+1] = (self.SS_W[t] + self.SS_W_in[t]) - self.SS_W_out[t] + BS_W_resuspend
                    # In case dredging placed any sediment here, begin with t+1 value here
                    self.BS_W[t+1] = 0 + self.BS_W[t+1]
            else:
                self.SS_W[t+1] = 0
                BS_W_resuspend = self.BS_W[t]
                # In case dredging placed any sediment here, begin with t+1 value here
                self.BS_W[t+1] = 0 + self.BS_W[t+1]

    def Import_External_Element_State(self, t):
        '''
        Imports data from other elements being simulated (e.g., junctions, reaches and reservoirs) for use in the
        simulation of this element.

        Input data are stored in the elem_xfer_input_dict (which originally came from other elements'
        elem_xfer_output_dict). This module simply takes this information and converts it into a specific use for
        this element. For this reason, each element will have their own unique implementation of this method. This
        method is executed at the beginning of every time step, before any processes are simulated.
        '''

        # Loop through elements sending data, and unpack it to provide use for this element.
        for keys in self.elem_xfer_input_dict:
            if keys in self.Element_Sub_Dict['Inflow Elements']:
                for sub_keys in self.elem_xfer_input_dict[keys]:
                    if sub_keys == "Upstream Junction Inflow Fraction":
                        self.fraction_Q_jct[t] = self.elem_xfer_input_dict[keys][sub_keys]
                    elif sub_keys == "Upstream Junction Inflow":
                        self.Q_jct[t] = self.elem_xfer_input_dict[keys][sub_keys]


class Diversion_Channel(Channel):
    '''
    Purpose: Defines Channel sub-class Diversion_Channel. A diversion is a channel or pipe that is a regulated (at
    inflow) conveyance channel, unregulated at outflow. This occurs at a dam, and may be a closed conduit.

    Note: This class is under development and currently provides limited functionality.
    '''

    def __init__(self, name, T, Input_Data_File, Element_Sub_Dict, stochastic_components=None):
        '''
        Purpose: Constructor for Diversion_Channel

        Inputs all defined in Channel class definition.
        '''
        if hasattr(Channel, '__init__'):
            Channel.__init__(self, name, T, Input_Data_File, Element_Sub_Dict, stochastic_components)  # If parent
            # class has constructor method, then call that first.