# -*- coding: utf-8 -*-

'''
This module defines the storage element superclass and methods.

All PySedSim simulation elements (reaches, reservoirs, and junctions) are subclasses of Storage Element. They all
share the characteristics that they store water and sediment mass and require methods to manage/route that mass in
each time step.
'''

# Import relevant libraries
from __future__ import division  # Ensures result of quotient of two integers will be a float, not an integer.
import numpy as np
from datetime import timedelta  # Used to add days/months/years to a datetime object

class Storage_Element:
    SPD = 86400  # SPD = Seconds Per Day

    def __init__(self, name, T, Element_Sub_Dict):
        '''
        Purpose: Constructor for Storage_Element class.

        This serves as a super class for all PySedSim simulation elements, including reservoirs, river channel
        segments, junctions, and diversion structures.

        :param name: User-assigned element name (string)
        :param T: Simulation Duration (days)
        :param Element_Sub_Dict: Dictionary for element that contains basic data (e.g., name, type, connectivity)
        :return: Nothing (sets self attributes)

        '''
        # Note: The order in which method calls appear below is important. For this reason, it is not advisable to
        # re-arrange the order in which these calls occur below.
        self.name = name
        self.T = T
        self.Element_Sub_Dict = Element_Sub_Dict
        self.current_date = self.Element_Sub_Dict["Simulation Start Date"] - timedelta(days = 1)
        self.Element_Sub_Dict["Daily Water Outflows"] = {}  # Key: holds name of outflow element. Item: outflow value.
        self.Element_Sub_Dict["Daily Sediment Outflows"] = {}  # Key: holds name of outflow element. Item: outflow value.
        self.Cum_Annual_SED_LOAD = 0  # Initialize to zero, but reset values later
        for i in self.Element_Sub_Dict["Outflow Elements"]:
            self.Element_Sub_Dict["Daily Water Outflows"][i] = 0  # Store names of outflow elements here in order
            self.Element_Sub_Dict["Daily Sediment Outflows"][i] = 0  # Store names of outflow elements here in order
        # Initialize as arrays critical reservoir attributes.
        Storage_Element.Array_Initialization_SE(self, T)

    def __repr__(self):
        '''
        Purpose: Added method that prints object information

        :return:
        '''

        return '[Element Info.: %s]' % (self.name)  # String to print

    def Master_Method_Caller(self, t):
        '''

        :param t: Current time period (integer, value range: [0,T])
        :return: Nothing (sets self attributes)
        '''
        # This method calls all other relevant methods in this class in the correct order.
        self.current_date += timedelta(days = 1)  # Add 1 day to current date

    def Array_Initialization_SE(self, T):
        '''
        Purpose: Initializes variable arrays, called once (as part of class constructor).

        :param T: Simulation Duration (days)
        :return: Nothing (sets self attributes)
        '''
        # In comments describing units below, PA=Period average (e.g., daily average), where EOP indicates end-of-period
        # value (e.g., end-of-day value)
        self.dt = 1  # time step = 1 day
        self.Q_in = np.zeros(T)  # Water inflow rate (PA, m^3/s)
        self.Q_out = np.zeros(T)  # Water outflow rate (PA, m^3/s)
        self.Q_out_unreg_for_param_calib = np.zeros(T)  # Cumulative unregulated outflow rate (PA, m^3/s)
        self.V_in = np.zeros(T+1)  # Water volume inflow (EOP, m^3)
        self.V_out = np.zeros(T+1)  # Water volume outflow (EOP, m^3)
        self.S = np.zeros(T+1)  # Water storage (EOP, m^3)
        self.BS_W = np.zeros(T+1)  # Settled (deposited) sediment mass (EOP, kg)
        self.SS_W_in = np.zeros(T)  # Suspended sediment load inflow during time period (PA, kg)
        self.SS_W_out = np.zeros(T)  # Suspended sediment load outflow during time period (PA, kg)
        self.SS_V_out = np.zeros(T)  # Suspended sediment volume inflow during time period (PA, m^3)
        self.SS_W_inc = np.zeros(T)  # Suspended sediment load incremental inflow during time period (PA, kg)
        self.BS_V = np.zeros(T+1)  # Settled (deposited) sediment volume (EOP, m^3)
        self.SS_W = np.zeros(T+1)  # Suspended sediment mass in element (EOP, kg)
        self.SS_V = np.zeros(T+1)  # Suspended sediment volume in element (EOP, m^3)
        self.SS_C = np.zeros(T)  # Suspended sediment concentration (EOP, kg/m^3)
        self.TS_W = np.zeros(T+1)  # Total sediment mass = SS_W + BS_W (EOP, kg)
        self.TS_V = np.zeros(T+1)  # total sediment volume (EOP, m^3)
        self.TS_surplus_deficit = np.zeros(T+1)  # Surplus/Deficit of sediment mass (kg) in element.
        self.Qbypass = 0  # bypass flow rate (PA, m^3/s); only used if bypass structure exists upstream of junction.
        self.SSWbypass = 0  # bypassed sediment mass (PA, kg); only used if bypass structure exists upstream of junction.
        self.QReservoir = 0  # reservoir inflow rate (PA, m^3/s); used if bypass structure exists upstream of junction.
        self.SedReservoir = 0  # Sediment load inflow (PA, kg) to reservoir if bypass structure upstream
        self.VinReservoir = 0  # Water volume inflow (PA, m^3) to reservoir if bypass structure upstream
        self.elem_xfer_input_dict = {}  # Stores data (at time t) sent by other system objects being simulated.
        self.elem_xfer_output_dict = {}  # Stores data (at time t) to send to other system objects being simulated.

    def Time_Zero_Initialization(self):
        '''
        Purpose: Sets current simulation date, and initializes time zero values for certain arrays in first time step
        if storage element type is either channel or reservoir.

        :return: Nothing (sets self attributes)
        '''
        self.current_date = self.Element_Sub_Dict["Simulation Start Date"] - timedelta(days = 1)

        # Initialize time zero values for certain arrays if element type is reservoir or channel.
        if self.__class__.__name__ == 'Reservoir' or self.__class__.__name__ == 'Channel':
            self.S[0] = self.Initial_Storage
            self.BS_W[0] = self.Initial_Sediment_Mass
            self.BS_V[0] = self.BS_W[0] / self.density_SS
            self.SS_W[0] = 0
            self.SS_V[0] = self.SS_W[0]/self.density_SS
            self.TS_W[0] = self.SS_W[0] + self.BS_W[0]
            self.TS_V[0] = self.SS_V[0] + self.BS_V[0]

    def Post_Mass_Balance_Calcs(self, t):
        '''
        Purpose: Perform basic mass balance calculations at end of each time step that are common to all storage
        elements, as well as some class-specific mass balance calculations.

        :param t: Current time period (integer, value range: [0,T])
        :return: Nothing (sets self attributes)
        '''
        self.BS_V[t+1] = self.BS_W[t+1] / self.density_SS
        self.SS_V[t+1] = self.SS_W[t+1] / self.density_SS
        self.TS_W[t+1] = self.SS_W[t+1] + self.BS_W[t+1]
        self.TS_V[t+1] = self.TS_W[t+1] / self.density_SS
        self.TS_surplus_deficit[t+1] = self.TS_W[t+1] - self.TS_W[0]
        if self.__class__.__name__ == 'Reservoir':
            self.capacity_total_reservoir[t+1] = self.capacity_active_reservoir[t+1] + self.capacity_dead_reservoir[t+1]
        # self.SS_V_out[t] = self.SS_W_out[t]*self.SPD*self.dt

    def Outflow_Element_Allocation(self, t):
        '''

        Purpose: This method allocates the element outflow to each of the downstream elements. Note that the Junction()
        class has its own Outflow_Element_Allocation() method to handle junction-specific circumstances.

        :param t: Current time period (integer, value range: [0,T])
        :return: Nothing (sets self attributes)
        '''

        if len(self.Element_Sub_Dict["Daily Water Outflows"]) == 1:
            # Allocate all outflow to next element's inflow. This is a reach or reservoir, as outflow goes into only
            # one downstream element (likely a junction, which can then have multiple outflow points). Method does not
            # execute for last element(s) in system that have no outflow elements stored in ["Daily Water Outflows"]
            # list.
            self.Element_Sub_Dict["Daily Water Outflows"][self.Element_Sub_Dict["Daily Water Outflows"].keys()[0]] = \
            self.Q_out[t]
            self.Element_Sub_Dict["Daily Sediment Outflows"][
                self.Element_Sub_Dict["Daily Sediment Outflows"].keys()[0]] = self.SS_W_out[t]

    def Element_Inflows(self, t, Flow_in, Sed_in):
        '''
        Purpose: Set element inflow array values for time period t.

        :param t: Current time period (integer, value range: [0,T])
        :param Flow_in: Inflow rate (m^3/s) during time period t
        :param Sed_in: Sediment inflow (kg) during time period t
        :return: Nothing (sets self attributes)
        '''
        self.Q_in[t] = Flow_in
        self.V_in[t] = self.Q_in[t] * self.SPD * self.dt
        self.SS_W_in[t] = Sed_in