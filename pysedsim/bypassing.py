# -*- coding: utf-8 -*-

'''
Module to define the Bypass Structure class and asssociated methods.

Module contains the reservoir sediment bypass structure class definition. A bypass structure serves to physically
divert some water/sediment into a bypass, while the remaining water/sediment enters the reservoir (e.g.,
comprised of weirs and/or check dams). It should have an outflow at the junction immediately upstream of the
reservoir.
'''

# Import relevant modules
from __future__ import division
import numpy as np
import data_processing
from storage_element import Storage_Element

class Bypass_Structure(Storage_Element):
    """An instance is a Bypass Structure."""
    def __init__(self, name, T, Input_Data_File, Element_Sub_Dict):
        self.Sed_Mgmt_Type = 'Bypassing'
        if hasattr(Storage_Element, '__init__'):
            Storage_Element.__init__(self, name, T, Element_Sub_Dict)  # Call parent constructor.
        Bypass_Structure.Array_Initialization(self, T)
        Bypass_Structure.Import_Data(self, T, Input_Data_File)

    def Master_Method_Caller(self, t, Flow_in, Sed_in):
        # This method calls all other relevant methods in this class in the correct order.
        Storage_Element.Master_Method_Caller(self, t)  # Increments date
        Storage_Element.Element_Inflows(self, t, Flow_in, Sed_in)
        Bypass_Structure.Mass_Balance(self, t)
        Bypass_Structure.Bypass_Flow_Allocation(self, t)
        Storage_Element.Outflow_Element_Allocation(self, t)  # This is a bypass, so flow will simply be allocated to single downstream junction.

    def Array_Initialization(self, T):
        self.SS_W_Bypass = np.zeros(T)  # Daily sediment flux (kg/day) to enter bypass channel. Must be sent to junction downstream of bypass structure first.
        self.Q_bypass = np.zeros(T) # Daily water flow rate (kg/day) to enter bypass channel. Must be sent to junction downstream of bypass structure first.
        self.Res_Q_in = np.zeros(T)  # Daily water flow rate (kg/day) to enter reservoir downstream. Must be sent to junction downstream of bypass structure first.
        self.Res_SS_W_in = np.zeros(T)
        self.Res_V_in = np.zeros(T)
        self.Q_bypass_target = np.zeros(T)  # Daily target flow to be maintained in bypass.
        self.Bypass_Flow_Success = np.zeros(T)  # Binary success variable. = 0 during t if did not meet Q_bypass_target[t], and 1 if did meet target.
        self.bypass_flow_resilience_tracker = np.zeros(T)  # Tracks days when a failure to achieve adequate fish pass inflow after the very same failure condition has occurred on a previous day.

    def Import_Data(self, T, Input_Data_File):
        if 'Bypassing' in Input_Data_File.sheetnames:
            [self.Bypass_Flow_Threshold, self.Bypass_Capacity, self.Bypass_Fraction, self.Min_Bypass_Flow, self.Min_Bypass_Flow_Fraction] = data_processing.Excel_Data_Import(self.name, Input_Data_File, 'Bypassing', 1, 5, max_distinct_data_types = None, data_name_offset = None)  # Optional worksheet
        else:
            self.error = 1
            print "Error: Sediment management of type %s does not have a corresponding and correctly named worksheet in the input file." % self.Sed_Mgmt_Type

    def Bypass_Flow_Allocation(self, t):
        # This method allocates the element outflow to each of the elements downstream of the bypass. Will be sent to downstream junction
        # to be allocated to whatever elements branch from there.
        self.Q_bypass[t] = min(max(self.Q_in[t] - self.Bypass_Flow_Threshold, self.Min_Bypass_Flow, self.Min_Bypass_Flow_Fraction * self.Q_in[t]), self.Q_in[t], self.Bypass_Capacity)
        self.Q_bypass_target[t] = max(self.Min_Bypass_Flow, self.Min_Bypass_Flow_Fraction * self.Q_in[t])  # Target is maximum of the two possible target types user may have defined.
        if self.Q_bypass[t] > 0:
            self.SS_W_Bypass[t] = (self.Q_bypass[t] / self.Q_in[t]) * self.SS_W_in[t] + (1 - self.Q_bypass[t] / self.Q_in[t]) * self.SS_W_in[t] * (self.Bypass_Fraction)
            self.Res_SS_W_in[t] = self.SS_W_in[t] - self.SS_W_Bypass[t]
            # Adjust sediment and water inflows now that mass is bypassing reservoir.
            self.Res_Q_in[t] = self.Q_in[t] - self.Q_bypass[t]
            self.Res_V_in[t] = self.Res_Q_in[t] * 86400 * self.dt
        Bypass_Structure.Bypass_Performance(self,t)  # Call performance evaluation method now that flows have been allocated to bypass.

    def Bypass_Performance(self, t):
        ''' This method defines time series variables that will be relevant to assessing  bypass performance.
        '''
        if self.Q_bypass[t] >= self.Q_bypass_target[t]:
            self.Bypass_Flow_Success[t] = 1  # Successfully met velocity target
        else:
            self.Bypass_Flow_Success[t] = 0  # Failed to meet velocity target
            if t > 0 and self.Bypass_Flow_Success[t-1] == 0:
                self.bypass_flow_resilience_tracker[t] += 1

    def Mass_Balance(self, t):
        self.Q_out[t] = self.Q_in[t]
        self.SS_W_out[t] = self.SS_W_in[t]