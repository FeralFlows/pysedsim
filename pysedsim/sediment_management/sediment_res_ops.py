# -*- coding: utf-8 -*-

'''

This module defines the Sediment_Res_Ops class and methods.

All reservoir sediment management types (i.e., flushing, sluicing, density current venting, bypassing and dredging),
which are to be contained in a reservoir's operating policy object, are sub-classes of this class, so they all share
the common behavior described here.

'''

from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0)
import logging
# from outlet import Outlet

class Sediment_Res_Ops:
    def __init__(self, element_name, T, Input_Data_File, Sediment_Mgmt_Type):
        self.Perform = 1 # If a class is a sub-type of this, then the user by definition is requesting that the procedure be considered in the simulation.
        Sediment_Res_Ops.Import_Data(self, element_name, T, Input_Data_File, Sediment_Mgmt_Type)

    def E_V_A_Adjustment(self):
        # Adjusts E-V-A curve. Can be done after trapped_load() function completes, or after various sediment management events have been completed.
        return

    def Import_Data(self, element_name, T, Input_Data_File, Sediment_Mgmt_Type):
        # See if flushing is scheduled to occur on this date. If so, import data and set Flushing = 1.
        self.Current = 0  # variable indicating whether management type is currently ongoing.
        self.error = 0  # Initialize sediment management data import error
        self.dt = 1  # time step = 1 day