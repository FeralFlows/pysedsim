# -*- coding: utf-8 -*-

'''

Module to define the Dredging class and asssociated methods.

'''

# This module contains the dredging class definition

# Import relevant modules
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0)
from pysedsim.sediment_management.sediment_res_ops import Sediment_Res_Ops
import numpy as np
from pysedsim.data_processing.data_processing import *
import logging

class Dredging(Sediment_Res_Ops):

    '''
    This class represents the dredging technique for removing sediment from a reservoir. Class instances (objects)
    become stored in reservoir instances as part of its operating policy dictionary. The primary function of a
    dredging object is to perform a daily check to determine whether a reservoir will have dredging performed on a
    given day, and adjust the reservoir's sediment mass balance if so.
    '''

    def __init__(self, element_name, T, Input_Data_File):
        self.Sed_Mgmt_Type = 'General Sediment Removal'
        if hasattr(Sediment_Res_Ops, '__init__'):
            Sediment_Res_Ops.__init__(self, element_name, T, Input_Data_File, self.Sed_Mgmt_Type)   # Call parent constructor.
        Dredging.Array_Initialization(self, T)
        Dredging.Import_Data(self, Input_Data_File, element_name)

    def Array_Initialization(self, T):
        self.Sediment_Load_Removed_Daily = np.zeros(T)

    def Import_Data(self, Input_Data_File, element_name):
        # See if sediment management activity is scheduled to occur on this date. If so, import data and set Dredging = 1.
        if 'General Sediment Removal' in Input_Data_File.sheetnames:
            self.Dredging_Data_Import = Excel_Data_Import(element_name, Input_Data_File, 'General Sediment Removal', 2, 5, max_distinct_data_types = None, data_name_offset = 2) # Optional worksheet
        else:
            self.error = 1
            logging.critical("Sediment management of type {0} does not have a corresponding and correctly named "
                             "worksheet in the input file.".format(self.Sed_Mgmt_Type))

        # User-specified dredging data will be stored in a dictionary called Specs
        self.Specs = {}  # Initialize Dictionary
        i = 0  # Initialize counter
        for date_item in self.Dredging_Data_Import[0]:
            self.Specs[date_item] = {}  # Sub-dictionary inside this dictionary will store all relevant specs.
            try:
                self.Specs[date_item]["Duration"] = self.Dredging_Data_Import[1][i]
            except IndexError:
                self.Specs[date_item]["Duration"] = self.Dredging_Data_Import[1][0]    # None specified, load user-specified value from first row. # User did not specify an appropriate value; Load the default value stored in the first row
            try:
                self.Specs[date_item]["Max Sediment Load Removed"] = 1000 * self.Dredging_Data_Import[2][i]  # Note conversion from metric tons to kg
            except IndexError:
                self.Specs[date_item]["Max Sediment Load Removed"] = 1000 * self.Dredging_Data_Import[2][0]    # None specified, load user-specified value from first row. # Note conversion from metric tons to kg
            self.Specs[date_item]["Max Sediment Load Removed Daily"] = self.Specs[date_item]["Max Sediment Load Removed"] / self.Specs[date_item]["Duration"]
            try:
                self.Specs[date_item]["Sediment Destination"] = self.Dredging_Data_Import[3][i]
            except IndexError:
                self.Specs[date_item]["Sediment Destination"] = self.Dredging_Data_Import[3][0]  # None specified, load user-specified value from first row.
            try:
                self.Specs[date_item]["Fraction removed from active storage"] = self.Dredging_Data_Import[4][0]
            except IndexError:
                 self.Specs[date_item]["Fraction removed from active storage"] = self.Dredging_Data_Import[4][i]  # None specified, load user-specified value from first row.
            i += 1  # Increment counter

    def daily_initiation_check(self, current_date, t, BS_W_input):
        # This module checks date during time step t, before other reservoirs calcs are done, to see if the particular sediment management type will be performed during time step t.
        # Import method inputs
        self.current_date = current_date
        try:
            # Imports all flushing data and makes flushing current, if there is a dictionary key with a name equal to today's date.
            self.Specs["Current Specs"] = self.Specs[self.current_date] # Create a new key and data in specs dictionary that sets data for current date to most current data.
            self.Current = 1
        except KeyError:
            pass  # Do nothing
        if self.Current == 1:
            BS_W = Dredging.SedMgmt_Main(self, t, BS_W_input)  # Call main mass balance routing as long as Self.Current = 1
        return BS_W

    def SedMgmt_Main(self, t, BS_W_input):
        # Main method for dredging mass balance. Removes sediment then adjusts mass balance.
        # This time period resulted in successful sediment removal; reduce number of remaining required sediment removal days
        BS_W = BS_W_input
        self.Specs["Current Specs"]["Duration"] -= 1
        self.Sediment_Load_Removed_Daily[t] = min(BS_W, self.Specs["Current Specs"]["Max Sediment Load Removed Daily"])
        BS_W -= self.Sediment_Load_Removed_Daily[t]  # This is simply updating the BS_W(t+1,n) that has already been updated in the sediment trapping section above.
        if self.Specs["Current Specs"]["Duration"] == 0:
            self.Current = 0
        return BS_W

    def mass_balance(self, constant, t=None):
        self.att = constant
        self.removed_load = 3
        return self.removed_load