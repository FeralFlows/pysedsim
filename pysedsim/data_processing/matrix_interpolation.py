# -*- coding: utf-8 -*-

'''

Module to define the Matrix_Interpolation() function for interpolating with two lists of data.

The module allows the user to specify two lists of values (each representing different variables), as well as a known
value of one of the two variables. The function interpolates to determine the value of the second variable given the
known value of the first variable.
'''

# This file contains a function that can interpolate within elevation-volume-area or flow-elevation tables to return values of interest.

# Import relevant libraries:
from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer. (e.g., 1/4 = 0.25, not 0).
import logging

def Matrix_Interpolation(name, input_matrix, provided, locate, input_value):
    # This is a generalized function that accepts a matrix and performs an interpolation in the matrix given a value in one column and a desired value in the next column over.
    # name: used only to produce error messages for the user, indicating for what reservoir (and in what worksheet) a possible error has occurred.
    # input_matrix: this is a matrix/table with at least two columns and at least two rows, within which the user
    # wishes to interpolate.
    # provided = value type provided by user: "elevation", "storage", "area", or "flow"
    # locate = value type user wishes to have returned via interpolation: "elevation", "storage", "area", or "flow"
    # value = the value (of type = provided) specified by the user
    # elevation: Enter an elevation value, for which a corresponding water storage value will be interpolated from the column directly to the right in input_matrix. . elevation = None should be entered if category does not apply.
    # water_storage: Enter a water_storage value, for which a corresponding elevation value will be interpolated from the column directly to the left in input_matrix. water_storage = None should be entered if category does not apply.
    # flow: Enter a flow value, for which a corresponding elevation value will be interpolated from the column directly to the left in input_matrix. flow = None should be entered if category does not apply.
    # This function returns a reservoir water storage from the Elevation-Volume curve as a function of a specified elevation

    # Example of typical function calls for this function:
    # 1. Matrix_Interpolation(self.name, self.E_V_A, elevation, storage, 20.5)
    # 2. self.STORAGE_TARGET[t+1] = Reservoir.Matrix_Interpolation(self.name, self.E_V_A, elevation, storage, self.ELEVATION_TARGET[t+1])
    # 3. self.ELEVATION_TARGET[t+1] = Reservoir.Matrix_Interpolation(self.name, self.E_V_A, storage, elevation, self.STORAGE_TARGET[t+1])
    # 4. self.ELEVATION_TARGET[t+1] = Reservoir.Matrix_Interpolation(self.name, self.Tailwater_Rating_Curve, flow, elevation, self.Incremental_Flow_Junction[t])

    # Depending on table type, set from which column the user-provided value comes, and from which column an interpolated value will be retrieved.
    if provided == "elevation":
        # Elevation provided.
        column_start = 0
        if locate == "storage":
            column_return = 1
        elif locate == "area":
            column_return = 2
        elif locate == "flow":
            column_return = 1  # Used only if elevation vs. outlet capacity is provided
        else:
            column_return = 1  # default for the function is to return storage if provided elevation.
        input_label = 'Elevation'  # For providing error messages to user
    elif provided == "storage":
        # Water storage provided.
        column_start = 1
        if locate == "elevation":
            column_return = 0
        elif locate == "area":
            column_return = 2
        else:
            column_return = 0  # default for the function is to return elevation if provided storage.
        input_label = 'Storage'
    elif provided == "area":
        # Reservoir water inundated area provided
        column_start = 2
        if locate == "elevation":
            column_return = 0
        elif locate == "storage":
            column_return = 1
        else:
            column_return = 0  # default for the function is to return elevation if provided area.
        input_label = 'Surface Area'
    elif provided == "flow":
        # Used for interpolating in: (1) tailwater table (flow-elevation table), (2) flow-fraction tables for junctions
        column_start = 0
        if locate == "elevation":
            column_return = 1
        elif locate == "fraction":
            column_return = 1
        else:
            column_return = 1  # default for the function is to return elevation if provided area.
        input_label = 'Flow'
    else:
        return

    # Main interpolation loop: loop through columns until provided value exceeds a value in the start column, and interpolate a corresponding value from the return column.
    for j in range(len(input_matrix[0])):
        if input_value >= 0:
            track_highest_value = input_matrix[column_return][j] #Store highest current value, in case we need to set interpolated_value = this value at the end.
        else:
            Interpolated_Value = input_matrix[column_return][j]
            logging.critical("for element {0}, provided {1} value {2} negative, which is physically "
                             "infeasible.".format(name, input_label, input_value))
            break
        if input_value > input_matrix[column_start][j]:
            # Skip to next row in E-V-A table, unless we are on last row of provided matrix, in which case store highest possible value and return that, along with an error message to user.
            if j == len(input_matrix[0])-1:
                Interpolated_Value = track_highest_value
#                print "Error: for element %s, provided %s value %s is higher than highest provided value in the following table %s" % (name, input_label, input_value, input_matrix)
                break
        elif input_value == input_matrix[column_start][j]:
            Interpolated_Value = input_matrix[column_return][j]
            break
        else:
            if j == 0:
                Interpolated_Value = input_matrix[column_return][j]
#                print "Error: for element %s, provided %s value %s is lower than lowest provided value in the following table %s" % (name, input_label, input_value, input_matrix)
                break
            else:
                Interpolated_Value = input_matrix[column_return][j-1] + (input_value - input_matrix[column_start][j-1]) * (input_matrix[column_return][j] - input_matrix[column_return][j-1]) / (input_matrix[column_start][j] - input_matrix[column_start][j-1])
                break
    return Interpolated_Value

def Return_Closest_Values(input_list, input_value):
    '''
    This method takes a list, sorts it, and finds the closest two list values to a provided value.

    Args:
        input_list: A 1-D list of unsorted values, e.g. [1,5,2,8].

    Returns:
        closest_value: a list of either 1 element (if an exact match in list exists) or 2 elements (the next smallest and largest elements)
    '''

    sorted_list = sorted(input_list)
    closest_value = []
    for i in range(len(sorted_list)):
        if input_value > sorted_list[i]:
            if i == len(sorted_list)-1:
                closest_value.append(sorted_list[i])
                break
        elif input_value == sorted_list[i]:
            closest_value.append(sorted_list[i])
            break
        else:
            if i == 0:
                closest_value.append(sorted_list[i])
                break
            else:
                closest_value.append(sorted_list[i])
                closest_value.append(sorted_list[i-1])
                break
    return closest_value