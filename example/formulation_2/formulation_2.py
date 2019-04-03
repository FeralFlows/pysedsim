# This module contains the pysedsim function calls required to execute Formulation 2 in Wild et al. (in review).
from pysedsim import PySedSim
import sys
import pysedsim.visualization.processing_reference_set

pys = PySedSim(file_name = 'formulation_2_serial.csv')
# Manually Set location of Borg library
#pys.execute_optimization(file_name = 'formulation_2_serial.csv', borg_path = 'Z:/twild/Documents/Publications/2018/Wild et al. (2018) - EMS/PySedSim/local_git_repo/pysedsim/optimization/Borg.dll')

[Borg_dict, DPS_dict] = pysedsim.visualization.processing_reference_set.Reference_Set(input_file_name = 'formulation_2_serial.csv')

# Create reference set file (.ref) from all random seeds (.set files).
#[Borg_dict, DPS_dict] = Reference_Set(input_file_name = 'formulation_2.csv')