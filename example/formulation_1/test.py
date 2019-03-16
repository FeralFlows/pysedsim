#import pysedsim
from pysedsim import PySedSim

# instantiate model
pys = PySedSim(file_name = 'formulation_1.csv')

# run intended configuration
pys.execute()