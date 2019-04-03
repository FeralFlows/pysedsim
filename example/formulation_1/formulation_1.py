# This script executes the PySedSim runs to produce the results associated with 
# Formulation I in Wild et al. (in review)

#import pysedsim
from pysedsim import PySedSim

# instantiate model
pys = PySedSim(file_name = 'formulation_1.csv')

# run a combined serial execution of multiple deterministic and stochastic simulations
pys.execute_simulation()