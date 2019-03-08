# -*- coding: utf-8 -*-

'''

This module is designed to run a Monte Carlo simulation with PySedSim on a cluster.

This module is designed to be called from a shell script, where this script will be run by every processor in the
submitted job. The module determines what portion of the total number of simulations that will be run should be
conducted by an individual processor.

'''

# Imports
from mpi4py import MPI
from pysedsim import PySedSim

# Initialize MPI environment
comm = MPI.COMM_WORLD

# Get the number of processors and the rank of processors
n_procs = comm.size  # Number of processors (e.g., 512)
rank = comm.rank  # Rank of individual processor from 0 to (nprocs-1) (e.g., 0, 1, 2, ..., 511)

# Use the processor rank to determine the chunk of work each processor will do. If a remainder exists,
# some processors will be required to do extra simulations.
number_of_simulations = 500  # If user is simulating multiple scenarios, this number needs to be the same in all scenarios.
count = number_of_simulations/n_procs  # Number of simulations to be executed per processor. This must be an integer.
remainder = number_of_simulations % n_procs
if rank < remainder:
    start = rank*(count+1)
    stop = start+count+1
else:
    start = count*rank + remainder
    stop = start+count

# Call PySedSim model on each processor. Number of calls on processor is defined above.
PySedSim(start_stop = [start, stop], rank = rank)  # Run PySedSim