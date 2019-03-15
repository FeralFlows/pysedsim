# This module contains the pysedsim function calls required to execute Case 3 in Wild et al. (in review).

from __future__ import division  # This ensures result of quotient of two integers will be a float, not an integer.
import sys
sys.path.append(r'\\essi12.umd.edu\documents\twild\Documents\PySedSim')

import processing_reference_set
import optimization_pysedsim

# Run optimization using the serial implementation of the Borg MOEA with the following function call:
optimization_pysedsim.Optimization_DPS_Borg(file_name = 'formulation_3_serial.csv')

# Alternatively, run optimization using the master-slave (parallelized) implementation of the Borg MOEA. This step must
#  be executed on a Linux cluster with the following function call:
optimization_pysedsim.Optimization_DPS_Borg(file_name = 'formulation_3_parallelized.csv')

# Create reference set file (.ref) from all random seeds (.set files).
[Borg_dict, DPS_dict] = processing_reference_set.Reference_Set(input_file_name = 'formulation_3.csv')

# Produce numerous 2D plots of four-objective tradeoffs
processing_reference_set.ref_set_plot_loop(file_name = 'formulation_3.csv')

# Find the energy maximizing policy in the reference set
# Based on the order of optimization objectives' listing in the "Optimization" worksheet of the input file,
# the objective order is as follows: 0) Annual Energy Production, 1) Dry season energy production, 2) Bypass flow,
# 3) Dry Season Larvae Flow Fraction, 4) Dry Season Larvae Flow Fraction, 5) Flushed load.

Objective_Range = {'Objective 0': {'Range': ['max', 'max'], 'Type': 'Max'},
                   'Objective 1': {'Range': ['max', 'max'], 'Type': 'Max'},
                   'Objective 2': {'Range': ['max', 'max'], 'Type': 'Max'},
                   'Objective 3': {'Range': ['max', 'max'], 'Type': 'Max'},
                   'Objective 4': {'Range': ['max', 'max'], 'Type': 'Max'},
                   'Objective 5': {'Range': ['max', 'max'], 'Type': 'Min'}}

formulation_3_ref_set_file_path = r'E:\Publications in prep\Paper ' \
                           r'2\Results\Output_Storage\Sambor_SA_formulation_3\sets\pysedsim_ref_set.ref'
formulation_3_objs_file_path = r'E:\Publications in prep\Paper ' \
                           r'2\Results\Output_Storage\Sambor_SA_formulation_3\sets\pysedsim_ref_set_no_vars.ref'

[rowlist, refset_objs, refset] = processing_reference_set.brush_pareto_policies(Objective_Range, objective_union='Yes',
                                                                                refset_file_path=formulation_3_ref_set_file_path,
                                                                                objs_file_path=formulation_3_objs_file_path,
                                                                                one_best_worst='Yes')

# Find a policy in the reference set that offers a "compromise"compromise policy that balanced


# NEW STUFF







# Probability plots

# Stochastically reevaluate policies identified through brushing above.
# 1. Entered values from policy list into "Reevaluation" sheet of input file, and specified preferences for variables
#  to be reevaluated and stored.
processing_reference_set.Policy_Reevaluation(file_name = 'formulation_3.csv')

# Deterministically reevaluate the policies over the historical record, using a new set of input files.
processing_reference_set.Policy_Reevaluation(file_name = 'formulation_3_deterministic_reeval.csv')