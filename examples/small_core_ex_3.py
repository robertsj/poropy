# examples/coretools/small_core_ex_3.py

import small_core

from poropy.coretools import OptimizerGA
from pypgapack import PGA
from mpi4py import MPI
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')


# Define the objective function.

def objective(self, k, p):
    """  User-defined objective
    """
    delta = 0
    if k < 1.1 :
        delta = k - 1.1
    return 1.0 * (1.5 - p) + 50.0 * delta
        
# Begin Optimization
# ==================

# Get the small reactor.  Pass rank for Laban initialization.
comm = MPI.COMM_WORLD   # Get the communicator.
rank = comm.Get_rank()  # Get my rank.  
reactor = small_core.make_small_core(rank)

# Create an optimizer.
opt = OptimizerGA(sys.argv, reactor, objective)

# Run it.
opt.execute()

best = opt.GetBestIndex(PGA.OLDPOP)     # Get the best string
bestpattern = opt.GetIntegerChromosome(best, PGA.OLDPOP)
print " best pattern ", bestpattern
reactor.shuffle(bestpattern)  
reactor.evaluate()
reactor.print_params()
reactor.print_pattern()
reactor.print_peaking()

opt.Destroy()  # Clean up PGAPack internals.
