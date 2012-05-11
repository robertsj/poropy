# examples/large_core_ex_2.py
#
# This demonstrates use of a genetic algorithm for 
# optimization of the large core example.

from poropy.coretools import OptimizerGA

try :
  from pypgapack import PGA
except ImportError:
  print "You need pypgapack for this example."
  exit()

try :
  from mpi4py import MPI
  HAVEMPI = True
except ImportError:
  HAVEMPI = False

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
import time

# ===============================
# Define the objective function
# ===============================

def objective(k, p):
  """  The objective function.
       
  Delta is *negative* when p exceeds the maximum.  
  A k less than 1.13 also produces a negative value, but by a 
  much smaller amount. 
  """
  delta = 0.0
  if k < 1.1 :
    delta = k - 1.1
  return 1.0 * (1.5 - p) + 50.0 * delta
           
# ===============================
# Get the reactor model
# ===============================
import large_core
  
# ===============================
# Begin Optimization
# ===============================
    
# MPI things.  

comm = MPI.COMM_WORLD   # Get the communicator.
rank = comm.Get_rank()  # Get my rank.  
size = comm.Get_size()  # Get number of processes
#time = MPI.Wtime()      # Start the clock.  

# Get the large reactor.  Pass rank for Laban initialization.
reactor = large_core.make_large_core(rank)


number_runs = 1

vals = np.zeros(number_runs)
kefs = np.zeros(number_runs)
peks = np.zeros(number_runs)
evas = np.zeros(number_runs)

for i in range(0, number_runs) :

    # Create an optimizer.
    opt = OptimizerGA(sys.argv, reactor, objective)

    # Change some defaults.
    opt.maximum_generations    = 20 # Maximum generations
    opt.population_size        = 100000  # Population size
    opt.number_replaced        = 90   # Number replaced each generation
    opt.seed                   = 123  # PGA random number seed
    opt.np_seed                = 123  # NumPy random number seed

    t1 = time.time()

    # Run the optimization.  
    opt.execute()

    etime = time.time() - t1

    if rank == 0 :
      print " elapsed time = ", etime, " seconds."
    print " rank = ", rank, " evals = ", reactor.number_evaluations
    

#    if rank == 0 :
#        best = opt.GetBestIndex(PGA.OLDPOP)     # Get the best string
#        bestpattern = opt.GetIntegerChromosome(best, PGA.OLDPOP)
#        print " best pattern ", bestpattern
#        reactor.shuffle(bestpattern)  
#        reactor.evaluate()
#        reactor.print_params()
#        reactor.print_pattern()
#        reactor.print_peaking()
#        vals[i]=opt.fun(reactor.evaluator.keff, reactor.evaluator.maxpeak)
#        kefs[i]=reactor.evaluator.keff
#        peks[i]=reactor.evaluator.maxpeak
#        evas[i]=opt.evals

#if rank > 0:
#    evals = np.array([opt.evals],dtype='i')
#    comm.Send([evals,MPI.INT], dest=0, tag=13)
#else :
#    evals = np.array([1],dtype='i')
#    for i in range(1, size) :
#        comm.Recv([evals,MPI.INT], source=i, tag=13)
#        opt.evals += evals[0]

#if rank == 0 :
#    print np.mean(vals), np.std(vals)
#    print np.mean(kefs), np.std(kefs)
#    print np.mean(peks), np.std(peks)
#    print np.mean(evas), np.std(evas)
#    print " ELAPSED TIME: ", MPI.Wtime() - time, " EVALS = ", opt.evals
#    #reactor.plot_pattern('burnup')
#    #print opt.best_eval
#    plt.plot( np.arange(0,num_iter+1), opt.best_eval, 'b',\
#              lw=2) # Plot the objective as a function of generations
#                    #   against the reference solution.
#    plt.title('Convergence of Objective')
#    plt.xlabel(' generation')
#    plt.ylabel(' objective ')
#    plt.grid(True)
#    plt.show()

opt.Destroy()  # Clean up PGAPack internals.
