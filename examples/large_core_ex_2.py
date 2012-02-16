# examples/coretools/large_core_ex_2.py

import large_core

from poropy.coretools import Optimizer
from pypgapack import PGA
from mpi4py import MPI
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

# Optimizer
# =========
#
# While Optimizer has some default definitions,
# we make an inherited class to set our own
# objective function.  Moreover, because we are
# fixing the central bundle to be twice burned,
# we need special initialization and mutation
# functions.

class OptimizeLargeCore(Optimizer) :
    """  Derive our own class from PGA.
    """
    def large_core_objective(self, p, pop) :
        """ Minimize peaking and maximize keff using weighted objective.
        
        Here, we seek to maximize keff for p < 1.50
        """
        pattern = self.GetIntegerChromosome(p, pop)
        self.reactor.shuffle(pattern) 
        keff, peak = self.reactor.evaluate()
        val = self.fun(keff,peak)
        self.evals += 1
        #print val, keff, peak, pattern
        del pattern
        return val
    
    def fun(self,k,p):
        """  Delta is *negative* when p exceeds the maximum.  

        A k less than 1.13 also produces a negative value, but by a 
        much smaller amount. 
        """
        delta = 0
        if p > 1.50 :
            delta = 1.50 - p
        return (k - 1.13) + 50.0 * delta
        

    def swap(self, p, pop, pm) :
        """  Random swap of bundles.
        
        This example allows swapping for all but the central
        element.  Moreover, no checking is done to ensure swaps
        are not done between identical or forbidden bundles.  In
        those cases, a swap simply doesn't happen.  This is one
        area where a lot more work can be done w/r to implementation.
        """
        #pass
        n = self.GetStringLength()
        pattern = self.GetIntegerChromosome(p, pop)
        # index is a random permutation of [1,2,... numberbundles-1]
        #   i.e. zero is excluded (and zero is the central bundle)
        index = np.random.permutation(n-1) + 1 
        i1 = index[0]
        i2 = index[1]
        tmp1 = pattern[i1]
        tmp2 = pattern[i2]
        pattern[i1] = tmp2  
        pattern[i2] = tmp1
        del pattern 
        return 1 # A positive value means something swapped.
            
    def init(self, p, pop) :
        """  Random initial states.  
        """
        n = self.GetStringLength()
        pattern = self.GetIntegerChromosome(p, pop)
        # perm is a random permutation of [0,1,... numberbundles-2]
        #   i.e. 48 is excluded.  The 49th bundle is by definition
        #   the least reactive, and so that be the default central.
        perm = np.random.permutation(n-1)
        pattern[0] = 48
        for i in range(1, n) :
            pattern[i] = perm[i-1]
        #print " pattern=", pattern   
        del pattern

# Begin Optimization
# ==================
    
# MPI things.  These should be commented out if 
# a serial version of pypgapack is used.
comm = MPI.COMM_WORLD   # Get the communicator.
rank = comm.Get_rank()  # Get my rank.  
size = comm.Get_size()  # Get number of processes
time = MPI.Wtime()      # Start the clock.  

# Get the large reactor.  Pass rank for Laban initialization.
reactor = large_core.make_large_core(rank)


number_runs = 1
vals = np.zeros(number_runs)
kefs = np.zeros(number_runs)
peks = np.zeros(number_runs)
evas = np.zeros(number_runs)

for i in range(0, number_runs) :
    # Create an optimizer.
    opt = OptimizeLargeCore(sys.argv, reactor)
    
    # Set initialization and swapping functions.
    opt.SetInitString(opt.init)
    opt.SetMutation(opt.swap)
    opt.SetNoDuplicatesFlag(PGA.TRUE) # Keep no duplicate patterns.
    # Set various PGA parameters 2
    opt.SetNumReplaceValue(190)
    opt.SetRandomSeed(i+11234)    # Set random seed for verification.  
    np.random.seed(i+2391)        # Do the same with Numpy.
    opt.SetPopSize(200)            # Large enough to see some success.
    num_iter = 200
    opt.SetMaxGAIterValue(num_iter)  # Small number for output.

    # Set various Optimizer parameters
    opt.set_track_best(True)     # Track best everything each generation
    opt.set_fixed_central(True)  # Fix the central bundle.
    
    # Run the optimization.  This implicitly performs both SetUp 
    #   and Run of PGA. 
    opt.run(opt.large_core_objective)

    if rank == 0 :
        best = opt.GetBestIndex(PGA.OLDPOP)     # Get the best string
        bestpattern = opt.GetIntegerChromosome(best, PGA.OLDPOP)
        print " best pattern ", bestpattern
        reactor.shuffle(bestpattern)  
        reactor.evaluate()
        reactor.print_params()
        reactor.print_pattern()
        reactor.print_peaking()
        vals[i]=opt.fun(reactor.evaluator.keff, reactor.evaluator.maxpeak)
        kefs[i]=reactor.evaluator.keff
        peks[i]=reactor.evaluator.maxpeak
        evas[i]=opt.evals
        #reactor.plot_pattern('burnup')
if rank > 0:
    evals = np.array([opt.evals],dtype='i')
    comm.Send([evals,MPI.INT], dest=0, tag=13)
else :
    evals = np.array([1],dtype='i')
    for i in range(1, size) :
        comm.Recv([evals,MPI.INT], source=i, tag=13)
        opt.evals += evals[0]

if rank == 0 :
    print np.mean(vals), np.std(vals)
    print np.mean(kefs), np.std(kefs)
    print np.mean(peks), np.std(peks)
    print np.mean(evas), np.std(evas)
    print " ELAPSED TIME: ", MPI.Wtime() - time, " EVALS = ", opt.evals
    #reactor.plot_pattern('burnup')
    #print opt.best_eval
    plt.plot( np.arange(0,num_iter+1), opt.best_eval, 'b',\
              lw=2) # Plot the objective as a function of generations
                    #   against the reference solution.
    plt.title('Convergence of Objective')
    plt.xlabel(' generation')
    plt.ylabel(' objective ')
    plt.grid(True)
    plt.show()

opt.Destroy()  # Clean up PGAPack internals.
