from __future__ import division

import os
import sys

import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from poropy.coretools import OptimizerGA
from pypgapack import PGA
from mpi4py import MPI

class PGAOpt(QThread):
    def __init__(self, model, parent):
        QThread.__init__(self,parent)
        # parent required to prevent early garbage collecting

        self.reactor = model.reactor


    def run(self):

        cwd = os.getcwd()
        os.chdir(os.path.join(sys.path[0],"tmpdir"))

        comm = MPI.COMM_WORLD   # Get the communicator.
        rank = comm.Get_rank()  # Get my rank. 
        time = MPI.Wtime()      # Start the clock.  

        number_runs = 1
        vals = np.zeros(number_runs)
        kefs = np.zeros(number_runs)
        peks = np.zeros(number_runs)
        evas = np.zeros(number_runs)

        for i in range(0, number_runs) :
            # Create an optimizer.
            opt = OptimizeSmallCore(sys.argv, self.reactor)
            
            # Set initialization and swapping functions.
            opt.SetInitString(opt.init)
            opt.SetMutation(opt.swap)
            
            # Set various PGA parameters 2
            opt.SetNumReplaceValue(290)
            opt.SetRandomSeed(i+1)       # Set random seed for verification.  
            np.random.seed(i+1)          # Do the same with Numpy.
            opt.SetPopSize(300)           # Large enough to see some success.
            opt.SetMaxGAIterValue(500)   # Small number for output.

            # Set various Optimizer parameters
            opt.set_track_best(True)     # Track best everything each generation
            opt.set_fixed_central(True)  # Fix the central bundle.
            
            # Run the optimization.  This implicitly performs both SetUp 
            #   and Run of PGA. 
            opt.run(opt.small_core_objective)

            if rank == 0 :
                best = opt.GetBestIndex(PGA.OLDPOP)     # Get the best string
                bestpattern = opt.GetIntegerChromosome(best, PGA.OLDPOP)
                print " best pattern ", bestpattern
                self.reactor.shuffle(bestpattern)
                self.reactor.evaluate()
                self.reactor.print_params()
                self.reactor.print_pattern()
                self.reactor.print_peaking()
                vals[i]=opt.fun(self.reactor.evaluator.keff, self.reactor.evaluator.maxpeak)
                kefs[i]=self.reactor.evaluator.keff
                peks[i]=self.reactor.evaluator.maxpeak
                evas[i]=opt.evals
                
        if rank == 0 :
            print np.mean(vals), np.std(vals)
            print np.mean(kefs), np.std(kefs)
            print np.mean(peks), np.std(peks)
            print np.mean(evas), np.std(evas)
            print " ELAPSED TIME: ", MPI.Wtime() - time

        opt.Destroy()  # Clean up PGAPack internals. 

        os.chdir(cwd)
        self.emit(SIGNAL("progressChanged(int)"),100)

class OptimizeSmallCore(OptimizerGA):
    """  Derive our own class from PGA.
    """
    def small_core_objective(self, p, pop) :
        """ Minimize peaking and maximize keff using weighted objective.
        
        Here, we seek to minimize p for k >= 1.135. Note that p=1.85 is
        about what we found manually, albeit with just one swap from the
        default pattern.
        """
        pattern = self.GetIntegerChromosome(p, pop)
        self.reactor.shuffle(pattern) 
        keff, peak = self.reactor.evaluate()
        val = self.fun(keff,peak)
        self.evals += 1
        del pattern
        return val
    
    def fun(self,k,p):
        delta = 0
        if k < 1.1 :
            delta = k - 1.1
        return 1.0 * (1.5 - p) + 50.0 * delta
        

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
        #   i.e. 17 is excluded.  The 17th bundle is by definition
        #   the least reactive, and so that be the default central.
        perm = np.random.permutation(n-1)
        pattern[0] = 17
        for i in range(1, n) :
            pattern[i] = perm[i-1]
        #print " pattern=", pattern   
        del pattern


