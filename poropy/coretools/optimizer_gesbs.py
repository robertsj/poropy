'''
Created on Feb 17, 2012

@author: robertsj
'''

from optimizer import Optimizer
import numpy as np

class OptimizerGESBS(Optimizer) :
    """  Optimizes using greedy exhaustive single binary sweeps.

    The basic algorithm is as follows.  For each bundle, swap with
    every other bundle.  Keep track of the single swap that improves
    the objective most.  For a quarter core of 50 assemblies, this
    means do 2500 swaps and record the best.  Then, do the swap, 
    and perform 2500 new swaps, keeping the best.  Do this until
    converged.

    As an alternative, a flag can be set to make the method
    greedier: in this case, as soon as a swap improves the objective,
    it is kept.
    """
    
    # Public Interface
    
    def __init__(self, argv, reactor, objective):
        """ Constructor.  
        """
        # Initialize the Optimizer object.
        Optimizer.__init__(self, argv, reactor, objective)

        
    def execute(self):
        """  Optimize the reactor for the objective.
        """     
        raise NotImplementedError 
    
    def end_of_iteration(self):
        """  Do something at the end of each iteration.
        """
        raise NotImplementedError 
