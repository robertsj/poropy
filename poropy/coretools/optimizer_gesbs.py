'''
Created on Feb 17, 2012

@author: robertsj
'''

from optimizer import Optimizer
import numpy as np

class OptimizerGESBS(Optimizer) :
    """  Optimizes using greedy exhaustive single binary sweeps.
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
