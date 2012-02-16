'''
Created on Feb 14, 2012

@author: robertsj
'''
from evaluator import Evaluator
from pyflare import *
import numpy as np


class Flare(Evaluator):
    """ Uses the FLARE model to evaluate loading patterns.
    """
    
    # Public Interface

    def __init__(self) :
        """ Constructor
        """
        self.keff = 1.0
        self.maxpeak = 1.0
        
    def setup(self, core) :     
        """  Set up some structures needed before solving.
        """
        # We need direct access to the core
        self.core = core
        # Validate a square quarter core. (Not applicable to 1/2 or 1/8)
        assert(len(self.core.stencil[0,:])==len(self.core.stencil[:,0]))
        
        # Set the geometry
        geometry.stencil = self.core.stencil
        geometry.number_bundles = len(self.core.pattern)
        geometry.stencil_dimension = len(self.core.stencil[0,:])
        geometry.initialize_geometry()
        geometry.delta = self.core.width
        geometry.build_geometry()
        
        # Set the group data.
        group_data.d1  = np.zeros(geometry.number_bundles)
        group_data.d2  = np.zeros(geometry.number_bundles)
        group_data.r1  = np.zeros(geometry.number_bundles)
        group_data.a2  = np.zeros(geometry.number_bundles)
        group_data.f1  = np.zeros(geometry.number_bundles)
        group_data.f2  = np.zeros(geometry.number_bundles)
        group_data.s12 = np.zeros(geometry.number_bundles)
        a = self.core.assemblies
        for i in range(0, geometry.number_bundles) :
            D1,D2,A1,A2,F1,F2,S12 = a[i].get_constants()
            group_data.d1[i]  = D1
            group_data.d2[i]  = D2
            group_data.r1[i]  = A1 + S12
            group_data.a2[i]  = A2
            group_data.f1[i]  = F1
            group_data.f2[i]  = F2
            group_data.s12[i] = S12
                      
        # Initialize everything                      
        # 0.9301625652404083	-0.014350080586832403	0.0431344280309891
        coefficients.initialize_coefficients()
        coefficients.set_model(0.9302, -0.0144, 0.0431)
        state.initialize_state()
        solver.initialize_solver()
        #solver.verbose = 1
        # Core size per dimension.
        self.dimension = len(self.core.stencil[0,:])
        # Peaking factor map
        self.peaking  = np.zeros((self.dimension, self.dimension))
        
    def evaluate(self) :
        """  Evaluate the current core.
        """
        # The core is a member variable, and so the updated
        #   one is always available once shuffle is called
        #   in Reactor.
        
        # Update the p0attern.
        geometry.set_pattern(self.core.pattern+1)
        # Solve the problem.
        solver.solve()
        # Get the solution.
        self.keff = state.get_keff()
        self.maxpeak = state.get_mppf()
        state.make_peaking_map()
        self.peaking = state.peaking_map
        # Return the evaluation parameters
        return self.keff, self.maxpeak
        
    def display(self) :
        """  Introduce myself.
        """
        print "    EVALUATOR:  FLARE"
   
