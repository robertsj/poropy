'''
Created on Feb 14, 2012

@author: robertsj
'''

import numpy as np

from evaluator import Evaluator

import pyflare

import sys

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
        pyflare.geometry.deallocate_geometry()
        pyflare.geometry.stencil = self.core.stencil
        pyflare.geometry.number_bundles = len(self.core.pattern)
        pyflare.geometry.stencil_dimension = len(self.core.stencil[0,:])
        pyflare.geometry.initialize_geometry()
        pyflare.geometry.delta = self.core.width
        pyflare.geometry.build_geometry()

        # Set the group data
        pyflare.group_data.deallocate_group_data()
        pyflare.group_data.allocate_group_data(pyflare.geometry.number_bundles)
        pyflare.group_data.d1  = np.zeros(pyflare.geometry.number_bundles)
        pyflare.group_data.d2  = np.zeros(pyflare.geometry.number_bundles)
        pyflare.group_data.r1  = np.zeros(pyflare.geometry.number_bundles)
        pyflare.group_data.a2  = np.zeros(pyflare.geometry.number_bundles)
        pyflare.group_data.f1  = np.zeros(pyflare.geometry.number_bundles)
        pyflare.group_data.f2  = np.zeros(pyflare.geometry.number_bundles)
        pyflare.group_data.s12 = np.zeros(pyflare.geometry.number_bundles)

        a = self.core.assemblies
        for i in range(0, pyflare.geometry.number_bundles) :
            D1,D2,A1,A2,F1,F2,S12 = a[i].get_constants()
            pyflare.group_data.d1[i]  = D1
            pyflare.group_data.d2[i]  = D2
            pyflare.group_data.r1[i]  = A1 + S12
            pyflare.group_data.a2[i]  = A2
            pyflare.group_data.f1[i]  = F1
            pyflare.group_data.f2[i]  = F2
            pyflare.group_data.s12[i] = S12
        
        # Initialize everything                      
        # 0.9301625652404083	-0.014350080586832403	0.0431344280309891
        pyflare.coefficients.initialize_coefficients()
        pyflare.coefficients.set_model(0.9302, -0.0144, 0.0431)
        pyflare.state.initialize_state()
        pyflare.solver.initialize_solver()
        #solver.verbose = 1
        # Core size per dimension.
        self.dimension = len(self.core.stencil[0,:])
        # Peaking factor map
        self.peaking = np.zeros(pyflare.geometry.number_bundles)
        self.peaking_map = np.zeros((self.dimension, self.dimension))
        
    def evaluate(self) :
        """  Evaluate the current core.
        """
        # The core is a member variable, and so the updated
        #   one is always available once shuffle is called
        #   in Reactor.
#        print "evaluating with flare"
        # Update the pattern.
        pyflare.geometry.set_pattern(self.core.pattern+1)
        # Solve the problem.
        pyflare.solver.solve()
        # Get the solution.
        self.keff = pyflare.state.get_keff()
        self.maxpeak = pyflare.state.get_mppf()
        pyflare.state.make_peaking_map()
        # Get the peaking map and a 1-d array of peaking
        self.peaking_map = pyflare.state.peaking_map
        self.peaking_map[:, 0] = self.peaking_map[0, :]
        self.peaking = pyflare.state.assembly_peaking

        # Return the evaluation parameters
        return self.keff, self.maxpeak
        
    def display(self) :
        """  Introduce myself.
        """
        print "    EVALUATOR:  FLARE"
   
