'''
Created on Feb 17, 2012

@author: robertsj
'''


import numpy as np

from poropy.coretools import Reactor, Assembly, Reflector, Laban, Flare

import plugin

import copy

class SmallCore(plugin.Plugin):
    """  Defines metadata for OptimizerGA
    """
    
    def __init__(self,model,parent):
      plugin.Plugin.__init__(self,model,parent)
      self.reactor = None

    def userInputs(self):
        
        # define the user inputs

        inputs =    (
#                        ("maximum_iterations",  {'type': "freetext", 'label': "maximum_iterations",
#                                                 'default': 123}),
#                        ("population_size",     {'type': "freetext", 'label': "population_size",
#                                                 'default': 123}),
#                        ("number_replaced",     {'type': "freetext", 'label': "number_replaced",
#                                                 'default': 123}),
#                        ("seed",                {'type': "freetext", 'label': "seed",
#                                                 'default': 123}),
#                        ("np_seed",             {'type': "freetext", 'label': "np_seed",
#                                                 'default': 123}),
#                        ("binary_sweep",        {'type': "freetext", 'label': "np_seed",
#                                                 'default': false}),                     
                    )
        return inputs

    def runActions(self):

        # define the buttons

        actions =   (
                       ("build",       {'label': "Build",
                                        'description': "Build the predefined small core geometry"}),
                    )
        return actions

    def title(self):
        return "Small Core Problem"

    def description(self):
        return """69 assembly core"""
    
    def build(self):
        """  Execute the parent class. 
        """
        if not self.reactor:
          self.reactor = self.make_small_core()
          
        self.model.set_reactor(self.reactor)
        
        self.close_plugin_diag()

    def make_small_core(self,rank=0) :
        """ This returns a Reactor object for the small benchmark.
        """
        
        # Problem
        # =======
        
        # This is a small 69 assembly core modeled as a quarter
        # core.  It has three roughly equal batches composed of fresh,
        # once-, and twice-burned fuel.  The materials actually
        # correspond a 17x17, 4.25 w/o, 92 IFBA PWR assembly at 
        # 0, 15, and 30 MWd/kg.  The assembly pitch is 20.5036, 
        # so this benchmark is not using the data realistically.
        # The default pattern is a
        
        # Stencil
        # =======
        
        # First, define the stencil.  Fuel regions are positive
        # integers, where each value corresponds to a "region". 
        # In an optimization sequence, the user can exclude regions
        # from the search space completely.  Currently, this isn't
        # used, but there is an option to fixed the central bundle 
        # directly.
        # Note the whole stencil must be specified, even if certain locations
        # are to be fixed by rotational symmetry.  This makes it easier
        # for the user visually (and makes implementation a bit easier).
        #"""
        stencil = np.array([[2, 1, 1, 1, 1, 0], \
                            [0, 1, 1, 1, 1, 0], \
                            [0, 1, 1, 1, 1, 0], \
                            [0, 1, 1, 1, 0, 0], \
                            [0, 1, 1, 0, 0, 0], \
                            [0, 0, 0, 0, 0, 0]] )
        # Note that 0 indicates reflector.  -1's can be placed as void,
        # but we skip that here.
        
        # Also, the regions are to be indexed in a more natural way than
        # is standard.  The stencil is indexed as would be any matrix.
        # Hence a fuel location [i,j] corresponds to the [i,j]th location
        # in the pattern (using 0-based indexing)  In other words, there is 
        # no alpha-numeric scheme.
        
        # Initial Loading Pattern
        # =======================
        
        # The pattern is a 1-D array of fuel identifiers.  These 
        # identifiers correspond to assemblies to be defined below.
        # Unlike above, fuel is not specified if it is constrained by
        # rotation symmetry.  The fuel indices as 
        # listed correspond to the their location in stencil using
        # a row-major storage.  We work with 1-D data for simplicity.
        #"""
        pattern = np.array([ 2, 0, 1, 0, 2,      \
                                1, 0, 1, 2,      \
                                0, 1, 0, 2,      \
                                1, 0, 1,         \
                                2, 2             ], dtype='i')

        # Assembly Definitions
        # ====================
        
        # Define the assemblies.  The pattern above has three unique
        # values, so we need to define 3 unique assemblies. 
        
        # For simple BOC cycle analysis, it would
        # be sufficient to keep only unique assemblies and then map them to 
        # core locations.  Once burnup is accounted for, assemblies become
        # unique, and having an individual assembly object for each physical
        # assembly makes more sense; that's what we do.  Note, the construction
        # used below assumes the unique assemblies are defined in an order
        # corresponding to their index in pattern.
        
        # Physical assemblies, as it were.
        assemblies = []
        
        # The unique types we have available.
        unique_assemblies = []
     
        # Assemblies are built with the following signature:
        #   ('type', enrichment, burnup, array([D1,D2,A1,A2,F1,F2,S12]))
     
        # Fresh
        unique_assemblies.append(Assembly('IFBA', 4.25, 0.0,  ))
                                         # np.array([1.4493e+00, 3.8070e-01, \
                                         #           9.9000e-03, 1.0420e-01, \
                                         #           7.9000e-03, 1.6920e-01, \
                                         #           1.5100e-02])))
        # Once burned                  
        unique_assemblies.append(Assembly('IFBA', 4.25, 15.0, ))
                                         # np.array([1.4479e+00, 3.7080e-01, \
                                         #           1.1000e-02, 1.2000e-01, \
                                         #           6.9000e-03, 1.7450e-01, \
                                         #           1.4800e-02])))
        # Twice burned                  
        unique_assemblies.append(Assembly('IFBA', 4.25, 30.0,)) 
                                         # np.array([1.4494e+00, 3.6760e-01, \
                                         #           1.1500e-02, 1.1910e-01, \
                                         #           6.0000e-03, 1.6250e-01, \
                                         #           1.4700e-02])))

        # Loop through and assign assemblies to each fuel location in the pattern.
        for i in range(0, len(pattern)) :
            assemblies.append(copy.copy(unique_assemblies[pattern[i]]))

        # Use the Biblis reflector material (the only one for now)
        reflector = Reflector('biblis')
        
        # Assembly dimension.
        width = 23.1226
        
        # Build the Reactor
        # =================
        
        physics = Flare()
#        physics = Laban()
        
        return Reactor(stencil, pattern, assemblies, reflector, width, physics)
        
