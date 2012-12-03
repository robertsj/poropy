'''
Created on Feb 17, 2012

@author: robertsj
'''


import numpy as np

from poropy.coretools import Reactor, Assembly, Reflector, Flare, Laban

import plugin


class Yamamoto(plugin.Plugin):
    """  Defines metadata for the large core setup
    """
    
    def __init__(self,model,parent=None):
        plugin.Plugin.__init__(self,model,parent=parent)
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
        return "Yamamoto's Benchmark"

    def description(self):
        return """"""
    
    def build(self):
        """  Execute the parent class. 
        """
        if not self.reactor:
          self.reactor = self.make_large_core()
          
        self.model.set_reactor(self.reactor)
        
        self.close_plugin_diag()
        


    def make_large_core(self,rank=0) :
        """ This returns a Reactor object for the small benchmark.
        """
        
        # Stencil
        # =======
        
        stencil = np.array([[2, 1, 1, 1, 1, 1, 1, 1, 0], \
                            [0, 1, 1, 1, 1, 1, 1, 1, 0], \
                            [0, 1, 1, 1, 1, 1, 1, 0, 0], \
                            [0, 1, 1, 1, 1, 1, 1, 0, 0], \
                            [0, 1, 1, 1, 1, 1, 0, 0, 0], \
                            [0, 1, 1, 1, 1, 0, 0, 0, 0], \
                            [0, 1, 1, 1, 0, 0, 0, 0, 0], \
                            [0, 1, 0, 0, 0, 0, 0, 0, 0], \
                            [0, 0, 0, 0, 0, 0, 0, 0, 0]] )
        
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

        # Burnups
        burnups = np.array([34.7, 34.7, 16.8, 19.5,  0.0, 28.9, 12.6, 12.6,  \
                                  32.7,  0.0, 23.8, 19.0, 18.4,  0.0,  0.0,  \
                                   0.0, 23.0, 18.9, 27.9,  0.0, 10.2,        \
                                  23.8, 18.9, 32.2,  0.0, 11.3,  0.0,        \
                                  19.0, 27.9,  0.0, 11.3,  0.0,              \
                                  18.4,  0.0, 11.3,  0.0,                    \
                                   0.0, 10.2,  0.0,                          \
                                   0.0], dtype='d')

        # Pattern
        pattern = np.array([ 2, 1, 0, 2, 0, 1, 2, 2, \
                                1, 2, 0, 1, 1, 0, 2, \
                                2, 0, 1, 1, 0, 0,    \
                                0, 1, 1, 1, 1, 0,    \
                                1, 1, 1, 1, 0,       \
                                1, 0, 1, 0,          \
                                0, 0, 0,             \
                                2                    ], dtype='i')
       

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
        
        # Assemblies are built with the following signature:
        #   ('type', enrichment, burnup, array([D1,D2,A1,A2,F1,F2,S12]))
     
        # Poison
        poison = np.array(['GAD',  'GAD',  'GAD',  'GAD',  'GAD', 'None', 'None', 'None',  \
                                   'GAD',  'GAD', 'None',  'GAD',  'GAD',  'GAD', 'None',  \
                                   'GAD', 'None',  'GAD', 'None',  'GAD', 'None',          \
                                  'None',  'GAD',  'GAD',  'GAD', 'None', 'None',          \
                                   'GAD', 'None',  'GAD',  'GAD', 'None',                  \
                                   'GAD',  'GAD', 'None', 'None',                          \
                                   'GAD', 'None', 'None',                                  \
                                  'None'], dtype='S')

        # Loop through and assign assemblies to each fuel location in the pattern.
        for i in range(0, len(burnups)) :
            assemblies.append(Assembly(poison[i], 4.10, burnups[i]))

        # Use the Biblis reflector material (the only one for now).  This is
        # only of use for Laban.
        reflector = Reflector('biblis')
        
        # Assembly dimension.
        width = 21.0000
        
        # Choose a physics model (Laban or Flare).  Note, Laban must be called 
        # with "rank" as a parameter.
        physics = Flare()

        # Build the Reactor
        # =================
        return Reactor(stencil, pattern, assemblies, reflector, width, physics)
        
