# examples/large_core.py -- large benchmark reactor

import numpy as np

from poropy.coretools import Reactor, Assembly, Reflector, Laban, Flare

def make_large_core(rank=0) :
    """ This returns a Reactor object for the small benchmark.
    """
    
    # Stencil
    # =======
    
    stencil = np.array([[2, 1, 1, 1, 1, 1, 1, 1, 0], \
                        [0, 1, 1, 1, 1, 1, 1, 1, 0], \
                        [0, 1, 1, 1, 1, 1, 1, 1, 0], \
                        [0, 1, 1, 1, 1, 1, 1, 1, 0], \
                        [0, 1, 1, 1, 1, 1, 1, 0, 0], \
                        [0, 1, 1, 1, 1, 1, 1, 0, 0], \
                        [0, 1, 1, 1, 1, 1, 0, 0, 0], \
                        [0, 1, 1, 1, 0, 0, 0, 0, 0], \
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

    # Checkerboard.
#    pattern = np.array([ 2, 0, 1, 0, 1, 0, 1, 2,  \
#                            1, 0, 1, 0, 1, 0, 2,  \
#                            0, 1, 0, 1, 0, 1, 2,  \
#                            1, 0, 1, 0, 1, 2, 2,  \
#                            0, 1, 0, 1, 0, 2,     \
#                            1, 0, 1, 0, 1, 2,     \
#                            0, 1, 2, 2, 2,        \
#                            2, 2, 2               ], dtype='i')

    # Ring of fire.
    pattern = np.array([ 2, 1, 0, 2, 0, 1, 2, 2, \
                            1, 2, 0, 1, 1, 0, 2, \
                            2, 0, 1, 1, 0, 0, 2, \
                            0, 1, 1, 1, 1, 0, 2, \
                            1, 1, 1, 1, 0, 0,    \
                            1, 0, 1, 0, 1, 2,    \
                            0, 0, 0, 0, 2,       \
                            2, 2, 2              ], dtype='i')
   

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
    unique_assemblies.append(Assembly('IFBA', 4.25, 0.0))
    # Once burned                  
    unique_assemblies.append(Assembly('IFBA', 4.25, 15.0))
    # Twice burned                  
    unique_assemblies.append(Assembly('IFBA', 4.25, 30.0))

    # Loop through and assign assemblies to each fuel location in the pattern.
    for i in range(0, len(pattern)) :
        assemblies.append(unique_assemblies[pattern[i]])

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
