# poropy/coretools/reactor.py

# things from coretools
from laban import Laban
from assembly import Assembly, Reflector
# others
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import rc


# Allows LaTeX in labels and uses nicer serif font.
rc('text', usetex=True)
rc('font', family='serif')

class Reactor(object):
    """ Represents the reactor
    """
    
    def __init__(self, stencil, pattern, assemblies, reflector, width, evaluator) :
        """ Constructor
        """
  
        # Create the core.
        self.core = Core(stencil, pattern, assemblies, reflector, width)  
        # Create the spent fuel pool.  Currently not implemented.
        self.pool = SpentFuelPool()
        # Create the evaluator.
        self.evaluator = evaluator
        # Set up the evaluator, i.e. the neutronics kernel
        self.evaluator.setup(self.core)
        
        # Some things that aren't used now but may be important
        #   for later additions.
        self.power_thermal  = 1000 # MWth
        self.power_electric =  300 # MWe
    
    def shuffle(self, pattern) :
        """ Update the pattern.  
        
        This applies *only* to within-core shuffling.  Accounting for different
        new fuel or swaps with the pool will need something more to handle the
        assemblies issue.
        """
        self.core.update_pattern(pattern)
        
    def swap(self, x, y) :
        """  Swap two bundles.
        """
        self.core.swap(x, y)
        
    def evaluate(self) :  
        """ Evaluate the current core.  Pattern must be up-to-date.
        """
        # currently returns only keff and maxpeaking
        return self.evaluator.evaluate()
    
    def number_bundles(self) :
        """  Return the number of fuel bundles.
        """
        return len(self.core.pattern)
    
    def display(self) :
        """  Print out all information (for debugging).
        """
        print " poropy - diagnostic output "
        print ""
        print "REACTOR:"
        print "    thermal power :" + '%5i' % (self.power_thermal) + " MWth"
        print "   electric power :" + '%5i' % (self.power_electric) + " MWel"
        self.core.display()
        self.pool.display()
        self.evaluator.display()
    
    def print_peaking(self) :
        """  Print out the peaking factor matrix.
        """
        # TODO(robertsj): Add some nice formatting.
        print ""
        print self.evaluator.peaking
        
    def print_pattern(self) :
        """  Print out the peaking factor matrix.
        """
        print ""
        print " current loading pattern "
        print " ----------------------- "
        self.core.print_pattern()   
        
    def print_params(self) :
        """  Print out optimization parameters.
        """
        print ""
        print " optimization parameters "
        print " ----------------------- "
        print "    keff   = ", self.evaluator.keff    
        print " maxpeak   = ", self.evaluator.maxpeak  
        
class Core(object):
    """ Represents the core.
    """
    
    def __init__(self, stencil, pattern, assemblies, reflector, width):
        """ Constructor
        """
        # Validate input.
        assert(len(pattern) == len(assemblies))
        assert(len(assemblies) > 0)
        assert(width > 0)

        # Assign the stencil
        self.stencil = stencil
        # Assign the pattern and assemblies
        self.pattern, self.assemblies = self.sort_assemblies(pattern, assemblies)
        # Make the fuel map
        self.make_fuel_map()
        # Assign the reflector
        self.reflector = reflector
        # Assign assembly dimension
        self.width = width
  
    def sort_assemblies(self, pattern, assemblies) :
        """ Sort the assemblies by reactivity.
        """ 
        # TODO(robertsj): Consider a cleaner approach for this sorting.
        
        # We build a 2-d array of [index,kinf] pairs.  Sorting this gives
        #   permuted index in the first entry.  The location of each
        #   original index will become the new pattern.  (Note that kinf
        #   is negated so we get descending order of reactivity. It seems
        #   argsort has no option for ascend/descend.
        pattern_length = len(pattern)
        index = np.zeros((pattern_length,2))
        for i in range(0, pattern_length) :
            index[i][0] = i
            index[i][1] = -assemblies[i].kinf()
        index=index[index[:,1].argsort(),0]
        
        # Define the sorted pattern and assemblies using the permuted 
        #   indices. Note that each pattern element will be unique, even 
        #   if a small number of unique assemblies  defined the pattern
        #   initially.
        sorted_pattern = np.zeros(len(pattern),dtype='i')
        sorted_assemblies = []
        for i in range(0, pattern_length) :
            sorted_pattern[i] = (np.where(index == i))[0][0]
            sorted_assemblies.append(assemblies[int(index[i])])

        return sorted_pattern, sorted_assemblies
            
    def make_fuel_map(self) :
        """  Map the fuel to (i,j) for easy swapping.
        """
        N = len(self.stencil[0,:])
        self.fuelmap = np.zeros((N, N),'i')
        self.fuelidx = np.ones((N, N),'i')*-1
        fuelindex = 0
        for i in range(0, N) :
            if i == 0 :
                jj = 0
            else :
                jj = 1
            for j in range(jj, N) :        
                if self.stencil[i, j] > 0 :    # a fuel
                    self.fuelmap[i, j] = self.pattern[fuelindex]
                    self.fuelidx[i, j] = fuelindex
                    fuelindex += 1
                elif self.stencil[i, j] == 0 : # a reflector
                    self.fuelmap[i, j] = -1
                else :                    # a void
                    self.fuelmap[i, j] = -2 
                    
        # Assign the ones constrained by symmetry.
        self.fuelmap[1:,0]=self.fuelmap[0,1:]     
          
    def update_pattern(self, pattern) :
        """  Update the pattern.
        """
        # Ensure we're not adding fuel (not for now anyway!).
        assert(len(pattern)==len(self.pattern))
        # Assign new pattern
        self.pattern[:] = pattern
        # and then update the map
        self.make_fuel_map()
        #TODO(robertsj): There should be checking.
     
    def swap(self, x, y) :
        """ Swap two bundles.
        """
        assert(len(x) == len(y))
        assert(len(x) == 2)
        assert(min(x) >= 0 and max(x) < len(self.fuelmap[0, :]))
        assert(min(y) >= 0 and max(y) < len(self.fuelmap[0, :]))
        if self.fuelmap[x[0], x[1]] < 0 :
            print " Bundle 1 isn't fuel!! "
            return
        if self.fuelmap[y[0], y[1]] < 0 :
            print " Bundle 2 isn't fuel!! "
            return

        tmp = self.fuelmap[x[0],x[1]]
        self.pattern[self.fuelidx[x[0],x[1]]] = self.pattern[self.fuelidx[y[0],y[1]]]
        self.pattern[self.fuelidx[y[0],y[1]]] = tmp
#        tmp = self.pattern[self.fuelmap[x[0],x[1]]]
#        self.pattern[self.fuelidx[x[0],x[1]]] = self.fuelidx[y[0],y[1]]
#        self.pattern[self.fuelidx[y[0],y[1]]] = tmp
        self.make_fuel_map()
        
    def inventory_size(self):
        """ Return the inventory size.
        """
        return len(self.assemblies)
 
    def print_pattern(self):
        """  Print the pattern in a nice 2-D format.
        """
        N = len(self.stencil[0,:]) 
        out = "    "
        for i in range(0, N-1) :
            out += '%4i' % (i)
        out += '\n' + "  --"
        for i in range(1, N) :
            out += "----"
        out += '\n' 
        for i in range(0, N-1) :
            out += '%2i' % (i)
            out += '| '
            if i == 0 :
                jj = 0
            else :
                jj = 1
                out+="  rs"
            for j in range(jj, N-1) :
                if self.fuelmap[i, j] >= 0 :
                    out += '%4i' % (self.fuelmap[i, j])
                elif self.stencil[i, j] == 0  :
                    out += ' ref'
                else :
                    out += '    '
            out += '\n'
        print out
        
    def display(self) :
        """  Print my information.
        """
        print "    CORE:           "
        print "        pattern:    ", self.pattern
        print "        ASSEMBLIES: "
        for i in range(0, len(self.assemblies)) :
            print "             assembly: ", i
            self.assemblies[i].display() 
        print "        REFLECTOR: "
        self.reflector.display()
 
class SpentFuelPool(object):
    """ Represents the spent fuel pool. *NOT IMPLEMENTED*
    """
    
    def __init__(self):
        """ Constructor.  
        """

    def display(self):
        """ Display my contents.
        """
        pass
