'''
Created on Feb 17, 2012

@author: robertsj
'''
from optimizer import Optimizer

try :
  from pypgapack import PGA
except ImportError :
  print "PGA could not be loaded; GA unavailable."
  raise

try :
  from mpi4py import MPI
except ImportError :
  print "MPI could not be loaded"
  pass

import numpy as np

class OptimizerGA(Optimizer, PGA) :
    """  Optimizes using a genetic algorithm.
    """
    
    # Public Interface
    
    def __init__(self, argv, reactor, objective):
        """ Constructor.  
        """
        
        # Initialize the Optimizer object.
        Optimizer.__init__(self, argv, reactor, objective)

        # Initialize the PGA object.
        PGA.__init__(self, argv, PGA.DATATYPE_INTEGER, self.reactor.number_bundles(), PGA.MAXIMIZE)
        
        # Set default operators.
        self.SetCrossover(self.htbx)            # Crossover
        self.SetEndOfGen(self.end_of_iteration) # End of generation info
        self.SetInitString(self.init)           # String initialization
        self.SetMutation(self.swap)             # Mutation via a single swap
        
        # Set default values.
        self.maximum_generations    = 100   # Maximum generations
        self.population_size        = 50    # Population size
        self.number_replaced        = 40    # Number replaced each generation
        self.seed                   = 123   # PGA random number seed
        self.np_seed                = 123   # NumPy random number seed
        self.binary_sweep           = False # Perform one sweep of binary exchanges
        
        # Optimizer-specific flags.
        self.track_best = False
        self.fixed_central = True
        # Counter for evaluations on each process.
        self.evals = 0 
        
    def execute(self):
        """  Optimize the reactor for the objective.
        """
        
        # Set PGA parameters
        self.SetMaxGAIterValue(self.maximum_generations)
        self.SetPopSize(self.population_size)
        self.SetNumReplaceValue(self.number_replaced)
        self.SetNoDuplicatesFlag(PGA.TRUE) 
        self.SetRandomSeed(self.seed)    
        np.random.seed(self.np_seed)       
        
        self.run(self.fun)
        
    def end_of_iteration(self):
        """  Do something at the end of each generation.
        
        In general, this a very customizable routine.  This is where 
        hill-climbing heuristics can be placed.  Additionally, tracking
        of objectives as a function of generation is easily done.  For this
        default implementation, the best keff and maxpeak are kept for 
        each generation.
        """
        best                = self.GetBestIndex(PGA.OLDPOP)
        bestpattern         = self.GetIntegerChromosome(best, PGA.OLDPOP)
        it                  = self.GetGAIterValue()
        # Note, we need to reshuffle and reevaluate for the best pattern.
        #   This is *probably* more efficient than keeping k and p for 
        #   all evaluations and then trying to find the values corresponding
        #   to the bestpattern.
        self.reactor.shuffle(bestpattern) 
        self.k, self.p = self.reactor.evaluate()

        # Tell mother what we've done.
        Optimizer.fire_signal(self)

        #print "iter = ",iter
        #print " it = ", it, " best = ", self.GetEvaluation(best, PGA.NEWPOP), " k p = ",k,p,bestpattern
        #print " *** ", self.GetEvaluation(best, PGA.OLDPOP)  
        if self.track_best :
            self.best_eval[it-1] = self.GetEvaluation(best, PGA.OLDPOP)                                                
            self.best_k[it-1] = self.k
            self.best_p[it-1] = self.p
            self.best_pattern[it-1,:] = bestpattern
        del bestpattern
        
    # Implementation
    
    def run(self, f):
        """  Optimize the reactor for the objective f
        """
        # Setup the optimization.
        PGA.SetUp(self)
        # Setup arrays to track best values.
        if self.track_best :
            max_iter = self.GetMaxGAIterValue()+1 # +1 since initial population is included
            num_fuel = self.GetStringLength()
            self.best_eval = np.zeros(max_iter)
            self.best_k = np.zeros(max_iter)
            self.best_p = np.zeros(max_iter)
            self.best_pattern = np.zeros((max_iter, num_fuel),dtype='i')
        PGA.Run(self, f)
    
    def fun(self, p, pop):
        """  This wraps the user-defined objective function.
        """
        pattern = self.GetIntegerChromosome(p, pop)
        self.reactor.shuffle(pattern) 
        k, p = self.reactor.evaluate()
        val = self.objective(k, p)
        return val
    
    def htbx(self, p1, p2, pop1, c1, c2, pop2) :
        """  Heuristic tie-breaking cross-over.
        
        Note, this implementation *assumes* that the assemblies
        are sorted by reactivity *before* optimization.  The 
        only geometric constraint considered is a fixed central
        bundle, and admittedly, the treatment isn't as clean as
        it could be.  Future work...
        """
        # Grab the city id's.
        paren1 = self.GetIntegerChromosome(p1, pop1)
        paren2 = self.GetIntegerChromosome(p2, pop1)
        child1 = self.GetIntegerChromosome(c1, pop2)
        child2 = self.GetIntegerChromosome(c2, pop2) 
        
        #print "P1 =", paren1
        #print "P2 =", paren2
         
        # Account for fixed central, which reduces
        # the number to crossover by 1 
        n = self.GetStringLength()
        size = n
        start = 0
        if self.fixed_central :
            size = n - 1
            start = 1
             
        # Copy the parents to temporary vector for manipulation.
        parent1 = np.zeros(size)
        parent2 = np.zeros(size)
        for i in range(start, n) :
            parent1[i-start] = paren1[i]
            parent2[i-start] = paren2[i]               
        
        # Code the parents using "position listing", which, because
        #   the pattern indices correspond to sorted assemblies, is
        #   equivalent to "reactivity listing".
        code1 = np.zeros(size)
        code2 = np.zeros(size)
        for i in range(0, size) :
            code1[parent1[i]] = i + 1
            code2[parent2[i]] = i + 1
        
        # Randomly choose two cross-over points.  The crossover is
        #   based on the one-dimensional pattern representation.  
        #   One could also work with 2-d approach, where a "chunk"
        #   or "chunks" are selected.  It's not clear what effect this
        #   would have, but it's worth investigating.
        perm = np.random.permutation(size)         
        point1 = np.min(perm[0:2])
        point2 = np.max(perm[0:2]) + 1 
        
        # Exchange all alleles between the two points.
        temp = np.zeros(point2 - point1)
        for i in range(point1, point2) :
            temp[i - point1] = parent1[i]
            parent1[i] = parent2[i]
            parent2[i] = temp[i - point1] 
        
        # Generate a cross-over map, a random ordering of the 0,1,...,n-1
        crossovermap = np.random.permutation(size) 
        
        # Multiply each allele of the strung by n and add the map.
        parent1 = parent1 * size + crossovermap
        parent2 = parent2 * size + crossovermap
        
        # Replace the lowest allele by 0, the next by 1, up to n-1.  Here,
        #   we sort the parents first, and then for each element, find
        #   where the increasing values are found in the original.  There
        #   is probably a simpler set of functions built in somewhere.
        sort1 = np.sort(parent1)
        sort2 = np.sort(parent2)
        for i in range(0, size) :
            index = np.where(parent1 == sort1[i])
            parent1[index[0][0]] = i
            index = np.where(parent2 == sort2[i])
            parent2[index[0][0]] = i 
            
        # Map the string back to elements.  These are the offspring.
        tempchild1 = np.zeros(size)
        tempchild2 = np.zeros(size)
        for i in range(0, size) :
            tempchild1[parent1[i]] = i
            tempchild2[parent2[i]] = i            
        for i in range(start, n) :
            child1[i] = tempchild1[i-start]
            child2[i] = tempchild2[i-start]   
             
        # Copy over the fixed central element.
        if self.fixed_central :
            child1[0] = paren1[0]
            child2[0] = paren2[0]
            
        #print "C1 =", child1
        #print "C2 =", child2
        
    def init(self, p, pop) :
        """  Random initial states.  
        
        Currently, only the central bundle may be fixed.
        """
        n = self.GetStringLength()        
        pattern = self.GetIntegerChromosome(p, pop)
        
        # perm is a random permutation of [0,1,... ].  If       
        if self.fixed_central:
            perm = np.random.permutation(n - 1)
            pattern[0] = n - 1
            for i in range(1, n) :
                pattern[i] = perm[i - 1]            
        else:
            perm = np.random.permutation(n - 1)
            for i in range(1, n + 1) :
                pattern[i] = perm[i - 1]     
        del pattern
        
    def swap(self, p, pop, pm) :
        """  Random swap of bundles.
        
        This example allows swapping for all but the central
        element.  Moreover, no checking is done to ensure swaps
        are not done between identical or forbidden bundles.  In
        those cases, a swap simply doesn't happen.  This is one
        area where a lot more work can be done w/r to implementation.
        """
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
                    
        
    def set_track_best(self, value=True) :
        """ Track the best evaluations.
        """
        self.track_best = value
