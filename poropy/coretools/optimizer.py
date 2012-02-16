# poropy/coretools/optimizer.py  --  optimization tools.

from pypgapack import PGA
from mpi4py import MPI
import numpy as np

class Optimizer(PGA) :
    """  Derive our own class from PGA.
    """
    
    def __init__(self, argv, reactor) :
        """ Constructor.
        """
        self.reactor = reactor
        self.evals = 0 # counter for evaluations on each process
        PGA.__init__(self, argv, PGA.DATATYPE_INTEGER, self.reactor.number_bundles(), PGA.MAXIMIZE)
        # Set default operators and other settings.
        self.SetCrossover(self.htbx)        # Crossover
        self.SetEndOfGen(self.eog)          # End of generation info
        self.SetNoDuplicatesFlag(PGA.TRUE)  # Keep no duplicate patterns.   
        # Optimizer-specific flags.
        self.track_best = False
        self.fixed_central = False
       
    def run(self, f) :
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
        
    def set_track_best(self, value=True) :
        """ Track the best evaluations.
        """
        self.track_best = value
        
    def set_fixed_central(self, value=True) :
        """ Fix the central bundle.
        
        By construction, the central bundle is the 
        first entry in the pattern.  This flag will
        indicate crossover should skip this.
        """
        self.fixed_central = value
        
    def objective(self, p, pop) :
        """ Minimize peaking and maximize keff using weighted objective.
        
        The default seeks to maximize k for p < 1.55.
        """
        pattern = self.GetIntegerChromosome(p, pop)
        self.reactor.shuffle(pattern) 
        k, p = self.reactor.evaluate()
        del pattern
        delta = 0
        if p > 1.55 :
            delta = p - 1.55
        val = 1.0 * k - 1.0 * delta
        self.evals += 1
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

    def eog(self):
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
        k, p = self.reactor.evaluate()
        #print "iter = ",iter
        #print " it = ", it, " best = ", self.GetEvaluation(best, PGA.NEWPOP), " k p = ",k,p,bestpattern
        #print " *** ", self.GetEvaluation(best, PGA.OLDPOP)  
        if self.track_best :
            self.best_eval[it-1] = self.GetEvaluation(best, PGA.OLDPOP)                                                
            self.best_k[it-1] = k
            self.best_p[it-1] = p
            self.best_pattern[it-1,:] = bestpattern
        del bestpattern
#        self.bestpeak[it-1] = p  
