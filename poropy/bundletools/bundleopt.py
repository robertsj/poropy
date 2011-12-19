'''
Created on Oct 31, 2011

@author: robertsj
'''

import numpy as np
import os, math, struct

class BundleOpt(object):
    '''
    Optimizes an assembly.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        
    def optimize(self, method) :
        """
        Perform optimization using the specified `method`.

        | We use a lattice of the form
        |
        | `1`                     
        | `1 1`
        | `1 1 1`
        | `1 1 1 1`
        | `1 1 1 1 1`
        | `1 1 1 2 2 1`
        | `1 1 1 2 2 1 1`
        | `1 1 1 1 1 1 1 1`
        | `1 1 1 1 1 1 1 1 1`
        | `1 1 1 1 1 1 1 1 1 1` 
        |
        | where "1" is fuel and "2" is water hole.  The fuel pins are numbered
        | as follows:
        |
        | `1`
        | `2 3`
        | `4 5 6` 

        and so on (skipping over water holes).  There are 51 pins, each of
        which has a fuel enrichment between 2.0 and 4.9% (by 0.1%) and a 
        possible gadolinia of weight percent 1.0 to 10.0% (by 1%).  The 
        density for fuel pins is assumed 10.5 g/cc, and for gad pins, 
        10.2 g/cc. 

        The following are used as helpful heuristics to cut down on some
        variables:
        * no gadolinia on outer edges
        * corner enrichments limited to 2.8\%
        * edge enrichments limited to 4.0\%

        Then, the unknowns are:
        * 51 fuel pin enrichments
        * 32 gad enrichments
        * 83 total
        """
#        # output case-by-case results
#        self.outfile    = open('data.out', 'w')
#        # size
#        n = 83
#        # initial guess
#        x = self.initial_guess('base.inp')     
#        # nlopt class (algorithm, number of variables)
#        #optG = nlopt.opt(nlopt.GN_CRS2_LM, n)        #Controlled Random Search
#        optG = nlopt.opt(nlopt.GN_ISRES, n)           #Improved Stochastic Ranking Evolution Strategy
#        # set lower and upper bounds on the variables
#        lb = []
#        ub = []
#        for i in range(0, 51) :
#            lb.append(2.0)
#            if i == 0 or i == 41 or i == 50 :
#                ub.append(3.0)
#            else :
#                ub.append(4.9)
#        for i in range(51, 83) :
#            lb.append(0.0)
#            ub.append(10.0)
#        # can be lists or arrays
#        optG.set_lower_bounds(lb)
#        optG.set_upper_bounds(ub)
#        # set to a minimization problem 
#        optG.set_max_objective(self.objective_function)
#        # set stopping criteria
#        #opt.set_xtol_rel(1e-15)
#        #opt.set_ftol_rel(1e-15)
#        optG.set_maxeval(1000)
#
#        # write initial
#        outI = " initial guess = \n " + str(x) + "\n" + "\n"     
#        self.outfile.write(outI)     
#
#        # solve
#        x    = optG.optimize(x)
#        maxf = optG.last_optimum_value()
#         
#        outG = " *** GLOBAL ***  \n" +                                 \
#               "      optimum = " + str(x) + "\n" +                    \
#               "minimum value = " + str(maxf) + "\n"                   \
#               "  result code = " + str(optG.last_optimize_result()) + \
#               "\n" + "\n" 
#        
#        # clean up by local
#        optL = nlopt.opt(nlopt.LN_COBYLA, n) 
#        # set lower and upper bounds on the variables
#        optL.set_lower_bounds(lb)
#        optL.set_upper_bounds(ub)
#        # set to a minimization problem 
#        optL.set_max_objective(self.objective)
#        # set stopping criteria
#        #opt2.set_xtol_rel(1e-8)
#        #opt2.set_ftol_rel(1e-15)
#        optL.set_maxeval(100)
#        # solve
#        x    = optL.optimize(x)
#        maxf = optL.last_optimum_value()
#        outL = " *** LOCAL ***  \n" +                                  \
#               "      optimum = " + str(x) + "\n" +                    \
#               "minimum value = " + str(maxf) + "\n"                   \
#               "  result code = " + str(optL.last_optimize_result()) + \
#               "\n" + "\n" 
#               
#        # print to file
#        out = outG + outL
#        self.outfile.write(out)
#        self.outfile.close()  
        return

    def optimize_nlopt(self) :
        """
        Wrapper for nlopt optimization.

        For this initial scoping study, we use the package nlopt, by S.
        Johnson here at MIT's math department.  This package is for 
        continuous real variable optimization with (non)linear bound
        constraints.  The problem statement calls for discrete enrichments,
        but we solve it continously and round the final results (and make
        a note of the difference)
        """       