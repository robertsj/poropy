'''
Created on Dec 11, 2011

@author: robertsj
'''
from evaluator import Evaluator
import numpy as np
import os


class Laban(Evaluator):
    """ Uses the LABAN-PEL code to evaluate loading patterns.
    """

    # Public Interface

    def __init__(self, rank=0) :
        """ Constructor
            
        Parameters
        ----------
        rank : integer
            MPI process rank.  
        """
        
        # LABAN options of potential interest
        self.order     = 4
        self.tolerance = 0.0001
        
        # Rank is needed for MPI runs.  This differentiates the inputs. 
        self.rank      = rank
        self.input     = "laban"+str(rank)+".inp"
        self.output    = "laban"+str(rank)+".out"    
        
        self.keff = 0
        self.maxpeak = 0
        
    def setup(self, core) :     
        """  Set up some structures needed before solving.
        """
        # We need direct access to the core
        self.core = core
        # Validate a square quarter core. (Not applicable to 1/2 or 1/8)
        assert(len(self.core.stencil[0,:])==len(self.core.stencil[:,0]))
        # Core size per dimension.
        self.dimension = len(self.core.stencil[0,:])
        # Assembly boundaries
        self.widths = np.zeros(self.dimension+1)
        self.widths[:] = self.core.width
        self.widths[0] = 0.5 * self.core.width
        # Subdivisions.  Not really used.
        self.subdivisions = np.ones(self.dimension,dtype='i')
        # Peaking factor map
        self.peaking_map  = np.zeros((self.dimension, self.dimension))
        self.peaking = np.zeros(len(self.core.assemblies))
        # Create the static top part of the LABAN-PEL input
        self.make_input_top()

    def evaluate(self) :
        """  Evaluate the current core.
        """
        # The core is a member variable, and so the updated
        #   one is always available once shuffle is called
        #   in Reactor.
        
        # Open the input file
        f = open(self.input, 'w')
        
        # Write
        f.write(self.input_top)
        self.make_input_map()
        f.write(self.input_map)
        
        self.make_input_materials()
        f.write(self.input_materials)
        
        # Close
        f.close()
        
        # Run LABAN-PEL
        self.run()
        
        # Read the output
        self.read()
        
        # Return the evaluation parameters
        return self.keff, self.maxpeak
    
    def display(self) :
        """  Introduce myself.
        """
        print "    EVALUATOR:  LABAN-PEL"
        print "        input: ", self.input
        print "       output: ", self.output

    
    # Implementation

    def make_input_top(self) :
        """  Create the string for the top invariant part of the input.
        """
        stencil_length = str(len(self.core.stencil[0,:]))
        self.input_top = \
                "LABANPEL -- poropy driver \n"                                +\
                "  2 " + stencil_length  + " " + stencil_length               +\
                " " + str(len(self.core.pattern)+1) + " 0    0    2  1.0\n"   +\
                " "+str(self.order)+" 0    0    0    0    0    0    0     0\n"+\
                "  1    0. 0.     0   100   10 " + str(self.tolerance)+" "    +\
                str(self.tolerance)+" "+str(self.tolerance) + "   1.0\n"      
        # horizontal subdivisions                                                      
        for i in range(0, self.dimension) :
            self.input_top += "  "+str(self.subdivisions[i])
        self.input_top += " \n"
        # vertical subdivisions
        for i in range(0, self.dimension) :
            self.input_top += "  "+str(self.subdivisions[i])
        self.input_top += " \n"            
        # horizontal mesh widths
        for i in range(0, self.dimension) :
            self.input_top += "  "+str(self.widths[i])
        self.input_top += " \n"            
        # vertical mesh widths
        for i in range(0, self.dimension) :
            self.input_top += "  "+str(self.widths[i])
        self.input_top += " \n"                        
 
    def make_input_map(self) :
        """  Print the map of assemblies.
        """

        self.input_map = ""
        stencil = self.core.stencil
        pattern = self.core.pattern
        reflect = len(pattern)+1 # reflector id, last material
        N = self.dimension
        coremap = np.zeros((N+2,N+2), dtype='i')
        
        # reflections and vacuum
        coremap[0,     1:N+1] = -1 
        coremap[1:N+1,     0] = -1
        coremap[N+1,   1:N+1] = -2
        coremap[1:N+1,   N+1] = -2
        
        fuelindex = 0
        
        for i in range(1, N+1) :
            for j in range(1, N+1) :
                if j == 1 and i > 1 :
                    pass
                else :
                    if stencil[i-1, j-1] > 0 :    # a fuel
                        coremap[i, j] = pattern[fuelindex]+1
                        fuelindex += 1
                    elif stencil[i-1, j-1] == 0 : # a reflector
                        coremap[i, j] = reflect
                    else :                    # a void
                        pass  
        # Copy elements such that rotational symmetry is enforced. 
        for j in range(2, N+1) :
            coremap[j, 1] = coremap[1, j]
        for i in range(0, N+2) :
            for j in range(0, N+2) :
                self.input_map +='%4i' % (coremap[i, j])
            self.input_map += '\n'
        
        
    def make_input_materials(self) :
        """  Print the materials definitions.
        """
        #       1    5    1  MATERIAL 1   (arbitrary line, i think)                                         
        #    1.4493e+00   9.9000e-03   7.9000e-03  1.  0.         0. 7.9000e-03  1.
        #    3.8070e-01   1.0420e-01   1.6920e-01  0   1.5100e-02 0. 1.6920e-01  1.
        self.input_materials = ""
        number_mats = len(self.core.pattern)+1
        a = self.core.assemblies
        for i in range(0, number_mats-1) :
            # Row 1: description.
            self.input_materials += "   " + str(i+1) + "    5    1  MATERIAL "   + \
               str(i+1) + "  ("   + \
               a[i].model   + ",  "                            + \
               str(a[i].enrichment) + " w/o,  "                        + \
               str(a[i].burnup)     + " MWd/kg)\n"  
            # Rows 2 and 3.
            D1,D2,A1,A2,F1,F2,S12 = a[i].get_constants()
            d = np.array([[D1,A1,F1,1.0,0.0,0.0,F1,1.0],[D2,A2,F2,0.0,S12,0.0,F2,1.0]])
            for j in range(0, 2) :
                for k in range(0, 8) :
                    self.input_materials +='%12.4e' %(d[j,k])
                self.input_materials += '\n'
            
        a = self.core.reflector
        # Row 1: description.
        self.input_materials += "   " + str(number_mats)  + "    5    1  MATERIAL "   + \
            str(number_mats) +  "  (REFLECTOR) \n"  
        # Rows 2 and 3.
        D1,D2,A1,A2,F1,F2,S12 = a.get_constants()
        d = np.array([[D1,A1,F1,1.0,0.0,0.0,F1,1.0],[D2,A2,F2,0.0,S12,0.0,F2,1.0]])
        for i in range(0, 2) :
            for j in range(0, 8) :
                self.input_materials +='%12.4e' %(d[i,j])
            self.input_materials += '\n'
        self.input_materials += "WHITE\n" + "BLACK\n" + "END\n"
        
 
   
    def read(self) : 
        """  Read a LABAN-PEL output file and load various data.
        """
        # Open the file.
        f     = open(self.output, 'r')
        lines = f.readlines()
        
        # Find the eigenvalue.
        count = 0
        while True :
            words = lines[count].split()
            if len(words) == 5 :
                if words[0] == "*" and words[1] == "K-EFF":
                    self.keff = float(words[3])
                    break
            count += 1
        
        # Find the peaking.
        a = 0 # Assembly index
        
        while True :
            words = lines[count].split()
            if len(words) == 8 :
                if words[0] == "NODE" and words[1] == "AVERAGE" and words[2] == "POWERS" :
                    count += 5 # Powers start 5 lines below title
                    for row in range(0, self.dimension) :
                        words = lines[count].split()
                        assert(len(words) >= self.dimension)
                        for col in range(0, self.dimension) :
                            self.peaking_map[row, col] = float(words[col+1])
                            if self.core.stencil[row, col] > 0:
                                #print " a=", a, " row=", row, " col=", col, len(self.peaking)
                                self.peaking[a] = self.peaking_map[row, col]
                                a += 1
                        count += 1
                    break
            count += 1             
        # Maximum peaking.
        self.maxpeak = np.max(self.peaking)

    def run(self) :
        """  Run LABAN-PEL (must set input first)
        """
        # currently, labanx reads from a preset file
        os.system('labanx '+str(self.rank)+" "+self.input+" "+self.output)     
        
