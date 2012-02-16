'''
Created on Dec 14, 2011

@author: robertsj
'''
from poropy.nucleardata import *

class Assembly(object):
    """  Represents a fuel assembly.

    Attributes
    ----------
        model : str
            Basic description of assembly, e.g. IFBA.
        enrichment : float
            Initial enrichment (weight percent)
        burnup : float
            Current burnup (MWd/kg)

    """
    
    def __init__(self,            \
                 model="IFBA",    \
                 enrichment=4.0,  \
                 burnup=0.0       ):
        """  Constructor.  

        By default, the built-in data module is used.  Users may 
        explicitly set constants for single statepoints via the
        set_constants command.
        """

        # Set characteristics used to define the group constants
        # within the data module.
        self.model = model
        self.enrichment = enrichment
        self.burnup = burnup 
       
        # Get the group data object.
        IFBA = 0
        WABA = 0
        GAD  = 0
        if (model == "IFBA"):
            IFBA = 1
        if (model == "WABA"):
            WABA = 1   
        if (model == "GAD"):
            GAD = 1  
        data = get_2g_parms(burnup, enrichment, IFBA, WABA, GAD)

        # Set the group data.
        self.D1   = data.DIFF1
        self.D2   = data.DIFF2
        self.A1   = data.ABS1
        self.A2   = data.ABS2
        self.F1   = data.NUFISS1
        self.F2   = data.NUFISS2
        self.S12  = data.REMOV1              # Note, this is strange but is
        self.R1   = data.REMOV1 + data.ABS1  # the way CASMO defines them.
        self.KINF = data.K_INF_NO_XE
        self.M2   = data.M2_XE

        # Peaking factor.
        self.peak = 0.0
        
    def set_constants(self, data=[1.4493e+00,3.8070e-01,9.9000e-03, \
                                  1.0420e-01,7.9000e-03,1.6920e-01, \
                                  1.5100e-02]                       ):
        """  Explicitly change the group constants.

        This is for use when a user wishes to bypass the data module for
        a run at a *single* state point, i.e. no burnup, etc.
        """
        self.D1   = data[0]
        self.D2   = data[1]
        self.A1   = data[2]
        self.A2   = data[3]   
        self.F1   = data[4]
        self.F2   = data[5]  
        self.S12  = data[6]
        self.R1   = self.A1 + self.S12
        self.KINF = (self.F1 + self.F2 * self.S12 / self.A2) / self.R1  
        self.M2   = self.D1/self.R1 + self.D2/self.A2
   
    def get_constants(self):    
        """ Grab individual group constants.
        """
        return self.D1, self.D2, self.A1,  self.A2,  \
               self.F1, self.F2, self.S12 
               
    def get_constants_list(self):    
        """ Grab group constants as a list.
        """
        return [self.D1, self.D2, self.A1,  self.A2,  \
                self.F1, self.F2, self.S12]   
    
    def kinf(self):
        """  Return kinf.
        """
        return self.KINF

    def peak(self):
        """  Return my peaking factor.
        """
        return self.peak

    def set_peak(self, p):
        """  Update my peaking factor.
        """
        self.peak = p
    
    def burn(self, burnup):
        """  Stub method for inserting a data model.
        """
        # Given the burnup, enrichment, and assembly type,
        #   the data model returns a new set of data.           
        
    def display(self):
        """ Display my contents.
        """
        print "                model: ", self.model
        print "           enrichment: ", self.enrichment
        print "               burnup: ", self.burnup        
        print "            constants: ", self.D1, self.D2, self.A1,  self.A2, \
                                self.F1, self.F2, self.S12, self.R1   
        print "                 kinf: ", self.kinf() 
        
class Reflector(Assembly):
    """ Represents a reflector.
    """
    
    def __init__(self, model):
        """ Constructor.  Default data corresponds to a fresh IFBA assembly.
        """
        if model == "biblis" :
            data = [1.3200e+00, 2.7720e-01, 2.6562e-03, \
                    7.1596e-02, 0.00000000, 0.00000000, \
                    2.3106e-02]   
        else :
            raise Exception("Reflector model not available")
        self.model = model
        # Default group data.
        self.D1  = data[0]
        self.D2  = data[1]
        self.A1  = data[2]
        self.A2  = data[3]   
        self.F1  = data[4]
        self.F2  = data[5]  
        self.S12 = data[6]   
        self.R1  = self.A1 + self.S12
        
    def display(self):
        """ Display my contents.
        """
        print "                model: ", self.model
        print "            constants: ", self.D1, self.D2, self.A1,  self.A2, \
                                   self.F1, self.F2, self.S12, self.R1   
        
    def burn(self, burnup):
        """ Stub method for inserting a data model.
        """
        print "Prepare to be extinguished."
        
class Void(Assembly):
    """ Represents a void.
    """
    
    def __init__(self):
        """ Constructor.  Default data corresponds to a fresh IFBA assembly.
        """
        pass

def assembly_compare(x, y) :
    """  Compare assemblies based on BOC kinf; sorted by descending reactivity.
    """
    if x.kinf() < y.kinf() :
        return 1
    elif x.kinf() == y.kinf() :
        return 0
    else :  #x.resultType < y.resultType
        return -1 
               
