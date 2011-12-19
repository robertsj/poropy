'''
Created on Dec 14, 2011

@author: robertsj
'''

class Assembly(object):
    """ Represents a fuel assembly.
    """
    
    def __init__(self,                                      \
                 model="IFBA",                              \
                 enrichment=4.0,                            \
                 burnup=1.0,                                \
                 data=[1.4493e+00,3.8070e-01,9.9000e-03,    \
                       1.0420e-01,7.9000e-03,1.6920e-01,    \
                       1.5100e-02]                          ):
        """  Constructor.  Default is a fresh 17x17 IFBA assembly.
        """
        # Identifier for potential use in data model, e.g. what type
        #   of assembly this is (WABA, IFBA, etc.)
        self.model = model
        # Initial enrichment w/o (also for use with model)
        self.enrichment = enrichment
        # Initial burnup MWd/kg (also for use with model)
        self.burnup = burnup        
        # Default group data.
        self.D1  = data[0]
        self.D2  = data[1]
        self.A1  = data[2]
        self.A2  = data[3]   
        self.F1  = data[4]
        self.F2  = data[5]  
        self.S12 = data[6]   
        self.R1  = self.A1 + self.S12
        
    def set_constants(self, data):    
        '''
        Explicitly change the group constants.
        '''
        self.D1  = data[0]
        self.D2  = data[1]
        self.A1  = data[2]
        self.A2  = data[3]   
        self.F1  = data[4]
        self.F2  = data[5]  
        self.S12 = data[6]
        self.R1  = self.A1 + self.S12
   
    def get_constants(self):    
        """ Explicitly change the group constants.
        """
        return self.D1, self.D2, self.A1,  self.A2,  \
               self.F1, self.F2, self.S12 
               
    def get_constants_list(self):    
        """ Explicitly change the group constants.
        
        This version returns a list, which might be convenient.
        """
        return [self.D1, self.D2, self.A1,  self.A2,  \
                self.F1, self.F2, self.S12]   
    
    def kinf(self) :
        """  Return kinf.
        """
        return (self.F1 + self.F2 * self.S12 / self.A2) / self.R1  
    
    def burn(self, burnup):
        '''
        Stub method for inserting a data model.
        '''
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
    """ Represents a reflector.
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
               
