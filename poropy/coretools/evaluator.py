'''
Created on Feb 14, 2012

@author: robertsj
'''
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

class Evaluator(object) :
    """ Abstract base class for evaluating patterns.
    
    The interface below is the public interface for all evaluators.
    """

    # Abstract Interface

    def __init__(self) :
        """ Constructor.  
        
        This is my approach to an abstract class in Python.  What it does is
        to prevent Evaluator from being instantiated directly.
        """
        if self.__class__ is Evaluator :
            raise NotImplementedError
        
    def setup(self, core) :     
        """  Set up some structures needed before solving.
        """     
        raise NotImplementedError 
    
    def evaluate(self) : 
        """ Evaluate the current core.
        """
        raise NotImplementedError   
    
        
    def display(self) :
        """  Introduce myself.
        """       
        raise NotImplementedError    
    
    # Concrete
    
    def plot_peaking(self) :
        """  Plot the power peaking factors.
        """
        plt.imshow(self.peaking, interpolation='nearest', cmap=plt.cm.hot)  
        plt.title('power peaking')
        plt.colorbar()
        plt.show() 