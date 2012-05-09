# poropy/coretools/optimizer.py  --  optimization tools.

import copy

from reactor import Reactor, Core, SpentFuelPool

try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    PYQT4 = True
    OBJECTBASE = QObject
except ImportError:
    PYQT4 = False
    OBJECTBASE = object

class Optimizer(OBJECTBASE):
    """  Abstract base class for optimization.
    """
    
    def __init__(self, argv, reactor, objective, parent=None):
        """ Constructor.  
        """
        if PYQT4:
            QObject.__init__(self, parent)
    
        # Set the reactor
        self.reactor = reactor
        # Set the objective function; must have signature f(k, p).
        self.objective = objective
        # Current evaluation of keff
        self.k = 1.0
        # Current maximum power peaking factor
        self.p = 0.0
        
    def set_fixed_central(self, value=True):
        """  Fix the central bundle.
        
        By construction, the central bundle is the  first entry in the 
        pattern. This flag will indicate crossover should skip this.
        """
        self.fixed_central = value      

        
    def fire_signal(self):
        """  Let mainWindow know a new evaluation is available.
        
        This *must* be called at the end of each inherited optimizer's
        end_of_iteration class.
        """
        # Have reactor update anything needed for displays.
        self.reactor.core.update_assembly_peaking()
        if PYQT4:
            self.emit(SIGNAL("patternUpdated"), self.k,
                                                self.p,
                                                self.objective(self.k,self.p),
                                                copy.copy(self.reactor.core.pattern))

    # Abstract Interface       
        
    def execute(self):
        """  Optimize the reactor for the objective.
        """     
        raise NotImplementedError   
        
    def end_of_iteration(self):
        """  Do something at the end of each iteration.
        
        Many optimization algorithms have relatively distinct 
        "iterations".  For genetic algorithms, this is the generation.
        For simulated annealing, this could be a temperature change.
        """
        raise NotImplementedError 

