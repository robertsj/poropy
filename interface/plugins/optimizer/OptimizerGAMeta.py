'''
Created on Feb 17, 2012

@author: robertsj
'''

import sys

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from poropy.coretools.optimizer_ga import OptimizerGA

import plugin


class OptimizerGAMeta(plugin.Plugin):
    """  Defines metadata for OptimizerGA
    """
    
    def __init__(self,model,parent=None):
        plugin.Plugin.__init__(self,model,parent=parent)

    def userInputs(self):
        
        # define the user inputs
        self.maximum_generations    = 100   # Maximum generations
        self.population_size        = 50    # Population size
        self.number_replaced        = 40    # Number replaced each generation
        self.seed                   = 123   # PGA random number seed
        self.np_seed                = 123   # NumPy random number seed
        self.binary_sweep           = False # Perform one sweep of binary exchanges
        inputs =    (
                        ("maximum_iterations",  {'type': "freetext", 'label': "Maximum generations",
                                                 'default': "100"}),
                        ("population_size",     {'type': "freetext", 'label': "Population size",
                                                 'default': "50"}),
                        ("number_replaced",     {'type': "freetext", 'label': "Number replaced",
                                                 'default': "40"}),
                        ("seed",                {'type': "freetext", 'label': "PGAPack RNG seed",
                                                 'default': "123"}),
                        ("np_seed",             {'type': "freetext", 'label': "NumPy RNG seed",
                                                 'default': "123"}),
                        ("binary_sweep",        {'type': "checkbox", 'label': "Binary sweep hill climber",}),                     
                    )
        return inputs

    def runActions(self):

        # define the buttons

        actions =   (
                       ("execute",     {'label': "Execute",
                                        'description': "Optimize the problem."}),
                    )
        return actions

    def title(self):
        return "Genetic Algorithm"

    def description(self):
        return """Optimize a loading pattern for a given objective using the genetic algorithm from PGAPack."""
    
    def execute(self):
        """  Execute the parent class. 
        """

        self.progressChanged(100)
        self.close_plugin_diag()
        
        argv = sys.argv
        reactor = self.model.reactor
        objective = self.model.objective

        opt = OptimizerGA(argv,reactor,objective)

        opt.maximum_generations    = int(self.maximum_iterations)
        opt.population_size        = int(self.population_size)
        opt.number_replaced        = int(self.number_replaced)
        opt.seed                   = int(self.seed)
        opt.np_seed                = int(self.np_seed)
        opt.binary_sweep           = self.binary_sweep
        
        self.connect(opt,SIGNAL("patternUpdated"),self.model.opt_pattern_updated)

        opt.execute()
        
#        self.close_plugin_diag()
