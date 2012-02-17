'''
Created on Feb 17, 2012

@author: robertsj
'''
from poropy.coretools.optimizer_ga import OptimizerGA

class OptimizerGAMeta(OptimizerGA):
    """  Defines metadata for OptimizerGA
    """
    
    def __init__(self, parent=None):
        control.ScriptControl.__init__(self, parent)

    def userInputs(self):
        
        # define the user inputs

        inputs =    (
                        ("maximum_iterations",  {'type': "freetext", 'label': "Maximum generations",
                                                 'default': 123}),
                        ("population_size",     {'type': "freetext", 'label': "Population size",
                                                 'default': 123}),
                        ("number_replaced",     {'type': "freetext", 'label': "Number replaced",
                                                 'default': 123}),
                        ("seed",                {'type': "freetext", 'label': "PGAPack RNG seed",
                                                 'default': 123}),
                        ("np_seed",             {'type': "freetext", 'label': "NumPy RNG seed",
                                                 'default': 123}),
                        ("binary_sweep",        {'type': "freetext", 'label': "Binary sweep hill climber",
                                                 'default': false}),                     
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
        
        # Deal with the inputs.  For example,
        #   opt.maximum_iterations = inputs.maximum_iterations
        # or whatever.
        
        OptimizerGA.execute()  
    
