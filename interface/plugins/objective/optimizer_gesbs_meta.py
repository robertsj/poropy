'''
Created on Feb 17, 2012

@author: robertsj
'''
from poropy.coretools.optimizer_gesbs import OptimizerGESBS

class OptimizerGESBSMeta(OptimizerGESBS):
    """  Defines metadata for OptimizerGESBS
    """
    
    def __init__(self, parent=None):
        control.ScriptControl.__init__(self, parent)

    def userInputs(self):
        
        # define the user inputs

        inputs =    (
                        ("maximum_iterations",  {'type': "freetext", 'label': "maximum_iterations",
                                                 'default': 123}),
                        ("seed",                {'type': "freetext", 'label': "seed",
                                                 'default': 123}),
                        ("np_seed",             {'type': "freetext", 'label': "np_seed",
                                                 'default': 123}),
                        ("binary_sweep",        {'type': "freetext", 'label': "np_seed",
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
        return "Greedy Exhaustive Single Binary Swap"

    def description(self):
        return """Optimize a loading pattern for a given objective using GESBS."""
    
    def execute(self):
        """  Execute the parent class. 
        """
        
        # Deal with the inputs.  For example,
        #   opt.maximum_iterations = inputs.maximum_iterations
        # or whatever.
        
        OptimizerGESBS.execute()  
    
