'''
Created on Feb 17, 2012

@author: robertsj
'''

class SmallCoreMeta(object):
    """  Defines metadata for OptimizerGA
    """
    
    def __init__(self, parent=None):
        control.ScriptControl.__init__(self, parent)

    def userInputs(self):
        
        # define the user inputs

        inputs =    (
#                        ("maximum_iterations",  {'type': "freetext", 'label': "maximum_iterations",
#                                                 'default': 123}),
#                        ("population_size",     {'type': "freetext", 'label': "population_size",
#                                                 'default': 123}),
#                        ("number_replaced",     {'type': "freetext", 'label': "number_replaced",
#                                                 'default': 123}),
#                        ("seed",                {'type': "freetext", 'label': "seed",
#                                                 'default': 123}),
#                        ("np_seed",             {'type': "freetext", 'label': "np_seed",
#                                                 'default': 123}),
#                        ("binary_sweep",        {'type': "freetext", 'label': "np_seed",
#                                                 'default': false}),                     
                    )
        return inputs

    def runActions(self):

        # define the buttons

        actions =   (
                       ("build",       {'label': "Build",
                                        'description': "Build the core."}),
                    )
        return actions

    def title(self):
        return "Small Core Problem"

    def description(self):
        return """69 assembly core."""
    
    def build(self):
        """  Execute the parent class. 
        """
        