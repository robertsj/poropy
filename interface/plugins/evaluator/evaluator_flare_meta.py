'''
Created on Feb 17, 2012

@author: robertsj
'''
from poropy.coretools.flare import Flare

class EvaluatorFlareMeta(Flare):
    """  Defines metadata for Flare
    """
    
    def __init__(self, parent=None):
        control.ScriptControl.__init__(self, parent)

    def userInputs(self):
        
        # define the user inputs
        inputs =    (
                        ("mixing_parameter",    {'type': "freetext", 'label': "mixing parameter",
                                                 'default': 0.98}),
                        ("alpha_1",             {'type': "freetext", 'label': "alpha 1",
                                                 'default': 0.02}),      
                        ("alpha_2",             {'type': "freetext", 'label': "alpha 2",
                                                 'default': 0.01}),                                     
                        ("tolerance",           {'type': "freetext", 'label': "tolerance",
                                                 'default': 0.0001}),               
                    )
        return inputs

    def runActions(self):

        # define the buttons
        actions =   (
                       # Can't think of anything.
                    )
        return actions

    def title(self):
        return "Flare"

    def description(self):
        return """Flare is a simple two dimensional, effectively one group model.  
               """
 