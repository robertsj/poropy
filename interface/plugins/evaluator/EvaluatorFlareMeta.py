'''
Created on Feb 17, 2012

@author: robertsj
'''
from poropy.coretools import Reactor, Flare

import plugin

class EvaluatorFlareMeta(plugin.Plugin):
    """  Defines metadata for Flare
    """
    
    def __init__(self,model,parent=None):
      plugin.Plugin.__init__(self,model,parent=parent)


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
                       ("set_evaluator",       {'label': "Choose",
                                                'description': "Set Flare as the Poropy evaluator."}),
                    )
        return actions

    def title(self):
        return "Flare"

    def description(self):
        return """Flare is a simple two dimensional, effectively one group model.  

               """
               
    def set_evaluator(self):
    

      physics = Flare()
      
      # use your userIntputs here to customize the Flare() instance
      
      self.model.set_evaluator(physics)
      
      self.close_plugin_diag()
      
      
