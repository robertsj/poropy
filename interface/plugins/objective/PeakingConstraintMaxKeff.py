'''
Created on Feb 17, 2012

@author: robertsj
'''

import plugin


class PeakingConstraintMaxKeff(plugin.Plugin):
    
    def __init__(self,model,parent=None):
      plugin.Plugin.__init__(self,model,parent=parent)


    def userInputs(self):
        
        # define the user inputs
        inputs =    (
                                   # don't need any  
                    )
        return inputs

    def runActions(self):

        # define the buttons
        actions =   (
                       ("set_objective",       {'label': "Choose",
                                                'description': "Set the objective."}),
                    )
        return actions

    def title(self):
        return "PeakingConstraintMaxKeff"

    def description(self):
        return """Constrained Peaking and Maximized Keff

               """
    def set_objective(self):
    
        self.model.set_objective(objective)
        
        self.close_plugin_diag()

def objective(k, p):
    """  Default objective function.
    """
    delta = 0
    if p > 1.5 :
        delta = 1.5 - p
    return 1.0 * (k - 1.13) + 50.0 * delta
