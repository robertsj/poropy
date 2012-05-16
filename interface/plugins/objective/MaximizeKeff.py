'''
Created on Feb 17, 2012

@author: robertsj
'''

import plugin


class MaximizeKeff(plugin.Plugin):
    
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
        return "Maximize keff"

    def description(self):
        return """Maximize keff.

               Sort of like maximizing energy (not quite, though).
               """
    def set_objective(self):
    
        self.model.set_objective(objective)
        
        self.close_plugin_diag()

def objective(k, p):
    """  Default objective function.
    """
    return k


