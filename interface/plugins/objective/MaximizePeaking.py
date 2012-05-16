'''
Created on Feb 17, 2012

@author: robertsj
'''

import plugin


class MaximizePeaking(plugin.Plugin):
    
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
        return "Maximize peaking"

    def description(self):
        return """Maximize peaking.

               Not recommended, but good for sanity checks!
               """
    def set_objective(self):
    
        self.model.set_objective(objective)
        
        self.close_plugin_diag()

def objective(k, p):
    """  Maximize peaking.
    """
    return p


