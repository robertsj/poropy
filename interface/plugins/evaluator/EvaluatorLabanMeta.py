'''
Created on Feb 17, 2012

@author: robertsj
'''
from poropy.coretools.laban import Laban

import plugin

class EvaluatorLabanMeta(plugin.Plugin):
    """  Defines metadata for Laban-Pel
    """
    
    def __init__(self,model,parent=None):
      plugin.Plugin.__init__(self,model,parent=parent)
    

    def userInputs(self):
        
        # define the user inputs
        inputs =    (
                        ("order",               {'type': "freetext", 'label': "response order",
                                                 'default': 2}),
                        ("tolerance",           {'type': "freetext", 'label': "tolerance",
                                                 'default': 0.0001}),    
                        ("maximum_iterations",  {'type': "freetext", 'label': "maximum iterations",
                                                 'default': 0.0001}),                                    
                    )
        return inputs

    def runActions(self):

        # define the buttons
        actions =   (
                       ("set_evaluator",       {'label': "Choose",
                                                'description': "Set Laban as the Poropy evaluator."}),
                    )
        return actions

    def title(self):
        return "LABAN-PEL"

    def description(self):
        return """LABAN-PEL is a 2-D multigroup response code.
                  See, for example:
                  E. Z. Muller, "LABAN-PEL: A Two-Dimensional, Multigroup Diffusion, 
                      High-Order Response Matrix Code," PEL-309 (June 1991).
                  E. Z. Muller, Z. J. Weiss, "Benchmarking with the Multigroup Diffusion 
                      High-Order Response Matrix Method," Ann. nucl. Energy, 18, 535-544 (1991).
                      
                  Note, in quarter core geometry, LABAN-PEL can only employ reflective
                  boundary conditions.    
               """
 

    def set_evaluator(self):

      physics = Laban()
      
      # use your userIntputs here to customize the Laban() instance
      
      self.model.set_evaluator(physics)
      
      self.close_plugin_diag()
      
