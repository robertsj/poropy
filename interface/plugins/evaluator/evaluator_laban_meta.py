'''
Created on Feb 17, 2012

@author: robertsj
'''
from poropy.coretools.laban import Laban

class EvaluatorLabanMeta(Laban):
    """  Defines metadata for Laban-Pel
    """
    
    def __init__(self, parent=None):
        control.ScriptControl.__init__(self, parent)

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
                       # Can't think of anything.
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
 