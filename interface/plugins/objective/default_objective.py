'''
Created on Feb 17, 2012

@author: robertsj
'''

def objective(self, k, p):
    """  Default objective function.
    """
    delta = 0
    if k < 1.1 :
        delta = k - 1.1
    return 1.0 * (1.5 - p) + 50.0 * delta