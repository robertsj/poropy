'''
Created on Oct 31, 2011

@author: robertsj
'''
import abc


class Bundle(object):
    '''
    An abstract base class for representing bundles.
    '''
    __metaclass__ = abs.ABCMeta

    def __init__(self):
        '''
        Constructor
        '''
        