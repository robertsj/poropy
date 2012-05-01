"""

Provides a container and access functions for the model, which is what actually
hold the problem parameters and carries out functions. All access to the model
should be done with access functions, not direct manipulation of attributes.
This ensures that whenever there is a change in the model in any way, the proper
signals are emmitted.

In the MVC framework

model - this file
views - the plugin control class and display widgets
controller - interface.pyw and the plugin files

controllers need only to interact with this model using the access functions,
and the views will update themselves appropriately

"""

import os
import sys
import warnings
import datetime
import copy

from PyQt4.QtCore import *
from PyQt4.QtGui import *

class Model(QWidget):
  """Class for all interactions with the Poropy model"""

  def __init__(self):
    QWidget.__init__(self)
    
    self.reactor = None
    self.evaluator = None # default will be set when the reactor is set
    self.objective = default_objective
    
    self.pattern_history = []
    self.current_pattern = 0

################################################################################
# Getters and Setters
################################################################################

  def set_reactor(self,reactor):
    self.reactor = reactor
    self.pattern_history = []
    self.emit(SIGNAL("reactorChanged()"))
  
  def get_reactor(self):
    return self.reactor
    
  def get_core(self):
    return self.reactor.core
    
  def set_evaluator(self,evaluator):

    if not self.reactor:
      warnings.warn("Can't set evaluator before setting up a core.")
      return
    
    self.evaluator = evaluator.display()
    print self.evaluator
    self.reactor.set_evaluator(evaluator)

    self.emit(SIGNAL("evaluatorChanged()"))
  
  def get_evaluator(self):
    return self.evaluator

  def set_objective(self, obj):
    self.objective = obj
    
  def get_objective(self):
    return self.objective

  def get_pattern_by_index(self,i):
    if i < 0 or i > len(self.pattern_history):
      raise Exception("Pattern history index out of range: {0}".format(i))
    return self.patten_history[i]
  
  def get_latest_pattern(self):
    if len(self.pattern_history) == 0:
      raise Exception("Trying to access latest pattern when no patterns saved")
    return len(self.pattern_history)-1,self.pattern_history[-1]
  
  def get_pattern_history(self):
    return self.pattern_history

  def get_current_pattern_index(self):
    return self.current_pattern

################################################################################
# Methods for interacting with the poropy Reactor model
################################################################################

  def swap_assemblies(self, pos1, pos2):
    self.reactor.swap(pos1,pos2)
    self.emit(SIGNAL("patternChanged()"))

  def evaluate_reactor(self):
    cwd = os.getcwd()
    os.chdir(os.path.join(sys.path[0],"tmpdir"))
    keff,maxpeak = self.reactor.evaluate()
    os.chdir(cwd)
    self.pattern_history.append({'timestamp': datetime.datetime.now(),
                                 'pattern': copy.copy(self.reactor.core.pattern),
                                 'keff': keff,
                                 'maxpeak': maxpeak,
                                 'objective': self.objective(keff,maxpeak)})
    self.current_pattern = len(self.pattern_history) - 1
    
    self.emit(SIGNAL("reactorEvaluated()"))

  def change_to_pattern(self, i):
    if i < 0 or i > len(self.pattern_history):
      raise Exception("Pattern history index out of range: {0}".format(i))

    self.reactor.shuffle(self.pattern_history[i]['pattern'])
    
    self.current_pattern = i
    
    self.emit(SIGNAL("patternChanged()"))
    
################################################################################
# Slots for the optimizers
################################################################################
    
  def opt_pattern_updated(self,k,p,obj,pattern):
    self.pattern_history.append({'timestamp': datetime.datetime.now(),
                                 'pattern': pattern,
                                 'keff': k,
                                 'maxpeak': p,
                                 'objective':obj})
    self.current_pattern = len(self.pattern_history) - 1
    
    self.emit(SIGNAL("reactorEvaluated()"))
    self.emit(SIGNAL("patternChanged()"))


def default_objective(k, p):
    """  Default objective function.
    """
    delta = 0
    if k < 1.1 :
        delta = k - 1.1
    return 1.0 * (1.5 - p) + 50.0 * delta

    
