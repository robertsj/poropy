#!/usr/bin/env python

from __future__ import division

import sys
import os
import platform
import copy

import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import widgets

from plugincontrol import PluginControl
from model import Model


__version__ = "0.1.0"

class MainWindow(QMainWindow):

    def __init__(self, parallel=False, parent=None):
        super(MainWindow, self).__init__(parent)

        self.showMaximized()

        self.setWindowTitle("Poropy PWR Core Optimization Interface")

        self.model = Model()

        # Add menu items

        self.menubar = QMenuBar(self)

        self.menuFile = QMenu("&File",self.menubar)
        self.actionBuild = QAction("Build &Core",self)
        self.connect(self.actionBuild, SIGNAL("triggered()"), self.build_core)
        self.actionExit = QAction("E&xit",self)
        self.connect(self.actionExit, SIGNAL("triggered()"), self.close)
        self.menuFile.addActions([self.actionBuild,self.actionExit])
        self.menuFile.insertSeparator(self.actionExit)

        self.menuSettings = QMenu("&Settings",self.menubar)
        self.actionEval = QAction("&Evaluator Settings",self)
        self.connect(self.actionEval, SIGNAL("triggered()"), self.evaluator_settings)
        self.actionObjective = QAction("&Objective Settings",self)
        self.connect(self.actionObjective, SIGNAL("triggered()"), self.objective_settings)
        self.actionCore = QAction("&Core Display Settings",self)
        self.connect(self.actionCore, SIGNAL("triggered()"), self.core_display_settings)
        self.menuSettings.addActions([self.actionEval,self.actionObjective,self.actionCore])

        self.menuTools = QMenu("&Tools",self.menubar)
        self.actionOptimize = QAction("&Run Optimization",self)
        self.connect(self.actionOptimize, SIGNAL("triggered()"), self.run_optimization)
        self.menuTools.addActions([self.actionOptimize])

        self.menuHelp = QMenu("&Help",self.menubar)
        self.actionAbout = QAction("&About",self)
        self.connect(self.actionAbout, SIGNAL("triggered()"), self.about)
        self.menuHelp.addActions([self.actionAbout])

        self.menubar.addActions([self.menuFile.menuAction(),
                                 self.menuSettings.menuAction(),
                                 self.menuTools.menuAction(),
                                 self.menuHelp.menuAction()])
        self.setMenuBar(self.menubar)


        # Instantiate widgets

        self.coreDisplay = widgets.CoreDisplay(None)
        self.plotObjective = widgets.PlotWidget()
        self.plotKeff = widgets.PlotWidget()
        self.plotPeaking = widgets.PlotWidget()
        self.allPatterns = widgets.PatternList()
        self.savedPatterns = widgets.PatternList()
        self.logView = widgets.LogWatcher("output.log")

        # Setup widget layouts

        rightLayout = QHBoxLayout()
        tabsPlot = QTabWidget()
        tabsPlot.addTab(self.plotKeff, "keff")
        tabsPlot.addTab(self.plotPeaking, "Max Peaking")
        tabsPlot.addTab(self.plotObjective, "Objective")

        tabsPatt = QTabWidget()
        tabsPatt.addTab(self.allPatterns, "Previous Patterns")
#        tabsPatt.addTab(self.savedPatterns, "Saved Patterns")

        outerHorSplit = QSplitter()
        leftVertSplit = QSplitter()
        rightVertSplit = QSplitter()
        leftVertSplit.setOrientation(Qt.Vertical)
        rightVertSplit.setOrientation(Qt.Vertical)

        self.coreDisplay.setMinimumSize(400,400)
        leftVertSplit.addWidget(self.coreDisplay)
        leftVertSplit.addWidget(self.logView)

        rightVertSplit.addWidget(tabsPlot)
        rightVertSplit.addWidget(tabsPatt)

        outerHorSplit.addWidget(leftVertSplit)
        outerHorSplit.addWidget(rightVertSplit)
        self.setCentralWidget(outerHorSplit)

        # this file is the controller, which handles the interplay between the
        # model and the views

        # connections to the views

        self.connect(self.coreDisplay,SIGNAL("assemblySwapped"),self.assembly_swap)
        self.connect(self.allPatterns,SIGNAL("patternChanged(QVariant)"),self.change_pattern)

        # connections to the model

        self.connect(self.model,SIGNAL("reactorChanged()"),self.new_reactor)
        self.connect(self.model,SIGNAL("reactorEvaluated()"),self.new_evaluation)
        self.connect(self.model,SIGNAL("patternChanged()"),self.pattern_changed)

        
################################################################################
#  Menu Slots
################################################################################

    def build_core(self):
      pluginDir = os.path.join(os.path.join(sys.path[0],"plugins"),"reactor")
      form = PluginControl(self.model,pluginDir,self)
      form.set_help("""<html>
<head>
<title>       </title>
<style type="text/css">
<!--
h1	{text-align:center;
	font-family:Arial, Helvetica, Sans-Serif;
	}

p	{text-indent:20px;
	}
-->
</style>
</head>
<body bgcolor = "#ffffcc" text = "#000000">
<h1>Hello, World!</h1>

<p>All of these descriptions are fed to a QTextEdit
(http://qt-project.org/doc/qt-4.8/qtextedit.html#setText), which supports some
simple HTML and CSS.</p>

</body>
</html>""")
      form.show()
    
    def evaluator_settings(self):
      pluginDir = os.path.join(os.path.join(sys.path[0],"plugins"),"evaluator")
      form = PluginControl(self.model,pluginDir,self)
      form.set_help("""Evaluator Settings Help""")
      form.show()
      
    def objective_settings(self):
      pluginDir = os.path.join(os.path.join(sys.path[0],"plugins"),"objective")
      form = PluginControl(self.model,pluginDir,self)
      form.set_help("""Evaluator Settings Help""")
      form.show()
    
    def core_display_settings(self):
    
      choices = ["Burnup","Enrichment","Power Peaking","K-Infinity"]
      choice,ok = QInputDialog.getItem(self,"Core Display Coloring",
                                       "Choose how to display core",
                                       choices,0,False)
      if ok:
        if choice == "Burnup":
          self.coreDisplay.set_coloring(widgets.CoreDisplay.COLOR_BURNUP)
        elif choice == "Enrichment":
          self.coreDisplay.set_coloring(widgets.CoreDisplay.COLOR_ENRICHMENT)
        elif choice == "Power Peaking":
          self.coreDisplay.set_coloring(widgets.CoreDisplay.COLOR_POWER)
        elif choice == "K-Infinity":
          self.coreDisplay.set_coloring(widgets.CoreDisplay.COLOR_KINF)
        self.coreDisplay.pattern_updated()


    def run_optimization(self):
      pluginDir = os.path.join(os.path.join(sys.path[0],"plugins"),"optimizer")
      form = PluginControl(self.model,pluginDir,self)
      form.set_help("""Optimizer Help""")
      form.show()


################################################################################
#  View Slots - for views that need to interact with the model
################################################################################

    def assembly_swap(self,toFrom):
        """Slot called on manual swap from the core display widget"""
        self.model.swap_assemblies(toFrom[0],toFrom[1])
        self.model.evaluate_reactor()


    def change_pattern(self, index):
        """Slot for when a pattern list item is clicked"""
        
        if isinstance(index,QVariant):
            index = index.toPyObject()

        self.model.change_to_pattern(index)
        

################################################################################
#  Model Slots - for updating views when the model changes
################################################################################

    def new_reactor(self):
      """Do a full reset on the gui with a new starting core pattern"""
      
      self.plotObjective.clear()
      self.plotPeaking.clear()
      self.plotKeff.clear()
      self.allPatterns.clear()
      self.savedPatterns.clear()
      
      self.coreDisplay.set_core(self.model.get_core())
      self.model.evaluate_reactor()

      self.allPatterns.resize()
      self.savedPatterns.resize()
      

    def new_evaluation(self):
      """Update the pattern list, the plots, and any printouts with the new data"""

      i,latest = self.model.get_latest_pattern()
      
      self.allPatterns.add_pattern( i,
                                    latest['timestamp'],
                                    latest['pattern'],
                                    latest['keff'],
                                    latest['maxpeak'],
                                    latest['objective'])
      
      self.plotPeaking.add_point(latest['maxpeak'])
      self.plotKeff.add_point(latest['keff'])
      self.plotObjective.add_point(latest['objective'])


    def pattern_changed(self):
        """Slot called whenever the current pattern of the Reactor is changed"""
        
        self.coreDisplay.pattern_updated()
        i = self.model.get_current_pattern_index()

        self.plotKeff.move_selection(i)
        self.plotPeaking.move_selection(i)
        self.plotObjective.move_selection(i)

    def about(self):
        QMessageBox.about(self, "About Poropy PWR Core Optimization Interface",
                          """<b>Py-Image PWR Core Optimization Interface</b> v %s
                          <p>Copyright &copy; 2012 Nick Horelik, Jeremy Roberts, 
                          All Rights Reserved.
                          <p>Python %s -- Qt %s -- PyQt %s on %s""" %
                          (__version__, platform.python_version(),
                           QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()
