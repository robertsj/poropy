from __future__ import division

import sys
import os
import platform

import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import widgets
from driver import PGAOpt
import small_core
import large_core



__version__ = "0.0.1"

class MainWindow(QMainWindow):

    def __init__(self, parallel=False, parent=None):
        super(MainWindow, self).__init__(parent)

        self.showMaximized()

        self.setWindowTitle("Poropy PWR Core Optimization Interface")

        #defaultReactor = small_core.make_small_core(0)
        defaultReactor = large_core.make_large_core(0)
        self.reactor = defaultReactor


        # Add menu items

        self.menubar = QMenuBar(self)

        self.menuFile = QMenu("&File",self.menubar)
        self.actionNewReactor = QAction("&New Reactor",self)
        self.connect(self.actionNewReactor, SIGNAL("triggered()"), self.new_reactor)
        self.actionLoadReactor = QAction("&Load Reactor",self)
        self.connect(self.actionLoadReactor, SIGNAL("triggered()"), self.load_reactor)
        self.actionExit = QAction("E&xit",self)
        self.connect(self.actionExit, SIGNAL("triggered()"), self.close)
        self.menuFile.addActions([self.actionNewReactor,
                                  self.actionLoadReactor,
                                  self.actionExit])
        self.menuFile.insertSeparator(self.actionExit)



        self.menuEdit = QMenu("&Edit",self.menubar)
        self.actionEditCoreLayout = QAction("Core L&ayout",self)
        self.connect(self.actionEditCoreLayout, SIGNAL("triggered()"), self.edit_core)
        self.actionEditAssemblies = QAction("&Assemblies",self)
        self.connect(self.actionEditAssemblies, SIGNAL("triggered()"), self.edit_assemblies)
        self.menuEdit.addActions([self.actionEditCoreLayout,
                                  self.actionEditAssemblies])

        self.menuTools = QMenu("&Tools",self.menubar)
        self.actionRunOptimization = QAction("&Run Optimization",self)
        self.connect(self.actionRunOptimization, SIGNAL("triggered()"), self.run_optimization)
        self.menuTools.addActions([self.actionRunOptimization])

        self.menuHelp = QMenu("&Help",self.menubar)
        self.actionAbout = QAction("&About",self)
        self.connect(self.actionAbout, SIGNAL("triggered()"), self.about)
        self.menuHelp.addActions([self.actionAbout])

        self.menubar.addActions([self.menuFile.menuAction(),
                                 self.menuEdit.menuAction(),
                                 self.menuTools.menuAction(),
                                 self.menuHelp.menuAction()])
        self.setMenuBar(self.menubar)


        # Instantiate widgets

        self.coreDisplay = widgets.CoreDisplay(defaultReactor.core)
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
        tabsPatt.addTab(self.savedPatterns, "Saved Patterns")

        outerHorSplit = QSplitter()
        leftVertSplit = QSplitter()
        rightVertSplit = QSplitter()
        leftVertSplit.setOrientation(Qt.Vertical)
        rightVertSplit.setOrientation(Qt.Vertical)

        leftVertSplit.addWidget(self.coreDisplay)
        leftVertSplit.addWidget(self.logView)

        rightVertSplit.addWidget(tabsPlot)
        rightVertSplit.addWidget(tabsPatt)

        outerHorSplit.addWidget(leftVertSplit)
        outerHorSplit.addWidget(rightVertSplit)
        self.setCentralWidget(outerHorSplit)

        # connections

        self.connect(self.coreDisplay,SIGNAL("assemblySwap"),self.assembly_swap)

        #self.connect(self.optimizer,SIGNAL("reactorEvaluated(float,float)"),self.reactor_evaluated)
        #self.connect(self.optimizer,SIGNAL("patternUpdated()"),self.coreDisplay.pattern_updated)

        self.connect(self.allPatterns,SIGNAL("patternChanged(QVariant)"),self.change_pattern)
        self.connect(self.savedPatterns,SIGNAL("patternChanged(QVariant)"),self.change_pattern)
        #self.connect(self.allPatterns,SIGNAL("saveItem(QTreeWidgetItem*)"),self.savedPatterns.add_item)
        #self.connect(self.allPatterns,SIGNAL("unsaveItem(QTreeWidgetItem*)"),self.savedPatterns.remove_item)

        keff,maxpeak = self.evaluate_reactor()
        self.allPatterns.add_pattern(self.reactor.core.pattern,keff,maxpeak)
        self.savedPatterns.add_pattern(self.reactor.core.pattern,keff,maxpeak)
        self.allPatterns.resize()
        self.savedPatterns.resize()


    def change_pattern(self, pattern):
        """slot for when you click on pattern list item"""
    
        if isinstance(pattern,QVariant):
            pattern = pattern.toPyObject()
        if str(pattern) != str(self.reactor.core.pattern):
            self.reactor.shuffle(pattern)
            keff,maxpeak = self.evaluate_reactor()
            self.update_plots(keff,maxpeak)
            self.coreDisplay.pattern_updated()


    def evaluate_reactor(self):
        cwd = os.getcwd()
        os.chdir(os.path.join(sys.path[0],"tmpdir"))
        keff,maxpeak = self.reactor.evaluate()
        os.chdir(cwd)
        self.update_plots(keff,maxpeak)
        return keff,maxpeak


    def reactor_evaluated(self,keff,maxpeak):
        """Slot for optimizers to call periodically to update plots and add the pattern"""
        self.update_plots(keff,maxpeak)
        self.allPatterns.add_pattern(self.reactor.core.pattern,keff,maxpeak)


    def assembly_swap(self,toFrom):
        """called on manual swap"""
        self.reactor.swap(toFrom[0],toFrom[1])
        keff,maxpeak = self.evaluate_reactor()
        self.coreDisplay.pattern_updated()
        self.allPatterns.add_pattern(self.reactor.core.pattern,keff,maxpeak)


    def update_plots(self,keff,maxpeak):
        self.plotKeff.add_point(keff)
        self.plotPeaking.add_point(maxpeak)


    def new_reactor(self):
        pass


    def load_reactor(self):
        pass


    def edit_core(self):
        pass


    def edit_assemblies(self):
        pass


    def run_optimization(self):

        cwd = os.getcwd()

        opt = PGAOpt(self.reactor,self)

        progress = QProgressDialog("Running task...","Cancel",0,100,self)
        progress.setWindowTitle("Processing")
        progress.setMinimumDuration(0)
        progress.setWindowModality(Qt.WindowModal)

        self.connect(opt,SIGNAL("progressChanged(int)"),progress,SLOT("setValue(int)"))
        # emit the "progressChanged(int)" signal from opt to have the bar actually track progress
        self.connect(progress,SIGNAL("canceled()"),opt,SLOT("terminate()"))
        # terminate is dangerous (thread may have been modifying data, it can't unlock mutexes, etc)
        # better way is to have the underlying loop check a flag before each iteration, and then 
        # here we connect to a method that flips that flag, allowing the loop to terminate gracefully 
        self.connect(opt,SIGNAL("finished()"),progress,SLOT("cancel"))
        progress.setValue(0)

        opt.start()
        os.chdir(cwd)


    def about(self):
        QMessageBox.about(self, "About Poropy PWR Core Optimization Interface",
                          """<b>Py-Image PWR Core Optimization Interface</b> v %s
                          <p>Copyright &copy; 2012 Nick Horelik, Jeremey Roberts, 
                          All Rights Reserved.
                          <p>Python %s -- Qt %s -- PyQt %s on %s""" %
                          (__version__, platform.python_version(),
                           QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))


