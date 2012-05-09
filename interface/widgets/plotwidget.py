from __future__ import division

import sys
import datetime

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import pyqtgraph as pg


class PlotWidget(QWidget):
    def __init__(self,parent=None):
        QWidget.__init__(self,parent)
        self.x = []
        self.y = []
        self.pw = pg.PlotWidget()
        self.plot = self.pw.plot()
        self.pw.getPlotItem().setWindowFrameMargins(10,10,10,10)
        l = QHBoxLayout()
        l.addWidget(self.pw)
        self.setLayout(l)
        
        self.rect = None

    def add_point(self,y,x=None):
        
        self.y.append(y)
        if x:
            self.x.append(x)
            cur_x = x
        else:
            try:
                self.x.append(self.x[-1]+1)
                cur_x = self.x[-1]
            except IndexError:
                self.x.append(0)
                cur_x = 0

        self.refresh()
        
        self.pw.removeItem(self.rect)
        self.rect = self.pw.plot(x=[cur_x], y=[y], pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')
        

    def refresh(self):
        self.plot.setData(y=self.y,x=self.x)
        if len(self.x) == 1:
          self.pw.setXRange(0,1,padding=0)
        else:
          self.pw.autoRange()
          
    def clear(self):
        self.x = []
        self.y = []
        self.refresh()
        
    def move_selection(self,i):

        y = self.y[i]
        
        self.pw.removeItem(self.rect)
        self.rect = self.pw.plot(x=[i], y=[y], pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')
        
        self.refresh()

