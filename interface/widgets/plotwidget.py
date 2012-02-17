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
        pw = pg.PlotWidget()
        self.plot = pw.plot()
        l = QHBoxLayout()
        l.addWidget(pw)
        self.setLayout(l)

    def add_point(self,y,x=None):
        
        self.y.append(y)
        if x:
            self.x.append(x)
        else:
            try:
                self.x.append(self.x[-1]+1)
            except IndexError:
                self.x.append(0)

        self.refresh()

    def refresh(self):
        self.plot.updateData(self.y,self.x)
