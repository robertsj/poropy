#!/usr/bin/python
# -*- coding: utf-8 -*-
## Add path to library (just for examples; you do not need this)
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from scipy import random
from PyQt4 import QtGui, QtCore
from pyqtgraph.PlotWidget import *
from pyqtgraph.graphicsItems import *
import pyqtgraph as pg

app = QtGui.QApplication([])

def mkPlot():
    global mw
    mw = QtGui.QMainWindow()
    pw = PlotWidget()
    mw.setCentralWidget(pw)
    mw.show()


def mkImages():
    global w, img
    w = pg.GraphicsWindow()
    img = []
    for i in range(10):
        im = pg.ImageItem(np.random.normal(size=(200,200)))
        im.setCacheMode(im.NoCache)
        img.append(im)
        
        
        if i%2 == 0:
            w.scene().addItem(img[-1])
            img[-1].translate(i*10, i*10)
        
mkImages()
