# -*- coding: utf-8 -*-

class Layout(QtGui.QGraphicsWidget):
    def __init__(self):
        QtGui.QGraphicsWidget.__init__(self)
        
    def add(self, item, row, col, rowspan=1, colspan=1):