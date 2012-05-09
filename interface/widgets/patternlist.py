from __future__ import division

import copy


from PyQt4.QtCore import *
from PyQt4.QtGui import *

class PatternList(QWidget):
    def __init__(self,parent=None):
        QWidget.__init__(self,parent)

        headers = ["Timestamp","Maxpeak","keff","Objective","Pattern"]
        self.numCols = 4

        self.display = QTreeWidget()
        self.display.setSortingEnabled(True)
        self.display.setColumnCount(self.numCols)
        self.display.setHeaderLabels(headers)
        l = QHBoxLayout()
        l.addWidget(self.display)
        self.setLayout(l)

        self.connect(self.display,SIGNAL("itemSelectionChanged()"),self.change_pattern)
        #self.connect(self.display,SIGNAL("itemClicked(QTreeWidgetItem*,int)"),self.checkbox_clicked)


    #def checkbox_clicked(self,item,col):
    #    if col != 4: return
    #    if item.checkState(col) == Qt.Checked:
    #        self.emit(SIGNAL("saveItem(QTreeWidgetItem*)"),item)
    #    else:
    #        self.emit(SIGNAL("unsaveItem(QTreeWidgetItem*)"),item)
    #
    #def add_item(self,item):
    #    self.display.addTopLevelItem(item)

    #def remove_item(self,item):
    #    self.display.removeItemWidget(item,0)


    def add_pattern(self,i,timestamp,pattern,keff,maxpeak,objective):
        strlst = QStringList()
        strlst.append(str(timestamp))
        strlst.append(str(maxpeak))
        strlst.append(str(keff))
        strlst.append(str(objective))
        strlst.append(str(pattern))
        item = QTreeWidgetItem(self.display,strlst)
        #item.setCheckState(4,Qt.Unchecked)
        item.setData(0,32,i)

    def resize(self):
        for c in range(self.numCols):
            self.display.resizeColumnToContents(c)

    def change_pattern(self):
        try:
          data = self.display.selectedItems()[0].data(0,32)
          self.emit(SIGNAL("patternChanged(QVariant)"),data)
        except IndexError: pass # selection was changed to nothing
        
    def clear(self):
        self.display.clear()
