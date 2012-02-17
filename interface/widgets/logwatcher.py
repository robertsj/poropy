from __future__ import division

from PyQt4.QtCore import *
from PyQt4.QtGui import *


class LogWatcher(QWidget):
    def __init__(self,logfile,parent=None):
        QWidget.__init__(self,parent)

        self.display = QTextEdit()
        self.display.setReadOnly(True)
        l = QHBoxLayout()
        l.addWidget(self.display)
        self.setLayout(l)

        timer = QTimer(self)
        QObject.connect(timer, SIGNAL("timeout()"), self.update_display)
        timer.start(500)

        self.fh = open(logfile,'r')
        self.lines = self.fh.readlines()
        for line in self.lines:
            self.display.append(line)

    def update_display(self):
        newlines = self.fh.readlines()
        newTextLines = newlines[len(newlines)-len(self.lines)-1:]
        for line in newTextLines:
            self.display.append(line)
        self.display.update()
        self.lines = newlines
