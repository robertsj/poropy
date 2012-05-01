from __future__ import division

import os

from PyQt4.QtCore import *
from PyQt4.QtGui import *


class LogWatcher(QWidget):
    def __init__(self,logfile,parent=None):
        QWidget.__init__(self,parent)

        self.file = logfile

        self.lines = []

        self.display = QTextEdit()
        self.display.setReadOnly(True)
        l = QHBoxLayout()
        l.addWidget(self.display)
        self.setLayout(l)

        timer = QTimer(self)
        QObject.connect(timer, SIGNAL("timeout()"), self.update_display)
        timer.start(1000)


    def update_display(self):

        if not os.path.isfile(self.file): return

        with open(self.file,'r') as fh:
            newlines = fh.readlines()
            if len(newlines) > len(self.lines):
                for line in newlines[len(self.lines)-len(newlines):]:
                    self.display.append(line[:-1])
                self.display.update()
                self.lines = newlines

