#!/usr/bin/env python

import sys, os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from dialogs.mainWindow import MainWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()
