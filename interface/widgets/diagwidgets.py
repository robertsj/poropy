"""

Provides custom QWidgets for use in automatically drawn plugin dialogs

"""

from __future__ import division

import os

from PyQt4.QtCore import *
from PyQt4.QtGui import *


class LinkText(QWidget):
    def __init__(self,name,parent=None):
        QWidget.__init__(self,parent)

        self.name = QString(name)
        self.line = QLineEdit()
    
        lay = QHBoxLayout()
        lay.addWidget(self.line)

        self.setLayout(lay)

        self.connect(self.line,SIGNAL("textChanged(QString)"),self.textChanged)
        self.connect(self.line,SIGNAL("textChanged(QString)"),self.textChanged2)

    def text(self):
        return self.line.text()

    @pyqtSignature("textChanged(QString,QString)")
    def textChanged(self,text):
        self.emit(SIGNAL("textChanged(QString,QString)"),self.name,text)

    @pyqtSignature("textChanged()")
    def textChanged2(self,text):
        self.emit(SIGNAL("textChanged()"))

    def setText(self,txt):
        self.line.setText(str(txt))


class FileChooser(QWidget):
    def __init__(self,parent=None,directory=False,save=False):
        QWidget.__init__(self,parent)

        self.dir = directory
        self.save = save

        self.line = QLineEdit()
        self.btn = QToolButton()

        layout = QHBoxLayout(self)
        layout.addWidget(self.line)
        layout.addWidget(self.btn)

        self.connect(self.btn,SIGNAL("clicked()"),self.chooser)
        self.connect(self.line,SIGNAL("textChanged(QString)"),self.pathChanged)

    def chooser(self):
        if not self.dir and not self.save:
            path = QString(QFileDialog.getOpenFileName(None, "Select file", "*"))
        elif self.save:
            path = QString(QFileDialog.getSaveFileName(None, "Select file", "*"))
        else:
            path = QString(QFileDialog.getExistingDirectory(None, "Select directory"))
        self.line.setText(path)

    @pyqtSignature("pathChanged(QString)")
    def pathChanged(self):
        self.emit(SIGNAL("pathChanged(QString)"),self.line.text())

    def getText(self):
        return self.line.text()

    def setText(self,txt):
        self.line.setText(str(txt))



class ProgressBarButton(QPushButton):
    def __init__(self,label,controller,func, parent=None):
        QPushButton.__init__(self, label, parent=parent)

        self.controller = controller
        self.func = func
        


    def progressBarRun(self):

        self.controller.cancel = False

        worker = Worker(self,self.controller,self.func)

        self.progress = QProgressDialog("Running task...","Cancel",0,100,self)
        self.progress.setWindowTitle("Processing")
        self.progress.setMinimumDuration(0)
        self.progress.setWindowModality(Qt.WindowModal)

        self.connect(self.controller,SIGNAL("progressChanged(int)"),self.updateProgress)
        self.connect(self.progress,SIGNAL("canceled()"),self.controller,SLOT("cancel()"))
        self.connect(worker,SIGNAL("error()"),self.error)
        self.connect(worker,SIGNAL("finished()"),self.progress,SIGNAL("cancel()"))

        self.progress.setValue(0)
        worker.start()

    def updateProgress(self,int):
        self.progress.setValue(int)   

    def error(self):
        QMessageBox.critical(self,"Error","An Exception was encountered when processing the plugin runAction.  Please check the terminal for the Python traceback.")
        self.progress.cancel()



class Worker(QThread):

    def __init__(self,parent,controller,func):
        '''Parent must be defined to stop early garbage collection'''
        QThread.__init__(self,parent)
        self.controller = controller
        self.func = func
        
    def run(self):
        try:
            getattr(self.controller,self.func)()  # where the run method in the plugin is called
            self.controller.progressChanged(100)
        except:
            self.emit(SIGNAL("error()"))
            raise
