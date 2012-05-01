from __future__ import division

import os
import sys
import warnings
import traceback

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import widgets

class PluginControl(QDialog):
    """Base control dialog for using plugins defined in plugin.py"""

    def __init__(self,model,pluginDir,parent=None):
        QDialog.__init__(self,parent)
      
        self.model = model
        self.pluginDir = pluginDir
        defaultHelp = """Define this help section with the set_help method of the
        PluginControl subclass.  This will be given to a QTextEdit (thus you 
        can use HTML, among other things)."""

        ### create the initial dialog

        self.setSizeGripEnabled(True)

        # create default widgets
        self.pluginList = QListWidget()
        self.pluginList.setMaximumWidth(200)
        self.pluginList.addItem("Help")
        desc = QGroupBox('Core Builder')
        desc.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Minimum)
        self.help = QTextEdit(defaultHelp)
        self.help.setReadOnly(True)
        dLayout = QVBoxLayout()
        dLayout.addWidget(self.help)
        desc.setLayout(dLayout)
        helpLayout = QVBoxLayout()
        helpLayout.addWidget(desc)
        defaultFrame = QFrame()
        defaultFrame.setLayout(helpLayout)
        sa = QScrollArea()
        sa.setMinimumSize(600,600)
        saLayout = QVBoxLayout()
        saLayout.addWidget(defaultFrame)
        sa.setLayout(saLayout)
        closeButton = QDialogButtonBox(QDialogButtonBox.Close)

        # create initial layout
        self.frameLayout = QStackedLayout()
        self.frameLayout.addWidget(sa)

        rightLayout = QVBoxLayout()
        rightLayout.addLayout(self.frameLayout)
        rightLayout.addWidget(closeButton)

        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.pluginList)
        mainLayout.addLayout(rightLayout)

        self.setLayout(mainLayout)

        self.connect(closeButton, SIGNAL("rejected()"), self.accept)
        self.connect(self.pluginList,SIGNAL("currentRowChanged(int)"),self.frameLayout,SLOT("setCurrentIndex(int)"))


        ### populate panels for each core plugin

        sys.path.append(self.pluginDir) # add this util directory to python path

        for root, dirs, files in os.walk(self.pluginDir):
            for f in files:
                if f.endswith(".py"):
                    mod_name,file_ext = os.path.splitext(f)
                    try:
                        py_mod = __import__(mod_name)
                    except Exception:
                        py_mod = None
                        traceback.print_exc()
                        warnings.warn("Unable to load module {0}".format(f))
                    if not mod_name in dir(py_mod):
                        warnings.warn("No control class named {0} found in {1}.  Is the filename exactly the same as the control class name?".format(mod_name,f))
                    else:
                        target = getattr(py_mod,mod_name)
                        try:
                            controller = target(self.model,parent=self.model) # where we instantiate the control class in the plugin
                            self.frameLayout.addWidget(self.get_control_frame(controller))
                            self.pluginList.addItem(controller.title())
                            # connect the controller to self to allow it to close the dialog
                            # since the parent of each controller is the model, this allows them
                            # to continue to do work in a threadsafe manner after closing the dialog,
                            # while retaining the progress bar
                            self.connect(controller,SIGNAL("closePluginDiag()"),self,SLOT("accept()"))
                        except Exception:
                            traceback.print_exc()
                            warnings.warn("Unable to process plugin in {0}".format(f))


    def set_help(self,helpString):
        self.help.setText(QString(helpString))

    def get_control_frame(self,controller):
        
        ## create frame

        frame = QFrame()
        frame.setMinimumSize(600,600)

        frame.controllerVals = {}
        frame.controller = controller

        inputs = controller.userInputs()
        actions = controller.runActions()
        description = controller.description()

        ## create description

        desc = QGroupBox('Description')
        dLayout = QVBoxLayout()
        dd = QTextEdit(description)
        dd.setReadOnly(True)
        dd.setMaximumHeight(75)
        dLayout.addWidget(dd)
        desc.setLayout(dLayout)

        ## create input fields

        ins = QGroupBox('Inputs')
        groupLayout = QVBoxLayout()
        for var,varDict in inputs:
            typ = varDict['type']
            itemLayout = QHBoxLayout()
            if typ == "openpath":
                itemLayout.addWidget(QLabel(varDict['label']))
                frame.controllerVals[var] = {}
                frame.controllerVals[var]['widget'] = widgets.FileChooser()
                frame.controllerVals[var]['ret'] = lambda widget: str(widget.getText())
                if 'default' in varDict:
                    frame.controllerVals[var]['widget'].setText(varDict['default'])
                itemLayout.addWidget(frame.controllerVals[var]['widget'])
            elif typ == "savepath":
                itemLayout.addWidget(QLabel(varDict['label']))
                frame.controllerVals[var] = {}
                frame.controllerVals[var]['widget'] = widgets.FileChooser(save=True)
                frame.controllerVals[var]['ret'] = lambda widget: str(widget.getText())
                if 'default' in varDict:
                    frame.controllerVals[var]['widget'].setText(varDict['default'])
                itemLayout.addWidget(frame.controllerVals[var]['widget'])
            elif typ == "dir":
                itemLayout.addWidget(QLabel(varDict['label']))
                frame.controllerVals[var] = {}
                frame.controllerVals[var]['widget'] = widgets.FileChooser(directory=True)
                frame.controllerVals[var]['ret'] = lambda widget: str(widget.getText())
                if 'default' in varDict:
                    frame.controllerVals[var]['widget'].setText(varDict['default'])
                itemLayout.addWidget(frame.controllerVals[var]['widget'])
            elif typ == "checkbox":
                itemLayout.addWidget(QLabel(varDict['label']))
                frame.controllerVals[var] = {}
                frame.controllerVals[var]['widget'] = QCheckBox()
                frame.controllerVals[var]['ret'] = lambda widget: widget.isChecked()
                itemLayout.addWidget(frame.controllerVals[var]['widget'])
            elif typ == "freetext":
                itemLayout.addWidget(QLabel(varDict['label']))
                frame.controllerVals[var] = {}
                frame.controllerVals[var]['widget'] = QLineEdit()
                frame.controllerVals[var]['ret'] = lambda widget: str(widget.text())
                if 'default' in varDict:
                    frame.controllerVals[var]['widget'].setText(str(varDict['default']))
                itemLayout.addWidget(frame.controllerVals[var]['widget'])
            elif typ == "linktext":
                itemLayout.addWidget(QLabel(varDict['label']))
                frame.controllerVals[var] = {}
                frame.controllerVals[var]['widget'] = widgets.LinkText(var)
                frame.controllerVals[var]['ret'] = lambda widget: str(widget.text())
                if 'default' in varDict:
                    frame.controllerVals[var]['widget'].setText(varDict['default'])
                itemLayout.addWidget(frame.controllerVals[var]['widget'])
            else:
                warnings.warn("Skipping unknown widget type: {0}".format(var))
            groupLayout.addLayout(itemLayout)
        ins.setLayout(groupLayout)

        ## create action buttons

        acts = QGroupBox('Actions')
        groupLayout = QVBoxLayout()
        for func,funcDict in actions:
            itemLayout = QHBoxLayout()
            btn = widgets.ProgressBarButton(funcDict['label'],controller,func,parent=self.model)
            funcDesc = QLabel(funcDict['description'])
            self.connect(btn,SIGNAL("pressed()"),self.resolve_frame)
            self.connect(btn,SIGNAL("released()"),btn.progressBarRun)
            itemLayout.addWidget(btn)
            itemLayout.addWidget(funcDesc)
            groupLayout.addLayout(itemLayout)    
        acts.setLayout(groupLayout)

        ## set frame layouts
        inFrame = QFrame()
        inFrame.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed)
        inLayout = QVBoxLayout()
        inLayout.addWidget(desc)
        inLayout.addWidget(ins)
        inLayout.addWidget(acts)
        inFrame.setLayout(inLayout)
        sa = QScrollArea()
        sa.setWidget(inFrame)
        sa.setWidgetResizable(True)
        outlayout = QVBoxLayout()
        outlayout.addWidget(sa)
        frame.setLayout(outlayout)

        return frame


    def resolve_frame(self):
        """Slot for all frame action buttons

            Sends input values to current frame's controller

        """

        frame = self.frameLayout.currentWidget()
        for var,varDict in frame.controllerVals.iteritems():
            setattr(frame.controller,var,varDict['ret'](varDict['widget']))
