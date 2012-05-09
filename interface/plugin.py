try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
except:
    class QObject:
        def __init__(self,parent):
            pass
    def pyqtSignature(s):
        return lambda x: True

class Plugin(QObject):
    """Base control class for adding functionality to the poropy gui"""

    def __init__(self,model,parent=None):
        QObject.__init__(self,parent)
        self.cancel = False
        self.model = model

    def userInputs(self):
        """Must return user inputs in a specific form

        See the included plugins for examples of the options

        """
        raise Exception('Must implement the userInputs method')

    def runActions(self):
        """Must return run actions in a specific form

        See the included plugins for examples of the options

        """
        raise Exception('Must implement the runActions method')

    def title(self):
        return 'unnamed'

    def description(self):
        return 'No description specified'

    @pyqtSignature("progressChanged(int)")
    def progressChanged(self,p):
        """Call this method from run action methods to implement the progress bar
        
            The implementing widget in widget.py (ProgressBarButton) assumes progress goes from 0 to 100        

        """
        self.emit(SIGNAL("progressChanged(int)"),p)

    @pyqtSignature("cancel()")
    def cancel(self):
        """Slot for setting the cancel flag

            If you watch for the value of self.cancel in your run actions at each step, you can break
            out of a loop to cancel a long run.

        """
        self.cancel = True
        
    def close_plugin_diag(self):
        self.emit(SIGNAL("closePluginDiag()"))

