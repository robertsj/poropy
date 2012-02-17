from __future__ import division

from PyQt4.QtCore import *
from PyQt4.QtGui import *

L = 100 # default pixel width of assembly displays

class CoreDisplay(QGraphicsView):
    def __init__(self,core,parent=None):
        QGraphicsView.__init__(self,parent)
        
        self.core = core
        self.needsRefresh = True

        self.mutex = QMutex()

        self.scene = QGraphicsScene(self)
        self.scene.maxZ = 0
        self.setScene(self.scene)
        self.connect(self.scene, SIGNAL("selectionChanged()"), self.refresh)

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.draw_core()

        timer = QTimer(self)
        QObject.connect(timer, SIGNAL("timeout()"), self.draw_core)
        timer.start(50)




    def pattern_updated(self):
        self.needsRefresh = True

    def resizeEvent(self, event):
        self.fitInView(self.scene.sceneRect(),Qt.KeepAspectRatio)

    def refresh(self):
        try:
            for item in self.scene.items():
                if isinstance(item,AssemblyDisplay):
                    item.refresh()
        except RuntimeError: pass


    def draw_core(self):
        if not self.needsRefresh: return
        
        self.mutex.lock()

        self.needsRefresh = False

        

        self.scene.clear()

        # map positions onto (row,col) for easier plotting
        mapped_pattern_ids = {}
        i = 0
        for r,row in enumerate(self.core.stencil):
            for c,assType in enumerate(row):
                if r != 0 and c == 0: # if on reflected surface
                    mapped_pattern_ids[r,c] = mapped_pattern_ids[0,r]
                else:
                    if assType != 0:
                        mapped_pattern_ids[r,c] = self.core.pattern[i]
                        i += 1
                    else:
                        mapped_pattern_ids[r,c] = -1 # reflector

        numRows = len(self.core.stencil)
        numCols = len(self.core.stencil[0])

        i = 0
        for r,row in enumerate(self.core.stencil):
            for c,assType in enumerate(row):
                patt_id = mapped_pattern_ids[r,c]
                if r != 0 and c == 0:    # if on reflected surface
                    a = self.core.assemblies[patt_id]
                    if patt_id == -1:
                        a.name = "R"
                        type_ = "reflector"
                    else:
                        a.name = patt_id
                        type_ = "reflected surface"
                    ass = AssemblyDisplay(r,c,a,type_,self,scene=self.scene)

                else:

                    if patt_id == -1:
                        a = self.core.reflector
                        a.name = "R"
                        ass = AssemblyDisplay(r,c,a,"reflector",self,scene=self.scene)
                    else:
                        a = self.core.assemblies[patt_id]
                        a.name = patt_id
                        ass = AssemblyDisplay(r,c,a,"fuel",self,scene=self.scene)
                        i += 1
                self.connect(ass.sigFire,SIGNAL("assemblySwap"),self.assembly_swap)
                

        m = 20 # pixel margin
        self.scene.setSceneRect(-m,-m,numCols*L+2*m,numRows*L+2*m)

        self.mutex.unlock()


    def assembly_swap(self,toFrom):
        self.emit(SIGNAL("assemblySwap"),toFrom)



class AssemblyDisplay(QGraphicsItem):
    def __init__(self,r,c,ass,type_,view,parent=None,scene=None):
        QGraphicsItem.__init__(self,parent,scene)

        self.loc = [r,c]
        self.assembly = ass
        self.type = type_
        self.view = view

        self.sigFire = SignalFire()

        bc = 50.0
        if self.type == "fuel":
            b = self.assembly.burnup
            if b <= bc: self.defaultColor = QColor(b/bc*255,(bc-b)/bc*255,0)
            else: self.defaultColor = Qt.lightGray
            self.setFlags(QGraphicsItem.ItemIsSelectable)
            self.setAcceptHoverEvents(True)
            self.setAcceptDrops(True)
        elif self.type == "reflector":
            self.defaultColor = Qt.cyan
        elif self.type == "reflected surface":
            b = self.assembly.burnup
            if b <= bc: self.defaultColor = QColor(b/bc*255,(bc-b)/bc*255,0)
            else: self.defaultColor = Qt.lightGray
        else:
            self.defaultColor = Qt.white

        self.draw_item()


    def draw_item(self):

        self.color = self.defaultColor
        self.path = QPainterPath()
        points = [(0,0),(0,L),(L,L),(L,0)]
        self.path.addPolygon(QPolygonF([QPointF(x,y) for x,y in points]))
        self.path.closeSubpath()
        self.make_labels()
        self.translate(self.loc[1]*L,self.loc[0]*L)

        s = self.get_scale()
        s = 1
        self.pixmap = QPixmap(L*s,L*s)
        self.pixpainter = QPainter(self.pixmap)
        rect = self.boundingRect()
        rect.translate(self.loc[1]*L,self.loc[0]*L)
        self.scene().render(self.pixpainter,QRectF(0,0,L*s,L*s),rect)


    def make_labels(self):
        id_ = QGraphicsTextItem(parent=self)

        f1 = self.assembly.name
        if self.type == "fuel":
            f2 = self.assembly.model
            f3 = self.assembly.burnup
        elif self.type == "reflector":
            f2 = ""
            f3 = ""
        elif self.type == "reflected surface":
            f2 = self.assembly.model
            f3 = self.assembly.burnup
        else:
            f2 = ""
            f3 = ""

        text = ""
        text += "<center>{0}<\center>".format(f1)
        text += "<center>{0}<\center>".format(f2)
        text += "<center>{0}<\center>".format(f3)
        id_.setHtml(text)
        id_.setTextWidth(L)


    def get_scale(self):
        sy = self.view.height()/self.scene().height()
        sx =  self.view.width()/self.scene().width()
        return min([sx,sy])

    def boundingRect(self):
        return self.path.boundingRect()

    def shape(self):
        return self.path

    def paint(self, painter, option, widget=None):
        painter.setPen(Qt.black)
        painter.setBrush(QBrush(self.color))
        painter.drawPath(self.path)

    def refresh(self):
        if self.isSelected():
            self.color = Qt.red
        else:
            self.color = self.defaultColor
        self.update()

    def hoverLeaveEvent(self, event):
        self.refresh()

    def hoverEnterEvent(self, event):
        self.color = Qt.gray
        self.update()

    def dragEnterEvent(self, event):

        if event.mimeData().assembly == self: return
        self.color = Qt.magenta
        self.update()

    def dragLeaveEvent(self, event):
        self.color = self.defaultColor
        self.update()

    def dropEvent(self, event):
        from_ = event.mimeData().assembly.loc
        to = self.loc
        if from_ != to:
            self.sigFire.fireSwap(to,from_)

    def mousePressEvent(self, event):
        if not self.type == "fuel": return

        for item in self.scene().selectedItems():
            item.setSelected(False)
        self.setSelected(True)

    def mouseMoveEvent(self, event):
        if not self.type == "fuel": return

        drag = QDrag(event.widget())
        mime = QMimeData()
        mime.assembly = self
        drag.setMimeData(mime)
        s = self.get_scale()
        drag.setPixmap(self.pixmap.scaledToWidth(L*s))
        drag.setHotSpot(QPoint(event.pos().x()*s,event.pos().y()*s))
        drag.start()


    def mouseReleaseEvent(self, event):
        if not self.type == "fuel": return


class SignalFire(QObject):
    def __init__(self,parent=None):
        QObject.__init__(self,parent)
    def fireSwap(self,to,from_):
        self.emit(SIGNAL("assemblySwap"),[from_,to])
