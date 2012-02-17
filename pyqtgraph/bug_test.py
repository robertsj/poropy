# -*- coding: utf-8 -*-
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.QtOpenGL import *

app = QApplication([])

## Create a graphicsview with GL viewport, show a pixmapitem
gv = QGraphicsView()
glw = QGLWidget()
gv.setViewport(glw)
s = QGraphicsScene()
gv.setScene(s)
gv.show()
data = "dhsuo987hifdpo2$fh834902rlohasf9afjhi34rhoerifaffjsiofp[w30a3fse"*100
img = QImage(data, 40, 40, QImage.Format_RGB32)
px = QPixmap(img)
pi = QGraphicsPixmapItem(px)
gv.scene().addItem(pi)

## create a second GL graphicsview with a proxywidget
#gv2 = QGraphicsView()
#glw2 = QGLWidget()
#gv2.setViewport(glw2)
#s2 = QGraphicsScene()
#gv2.setScene(s2)
#item = QGraphicsWidget()
#btn = QToolButton()
#btn.setText("?")
#proxy = QGraphicsProxyWidget(item)
#proxy.setWidget(btn)
#s2.addItem(item)
#gv2.show()
    
#app.processEvents()

### update the first view
#gv.translate(1,0)

### delete the second view; the pixmap in the first view goes blank.
#del gv2

