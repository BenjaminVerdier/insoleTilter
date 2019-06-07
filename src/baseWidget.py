import copy

#3rd party libs
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import trimesh as tm
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl

#custom files
from debugThingies import *

class baseWidget(QWidget):

    def __init__(self, nextTab = -1):
        super(baseWidget, self).__init__()
        self.nextTab = nextTab
        self.initBaseUI()

    def initBaseUI(self):
        #main widget definition
        self.mainLayout = QVBoxLayout()
        self.setLayout(self.mainLayout)

        #3D view
        self.view = gl.GLViewWidget()
        zgrid = gl.GLGridItem()
        zgrid.scale(20, 20, 20)
        axis = gl.GLAxisItem()
        axis.setSize(x=2, y=1, z=1)
        axis.scale(100,100,100)
        self.view.addItem(zgrid)
        self.view.addItem(axis)
        self.view.setCameraPosition(distance=288,elevation=32,azimuth=324)
        self.meshItem = None
        self.recompute = True

        self.mainLayout.addWidget(self.view)

    def initCustomUI(self):
        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)

        self.initDoneButtonAndShader()

    def initDoneButtonAndShader(self, vbxlayout):
        if self.nextTab > 0:
            doneButton = QPushButton("Done")
            doneButton.clicked.connect(self.doneBtnClicked)
            vbxlayout.addWidget(doneButton)
        #shader choice
        shaderLabel = QLabel("Choose display shader:")
        vbxlayout.addWidget(shaderLabel)
        shaderCB = QComboBox(self)
        shaderCB.addItem("balloon")
        shaderCB.addItem("normalColor")
        shaderCB.addItem("viewNormalColor")
        shaderCB.addItem("shaded")
        shaderCB.setCurrentIndex(1)
        shaderCB.activated.connect(self.shaderSelect)
        vbxlayout.addWidget(shaderCB)

        self.shader = "normalColor"
        self.glOptions = "opaque"

    def reshade(self):
        if not self.meshItem is None:
            self.meshItem.setShader(self.shader)
            self.meshItem.setGLOptions(self.glOptions)
        pass

    def start(self, transmittedData):
        pass

    def resetNextTabs(self):
        if self.nextTab > 0:
            self.parent().parent().widget(self.nextTab).resetNextTabs()
            self.parent().parent().setTabEnabled(self.nextTab, False)


    @pyqtSlot()
    def doneBtnClicked(self):
        self.parent().parent().setTabEnabled(self.nextTab, True)
        self.parent().parent().widget(self.nextTab).start(self.toTransmit)
        self.parent().parent().setCurrentIndex(self.nextTab)

    @pyqtSlot()
    def shaderSelect(self):
        cb = self.sender()
        self.shader = cb.currentText()
        if self.shader == "balloon":
            self.glOptions = "additive"
        else:
            self.glOptions = "opaque"
        self.reshade()
