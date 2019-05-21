#default libs
import sys
import time
from enum import Enum, auto
import math
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

class Rotation(Enum):
    FULL = auto()
    STRETCH = auto()

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()

    def initUI(self):
        mw = MainWidget()
        self.setCentralWidget(mw)
        self.setWindowTitle('Insole Tilter')
        self.resize(1300,1200)

        self.show()

class MainWidget(QWidget):

    def __init__(self):
        super(MainWidget, self).__init__()
        self.initUI()

    def initUI(self):
        #main widget definition
        mainLayout = QVBoxLayout()
        self.setLayout(mainLayout)

        #3D view
        self.view = gl.GLViewWidget()
        zgrid = gl.GLGridItem()
        zgrid.scale(20, 20, 20)
        self.view.addItem(zgrid)
        self.meshItem = gl.GLMeshItem()
        self.view.addItem(self.meshItem)
        self.view.setCameraPosition(distance=288,elevation=32,azimuth=324)

        mainLayout.addWidget(self.view)

        #Buttons and stuff
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(250)
        mainLayout.addWidget(paramWidget)
        #stl load
        self.mesh = None
        loadMeshBtn = QPushButton("Load mesh")
        paramLayout.addWidget(loadMeshBtn)
        loadMeshBtn.clicked.connect(self.loadMesh)

        #Rotation
        self.rot = Rotation.FULL
        rotChoiceLabel = QLabel("Choose the rotation method:")
        paramLayout.addWidget(rotChoiceLabel)
        rotLayout = QHBoxLayout()
        rb = QRadioButton("Full rotation")
        rb.setChecked(True)
        rb.rot = Rotation.FULL
        rb.toggled.connect(self.toggleRot)
        rotLayout.addWidget(rb)
        rb = QRadioButton("Stretched rotation")
        rb.rot = Rotation.STRETCH
        rb.toggled.connect(self.toggleRot)
        rotLayout.addWidget(rb)
        paramLayout.addLayout(rotLayout)

        self.rotAngle = 0
        rotAngleLabel = QLabel("Choose rotation angle:")
        paramLayout.addWidget(rotAngleLabel)
        rotAngleSb = QSpinBox()
        rotAngleSb.setRange(-180, 180)
        rotAngleSb.setSingleStep(1)
        rotAngleSb.setValue(0)
        rotAngleSb.valueChanged.connect(self.setRotAngle)
        paramLayout.addWidget(rotAngleSb)

        self.loadMesh()

    #Utility functions
    def displayMesh(self):
        glmesh = gl.MeshData(vertexes=self.dpMesh.vertices, faces=self.dpMesh.faces)
        self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh)
        self.view.addItem(self.meshItem)

    def rotate(self):
        agl = math.radians(self.rotAngle)
        rotMatrix = np.array([[math.cos(agl), 0, math.sin(agl)],[0, 1, 0],[-math.sin(agl), 0, math.cos(agl)]])
        if self.rot == Rotation.FULL:
            for i in range(len(self.mesh.vertices)):
                v = np.array(self.mesh.vertices[i])
                self.dpMesh.vertices[i] = np.matmul(rotMatrix,v)
        else:
            for i in range(len(self.mesh.vertices)):
                v = np.array(self.mesh.vertices[i])
                self.dpMesh.vertices[i][2] = np.matmul(rotMatrix,v)[2]
        self.displayMesh()

    def resetdpMeshValues(self):
        self.dpMesh = copy.deepcopy(self.mesh)

    #Slot functions
    @pyqtSlot()
    def loadMesh(self):
        #self.meshPath,_ = QFileDialog.getOpenFileName(self, 'Open mesh', '../resources', '*.stl')
        #f self.meshPath == '':
        #    print("No file chosen")
        #    return
        #self.mesh = tm.load(self.meshPath)
        self.mesh = tm.load("../resources/basic_sole.stl")

        #We process the mesh so that it is centered and lays flat on the xy plane.
        #We then copy it t a display mesh that will sustain all the tranformation we apply

        self.mesh.vertices -= self.mesh.centroid
        offset = np.amin(self.mesh.vertices,axis=0)[2]
        self.mesh.vertices -= [0,0,offset]
        self.dpMesh = copy.deepcopy(self.mesh)
        self.displayMesh()


    @pyqtSlot()
    def toggleRot(self):
        rb = self.sender()
        if rb.isChecked():
            self.rot = rb.rot
            self.resetdpMeshValues()
            self.rotate()

    @pyqtSlot()
    def setRotAngle(self):
        sb = self.sender()
        self.rotAngle = sb.value()
        self.rotate()

#residue of development, can be handy
    @pyqtSlot()
    def displayCameraInfo(self):
        dbg(self.view.opts['distance'])
        dbg(self.view.opts['elevation'])
        dbg(self.view.opts['azimuth'])
