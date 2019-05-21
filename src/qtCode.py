#default libs
import sys
import time

#3rd party libs
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import trimesh as tm
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl

#custom files
from debugThingies import *

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
        mainLayout = QVBoxLayout()
        self.setLayout(mainLayout)


        self.view = gl.GLViewWidget()
        zgrid = gl.GLGridItem()
        zgrid.scale(20, 20, 20)
        self.view.addItem(zgrid)
        self.meshItem = gl.GLMeshItem()
        self.view.addItem(self.meshItem)
        self.view.setCameraPosition(distance=288,elevation=32,azimuth=324)


        mainLayout.addWidget(self.view)

        cameraBtn = QPushButton("camera")
        cameraBtn.clicked.connect(self.displayCameraInfo)
        mainLayout.addWidget(cameraBtn)


        self.loadSTL()

    @pyqtSlot()
    def loadSTL(self):
        #self.meshPath,_ = QFileDialog.getOpenFileName(self, 'Open mesh', '../resources', '*.stl')
        #f self.meshPath == '':
        #    print("No file chosen")
        #    return
        #self.mesh = tm.load(self.meshPath)
        self.mesh = tm.load("../resources/basic_sole.stl")
        self.mesh.vertices -= self.mesh.centroid
        offset = np.amin(self.mesh.vertices,axis=0)[2]
        self.mesh.vertices -= [0,0,offset]
        self.displayMesh()

    def displayMesh(self):
        glmesh = gl.MeshData(vertexes=self.mesh.vertices, faces=self.mesh.faces)
        self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh)
        self.view.addItem(self.meshItem)

    @pyqtSlot()
    def displayCameraInfo(self):
        dbg(self.view.opts['distance'])
        dbg(self.view.opts['elevation'])
        dbg(self.view.opts['azimuth'])
