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
from scipy.spatial import ConvexHull, convex_hull_plot_2d

#custom files
from debugThingies import *

class Rotation(Enum):
    FULL = auto()
    STRETCH = auto()

class Display(Enum):
    SOLE = auto()
    INSOLE = auto()
    MOLD = auto()

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


        tabs = QTabWidget()

        #2D view
        self.canvas = pg.GraphicsLayoutWidget()
        self.plot = self.canvas.addPlot()
        self.pt = self.plot.plot(pen='w')

        #3D view
        self.view = gl.GLViewWidget()
        zgrid = gl.GLGridItem()
        zgrid.scale(20, 20, 20)
        self.view.addItem(zgrid)
        self.meshItem = gl.GLMeshItem()
        self.view.addItem(self.meshItem)
        self.view.setCameraPosition(distance=288,elevation=32,azimuth=324)

        #Adding tabs
        tabs.addTab(self.view, "3D")
        tabs.addTab(self.canvas, "2D")

        mainLayout.addWidget(tabs)

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
        rotGroup = QButtonGroup(paramWidget)
        rb = QRadioButton("Full rotation")
        rb.setChecked(True)
        rb.rot = Rotation.FULL
        rb.toggled.connect(self.toggleRot)
        rotGroup.addButton(rb)
        rotLayout.addWidget(rb)
        rb = QRadioButton("Stretched rotation")
        rb.rot = Rotation.STRETCH
        rb.toggled.connect(self.toggleRot)
        rotGroup.addButton(rb)
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

        #Choice of display between sole, insole and mold
        self.disp = Display.SOLE
        displayChoiceLabel = QLabel("Choose the model to display:")
        paramLayout.addWidget(displayChoiceLabel)
        dispLayout = QHBoxLayout()
        dispGroup = QButtonGroup(paramWidget)
        rb2 = QRadioButton("Sole")
        rb2.setChecked(True)
        rb2.disp = Display.SOLE
        rb2.toggled.connect(self.toggleDisp)
        dispGroup.addButton(rb2)
        dispLayout.addWidget(rb2)
        rb2 = QRadioButton("Insole")
        rb2.disp = Display.INSOLE
        rb2.toggled.connect(self.toggleDisp)
        dispGroup.addButton(rb2)
        dispLayout.addWidget(rb2)
        rb2 = QRadioButton("Mold")
        rb2.disp = Display.MOLD
        rb2.toggled.connect(self.toggleDisp)
        dispGroup.addButton(rb2)
        dispLayout.addWidget(rb2)
        paramLayout.addLayout(dispLayout)

        self.loadMesh()

    #Utility functions
    def displayMesh(self):
        self.resetMeshesValues()
        self.rotate()
        self.makeInsoleAndMold()
        mesh = None
        if self.disp == Display.SOLE:
            mesh = self.soleMesh
        if self.disp == Display.INSOLE:
            mesh = self.insoleMesh
        if self.disp == Display.MOLD:
            mesh = self.moldMesh
        glmesh = gl.MeshData(vertexes=mesh.vertices, faces=mesh.faces)
        self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh)
        self.view.addItem(self.meshItem)

    def rotate(self):
        agl = math.radians(self.rotAngle)
        rotMatrix = np.array([[math.cos(agl), 0, math.sin(agl)],[0, 1, 0],[-math.sin(agl), 0, math.cos(agl)]])
        if self.rot == Rotation.FULL:
            for i in range(len(self.mesh.vertices)):
                v = np.array(self.mesh.vertices[i])
                self.soleMesh.vertices[i] = np.matmul(rotMatrix,v)
        else:
            for i in range(len(self.mesh.vertices)):
                v = np.array(self.mesh.vertices[i])
                self.soleMesh.vertices[i][2] = np.matmul(rotMatrix,v)[2]

    def resetMeshesValues(self):
        self.soleMesh = copy.deepcopy(self.mesh)
        self.insoleMesh = copy.deepcopy(self.mesh)
        self.moldMesh = copy.deepcopy(self.mesh)

    def getOuterVerticesIndexes(self) -> list:
        #Convex hull does not work because the 2d projection is not convex.
        #So we check the edges that only belong to one polygon
        numUniqueEdges = len(self.insoleMesh.edges_unique)
        borderEdges = np.zeros(numUniqueEdges)
        for face in self.insoleMesh.faces_unique_edges:
            borderEdges[face[0]] += 1
            borderEdges[face[1]] += 1
            borderEdges[face[2]] += 1
        edges = np.where(borderEdges < 2)[0]
        #Now that we have the edges, we have all the vertices indices.
        #We need to order them tho.
        #Since every vertex is part of two edges, we take the first edges, take the next vertex and so on
        neighbors = self.insoleMesh.vertex_neighbors
        #the casting to list might be useless
        indices_unordered = [list(x) for x in self.insoleMesh.edges_unique[edges]]
        indices_ordered = copy.copy(indices_unordered[0])
        for i in range(len(edges) - 2):
            eds = [x for x in indices_unordered if indices_ordered[-1] in x]
            for i in range(2):
                for j in range(2):
                    if not eds[i][j] in indices_ordered:
                        indices_ordered.append(eds[i][j])
        return indices_ordered

    def makeInsoleAndMold(self):
        #----------------------------Insole
        zPos = np.amin(self.soleMesh.vertices,axis=0)[2] - 2
        indices = self.getOuterVerticesIndexes()
        n = len(self.mesh.vertices)
        #New vertices and faces vectors for new mesh
        verts = [list(x) for x in self.soleMesh.vertices]
        faces = [list(x) for x in self.soleMesh.faces]
        #Create first extruded vertex
        newVertex = copy.copy(verts[indices[0]])
        newVertex[2] = zPos
        verts.append(newVertex)
        for i in range(1,len(indices)):
            #Create new extruded vertex
            newVertex = copy.copy(verts[indices[i]])
            newVertex[2] = zPos
            verts.append(newVertex)
            #Create new extruded faces
            faces.append([indices[i-1], n + i, n + i - 1])
            faces.append([indices[i], n + i, indices[i-1]])
        #Create last extruded faces
        faces.append([indices[-1], n, n + len(indices) - 1])
        faces.append([indices[0], n, indices[-1]])
        #Create bottom faces, centroid is in 0,0,0 so we take 0,0,zPos
        verts.append([0,0,zPos])
        lastIndex = len(verts) -1
        n2 = len(faces)
        bottomVertices = verts[n:]
        for i in range(n, lastIndex):
            faces.append([i, i+1, lastIndex])
        faces.append([lastIndex - 1, n, lastIndex])
        bottomFaces = faces[n2:]
        self.insoleMesh = tm.Trimesh(vertices=verts, faces=faces)
        self.insoleMesh.export("../resources/truc.stl")

        #----------------------------Mold
        #first step: extract bottom of insole
        bottomFaces = list(map(lambda x: list(map(lambda y: y-n,x)), bottomFaces))
        scaleMat = np.array([[1.2,0,0],[0,1.2,0],[0,0,1]])
        for i in range(len(bottomVertices)):
            bottomVertices[i] = np.matmul(scaleMat,np.array(bottomVertices[i]))

        zPos2 = np.amax(self.soleMesh.vertices,axis=0)[2] + 5
        nBottom = len(bottomVertices)
        #Create first extruded vertex
        newVertex = copy.copy(bottomVertices[0])
        newVertex[2] = zPos2
        bottomVertices.append(newVertex)
        for i in range(1,nBottom-1):
            #Create new extruded vertex
            newVertex = copy.copy(bottomVertices[i])
            newVertex[2] = zPos2
            bottomVertices.append(newVertex)
            #Create new extruded faces
            bottomFaces.append([nBottom + i - 1, i, i - 1])
            bottomFaces.append([nBottom + i, i, nBottom + i - 1])
        #Create last extruded faces
        bottomFaces.append([2*nBottom - 2, 0, nBottom - 2])
        bottomFaces.append([nBottom, 0, 2*nBottom - 2])

        #Top faces
        bottomVertices.append([0,0,zPos2])
        lastIndex2 = len(bottomVertices) -1
        for i in range(nBottom, lastIndex2):
            bottomFaces.append([i, i+1, lastIndex2])
        bottomFaces.append([lastIndex2 - 1, nBottom, lastIndex2])

        for v in bottomVertices:
            v[2] += .5

        bloc = tm.Trimesh(vertices=bottomVertices, faces=bottomFaces)

        self.moldMesh = bloc.difference(self.insoleMesh)
        bloc.export("../resources/machin.stl")
        self.soleMesh.export("../resources/sole.stl")
        self.moldMesh.export("../resources/bidule.stl")


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
        self.soleMesh = copy.deepcopy(self.mesh)
        self.insoleMesh = copy.deepcopy(self.mesh)
        self.moldMesh = copy.deepcopy(self.mesh)
        self.displayMesh()

    @pyqtSlot()
    def toggleRot(self):
        rb = self.sender()
        if rb.isChecked():
            self.rot = rb.rot
            self.displayMesh()

    @pyqtSlot()
    def setRotAngle(self):
        sb = self.sender()
        self.rotAngle = sb.value()
        self.displayMesh()

    @pyqtSlot()
    def toggleDisp(self):
        rb = self.sender()
        if rb.isChecked():
            self.disp = rb.disp
            self.displayMesh()

#residue of development, can be handy
    @pyqtSlot()
    def displayCameraInfo(self):
        dbg(self.view.opts['distance'])
        dbg(self.view.opts['elevation'])
        dbg(self.view.opts['azimuth'])
