#default libs
import sys
import time
from enum import Enum, auto
import math
import copy
import os

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

        #3D view
        self.view = gl.GLViewWidget()
        zgrid = gl.GLGridItem()
        zgrid.scale(20, 20, 20)
        axis = gl.GLAxisItem()
        axis.setSize(x=2, y=1, z=1)
        axis.scale(100,100,100)
        self.view.addItem(zgrid)
        self.view.addItem(axis)
        self.meshItem = gl.GLMeshItem()
        self.view.addItem(self.meshItem)
        self.view.setCameraPosition(distance=288,elevation=32,azimuth=324)

        mainLayout.addWidget(self.view)

        #Buttons and stuff
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(500)
        mainLayout.addWidget(paramWidget)
        #stl load
        self.mesh = None
        loadMeshBtn = QPushButton("Load mesh")
        paramLayout.addWidget(loadMeshBtn)
        loadMeshBtn.clicked.connect(self.loadMesh)

        #Z rotate
        self.zRotAngle = 0
        zRotAngleLabel = QLabel("Choose rotation angle around z:")
        paramLayout.addWidget(zRotAngleLabel)
        zRotAngleSb = QSpinBox()
        zRotAngleSb.setRange(-180, 180)
        zRotAngleSb.setSingleStep(1)
        zRotAngleSb.setValue(0)
        zRotAngleSb.valueChanged.connect(self.setZRotAngle)
        paramLayout.addWidget(zRotAngleSb)

        #Rotation
        self.linearDescentPortion = .5
        linDescentLabel = QLabel("Choose portion of linearly descending angle:")
        paramLayout.addWidget(linDescentLabel)
        linDescentSb = QDoubleSpinBox()
        linDescentSb.setRange(0,1)
        linDescentSb.setSingleStep(.01)
        linDescentSb.setValue(0.5)
        linDescentSb.setDecimals(2)
        linDescentSb.valueChanged.connect(self.setLinDescent)
        paramLayout.addWidget(linDescentSb)

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
        modelCkBx = QCheckBox("Compute Insole and Mold model")
        modelCkBx.stateChanged.connect(self.toggleModelCompute)
        paramLayout.addWidget(modelCkBx)
        self.disp = Display.SOLE
        displayChoiceLabel = QLabel("Choose the model to display:")
        paramLayout.addWidget(displayChoiceLabel)
        dispLayout = QHBoxLayout()
        dispGroup = QButtonGroup(paramWidget)
        self.rbSole = QRadioButton("Sole")
        self.rbSole.setChecked(True)
        self.rbSole.disp = Display.SOLE
        self.rbSole.toggled.connect(self.toggleDisp)
        dispGroup.addButton(self.rbSole)
        dispLayout.addWidget(self.rbSole)
        self.rbInsole = QRadioButton("Insole")
        self.rbInsole.disp = Display.INSOLE
        self.rbInsole.toggled.connect(self.toggleDisp)
        self.rbInsole.setEnabled(False)
        dispGroup.addButton(self.rbInsole)
        dispLayout.addWidget(self.rbInsole)
        self.rbMold = QRadioButton("Mold")
        self.rbMold.disp = Display.MOLD
        self.rbMold.toggled.connect(self.toggleDisp)
        self.rbMold.setEnabled(False)
        dispGroup.addButton(self.rbMold)
        dispLayout.addWidget(self.rbMold)
        paramLayout.addLayout(dispLayout)

        #export buttons
        expLayout = QHBoxLayout()
        expBtn = QPushButton("Export Sole Model")
        expBtn.model = Display.SOLE
        expBtn.clicked.connect(self.exportModel)
        expLayout.addWidget(expBtn)
        self.expBtnInsole = QPushButton("Export Insole Model")
        self.expBtnInsole.model = Display.INSOLE
        self.expBtnInsole.clicked.connect(self.exportModel)
        self.expBtnInsole.setEnabled(False)
        expLayout.addWidget(self.expBtnInsole)
        self.expBtnMold = QPushButton("Export Mold Model")
        self.expBtnMold.model = Display.MOLD
        self.expBtnMold.clicked.connect(self.exportModel)
        self.expBtnMold.setEnabled(False)
        expLayout.addWidget(self.expBtnMold)
        paramLayout.addLayout(expLayout)

        #shader choice
        shaderLabel = QLabel("Choose display shader:")
        paramLayout.addWidget(shaderLabel)
        shaderCB = QComboBox(self)
        shaderCB.addItem("balloon")
        shaderCB.addItem("normalColor")
        shaderCB.addItem("viewNormalColor")
        shaderCB.addItem("shaded")
        shaderCB.setCurrentIndex(1)
        shaderCB.activated.connect(self.shaderSelect)
        paramLayout.addWidget(shaderCB)

        #Values
        self.recompute = True
        self.shader = "normalColor"
        self.glOptions = "opaque"
        self.computeInsoleAndMold = False

        self.loadBasicMesh()

    #Utility functions
    def displayMesh(self):
        if self.recompute:
            self.resetMeshesValues()
            self.rotateZ()
            self.rotate()
            if self.computeInsoleAndMold:
                self.makeInsoleAndMold()
            self.recompute = False
        mesh = None
        if self.disp == Display.SOLE:
            mesh = self.soleMesh
        if self.disp == Display.INSOLE:
            mesh = self.insoleMesh
        if self.disp == Display.MOLD:
            mesh = self.moldMesh
        glmesh = gl.MeshData(vertexes=mesh.vertices, faces=mesh.faces)
        self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.meshItem)

    def computeAngles(self) -> list:
        agl = math.radians(self.rotAngle)
        angles = [0]*len(self.soleMesh.vertices)
        minX = np.amin(self.soleMesh.vertices,axis=0)[0]
        maxX = np.amax(self.soleMesh.vertices,axis=0)[0]
        mid = minX + (maxX - minX) * (1-self.linearDescentPortion)
        for i in range(len(self.soleMesh.vertices)):
            if self.soleMesh.vertices[i][0] < mid:
                angles[i] = agl
            else:
                angles[i] = agl * (1 - (self.soleMesh.vertices[i][0] - mid)/(maxX - mid))
        return angles

    def rotate(self):
        agls = self.computeAngles()
        if self.rot == Rotation.FULL:
            for i in range(len(self.soleMesh.vertices)):
                rotMatrix = np.array([[1,0,0],[0,math.cos(agls[i]),-math.sin(agls[i])],[0,math.sin(agls[i]),math.cos(agls[i])]])
                v = np.array(self.soleMesh.vertices[i])
                self.soleMesh.vertices[i] = np.matmul(rotMatrix,v)
        else:
            for i in range(len(self.soleMesh.vertices)):
                rotMatrix = np.array([[1,0,0],[0,math.cos(agls[i]),-math.sin(agls[i])],[0,math.sin(agls[i]),math.cos(agls[i])]])
                v = np.array(self.soleMesh.vertices[i])
                self.soleMesh.vertices[i][2] = np.matmul(rotMatrix,v)[2]
        offset = np.amin(self.soleMesh.vertices,axis=0)[2] - 2
        self.soleMesh.vertices -= [0,0,offset]

    def rotateZ(self):
        agl = math.radians(self.zRotAngle)
        c = math.cos(agl)
        s = math.sin(agl)
        rotMatrix = np.array([[c,-s,0],[s,c,0],[0,0,1]])
        for i in range(len(self.mesh.vertices)):
            v = np.array(self.mesh.vertices[i])
            self.soleMesh.vertices[i] = np.matmul(rotMatrix,v)

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
        #We flip the mold for better viewing
        zOffset = np.amax(self.soleMesh.vertices,axis=0)[2]
        rotMatrix = np.array([[-1, 0, 0],[0, 1, 0],[0, 0, -1]])
        for i in range(len(self.moldMesh.vertices)):
            v = np.array(self.moldMesh.vertices[i])
            self.moldMesh.vertices[i] = np.matmul(rotMatrix,v)
            self.moldMesh.vertices[i][2] += zOffset

    def loadBasicMesh(self):
        path = str(os.path.dirname(os.path.realpath(__file__)))
        dbg(path)
        self.mesh = tm.load(path+"\\..\\resources\\basic_sole.stl")
        #We process the mesh so that it is centered and lays flat on the xy plane.
        #We then copy it t a display mesh that will sustain all the tranformation we apply
        self.mesh.vertices -= self.mesh.centroid
        transforms, probs = self.mesh.compute_stable_poses()
        self.mesh.apply_transform(transforms[np.argmax(probs)])
        offset = np.amin(self.mesh.vertices,axis=0)[2]
        self.mesh.vertices -= [0,0,offset]
        self.soleMesh = copy.deepcopy(self.mesh)
        self.insoleMesh = copy.deepcopy(self.mesh)
        self.moldMesh = copy.deepcopy(self.mesh)
        self.recompute = True
        self.displayMesh()


    #Slot functions
    @pyqtSlot()
    def loadMesh(self):
        path = str(os.path.dirname(os.path.realpath(__file__)))
        self.meshPath,_ = QFileDialog.getOpenFileName(self, 'Open mesh', path + "\\..\\resources", '*.stl')
        if self.meshPath == '':
            print("No file chosen")
            return
        self.mesh = tm.load(self.meshPath)

        #We process the mesh so that it is centered and lays flat on the xy plane.
        #We then copy it t a display mesh that will sustain all the tranformation we apply

        self.mesh.vertices -= self.mesh.centroid
        transforms, probs = self.mesh.compute_stable_poses()
        self.mesh.apply_transform(transforms[np.argmax(probs)])
        offset = np.amin(self.mesh.vertices,axis=0)[2]
        self.mesh.vertices -= [0,0,offset]
        self.soleMesh = copy.deepcopy(self.mesh)
        self.insoleMesh = copy.deepcopy(self.mesh)
        self.moldMesh = copy.deepcopy(self.mesh)
        self.recompute = True
        self.displayMesh()

    @pyqtSlot()
    def toggleRot(self):
        rb = self.sender()
        if rb.isChecked():
            self.rot = rb.rot
            self.recompute = True
            self.displayMesh()

    @pyqtSlot()
    def setRotAngle(self):
        sb = self.sender()
        self.rotAngle = sb.value()
        self.recompute = True
        self.displayMesh()

    @pyqtSlot()
    def setZRotAngle(self):
        sb = self.sender()
        self.zRotAngle = sb.value()
        self.recompute = True
        self.displayMesh()

    @pyqtSlot()
    def toggleDisp(self):
        rb = self.sender()
        if rb.isChecked():
            self.disp = rb.disp
            self.displayMesh()


    @pyqtSlot()
    def exportModel(self):
        btn = self.sender()
        mesh = None
        if btn.model == Display.SOLE:
            mesh = self.soleMesh
        if btn.model == Display.INSOLE:
            mesh = self.insoleMesh
        if btn.model == Display.MOLD:
            mesh = self.moldMesh
        path = QFileDialog.getSaveFileName(self, "Save Model", "../resources", "STL(*.stl)")[0]
        if path == "":
            print("No file chosen")
            return
        mesh.export(path)

    @pyqtSlot()
    def shaderSelect(self):
        cb = self.sender()
        self.shader = cb.currentText()
        if self.shader == "balloon":
            self.glOptions = "additive"
        else:
            self.glOptions = "opaque"
        self.meshItem.setShader(self.shader)
        self.meshItem.setGLOptions(self.glOptions)

    @pyqtSlot()
    def toggleModelCompute(self):
        cb = self.sender()
        self.computeInsoleAndMold = cb.isChecked()
        if self.computeInsoleAndMold:
            self.recompute = True
            self.displayMesh()
        else:
            self.rbSole.setChecked(True)
        self.rbInsole.setEnabled(self.computeInsoleAndMold)
        self.rbMold.setEnabled(self.computeInsoleAndMold)
        self.expBtnInsole.setEnabled(self.computeInsoleAndMold)
        self.expBtnMold.setEnabled(self.computeInsoleAndMold)

    @pyqtSlot()
    def setLinDescent(self):
        sb = self.sender()
        self.linearDescentPortion = sb.value()
        self.recompute = True
        self.displayMesh()

#residue of development, can be handy
    @pyqtSlot()
    def displayCameraInfo(self):
        dbg(self.view.opts['distance'])
        dbg(self.view.opts['elevation'])
        dbg(self.view.opts['azimuth'])
