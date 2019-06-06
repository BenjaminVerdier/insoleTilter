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
        self.resize(1300,1800)

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
        self.planeXMesh = None
        self.planeZMesh = None
        self.view.setCameraPosition(distance=288,elevation=32,azimuth=324)

        mainLayout.addWidget(self.view)

        #Buttons and stuff
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(800)
        mainLayout.addWidget(paramWidget)
        #mesh load and invertion
        meshLoadLayout = QHBoxLayout()

        self.mesh = None
        loadMeshBtn = QPushButton("Load mesh")
        meshLoadLayout.addWidget(loadMeshBtn)
        loadMeshBtn.clicked.connect(self.loadMesh)

        self.invertNormals = False
        NrmInvCkBx = QCheckBox("Invert normals (top should look purple-ish in normalColor shader)")
        NrmInvCkBx.stateChanged.connect(self.toggleInvertNormals)
        meshLoadLayout.addWidget(NrmInvCkBx)
        paramLayout.addLayout(meshLoadLayout)

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

        #Cuts locations
        cutsLayout = QHBoxLayout()

        cutXLayout = QVBoxLayout()
        self.doCutX = False
        cutCkBx = QCheckBox("Cut part of the model on X axis")
        cutCkBx.stateChanged.connect(self.toggleCutX)
        cutXLayout.addWidget(cutCkBx)
        self.cutXLocation = .9
        cutLabel = QLabel("Choose where the cut should be located:")
        cutXLayout.addWidget(cutLabel)
        cutSb = QDoubleSpinBox()
        cutSb.setRange(0,1)
        cutSb.setSingleStep(.01)
        cutSb.setValue(0.9)
        cutSb.setDecimals(2)
        cutSb.valueChanged.connect(self.setCutXLocation)
        cutXLayout.addWidget(cutSb)
        self.showPlaneCutX = False
        cutPlaneCkBx = QCheckBox("Show cut location")
        cutPlaneCkBx.stateChanged.connect(self.toggleShowPlaneCutX)
        cutXLayout.addWidget(cutPlaneCkBx)

        cutZLayout = QVBoxLayout()
        self.doCutZ = False
        cutCkBx = QCheckBox("Cut part of the model on Z axis")
        cutCkBx.stateChanged.connect(self.toggleCutZ)
        cutZLayout.addWidget(cutCkBx)
        self.cutZLocation = .5
        cutLabel = QLabel("Choose where the cut should be located:")
        cutZLayout.addWidget(cutLabel)
        cutSb = QDoubleSpinBox()
        cutSb.setRange(0,1)
        cutSb.setSingleStep(.01)
        cutSb.setValue(0.5)
        cutSb.setDecimals(2)
        cutSb.valueChanged.connect(self.setCutZLocation)
        cutZLayout.addWidget(cutSb)
        self.showPlaneCutZ = False
        cutPlaneCkBx = QCheckBox("Show cut location")
        cutPlaneCkBx.stateChanged.connect(self.toggleShowPlaneCutZ)
        cutZLayout.addWidget(cutPlaneCkBx)

        cutsLayout.addLayout(cutXLayout)
        cutsLayout.addLayout(cutZLayout)
        paramLayout.addLayout(cutsLayout)

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
            if self.doCutZ:
                self.cutTopOfModel()
            if self.doCutX:
                self.cutSoleX()
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
        self.moldMesh = copy.deepcopy(self.mesh)

    def getOuterVerticesIndexes(self) -> list:
        #Convex hull does not work because the 2d projection is not convex.
        #So we check the edges that only belong to one polygon
        numUniqueEdges = len(self.soleMesh.edges_unique)
        borderEdges = np.zeros(numUniqueEdges)
        for face in self.soleMesh.faces_unique_edges:
            borderEdges[face[0]] += 1
            borderEdges[face[1]] += 1
            borderEdges[face[2]] += 1
        edges = np.where(borderEdges < 2)[0]
        #Now that we have the edges, we have all the vertices indices.
        #We need to order them tho.
        #Since every vertex is part of two edges, we take the first edges, take the next vertex and so on
        neighbors = self.soleMesh.vertex_neighbors
        #the casting to list might be useless
        indices_unordered = [list(x) for x in self.soleMesh.edges_unique[edges]]
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
        n = len(self.soleMesh.vertices)
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
            faces.append([n + i, indices[i-1], n + i - 1])
            faces.append([n + i, indices[i], indices[i-1]])
        #Create last extruded faces
        faces.append([n, indices[-1], n + len(indices) - 1])
        faces.append([n, indices[0], indices[-1]])
        #Create bottom faces, centroid is in 0,0,0 so we take 0,0,zPos
        verts.append([0,0,zPos])
        lastIndex = len(verts) -1
        n2 = len(faces)
        bottomVertices = verts[n:]
        for i in range(n, lastIndex):
            faces.append([i+1, i, lastIndex])
        faces.append([n, lastIndex - 1, lastIndex])
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
            bottomFaces.append([i, nBottom + i - 1, i - 1])
            bottomFaces.append([i, nBottom + i, nBottom + i - 1])
        #Create last extruded faces
        bottomFaces.append([0, 2*nBottom - 2, nBottom - 2])
        bottomFaces.append([0, nBottom, 2*nBottom - 2])

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

    def cutSoleXDeprecated(self):
        #let's copy our vertices and faces
        verts = [np.array(x) for x in self.soleMesh.vertices]
        faces = [list(x) for x in self.soleMesh.faces]
        #indices of vertices to delete
        indDel = []
        minX = np.amin(verts,axis=0)[0]
        maxX = np.amax(verts,axis=0)[0]
        mid = minX + (maxX - minX) * self.cutLocation
        for i in range(len(verts)):
            if verts[i][0] > mid:
                indDel.append(i)
        #indices of faces to delete
        facesToDel = []
        for i in range(len(faces)):
            if len(list(set(faces[i]) & set(indDel))) > 0:
                facesToDel.append(i)
        #indices of orphan vertices
        orphanVerts = []
        for i in facesToDel:
            for j in faces[i]:
                if not j in indDel:
                    orphanVerts.append(j)
        orphanVerts = list(set(orphanVerts))
        #ordering by y:
        orphAndY = [(i, verts[i][1]) for i in orphanVerts]
        sorted_by_second = sorted(orphAndY, key=lambda tup: tup[1], reverse = True)
        orphanVertsOrdered = [x[0] for x in sorted_by_second]
        #Computation to straighten the edge
        newPtsIndexes = dict()
        for i in range(len(orphanVertsOrdered)-1):
            orph = orphanVertsOrdered[i]
            facesWithOrph = []
            for f in facesToDel:
                if orph in faces[f]:
                    facesWithOrph.append(faces[f])
            for f in facesWithOrph:
                if orphanVertsOrdered[i-1] in f:
                    #We did it when we were doing i-1
                    continue
                if orphanVertsOrdered[i+1] in f:
                    #Only 1 point is missing, we get its index
                    otherPoint = -1
                    for pt in f:
                        if not pt == orph and not pt == orphanVertsOrdered[i+1]:
                            otherPoint = pt
                            break
                    #we create two new points: one between orph and otherpoint, one between orph+1 and otherPoint
                    orphVal = verts[orph]
                    orphPlus1Val = verts[orphanVertsOrdered[i+1]]
                    otherPtVal = verts[otherPoint]
                    if not (orph,otherPoint) in newPtsIndexes:
                        k1 = (mid - orphVal[0])/(otherPtVal[0] - orphVal[0])
                        newPt1 = orphVal + k1 * (otherPtVal - orphVal)
                        newPtsIndexes[(orph,otherPoint)] = len(verts)
                        verts.append(newPt1)
                    if not (orphanVertsOrdered[i+1],otherPoint) in newPtsIndexes:
                        k2 = (mid - orphPlus1Val[0])/(otherPtVal[0] - orphPlus1Val[0])
                        newPt2 = orphPlus1Val + k2 * (otherPtVal - orphPlus1Val)
                        newPtsIndexes[(orphanVertsOrdered[i+1],otherPoint)] = len(verts)
                        verts.append(newPt2)
                    #We create the corresponding faces
                    faces.append([orphanVertsOrdered[i+1], newPtsIndexes[(orph,otherPoint)], orph])
                    faces.append([newPtsIndexes[(orphanVertsOrdered[i+1],otherPoint)], newPtsIndexes[(orph,otherPoint)], orphanVertsOrdered[i+1]])
                else:
                    #Two points are missing, we get their indices
                    otherPoint1 = -1
                    otherPoint2 = -1
                    for pt in f:
                        if not pt == orph:
                            if otherPoint1 == -1:
                                otherPoint1 = pt
                            else:
                                otherPoint2 = pt
                    #We make sure otherPoint1 in the leftmost (biggest y)
                    if verts[otherPoint1][1] < verts[otherPoint2][1]:
                        otherPoint1, otherPoint2 = otherPoint2, otherPoint1
                    #We create two new points
                    orphVal = verts[orph]
                    otherPt1Val = verts[otherPoint1]
                    otherPt2Val = verts[otherPoint2]
                    if not (orph,otherPoint1) in newPtsIndexes:
                        k1 = (mid - orphVal[0])/(otherPt1Val[0] - orphVal[0])
                        newPt1 = orphVal + k1 * (otherPt1Val - orphVal)
                        newPtsIndexes[(orph,otherPoint1)] = len(verts)
                        verts.append(newPt1)
                    if not (orph,otherPoint2) in newPtsIndexes:
                        k2 = (mid - orphVal[0])/(otherPt2Val[0] - orphVal[0])
                        newPt2 = orphVal + k2 * (otherPt2Val - orphVal)
                        newPtsIndexes[(orph,otherPoint2)] = len(verts)
                        verts.append(newPt2)
                    #We create the corresponding face
                    faces.append([newPtsIndexes[(orph,otherPoint2)], newPtsIndexes[(orph,otherPoint1)], orph])

        #We might need one last face
        lastPoint = orphanVertsOrdered[-1]
        facesWithOrph = []
        for f in facesToDel:
            if lastPoint in faces[f]:
                facesWithOrph.append(faces[f])
        for f in facesWithOrph:
            if orphanVertsOrdered[-2] in f:
                #We already did it
                continue
            else:
                #Two points are missing, we get their indices
                otherPoint1 = -1
                otherPoint2 = -1
                for pt in f:
                    if not pt == lastPoint:
                        if otherPoint1 == -1:
                            otherPoint1 = pt
                        else:
                            otherPoint2 = pt
                #We make sure otherPoint1 in the leftmost (biggest y)
                if verts[otherPoint1][1] < verts[otherPoint2][1]:
                    otherPoint1, otherPoint2 = otherPoint2, otherPoint1
                #We create two new points
                orphVal = verts[lastPoint]
                otherPt1Val = verts[otherPoint1]
                otherPt2Val = verts[otherPoint2]
                if not (lastPoint,otherPoint1) in newPtsIndexes:
                    k1 = (mid - orphVal[0])/(otherPt1Val[0] - orphVal[0])
                    newPt1 = orphVal + k1 * (otherPt1Val - orphVal)
                    newPtsIndexes[(lastPoint,otherPoint1)] = len(verts)
                    verts.append(newPt1)
                if not (lastPoint,otherPoint2) in newPtsIndexes:
                    k2 = (mid - orphVal[0])/(otherPt2Val[0] - orphVal[0])
                    newPt2 = orphVal + k2 * (otherPt2Val - orphVal)
                    newPtsIndexes[(lastPoint,otherPoint2)] = len(verts)
                    verts.append(newPt2)
                #We create the corresponding face
                faces.append([newPtsIndexes[(lastPoint,otherPoint2)], newPtsIndexes[(lastPoint,otherPoint1)], lastPoint])

        #Remove vertices, modify faces
        for i in sorted(facesToDel, reverse=True):
            del faces[i]
        for face in faces:
            for k in range(3):
                face[k] -= sum(1 for i in indDel if i < face[k])
        for i in sorted(indDel, reverse=True):
            del verts[i]

        self.soleMesh = tm.Trimesh(vertices=verts, faces=faces)

    def displayPlaneCutX(self):
        zpos = np.amax(self.soleMesh.vertices,axis=0)[2]
        zneg = np.amin(self.soleMesh.vertices,axis=0)[2]
        xpos = np.amax(self.soleMesh.vertices,axis=0)[0]
        ypos = np.amax(self.soleMesh.vertices,axis=0)[1]
        yneg = np.amin(self.soleMesh.vertices,axis=0)[1]
        minX = np.amin(self.soleMesh.vertices,axis=0)[0]
        xneg = minX + (xpos - minX) * self.cutXLocation
        #cube = tm.creation.box(extents=(maxLength, maxWidth, maxHeight))
        planeVerts = np.array([[xneg,yneg,zneg],[xneg,ypos,zneg],[xneg,yneg,zpos],[xneg,ypos,zpos]])
        planeFaces = np.array([[0,1,2],[3,2,1]])

        glmesh = gl.MeshData(vertexes=planeVerts, faces=planeFaces)
        if self.planeXMesh:
            self.view.removeItem(self.planeXMesh)
        self.planeXMesh = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.planeXMesh)

    def cutSoleX(self):
        zpos = np.amax(self.soleMesh.vertices,axis=0)[2]
        zneg = np.amin(self.soleMesh.vertices,axis=0)[2]
        xpos = np.amax(self.soleMesh.vertices,axis=0)[0]
        ypos = np.amax(self.soleMesh.vertices,axis=0)[1]
        yneg = np.amin(self.soleMesh.vertices,axis=0)[1]
        minX = np.amin(self.soleMesh.vertices,axis=0)[0]
        xneg = minX + (xpos - minX) * self.cutXLocation
        #cube = tm.creation.box(extents=(maxLength, maxWidth, maxHeight))
        cubeVerts = [[xpos,ypos,zneg],[xpos,yneg,zneg],[xneg,yneg,zneg],[xneg,ypos,zneg],[xpos,ypos,zpos],[xpos,yneg,zpos],[xneg,yneg,zpos],[xneg,ypos,zpos]]
        cubeFaces = [[0,4,1],[4,5,1],[1,5,2],[5,6,2],[2,6,3],[6,7,3],[3,7,0],[7,4,0],[7,4,5],[7,5,6],[3,0,1],[3,1,2]]
        cube = tm.Trimesh(vertices=cubeVerts, faces=cubeFaces)
        self.soleMesh = self.soleMesh.difference(cube)


    def displayPlaneCutZ(self):
        z = np.amax(self.soleMesh.vertices,axis=0)[2]
        xpos = np.amax(self.soleMesh.vertices,axis=0)[0]
        ypos = np.amax(self.soleMesh.vertices,axis=0)[1]
        xneg = np.amin(self.soleMesh.vertices,axis=0)[0]
        yneg = np.amin(self.soleMesh.vertices,axis=0)[1]
        zmin = np.amin(self.soleMesh.vertices,axis=0)[2]

        zToCut = zmin + (z - zmin) * self.cutZLocation
        #cube = tm.creation.box(extents=(maxLength, maxWidth, maxHeight))
        planeVerts = np.array([[xpos,ypos,zToCut],[xpos,yneg,zToCut],[xneg,yneg,zToCut],[xneg,ypos,zToCut]])
        planeFaces = np.array([[0,1,2],[2,3,0]])

        glmesh = gl.MeshData(vertexes=planeVerts, faces=planeFaces)
        if self.planeZMesh:
            self.view.removeItem(self.planeZMesh)
        self.planeZMesh = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.planeZMesh)

    def cutTopOfModel(self):
        verts = [np.array(x) for x in self.soleMesh.vertices]
        faces = [list(x) for x in self.soleMesh.faces]
        z = np.amax(verts,axis=0)[2]
        xpos = np.amax(verts,axis=0)[0]
        ypos = np.amax(verts,axis=0)[1]
        xneg = np.amin(verts,axis=0)[0]
        yneg = np.amin(verts,axis=0)[1]
        zmin = np.amin(verts,axis=0)[2]

        zToCut = zmin + (z - zmin) * self.cutZLocation
        #cube = tm.creation.box(extents=(maxLength, maxWidth, maxHeight))
        cubeVerts = [[xpos,ypos,zToCut],[xpos,yneg,zToCut],[xneg,yneg,zToCut],[xneg,ypos,zToCut],[xpos,ypos,z],[xpos,yneg,z],[xneg,yneg,z],[xneg,ypos,z]]
        cubeFaces = [[0,4,1],[4,5,1],[1,5,2],[5,6,2],[2,6,3],[6,7,3],[3,7,0],[7,4,0],[7,4,5],[7,5,6],[3,0,1],[3,1,2]]
        cube = tm.Trimesh(vertices=cubeVerts, faces=cubeFaces)
        newSole = self.soleMesh.difference(cube)
        verts = [np.array(x) for x in newSole.vertices]
        faces = [list(x) for x in newSole.faces]
        indicesOfMaxHeight = []
        for i in range(len(verts)):
            if abs(verts[i][2] - zToCut) < 0.0001:
                indicesOfMaxHeight.append(i)
        facesToDel = []
        for i in range(len(faces)):
            if all(v in indicesOfMaxHeight for v in faces[i]):
                facesToDel.append(i)

        for i in sorted(facesToDel, reverse=True):
            del faces[i]

        self.soleMesh = tm.Trimesh(vertices=verts, faces=faces)


    #Slot functions
    @pyqtSlot()
    def loadMesh(self):
        path = str(os.path.dirname(os.path.realpath(__file__)))
        self.meshPath,_ = QFileDialog.getOpenFileName(self, 'Open mesh', path + "\\..\\resources", 'Mesh (*.stl *.obj)')
        if self.meshPath == '':
            print("No file chosen")
            return
        self.mesh = tm.load(self.meshPath)

        if self.invertNormals:
            self.mesh.invert()

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

    #TODO: refactor 6 following methods below to only have 3
    @pyqtSlot()
    def setCutXLocation(self):
        sb = self.sender()
        self.cutXLocation = sb.value()
        self.recompute = True
        self.displayMesh()


    @pyqtSlot()
    def setCutZLocation(self):
        sb = self.sender()
        self.cutZLocation = sb.value()
        self.recompute = True
        self.displayMesh()

    @pyqtSlot()
    def toggleCutX(self):
        cb = self.sender()
        self.doCutX = cb.isChecked()
        self.recompute = True
        self.displayMesh()

    @pyqtSlot()
    def toggleCutZ(self):
        cb = self.sender()
        self.doCutZ = cb.isChecked()
        self.recompute = True
        self.displayMesh()

    pyqtSlot()
    def toggleShowPlaneCutX(self):
        cb = self.sender()
        if cb.isChecked():
            self.displayPlaneCutX()
        else:
            if self.planeXMesh:
                self.view.removeItem(self.planeXMesh)
                self.planeXMesh = None

    pyqtSlot()
    def toggleShowPlaneCutZ(self):
        cb = self.sender()
        if cb.isChecked():
            self.displayPlaneCutZ()
        else:
            if self.planeZMesh:
                self.view.removeItem(self.planeZMesh)
                self.planeZMesh = None

    @pyqtSlot()
    def toggleInvertNormals(self):
        cb = self.sender()
        self.invertNormals = cb.isChecked()
        self.mesh.invert()
        self.soleMesh.invert()
        self.insoleMesh.invert()
        self.moldMesh.invert()
        self.displayMesh()



#residue of development, can be handy
    @pyqtSlot()
    def displayCameraInfo(self):
        dbg(self.view.opts['distance'])
        dbg(self.view.opts['elevation'])
        dbg(self.view.opts['azimuth'])
