from baseWidget import *

from enum import Enum, auto
import sys
import math

class Display(Enum):
    SOLE = auto()
    INSOLE = auto()
    MOLD = auto()

class insoleMoldExportWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(insoleMoldExportWidget, self).__init__(nextTab)

    def initCustomUI(self):

        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)


        #Choice of display between sole, insole and mold
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
        dispGroup.addButton(self.rbInsole)
        dispLayout.addWidget(self.rbInsole)
        self.rbMold = QRadioButton("Mold")
        self.rbMold.disp = Display.MOLD
        self.rbMold.toggled.connect(self.toggleDisp)
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
        expLayout.addWidget(self.expBtnInsole)
        self.expBtnMold = QPushButton("Export Mold Model")
        self.expBtnMold.model = Display.MOLD
        self.expBtnMold.clicked.connect(self.exportModel)
        expLayout.addWidget(self.expBtnMold)
        paramLayout.addLayout(expLayout)

        self.initDoneButtonAndShader(paramLayout)

    def start(self, transmittedData):
        self.soleMesh = transmittedData
        self.disp = Display.SOLE
        self.rbSole.setChecked(True)
        self.makeInsoleAndMold()
        self.displayMesh()

    def displayMesh(self):
        mesh = None
        if self.disp == Display.SOLE:
            mesh = self.soleMesh
        if self.disp == Display.INSOLE:
            mesh = self.insoleMesh
        if self.disp == Display.MOLD:
            mesh = self.moldMesh
        glmesh = gl.MeshData(vertexes=mesh.vertices, faces=mesh.faces)
        if not self.meshItem is None:
            self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.meshItem)


    def getOuterVerticesIndexes(self, mesh = None) -> list:
        #Convex hull does not work because the 2d projection is not convex.
        #So we check the edges that only belong to one polygon
        if mesh == None:
            mesh = self.soleMesh
        numUniqueEdges = len(mesh.edges_unique)
        borderEdges = np.zeros(numUniqueEdges)
        for face in mesh.faces_unique_edges:
            borderEdges[face[0]] += 1
            borderEdges[face[1]] += 1
            borderEdges[face[2]] += 1
        edges = np.where(borderEdges < 2)[0]
        #Now that we have the edges, we have all the vertices indices.
        #We need to order them tho.
        #Since every vertex is part of two edges, we take the first edges, take the next vertex and so on
        neighbors = mesh.vertex_neighbors
        #the casting to list might be useless
        indices_unordered = [list(x) for x in mesh.edges_unique[edges]]
        indices_ordered = copy.copy(indices_unordered[0])
        for i in range(len(edges) - 2):
            eds = [x for x in indices_unordered if indices_ordered[-1] in x]
            for i in range(2):
                for j in range(2):
                    if not eds[i][j] in indices_ordered:
                        indices_ordered.append(eds[i][j])
        return indices_ordered

    def makeInsole(self):
        indices = self.getOuterVerticesIndexes()
        n = len(self.soleMesh.vertices)

        verts = [list(x) for x in self.soleMesh.vertices]
        faces = [list(x) for x in self.soleMesh.faces]

        #We make all faces' normals face up-ish
        normals = self.soleMesh.face_normals
        for i in range(len(normals)):
            if normals[i][2] < 0:
                faces[i] = [faces[i][0],faces[i][2],faces[i][1]]

        #We make all extruded vertices at z = 0 and we invert faces so the extruded faces are facing down

        newVerts = [[v[0],v[1],0] for v in verts]
        newFaces = [[f[0]+n,f[2]+n,f[1]+n] for f in faces]

        sideFaces = 2*len(indices)*[[0,0,0]]
        for k in range(-1,len(indices)-1):
            i1 = indices[k]
            i2 = indices[k+1]
            sideFaces[2*k] = [i1,n+i1,i2]
            sideFaces[1 + 2*k] = [n+i1,n+i2,i2]

        insoleVerts = verts + newVerts
        insoleFaces = faces + newFaces + sideFaces

        self.insoleMesh = tm.Trimesh(vertices=insoleVerts, faces=insoleFaces)

    def makeMold(self):
        indices = self.getOuterVerticesIndexes()
        n = len(self.soleMesh.vertices)
        n2 = len(indices)

        offset = 10
        verts = [np.array(x) for x in self.soleMesh.vertices]
        zOffsetPlus = 5
        zOffsetMinus = 5
        maxZ = np.amax(verts,axis=0)[2]
        minZ = np.amin(verts,axis=0)[2]

        verts = [list(x) for x in self.soleMesh.vertices]
        faces = [list(x) for x in self.soleMesh.faces]

        #We make all faces' normals face up-ish
        normals = self.soleMesh.face_normals
        for i in range(len(normals)):
            if normals[i][2] > 0:
                faces[i] = [faces[i][0],faces[i][2],faces[i][1]]

        #We compute all required vertices

        topVerticesInner = [[verts[i][0],verts[i][1], minZ - zOffsetMinus] for i in indices]

        topVerticesOuter =  copy.copy(topVerticesInner)

        for i in range(n2):
            curV = topVerticesOuter[i]
            length = math.sqrt(curV[0]*curV[0] + curV[1]*curV[1])
            offsetVec = [offset*curV[0]/length, offset*curV[1]/length]
            topVerticesOuter[i] = [curV[0] + offsetVec[0], curV[1] + offsetVec[1],minZ - zOffsetMinus]

        bottomVerticesOuter = [[v[0],v[1], maxZ + zOffsetPlus] for v in topVerticesOuter]

        bottomVerticesInner = [[v[0],v[1], maxZ + zOffsetPlus] for v in verts]

        #We stitch everything together
        sideFacesInner = 2*len(indices)*[[0,0,0]]
        for k in range(n2-1):
            i1 = indices[k]
            i2 = indices[k+1]
            sideFacesInner[2*k] = [i1,n+k,i2]
            sideFacesInner[1 + 2*k] = [n+k,n+k+1,i2]
        i1 = indices[-1]
        i2 = indices[0]
        sideFacesInner[-2] = [i1,n+n2-1,i2]
        sideFacesInner[-1] = [n+n2-1,n,i2]

        topFaces = 2*len(indices)*[[0,0,0]]

        for k in range(n2-1):
            topFaces[2*k] = [n+k,n+n2+k,n+k+1]
            topFaces[1 + 2*k] = [n+n2+k,n+n2+k+1,n+k+1]
        topFaces[-2] = [n+n2-1,n+2*n2-1,n]
        topFaces[-1] = [n+2*n2-1,n+n2,n]

        sideFacesOuter = 2*len(indices)*[[0,0,0]]

        for k in range(n2-1):
            sideFacesOuter[2*k] = [n+n2+k,n+2*n2+k,n+n2+k+1]
            sideFacesOuter[1 + 2*k] = [n+2*n2+k,n+2*n2+k+1,n+n2+k+1]
        sideFacesOuter[-2] = [n+2*n2-1,n+3*n2-1,n+n2]
        sideFacesOuter[-1] = [n+3*n2-1,n+2*n2,n+n2]

        botFacesOuter = 2*len(indices)*[[0,0,0]]
        for k in range(n2-1):
            i1 = indices[k]
            i2 = indices[k+1]
            botFacesOuter[2*k] = [n+2*n2+k,n+3*n2+i1,n+2*n2+k+1]
            botFacesOuter[1 + 2*k] = [n+3*n2+i1,n+3*n2+i2,n+2*n2+k+1]

        i1 = indices[-1]
        i2 = indices[0]
        botFacesOuter[2*k] = [n+3*n2-1,n+3*n2+i1,n+2*n2]
        botFacesOuter[1 + 2*k] = [n+3*n2+i1,n+3*n2+i2,n+2*n2]

        botFacesInner = [[f[0]+n+3*n2,f[1]+n+3*n2,f[2]+n+3*n2] for f in faces]

        moldVerts = verts + topVerticesInner + topVerticesOuter + bottomVerticesOuter + bottomVerticesInner

        moldFaces = faces + sideFacesInner + topFaces + sideFacesOuter + botFacesOuter + botFacesInner

        self.moldMesh = tm.Trimesh(vertices=moldVerts, faces=moldFaces)


    def makeInsoleAndMold(self):
        self.makeInsole()

        self.makeMold()

        #We flip the mold for better viewing
        zOffset = np.amax(self.soleMesh.vertices,axis=0)[2]
        rotMatrix = np.array([[-1, 0, 0],[0, 1, 0],[0, 0, -1]])
        for i in range(len(self.moldMesh.vertices)):
            v = np.array(self.moldMesh.vertices[i])
            self.moldMesh.vertices[i] = np.matmul(rotMatrix,v)
            self.moldMesh.vertices[i][2] += zOffset
        tm.repair.fix_normals(self.insoleMesh)
        tm.repair.fix_normals(self.moldMesh)

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
        path = QFileDialog.getSaveFileName(self, "Save Model", "../saved", "STL(*.stl)")[0]
        if path == "":
            print("No file chosen")
            return
        mesh.export(path)
