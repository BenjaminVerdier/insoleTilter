from baseWidget import *

from enum import Enum, auto
import sys

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
