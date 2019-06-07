from baseWidget import *

class cutZWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(cutZWidget, self).__init__(nextTab)

    def initCustomUI(self):

        self.planeZMesh = None
        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)

        self.doCutZ = False
        cutCkBx = QCheckBox("Cut part of the model on Z axis")
        cutCkBx.stateChanged.connect(self.toggleCutZ)
        paramLayout.addWidget(cutCkBx)
        self.cutZLocation = .5
        cutLabel = QLabel("Choose where the cut should be located:")
        paramLayout.addWidget(cutLabel)
        cutSb = QDoubleSpinBox()
        cutSb.setRange(0,1)
        cutSb.setSingleStep(.01)
        cutSb.setValue(0.5)
        cutSb.setDecimals(2)
        cutSb.valueChanged.connect(self.setCutZLocation)
        paramLayout.addWidget(cutSb)
        self.showPlaneCutZ = False
        cutPlaneCkBx = QCheckBox("Show cut location")
        cutPlaneCkBx.stateChanged.connect(self.toggleShowPlaneCutZ)
        paramLayout.addWidget(cutPlaneCkBx)

        self.initDoneButtonAndShader(paramLayout)
        pass

    def start(self, transmittedData):
        self.soleMesh = transmittedData
        self.soleMeshWithZCut = copy.deepcopy(self.soleMesh)
        self.recompute = True
        self.displayMesh()

    def displayMesh(self):
        if self.recompute:
            self.resetNextTabs()
            mesh = None
            if self.doCutZ:
                self.cutZ()
                mesh = self.soleMeshWithZCut
            else:
                mesh = self.soleMesh
            self.recompute = False
            self.toTransmit = mesh
        glmesh = gl.MeshData(vertexes=mesh.vertices, faces=mesh.faces)
        if not self.meshItem is None:
            self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.meshItem)


    def displayPlaneCutZ(self):
        if not self.displayZPlane:
            return
        z = np.amax(self.soleMesh.vertices,axis=0)[2]
        xpos = np.amax(self.soleMesh.vertices,axis=0)[0]
        ypos = np.amax(self.soleMesh.vertices,axis=0)[1]
        xneg = np.amin(self.soleMesh.vertices,axis=0)[0]
        yneg = np.amin(self.soleMesh.vertices,axis=0)[1]
        zmin = np.amin(self.soleMesh.vertices,axis=0)[2]

        zToCut = zmin + (z - zmin) * self.cutZLocation
        planeVerts = np.array([[xpos,ypos,zToCut],[xpos,yneg,zToCut],[xneg,yneg,zToCut],[xneg,ypos,zToCut]])
        planeFaces = np.array([[0,1,2],[2,3,0]])

        glmesh = gl.MeshData(vertexes=planeVerts, faces=planeFaces)
        if self.planeZMesh:
            self.view.removeItem(self.planeZMesh)
        self.planeZMesh = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.planeZMesh)


    def cutZ(self):
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

        self.soleMeshWithZCut = tm.Trimesh(vertices=verts, faces=faces)

    @pyqtSlot()
    def setCutZLocation(self):
        sb = self.sender()
        self.cutZLocation = sb.value()
        self.recompute = True
        self.displayMesh()
        self.displayPlaneCutZ()

    @pyqtSlot()
    def toggleCutZ(self):
        cb = self.sender()
        self.doCutZ = cb.isChecked()
        self.recompute = True
        self.displayMesh()

    pyqtSlot()
    def toggleShowPlaneCutZ(self):
        cb = self.sender()

        if cb.isChecked():
            self.displayZPlane = True
            self.displayPlaneCutZ()
        else:
            self.displayZPlane = False
            if self.planeZMesh:
                self.view.removeItem(self.planeZMesh)
                self.planeZMesh = None
