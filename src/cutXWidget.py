from baseWidget import *

class cutXWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(cutXWidget, self).__init__(nextTab)

    def initCustomUI(self):

        self.planeXMesh = None
        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)

        self.doCutX = False
        cutCkBx = QCheckBox("Cut part of the model on X axis")
        cutCkBx.stateChanged.connect(self.toggleCutX)
        paramLayout.addWidget(cutCkBx)
        self.cutXLocation = .5
        cutLabel = QLabel("Choose where the cut should be located:")
        paramLayout.addWidget(cutLabel)
        cutSb = QDoubleSpinBox()
        cutSb.setRange(0,1)
        cutSb.setSingleStep(.01)
        cutSb.setValue(0.5)
        cutSb.setDecimals(2)
        cutSb.valueChanged.connect(self.setCutXLocation)
        paramLayout.addWidget(cutSb)
        self.showPlaneCutX = False
        cutPlaneCkBx = QCheckBox("Show cut location")
        cutPlaneCkBx.stateChanged.connect(self.toggleShowPlaneCutX)
        paramLayout.addWidget(cutPlaneCkBx)

        self.initDoneButtonAndShader(paramLayout)
        pass

    def start(self, transmittedData):
        self.soleMesh = transmittedData
        self.soleMeshWithXCut = copy.deepcopy(self.soleMesh)
        self.recompute = True
        self.displayMesh()

    def displayMesh(self):
        if self.recompute:
            self.resetNextTabs()
            mesh = None
            if self.doCutX:
                self.cutX()
                mesh = self.soleMeshWithXCut
            else:
                mesh = self.soleMesh
            self.recompute = False
            self.toTransmit = mesh
        glmesh = gl.MeshData(vertexes=mesh.vertices, faces=mesh.faces)
        if not self.meshItem is None:
            self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.meshItem)


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


    def cutX(self):
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
        self.soleMeshWithXCut = self.soleMesh.difference(cube)

    @pyqtSlot()
    def setCutXLocation(self):
        sb = self.sender()
        self.cutXLocation = sb.value()
        self.recompute = True
        self.displayMesh()
        self.displayPlaneCutX()

    @pyqtSlot()
    def toggleCutX(self):
        cb = self.sender()
        self.doCutX = cb.isChecked()
        self.recompute = True
        self.displayMesh()

    pyqtSlot()
    def toggleShowPlaneCutX(self):
        cb = self.sender()

        if cb.isChecked():
            self.displayXPlane = True
            self.displayPlaneCutX()
        else:
            self.displayXPlane = False
            if self.planeXMesh:
                self.view.removeItem(self.planeXMesh)
                self.planeXMesh = None
