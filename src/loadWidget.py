import math
import os

from baseWidget import *

class loadWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(loadWidget, self).__init__(nextTab)

    def initCustomUI(self):
        self.toTransmit = None

        #Adding straight blue line on X axis to help with placement
        self.view.addItem(gl.GLLinePlotItem(pos=np.array([[0,0,0],[500,0,0]]),width=50,color=[0,0,1,1]))
        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)

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

        self.initDoneButtonAndShader(paramLayout)

        self.loadBasicMesh()


    def loadBasicMesh(self):
        path = str(os.path.dirname(os.path.realpath(__file__)))
        self.mesh = tm.load(path+"\\..\\resources\\scan.obj")
        #We process the mesh so that it is centered and lays flat on the xy plane.
        #We then copy it t a display mesh that will sustain all the tranformation we apply
        self.mesh.vertices -= self.mesh.centroid
        transforms, probs = self.mesh.compute_stable_poses()
        self.mesh.apply_transform(transforms[np.argmax(probs)])
        offset = np.amin(self.mesh.vertices,axis=0)[2]
        self.mesh.vertices -= [0,0,offset]
        self.soleMesh = copy.deepcopy(self.mesh)
        self.recompute = True
        self.displayMesh()

    def displayMesh(self):
        if self.recompute:
            self.resetNextTabs()
            self.rotateZ()
            self.recompute = False
            self.toTransmit = self.soleMesh
        glmesh = gl.MeshData(vertexes=self.soleMesh.vertices, faces=self.soleMesh.faces)
        if not self.meshItem is None:
            self.view.removeItem(self.meshItem)
        self.meshItem = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.meshItem)

    def rotateZ(self):
        agl = math.radians(self.zRotAngle)
        c = math.cos(agl)
        s = math.sin(agl)
        rotMatrix = np.array([[c,-s,0],[s,c,0],[0,0,1]])
        for i in range(len(self.mesh.vertices)):
            v = np.array(self.mesh.vertices[i])
            self.soleMesh.vertices[i] = np.matmul(rotMatrix,v)

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
        self.displayMesh()

    @pyqtSlot()
    def toggleInvertNormals(self):
        cb = self.sender()
        self.invertNormals = cb.isChecked()
        self.mesh.invert()
        self.soleMesh.invert()
        self.displayMesh()

    @pyqtSlot()
    def setZRotAngle(self):
        sb = self.sender()
        self.zRotAngle = sb.value()
        self.recompute = True
        self.displayMesh()
