from baseWidget import *

from enum import Enum, auto
import math


class Rotation(Enum):
    FULL = auto()
    STRETCH = auto()

class rotateXWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(rotateXWidget, self).__init__(nextTab)

    def initCustomUI(self):

        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)


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

        self.initDoneButtonAndShader(paramLayout)

    def start(self, transmittedData):
        self.soleMesh = transmittedData
        self.rotatedMesh = copy.deepcopy(self.soleMesh)
        self.recompute = True
        self.displayMesh()

    def displayMesh(self):
        if self.recompute:
            self.resetNextTabs()
            self.rotate()
            self.recompute = False
            self.toTransmit = self.rotatedMesh
        glmesh = gl.MeshData(vertexes=self.rotatedMesh.vertices, faces=self.rotatedMesh.faces)
        if not self.meshItem is None:
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
                self.rotatedMesh.vertices[i] = np.matmul(rotMatrix,v)
        else:
            for i in range(len(self.soleMesh.vertices)):
                rotMatrix = np.array([[1,0,0],[0,math.cos(agls[i]),-math.sin(agls[i])],[0,math.sin(agls[i]),math.cos(agls[i])]])
                v = np.array(self.soleMesh.vertices[i])
                self.rotatedMesh.vertices[i][2] = np.matmul(rotMatrix,v)[2]
        offset = np.amin(self.rotatedMesh.vertices,axis=0)[2] - 2
        self.rotatedMesh.vertices -= [0,0,offset]


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
    def setLinDescent(self):
        sb = self.sender()
        self.linearDescentPortion = sb.value()
        self.recompute = True
        self.displayMesh()
