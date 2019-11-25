from baseWidget import *

import splipy as sp
import math


class cutZWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(cutZWidget, self).__init__(nextTab)

    def initCustomUI(self):

        self.planeZMesh = None
        self.curvePlot = None
        self.ctrlPtsPlot = None
        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(300)
        self.mainLayout.addWidget(paramWidget)

        cutLabel = QLabel("Choose where the cut should be located:")
        paramLayout.addWidget(cutLabel)
        #We start with 8 control points
        ctPtsLayout = QHBoxLayout()
        self.controlpoints = 8*[[0,0,0]]
        for i in range(8):
            name = "Ctrl Pt "+ str(i+1)
            ctpLabel = QLabel(name)
            ctPtsLayout.addWidget(ctpLabel)
            btnUp = QPushButton("^")
            btnUp.setMinimumSize(5,5)
            btnUp.resize(25,25)
            btnUp.ctpIndex = i
            btnUp.val = 5
            btnUp.pressed.connect(self.ctrlPtsMoved)
            ctPtsLayout.addWidget(btnUp)
            btnDown = QPushButton("v")
            btnDown.setMinimumSize(5,5)
            btnDown.resize(25,25)
            btnDown.ctpIndex = i
            btnDown.val = -5
            btnDown.pressed.connect(self.ctrlPtsMoved)
            ctPtsLayout.addWidget(btnDown)
        paramLayout.addLayout(ctPtsLayout)

        self.doCutZ = False

        #Btn to recompute spline
        splineBtn = QPushButton("Recompute Spline")
        splineBtn.pressed.connect(self.recomputeSpline)
        paramLayout.addWidget(splineBtn)

        self.initDoneButtonAndShader(paramLayout)


    def start(self, transmittedData):
        self.soleMesh = transmittedData
        self.soleMeshWithZCut = copy.deepcopy(self.soleMesh)
        self.recompute = True
        self.displayMesh()
        self.initSpline()

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
        z = np.amax(verts,axis=0)[2]#Add 10?
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
        #tm.repair.fix_normals(self.soleMeshWithZCut)

    def initSpline(self):
        verts = [np.array(x) for x in self.soleMesh.vertices]
        self.z = np.amax(verts,axis=0)[2]
        xpos = np.amax(verts,axis=0)[0]
        def circle_sampler(n):
            angle_increment = (3 * math.pi) / (2*n)
            points = []
            for i in range(n):
                angle = (math.pi/4) + i*angle_increment
                points += [[xpos*math.cos(angle),xpos*math.sin(angle),self.z]]
            return points

        controlpoints = circle_sampler(8)
        contour = []
        for pt in controlpoints:
            while not self.soleMesh.ray.intersects_any([[0,0,self.z]],[[pt[0],pt[1],pt[2]-self.z]])[0]:
                pt[2] -= 1
            contour += [self.soleMesh.ray.intersects_location([[0,0,self.z]],[[pt[0],pt[1],pt[2]-self.z]])[0][0]]

        self.controlpoints = contour
        self.displayControlPoints()
        self.doSpline()

    def displayControlPoints(self):
        if self.ctrlPtsPlot:
            self.view.removeItem(self.ctrlPtsPlot)
        self.ctrlPtsPlot = gl.GLScatterPlotItem(pos=np.array(self.controlpoints),size=10,color=[0,0,1,1])
        self.view.addItem(self.ctrlPtsPlot)

    def doSpline(self):

        def get_spline_points(cont_pts, order = 4, sampling = 150):
            n_control_points = len(cont_pts)
            kn = [0] * order + list(range(1, n_control_points - order + 1)) + [n_control_points - order + 1] * order
            basis = sp.BSplineBasis(order=order, knots=kn)
            curve = sp.Curve(basis, cont_pts)
            t = np.linspace(0,n_control_points - order + 1,sampling)
            x = curve.evaluate(t)
            return x

        spline_pts_not_projected = get_spline_points(self.controlpoints, sampling=300)
        #spline_pts_projected = [[0,0,z]]
        spline_pts_projected = []
        for pt in spline_pts_not_projected:
            initZ = pt[2]
            while not self.soleMesh.ray.intersects_any([[0,0,self.z]],[[pt[0],pt[1],pt[2]-self.z]])[0]:
                pt[2] -= 1
            pjt = self.soleMesh.ray.intersects_location([[0,0,self.z]],[[pt[0],pt[1],pt[2]-self.z]])[0][0]
            spline_pts_projected += [[pjt[0],pjt[1],initZ]]
        #After that, need to project x to insole, so same as before with the rays, but without the while loop (hopefully)

        if self.curvePlot:
            self.view.removeItem(self.curvePlot)
        self.curvePlot = gl.GLLinePlotItem(pos=np.array(spline_pts_projected),width=5,color=[1,0,0,1])
        self.view.addItem(self.curvePlot)

        #faces = []
        #for i in range(1,len(spline_pts_projected)-1):
        #    faces += [[i,i+1,0]]
        #faces += [[len(spline_pts_projected)-1,1,0]]

        #glmesh = gl.MeshData(vertexes=np.array(spline_pts_projected), faces=np.array(faces))
        #if self.planeZMesh:
        #    self.view.removeItem(self.planeZMesh)
        #self.planeZMesh = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        #self.view.addItem(self.planeZMesh)


    @pyqtSlot()
    def ctrlPtsMoved(self):
        btn = self.sender()
        ctp = self.controlpoints[btn.ctpIndex].copy()
        a = math.sqrt(ctp[1]*ctp[1] + ctp[0]*ctp[0])
        b = self.z - ctp[2]
        alpha = math.atan(b/a)
        ctp[2] -= a * math.tan(alpha-math.radians(btn.val)) - b

        if self.soleMesh.ray.intersects_any([[0,0,self.z]],[[ctp[0],ctp[1],ctp[2]-self.z]])[0]:
            ctp = self.soleMesh.ray.intersects_location([[0,0,self.z]],[[ctp[0],ctp[1],ctp[2]-self.z]])[0][0]
        self.controlpoints[btn.ctpIndex] = ctp
        self.displayControlPoints()


    @pyqtSlot()
    def recomputeSpline(self):
        self.doSpline()
