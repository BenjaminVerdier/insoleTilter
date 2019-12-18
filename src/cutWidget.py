from baseWidget import *

import splipy as sp
import math


class cutWidget(baseWidget):

    def __init__(self, nextTab = -1):
        super(cutWidget, self).__init__(nextTab)
        self.horValChange = 5
        self.vertValChange = math.radians(5)

    def initCustomUI(self):

        self.planeZMesh = None
        self.curvePlot = None
        self.ctrlPtsPlot = None
        self.curCtrlPtSphere = None
        #Buttons and stuff, to be modified in child class
        paramWidget = QWidget()
        paramLayout = QVBoxLayout()
        paramWidget.setLayout(paramLayout)
        paramWidget.setMaximumHeight(425)
        self.mainLayout.addWidget(paramWidget)


        #Btn to move cam
        camLayout = QHBoxLayout()
        camBtn = QPushButton("Cam Left")
        camBtn.pressed.connect(self.moveCamLeft)
        camLayout.addWidget(camBtn)

        camBtn = QPushButton("Cam Up")
        camBtn.pressed.connect(self.moveCamUp)
        camLayout.addWidget(camBtn)

        camBtn = QPushButton("Cam Right")
        camBtn.pressed.connect(self.moveCamRight)
        camLayout.addWidget(camBtn)

        paramLayout.addLayout(camLayout)


        cutLabel = QLabel("Choose where the cut should be located:")
        paramLayout.addWidget(cutLabel)

        #Adding/Removing control points
        addRemLayout = QHBoxLayout()
        self.remCtptBtn = QPushButton("-")
        self.remCtptBtn.pressed.connect(self.remCtpt)
        addRemLayout.addWidget(self.remCtptBtn)
        self.addCtptBtn = QPushButton("+")
        self.addCtptBtn.pressed.connect(self.addCtpt)
        addRemLayout.addWidget(self.addCtptBtn)
        paramLayout.addLayout(addRemLayout)

        #We start with 8 control points
        ctPtsLayout = QHBoxLayout()
        self.currentControlPoint = 0
        self.minControlPoints = 4

        btnLayout = QVBoxLayout()
        self.ctptLabel = QLabel("Control Point : 1")
        ctPtsLayout.addWidget(self.ctptLabel)


        self.btnPrev = QPushButton("<<")
        self.btnPrev.setEnabled(False)
        self.btnPrev.pressed.connect(self.prevButton)
        ctPtsLayout.addWidget(self.btnPrev)

        self.btnLeft = QPushButton("<")
        self.btnLeft.setMinimumSize(5,5)
        self.btnLeft.resize(25,25)
        self.btnLeft.val = -self.horValChange
        self.btnLeft.pressed.connect(self.edgesMoved)
        ctPtsLayout.addWidget(self.btnLeft)

        btnUp = QPushButton("^")
        btnUp.setMinimumSize(5,5)
        btnUp.resize(25,25)
        btnUp.val = self.vertValChange
        btnUp.pressed.connect(self.ctrlPtsMoved)
        btnUp.setMinimumHeight(30)
        btnLayout.addWidget(btnUp)

        btnDown = QPushButton("v")
        btnDown.setMinimumSize(5,5)
        btnDown.resize(25,25)
        btnDown.val = -self.vertValChange
        btnDown.pressed.connect(self.ctrlPtsMoved)
        btnDown.setMinimumHeight(30)
        btnLayout.addWidget(btnDown)

        ctPtsLayout.addLayout(btnLayout)

        self.btnRight = QPushButton(">")
        self.btnRight.setMinimumSize(5,5)
        self.btnRight.resize(25,25)
        self.btnRight.val = self.horValChange
        self.btnRight.pressed.connect(self.edgesMoved)
        ctPtsLayout.addWidget(self.btnRight)

        self.btnNext = QPushButton(">>")
        self.btnNext.pressed.connect(self.nextButton)
        ctPtsLayout.addWidget(self.btnNext)

        paramLayout.addLayout(ctPtsLayout)

        self.doCutZ = False

        #Btn to recompute spline
        splineLayout = QHBoxLayout()
        splineBtn = QPushButton("Compute Full Spline")
        splineBtn.pressed.connect(self.recomputeSpline)
        splineLayout.addWidget(splineBtn)

        splineBtn = QPushButton("Compute Simple Spline")
        splineBtn.pressed.connect(self.recomputeSplineSimple)
        splineLayout.addWidget(splineBtn)

        paramLayout.addLayout(splineLayout)


        #Btn to cut
        cutBtn = QPushButton("Cut")
        cutBtn.pressed.connect(self.doCut)
        paramLayout.addWidget(cutBtn)

        #Btn to undo cut
        cutBtn = QPushButton("Undo Cut")
        cutBtn.pressed.connect(self.undoCut)
        paramLayout.addWidget(cutBtn)

        self.initDoneButtonAndShader(paramLayout)


    def start(self, transmittedData):
        self.soleMesh = transmittedData
        self.soleMeshWithZCut = copy.deepcopy(self.soleMesh)
        self.recompute = True
        verts = [np.array(x) for x in self.soleMesh.vertices]
        minX = np.amin(verts,axis=0)[0]
        maxX = np.amax(verts,axis=0)[0]
        self.footLength = maxX - minX
        self.first_point_x = minX + 0.75*self.footLength
        self.last_point_x = self.first_point_x - 1.51*self.footLength #there's a bug when they're totally symmetric for an odd number of points
        self.z = np.amax(verts,axis=0)[2]
        self.controlpoints = 8*[[0,0,self.z/3]]
        self.controlpoints[3] = [0,0,self.z/4]
        self.controlpoints[4] = [0,0,self.z/4]
        self.spline = None

        self.displayMesh()
        self.initSpline()

    def displayMesh(self):
        if self.recompute:
            self.resetNextTabs()
            mesh = None
            if self.doCutZ:
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
        self.displaySpline()
        self.displayControlPoints()


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
        minX = np.amin(verts,axis=0)[0]

        minY = np.amin(verts,axis=0)[1]
        maxY = np.amax(verts,axis=0)[1]
        def point_sampler(n):
            pos_increment = (self.first_point_x-self.last_point_x)/(n-1)
            points = []
            for i in range(n):
                x_position = self.first_point_x - i*pos_increment
                if x_position - minX > 0:
                    points += [[x_position,maxY,self.controlpoints[i][2]]]
                else:
                    points += [[minX - (x_position - minX),minY,self.controlpoints[i][2]]]
            return points

        controlpoints = point_sampler(len(self.controlpoints))
        contour = []
        for pt in controlpoints:
            while not self.soleMesh.ray.intersects_any([[pt[0],0,pt[2]]],[[0,pt[1],0]])[0]:
                pt[2] -= 1
            contour += [self.soleMesh.ray.intersects_location([[pt[0],0,pt[2]]],[[0,pt[1],0]])[0][0]]

        self.controlpoints = contour
        self.displayControlPoints()

    def displayControlPoints(self):
        if self.ctrlPtsPlot:
            self.view.removeItem(self.ctrlPtsPlot)
        col = np.ones((len(self.controlpoints),4))
        col[self.currentControlPoint] = [0,0,1.,1.]
        self.ctrlPtsPlot = gl.GLScatterPlotItem(pos=np.array(self.controlpoints),size=10,color=col)
        self.view.addItem(self.ctrlPtsPlot)

        if self.curCtrlPtSphere:
            self.view.removeItem(self.curCtrlPtSphere)
        glmesh = gl.MeshData.sphere(rows=10, cols=20)
        colors = np.ones((glmesh.faceCount(), 4), dtype=float)
        colors[:,0] = 0
        colors[:,1] = 0
        glmesh.setFaceColors(colors)
        self.curCtrlPtSphere = gl.GLMeshItem(meshdata=glmesh)
        self.curCtrlPtSphere.translate(*self.controlpoints[self.currentControlPoint])
        self.view.addItem(self.curCtrlPtSphere)

    def doSpline(self,samples = 300):

        def get_spline_points(cont_pts, order = 4, sampling = 150):
            n_control_points = len(cont_pts)
            kn = [0] * order + list(range(1, n_control_points - order + 1)) + [n_control_points - order + 1] * order
            basis = sp.BSplineBasis(order=order, knots=kn)
            curve = sp.Curve(basis, cont_pts)
            t = np.linspace(0,n_control_points - order + 1,sampling)
            x = curve.evaluate(t)
            return x

        spline_pts_not_projected = get_spline_points(self.controlpoints, sampling=samples)
        #spline_pts_projected = [[0,0,z]]
        spline_pts_projected = []
        for pt in spline_pts_not_projected:
            initZ = pt[2]
            while not self.soleMesh.ray.intersects_any([[0,0,self.z]],[[pt[0],pt[1],pt[2]-self.z]])[0]:
                pt[2] -= 1
            pjt = self.soleMesh.ray.intersects_location([[0,0,self.z]],[[pt[0],pt[1],pt[2]-self.z]])[0][0]
            spline_pts_projected += [[pjt[0],pjt[1],initZ]]
        #After that, need to project x to insole, so same as before with the rays, but without the while loop (hopefully)

        self.spline = spline_pts_projected

        self.displaySpline()

    def displaySpline(self):
        if self.spline == None:
            return
        if self.curvePlot:
            self.view.removeItem(self.curvePlot)
        self.curvePlot = gl.GLLinePlotItem(pos=np.array(self.spline),width=5,color=[1,0,0,1])
        self.view.addItem(self.curvePlot)


    def cutTop(self):
        #Top part
        partToCut = [[0,0,self.z]]
        partToCut += [[4*self.spline[0][0]+200, 4*self.spline[0][1],4*self.spline[0][2] - 3*self.z]]
        n = len(self.spline) -1
        for pt in self.spline:
            partToCut += [[4*pt[0], 4*pt[1],4*pt[2] - 3*self.z]]
        partToCut += [[4*self.spline[n][0]+200, 4*self.spline[n][1],4*self.spline[n][2] - 3*self.z]]
        n = len(partToCut)
        faces = []
        for i in range(1,n-1):
            faces += [[i,i+1,0]]
        faces += [[n-1,1,0]]
        for i in range(len(partToCut)):
            partToCut += [[partToCut[i][0], partToCut[i][1],partToCut[i][2] + 10*self.z]]
        for i in range(n+1,2*n-1):
            faces += [[i,i+1,n]]
        faces += [[2*n-1,n+1,n]]
        for i in range(1,n-1):
            faces += [[i,i+n,i+1]]
            faces += [[i+n+1,i+1,i+n]]
        faces += [[n-1,2*n-1,1]]
        faces += [[n+1,1,2*n-1]]
        """
        glmesh = gl.MeshData(vertexes=np.array(partToCut), faces=np.array(faces))

        if self.cubeMesh:
            self.view.removeItem(self.cubeMesh)
        self.cubeMesh = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.cubeMesh)
        """
        rem = tm.Trimesh(vertices=partToCut, faces=faces)
        self.soleMeshWithZCut = self.soleMesh.difference(rem)
        verts = list(self.soleMeshWithZCut.vertices)
        faces = list(self.soleMeshWithZCut.faces)
        vertToRemove = -1
        for i in range(len(verts)):
            if abs(verts[i][0]) < 1 and abs(verts[i][1]) < 1 and abs(verts[i][2]-self.z) < 1:
                vertToRemove = i
                break
        if vertToRemove > -1:
            del verts[vertToRemove]
            facesToDel = []
            for i in range(len(faces)):
                if vertToRemove in faces[i]:
                    facesToDel.append(i)

            for i in sorted(facesToDel, reverse=True):
                del faces[i]

            for f in self.soleMeshWithZCut.faces:
                if f[0] > vertToRemove:
                    f[0] -= 1
                if f[1] > vertToRemove:
                    f[1] -= 1
                if f[2] > vertToRemove:
                    f[2] -= 1

            self.soleMeshWithZCut = tm.Trimesh(verts,faces)

    def cutFront(self):
        #Front part
        cubeVerts = []

        ctp = self.controlpoints[0].copy()
        a = math.sqrt(ctp[1]*ctp[1] + ctp[0]*ctp[0])
        b = self.z - ctp[2]
        alpha = math.atan(b/a)
        ctp[2] -= a * math.tan(alpha-self.vertValChange) - b

        if self.soleMesh.ray.intersects_any([[0,0,self.z]],[[ctp[0],ctp[1],ctp[2]-self.z]])[0]:
            ctp = self.soleMesh.ray.intersects_location([[0,0,self.z]],[[ctp[0],ctp[1],ctp[2]-self.z]])[0][0]
        pt1 = np.array(ctp)

        ctp = self.controlpoints[len(self.controlpoints)-1].copy()
        a = math.sqrt(ctp[1]*ctp[1] + ctp[0]*ctp[0])
        b = self.z - ctp[2]
        alpha = math.atan(b/a)
        ctp[2] -= a * math.tan(alpha-self.vertValChange) - b

        if self.soleMesh.ray.intersects_any([[0,0,self.z]],[[ctp[0],ctp[1],ctp[2]-self.z]])[0]:
            ctp = self.soleMesh.ray.intersects_location([[0,0,self.z]],[[ctp[0],ctp[1],ctp[2]-self.z]])[0][0]
        pt2 = np.array(ctp)

        seg = pt1 - pt2
        left_side = seg + pt1
        right_side = pt2 - seg
        cubeVerts += [[right_side[0],right_side[1], - 100]]
        cubeVerts += [[right_side[0],right_side[1], self.z + 100]]
        cubeVerts += [[left_side[0],left_side[1], - 100]]
        cubeVerts += [[left_side[0],left_side[1], self.z + 100]]
        cubeVerts += [[right_side[0]+200,right_side[1], - 100]]
        cubeVerts += [[right_side[0]+200,right_side[1], self.z + 100]]
        cubeVerts += [[left_side[0]+200,left_side[1], - 100]]
        cubeVerts += [[left_side[0]+200,left_side[1], self.z + 100]]

        cubeFaces = [[0,1,2],[3,2,1],[4,5,6],[7,6,5],[0,1,4],[5,4,1],[2,3,6],[7,6,3],[1,3,5],[7,5,3],[0,2,4],[6,4,2]]
        cube = tm.Trimesh(cubeVerts,cubeFaces)
        self.soleMeshWithZCut = self.soleMeshWithZCut.difference(cube)
        """
        glmesh = gl.MeshData(vertexes=np.array(cubeVerts), faces=np.array(cubeFaces))

        if self.cubeMesh:
            self.view.removeItem(self.cubeMesh)
        self.cubeMesh = gl.GLMeshItem(meshdata=glmesh, shader=self.shader, glOptions=self.glOptions)
        self.view.addItem(self.cubeMesh)
        """

    def UpdateBtnLabel(self):
        t = "Control Point : " + str(self.currentControlPoint +1)
        self.ctptLabel.setText(t)

    def UpdateLeftRightBtns(self):
        if self.currentControlPoint == 0 or self.currentControlPoint == len(self.controlpoints) - 1:
            self.btnLeft.setEnabled(True)
            self.btnRight.setEnabled(True)
        else:
            self.btnLeft.setEnabled(False)
            self.btnRight.setEnabled(False)

        if self.currentControlPoint == len(self.controlpoints) - 1:
            self.btnNext.setEnabled(False)
        else:
            self.btnNext.setEnabled(True)

        if self.currentControlPoint == 0:
            self.btnPrev.setEnabled(False)
        else:
            self.btnPrev.setEnabled(True)

    @pyqtSlot()
    def ctrlPtsMoved(self):
        btn = self.sender()
        ctp = self.controlpoints[self.currentControlPoint].copy()
        a = math.sqrt(ctp[1]*ctp[1] + ctp[0]*ctp[0])
        b = self.z - ctp[2]
        alpha = math.atan(b/a)
        ctp[2] -= a * math.tan(alpha-btn.val) - b

        if self.soleMesh.ray.intersects_any([[ctp[0],0,self.z]],[[0,ctp[1],ctp[2]-self.z]])[0]:
            ctp = self.soleMesh.ray.intersects_location([[ctp[0],0,self.z]],[[0,ctp[1],ctp[2]-self.z]])[0][0]
        self.controlpoints[self.currentControlPoint] = ctp
        self.displayControlPoints()

    @pyqtSlot()
    def edgesMoved(self):
        btn = self.sender()
        if self.currentControlPoint == 0:
            self.first_point_x += btn.val
        else:
            self.last_point_x += -btn.val
        self.initSpline()


    @pyqtSlot()
    def doCut(self):
        self.doCutZ = True
        self.recompute = True
        self.cutTop()
        self.cutFront()
        self.displayMesh()


    @pyqtSlot()
    def undoCut(self):
        self.doCutZ = False
        self.recompute = True
        self.displayMesh()

    @pyqtSlot()
    def recomputeSpline(self):
        self.doSpline()


    @pyqtSlot()
    def recomputeSplineSimple(self):
        self.doSpline(25)


    @pyqtSlot()
    def moveCamRight(self):
        verts = [np.array(x) for x in self.soleMesh.vertices]
        maxz = np.amax(verts,axis=0)[2]
        minz = np.amax(verts,axis=0)[2]
        maxx = np.amax(verts,axis=0)[0]
        self.view.opts['center'].setZ(maxz/2)
        self.view.setCameraPosition(distance=2*maxx,elevation = 0, azimuth = 270)


    @pyqtSlot()
    def moveCamLeft(self):
        verts = [np.array(x) for x in self.soleMesh.vertices]
        maxz = np.amax(verts,axis=0)[2]
        minz = np.amax(verts,axis=0)[2]
        maxx = np.amax(verts,axis=0)[0]
        self.view.opts['center'].setZ(maxz/2)
        self.view.setCameraPosition(distance=2*maxx,elevation = 0, azimuth = 90)


    @pyqtSlot()
    def moveCamUp(self):
        verts = [np.array(x) for x in self.soleMesh.vertices]
        maxz = np.amax(verts,axis=0)[2]
        minz = np.amax(verts,axis=0)[2]
        maxx = np.amax(verts,axis=0)[0]
        self.view.opts['center'].setZ(maxz/2)
        self.view.setCameraPosition(distance=2*maxx,elevation = 90, azimuth = 180)


    @pyqtSlot()
    def nextButton(self):
        self.currentControlPoint += 1
        if self.currentControlPoint == 1:
            self.btnPrev.setEnabled(True)
        self.UpdateBtnLabel()
        self.UpdateLeftRightBtns()
        self.displayControlPoints()


    @pyqtSlot()
    def prevButton(self):
        self.currentControlPoint -= 1
        if self.currentControlPoint == len(self.controlpoints) - 2:
            self.btnNext.setEnabled(True)
        self.UpdateBtnLabel()
        self.UpdateLeftRightBtns()
        self.displayControlPoints()

    @pyqtSlot()
    def remCtpt(self):
        if self.currentControlPoint == 0:
            del self.controlpoints[1]
        elif self.currentControlPoint == len(self.controlpoints)-1:
            del self.controlpoints[-2]
            self.currentControlPoint -= 1
        else:
            del self.controlpoints[self.currentControlPoint]
            self.currentControlPoint -= 1
        if len(self.controlpoints) == self.minControlPoints:
            self.remCtptBtn.setEnabled(False)
        self.UpdateBtnLabel()
        self.UpdateLeftRightBtns()
        self.initSpline()

    @pyqtSlot()
    def addCtpt(self):
        if self.currentControlPoint == len(self.controlpoints)-1:
            newPt = self.controlpoints[-2].copy()
            newPt[2] = (newPt[2] + self.controlpoints[-1][2])/2
            self.controlpoints = self.controlpoints[:-1] + [newPt] + [self.controlpoints[-1]]
        else:
            newPt = self.controlpoints[self.currentControlPoint].copy()
            newPt[2] = (newPt[2] + self.controlpoints[self.currentControlPoint+1][2])/2
            self.controlpoints = self.controlpoints[:-1] + [newPt] + [self.controlpoints[-1]]
        self.currentControlPoint += 1
        if len(self.controlpoints) == self.minControlPoints + 1:
            self.remCtptBtn.setEnabled(True)
        self.UpdateBtnLabel()
        self.UpdateLeftRightBtns()
        self.initSpline()
