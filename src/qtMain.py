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
from loadWidget import *
from cutZWidget import *
from cutXWidget import *
from rotateXWidget import *
from insoleMoldExportWidget import *

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()

    def initUI(self):
        mw = MainTabsWidget()
        self.setCentralWidget(mw)
        self.setWindowTitle('Insole Tilter')
        self.resize(1300,1800)

        self.show()

class MainTabsWidget(QTabWidget):

    def __init__(self):
        super(MainTabsWidget, self).__init__()
        self.initUI()

    def initUI(self):
        loadWid = loadWidget(1)
        self.addTab(loadWid, "Load + Orient")
        cutZWid = cutZWidget(2)
        self.addTab(cutZWid, "Cut")
        self.setTabEnabled(1, False)
        rotXWid = rotateXWidget(3)
        self.addTab(rotXWid, "Rotate around X axis")
        self.setTabEnabled(2, False)
        stlsWid = insoleMoldExportWidget()
        self.addTab(stlsWid, "Compute and export insole and mold")
        self.setTabEnabled(3, False)

        loadWid.initCustomUI()
        cutZWid.initCustomUI()
        rotXWid.initCustomUI()
        stlsWid.initCustomUI()
