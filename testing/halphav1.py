# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'halpha.v1.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1266, 550)
        # this is the main window
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")

        ###############################################
        # left panel that contains the mosaic image
        # QFrame with verticalLayout
        ###############################################
        self.mainframe = QtWidgets.QFrame(self.centralwidget)
        self.mainframe.setMinimumSize(QtCore.QSize(512, 0))
        self.mainframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.mainframe.setFrameShadow(QtWidgets.QFrame.Raised)
        self.mainframe.setObjectName("mainframe")

        # left panel, vertical layout
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.mainframe)
        self.verticalLayout_2.setObjectName("verticalLayout_2")

        # place holder for ginga CanvasView
        #self.widget = QtWidgets.QWidget(self.mainframe)
        #self.widget.setMinimumSize(QtCore.QSize(480, 400))
        #self.widget.setMaximumSize(QtCore.QSize(16777215, 542))
        #self.widget.setObjectName("widget")
        #self.verticalLayout_2.addWidget(self.widget, 0, QtCore.Qt.AlignLeft)

        # button bar on left frame
        # widget with Horizontal Box layout
        self.button_bar = QtWidgets.QWidget(self.mainframe)
        self.button_bar.setObjectName("button_bar")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.button_bar)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.wopen = QtWidgets.QPushButton(self.button_bar)
        self.wopen.setObjectName("wopen")
        self.horizontalLayout_2.addWidget(self.wopen)
        self.wmark = QtWidgets.QPushButton(self.button_bar)
        self.wmark.setObjectName("wmark")
        self.horizontalLayout_2.addWidget(self.wmark)
        self.wgalid = QtWidgets.QComboBox(self.button_bar)
        self.wgalid.setObjectName("wgalid")
        self.horizontalLayout_2.addWidget(self.wgalid)
        self.wclear = QtWidgets.QPushButton(self.button_bar)
        self.wclear.setObjectName("wclear")
        self.horizontalLayout_2.addWidget(self.wclear)
        self.verticalLayout_2.addWidget(self.button_bar, 0, QtCore.Qt.AlignBottom)
        self.horizontalLayout.addWidget(self.mainframe)

        ###############################################
        # middle panel
        # groupbox with GridLayout
        ###############################################
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.profiles_r = QtWidgets.QGraphicsView(self.groupBox)
        self.profiles_r.setObjectName("profiles_r")
        self.gridLayout.addWidget(self.profiles_r, 6, 0, 1, 1)
        self.cutoutsizeslider = QtWidgets.QSlider(self.groupBox)
        self.cutoutsizeslider.setOrientation(QtCore.Qt.Horizontal)
        self.cutoutsizeslider.setObjectName("cutoutsizeslider")
        self.gridLayout.addWidget(self.cutoutsizeslider, 3, 0, 1, 1)
        self.cshalpha = QtWidgets.QGraphicsView(self.groupBox)
        self.cshalpha.setObjectName("cshalpha")
        self.gridLayout.addWidget(self.cshalpha, 2, 1, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 1, 1, 1)
        self.filterratioslider = QtWidgets.QSlider(self.groupBox)
        self.filterratioslider.setOrientation(QtCore.Qt.Horizontal)
        self.filterratioslider.setObjectName("filterratioslider")
        self.gridLayout.addWidget(self.filterratioslider, 3, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.groupBox)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 4, 1, 1, 1)
        self.profiles_ha = QtWidgets.QGraphicsView(self.groupBox)
        self.profiles_ha.setObjectName("profiles_ha")
        self.gridLayout.addWidget(self.profiles_ha, 6, 1, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.rcutout = QtWidgets.QGraphicsView(self.groupBox)
        self.rcutout.setObjectName("rcutout")
        self.gridLayout.addWidget(self.rcutout, 2, 0, 1, 1)
        self.horizontalLayout.addWidget(self.groupBox)

        ###############################################
        # right panel
        # groupbox with verticalLayout
        ###############################################
        # proceeding widgets are added to the verticalLayout
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_3 = QtWidgets.QLabel(self.groupBox_2)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)

        # place holder for mask image
        self.maskimage = QtWidgets.QGraphicsView(self.groupBox_2)
        self.maskimage.setObjectName("maskimage")
        self.verticalLayout.addWidget(self.maskimage)

        # button for editing the mask
        self.emask = QtWidgets.QPushButton(self.groupBox_2)
        self.emask.setObjectName("emask")
        self.verticalLayout.addWidget(self.emask)

        # GraphicsView - placeholder for
        # normalized profiles
        self.profiles_norm = QtWidgets.QGraphicsView(self.groupBox_2)
        self.profiles_norm.setObjectName("profiles_norm")
        self.verticalLayout.addWidget(self.profiles_norm)

        # why is this added to horizontalLayout
        self.horizontalLayout.addWidget(self.groupBox_2)

        # add centralwidget to the main window
        MainWindow.setCentralWidget(self.centralwidget)

        ###############################################
        # set up menu bars
        # QMenuBar
        ###############################################
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1266, 22))
        self.menubar.setObjectName("menubar")

        # File menu
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")

        # session menu - to make sure everthing is saved
        self.menuSession = QtWidgets.QMenu(self.menubar)
        self.menuSession.setObjectName("menuSession")

        # NSA menu - to provide path to NSA catalog
        self.menuNSA_Catalog = QtWidgets.QMenu(self.menubar)
        self.menuNSA_Catalog.setObjectName("menuNSA_Catalog")
        MainWindow.setMenuBar(self.menubar)

        
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionSave = QtWidgets.QAction(MainWindow)
        self.actionSave.setObjectName("actionSave")
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionsave = QtWidgets.QAction(MainWindow)
        self.actionsave.setObjectName("actionsave")

        # add options to File menu
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)

        # add options to session menu
        self.menuSession.addAction(self.actionsave)
        self.menuSession.addSeparator()

        # add to menubar
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuNSA_Catalog.menuAction())
        self.menubar.addAction(self.menuSession.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.wopen.setText(_translate("MainWindow", "Open File"))
        self.wmark.setText(_translate("MainWindow", "Mark Galaxies"))
        self.wclear.setText(_translate("MainWindow", "Clear Canvas"))
        self.groupBox.setTitle(_translate("MainWindow", "panel2"))
        self.label.setText(_translate("MainWindow", "R-band "))
        self.label_2.setText(_translate("MainWindow", "CS H-alpha"))
        self.label_5.setText(_translate("MainWindow", "adjust filter ratio"))
        self.label_4.setText(_translate("MainWindow", "adjust cutout size"))
        self.groupBox_2.setTitle(_translate("MainWindow", "panel3"))
        self.label_3.setText(_translate("MainWindow", "Mask"))
        self.emask.setText(_translate("MainWindow", "Edit Mask"))
        self.menuFile.setTitle(_translate("MainWindow", "Image Path"))
        self.menuSession.setTitle(_translate("MainWindow", "Session"))
        self.menuNSA_Catalog.setTitle(_translate("MainWindow", "NSA_Catalog"))
        self.actionOpen.setText(_translate("MainWindow", "R-band coadd"))
        self.actionSave.setText(_translate("MainWindow", "Halpha coadd"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionsave.setText(_translate("MainWindow", "save"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
