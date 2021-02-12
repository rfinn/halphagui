# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'maskGui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_maskWindow(object):
    def setupUi(self, maskWindow):
        maskWindow.setObjectName("maskWindow")
        maskWindow.resize(665, 460)
        self.centralwidget = QtWidgets.QWidget(maskWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.mainframe = QtWidgets.QFrame(self.centralwidget)
        self.mainframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.mainframe.setFrameShadow(QtWidgets.QFrame.Raised)
        self.mainframe.setObjectName("mainframe")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.mainframe)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.widget = QtWidgets.QWidget(self.mainframe)
        self.widget.setObjectName("widget")
        self.gridLayout_2.addWidget(self.widget, 0, 1, 1, 1)
        self.widget_2 = QtWidgets.QWidget(self.mainframe)
        self.widget_2.setObjectName("widget_2")
        self.gridLayout_2.addWidget(self.widget_2, 2, 1, 1, 1)
        self.buttonWidget = QtWidgets.QWidget(self.mainframe)
        self.buttonWidget.setMouseTracking(True)
        self.buttonWidget.setObjectName("buttonWidget")
        self.maskButtonLayout = QtWidgets.QGridLayout(self.buttonWidget)
        self.maskButtonLayout.setObjectName("maskButtonLayout")
        self.mcenterButton = QtWidgets.QPushButton(self.buttonWidget)
        self.mcenterButton.setObjectName("mcenterButton")
        self.maskButtonLayout.addWidget(self.mcenterButton, 0, 3, 1, 1)
        self.mremoveButton = QtWidgets.QPushButton(self.buttonWidget)
        self.mremoveButton.setObjectName("mremoveButton")
        self.maskButtonLayout.addWidget(self.mremoveButton, 0, 0, 1, 1)
        self.mquitButton = QtWidgets.QPushButton(self.buttonWidget)
        self.mquitButton.setObjectName("mquitButton")
        self.maskButtonLayout.addWidget(self.mquitButton, 5, 1, 1, 1)
        self.msethresholdButton = QtWidgets.QPushButton(self.buttonWidget)
        self.msethresholdButton.setObjectName("msethresholdButton")
        self.maskButtonLayout.addWidget(self.msethresholdButton, 4, 1, 1, 1)
        self.mboxsizeButton = QtWidgets.QPushButton(self.buttonWidget)
        self.mboxsizeButton.setObjectName("mboxsizeButton")
        self.maskButtonLayout.addWidget(self.mboxsizeButton, 0, 2, 1, 1)
        self.msesnrButton = QtWidgets.QPushButton(self.buttonWidget)
        self.msesnrButton.setObjectName("msesnrButton")
        self.maskButtonLayout.addWidget(self.msesnrButton, 4, 0, 1, 1)
        self.maddButton = QtWidgets.QPushButton(self.buttonWidget)
        self.maddButton.setObjectName("maddButton")
        self.maskButtonLayout.addWidget(self.maddButton, 0, 1, 1, 1)
        self.msaveButton = QtWidgets.QPushButton(self.buttonWidget)
        self.msaveButton.setObjectName("msaveButton")
        self.maskButtonLayout.addWidget(self.msaveButton, 5, 0, 1, 1)
        self.gridLayout_2.addWidget(self.buttonWidget, 2, 0, 1, 1)
        self.widget_3 = QtWidgets.QWidget(self.mainframe)
        self.widget_3.setObjectName("widget_3")
        self.gridLayout_2.addWidget(self.widget_3, 1, 1, 1, 1)
        self.cutouts = QtWidgets.QFrame(self.mainframe)
        self.cutouts.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.cutouts.setFrameShadow(QtWidgets.QFrame.Raised)
        self.cutouts.setObjectName("cutouts")
        self.cutoutsLayout = QtWidgets.QGridLayout(self.cutouts)
        self.cutoutsLayout.setObjectName("cutoutsLayout")
        self.dummyWidget = QtWidgets.QWidget(self.cutouts)
        self.dummyWidget.setObjectName("dummyWidget")
        self.cutoutsLayout.addWidget(self.dummyWidget, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.cutouts, 0, 0, 2, 1)
        self.gridLayout.addWidget(self.mainframe, 0, 1, 1, 1)
        maskWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(maskWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 665, 22))
        self.menubar.setObjectName("menubar")
        maskWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(maskWindow)
        self.statusbar.setObjectName("statusbar")
        maskWindow.setStatusBar(self.statusbar)

        self.retranslateUi(maskWindow)
        QtCore.QMetaObject.connectSlotsByName(maskWindow)

    def retranslateUi(self, maskWindow):
        _translate = QtCore.QCoreApplication.translate
        maskWindow.setWindowTitle(_translate("maskWindow", "MainWindow"))
        self.mcenterButton.setText(_translate("maskWindow", "Off Center Object"))
        self.mremoveButton.setText(_translate("maskWindow", "Remove Object"))
        self.mquitButton.setText(_translate("maskWindow", "Quit"))
        self.msethresholdButton.setText(_translate("maskWindow", "SE threshold"))
        self.mboxsizeButton.setText(_translate("maskWindow", "Change Box Size"))
        self.msesnrButton.setText(_translate("maskWindow", "SE SNR"))
        self.maddButton.setText(_translate("maskWindow", "Add Mask Pixels"))
        self.msaveButton.setText(_translate("maskWindow", "Save Mask"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    maskWindow = QtWidgets.QMainWindow()
    ui = Ui_maskWindow()
    ui.setupUi(maskWindow)
    maskWindow.show()
    sys.exit(app.exec_())

