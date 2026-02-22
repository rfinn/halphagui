# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'maskWidgetv2.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(665, 555)
        self.gridLayout = QtWidgets.QGridLayout(Form)
        self.gridLayout.setObjectName("gridLayout")
        self.mainframe = QtWidgets.QFrame(Form)
        self.mainframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.mainframe.setFrameShadow(QtWidgets.QFrame.Raised)
        self.mainframe.setObjectName("mainframe")
        self.cutoutsLayout = QtWidgets.QGridLayout(self.mainframe)
        self.cutoutsLayout.setObjectName("cutoutsLayout")
        self.buttonWidget = QtWidgets.QWidget(self.mainframe)
        self.buttonWidget.setMouseTracking(True)
        self.buttonWidget.setObjectName("buttonWidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.buttonWidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_3 = QtWidgets.QLabel(self.buttonWidget)
        self.label_3.setObjectName("label_3")
        self.gridLayout_2.addWidget(self.label_3, 0, 0, 1, 1)
        self.boxSizeLineEdit = QtWidgets.QLineEdit(self.buttonWidget)
        self.boxSizeLineEdit.setObjectName("boxSizeLineEdit")
        self.gridLayout_2.addWidget(self.boxSizeLineEdit, 0, 1, 1, 1)
        self.mhelpButton = QtWidgets.QPushButton(self.buttonWidget)
        self.mhelpButton.setObjectName("mhelpButton")
        self.gridLayout_2.addWidget(self.mhelpButton, 0, 2, 1, 1)
        self.label = QtWidgets.QLabel(self.buttonWidget)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 1, 0, 1, 1)
        self.seThresholdLineEdit = QtWidgets.QLineEdit(self.buttonWidget)
        self.seThresholdLineEdit.setObjectName("seThresholdLineEdit")
        self.gridLayout_2.addWidget(self.seThresholdLineEdit, 1, 1, 1, 1)
        self.mquitButton = QtWidgets.QPushButton(self.buttonWidget)
        self.mquitButton.setObjectName("mquitButton")
        self.gridLayout_2.addWidget(self.mquitButton, 1, 2, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.buttonWidget)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 2, 0, 1, 1)
        self.seSNRLineEdit = QtWidgets.QLineEdit(self.buttonWidget)
        self.seSNRLineEdit.setObjectName("seSNRLineEdit")
        self.gridLayout_2.addWidget(self.seSNRLineEdit, 2, 1, 1, 1)
        self.cutoutsLayout.addWidget(self.buttonWidget, 0, 0, 1, 1)
        self.gridLayout.addWidget(self.mainframe, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label_3.setText(_translate("Form", "Box Size"))
        self.mhelpButton.setText(_translate("Form", "Help          "))
        self.label.setText(_translate("Form", "SE threshold (0=lots, 1=no deblend, def 0.05)"))
        self.mquitButton.setText(_translate("Form", "Quit"))
        self.label_2.setText(_translate("Form", "SE SNR (default=2)"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

