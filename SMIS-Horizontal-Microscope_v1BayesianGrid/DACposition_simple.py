# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DACposition_simple.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.setWindowModality(QtCore.Qt.WindowModal)
        Form.resize(360, 220)
        Form.setMinimumSize(QtCore.QSize(360, 220))
        Form.setMaximumSize(QtCore.QSize(360, 220))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8("Logo/3107b825-90be-4ae6-8b9f-d4a29a72419f.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Form.setWindowIcon(icon)
        self.groupBox = QtGui.QGroupBox(Form)
        self.groupBox.setGeometry(QtCore.QRect(10, 10, 341, 201))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.formLayoutWidget = QtGui.QWidget(self.groupBox)
        self.formLayoutWidget.setGeometry(QtCore.QRect(10, 20, 321, 71))
        self.formLayoutWidget.setObjectName(_fromUtf8("formLayoutWidget"))
        self.formLayout = QtGui.QFormLayout(self.formLayoutWidget)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(self.formLayoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label)
        self.radioCamera = QtGui.QRadioButton(self.formLayoutWidget)
        self.radioCamera.setObjectName(_fromUtf8("radioCamera"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.radioCamera)
        self.radioSch = QtGui.QRadioButton(self.formLayoutWidget)
        self.radioSch.setObjectName(_fromUtf8("radioSch"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.radioSch)
        self.line = QtGui.QFrame(self.formLayoutWidget)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.line)
        self.formLayoutWidget_2 = QtGui.QWidget(self.groupBox)
        self.formLayoutWidget_2.setGeometry(QtCore.QRect(10, 90, 321, 61))
        self.formLayoutWidget_2.setObjectName(_fromUtf8("formLayoutWidget_2"))
        self.formLayout_2 = QtGui.QFormLayout(self.formLayoutWidget_2)
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.label_2 = QtGui.QLabel(self.formLayoutWidget_2)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_2)
        self.radioIR = QtGui.QRadioButton(self.formLayoutWidget_2)
        self.radioIR.setObjectName(_fromUtf8("radioIR"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.radioIR)
        self.radioRamanFL = QtGui.QRadioButton(self.formLayoutWidget_2)
        self.radioRamanFL.setObjectName(_fromUtf8("radioRamanFL"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.FieldRole, self.radioRamanFL)
        self.btnSend = QtGui.QPushButton(self.groupBox)
        self.btnSend.setGeometry(QtCore.QRect(180, 160, 151, 31))
        self.btnSend.setObjectName(_fromUtf8("btnSend"))

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Microscope Setup", None))
        self.groupBox.setTitle(_translate("Form", "Set Current Configuration", None))
        self.label.setText(_translate("Form", "DAC position", None))
        self.radioCamera.setText(_translate("Form", "Camera", None))
        self.radioSch.setText(_translate("Form", "Schwarzschild", None))
        self.label_2.setText(_translate("Form", "Microscope mode", None))
        self.radioIR.setText(_translate("Form", "Infrared", None))
        self.radioRamanFL.setText(_translate("Form", "Raman / Fluorescence", None))
        self.btnSend.setText(_translate("Form", "Send Configuration", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Form = QtGui.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

