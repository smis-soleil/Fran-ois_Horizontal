# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DACposition.ui'
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
        Form.resize(880, 480)
        self.label = QtGui.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(0, 10, 880, 31))
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName(_fromUtf8("label"))
        self.widget = QtGui.QWidget(Form)
        self.widget.setGeometry(QtCore.QRect(20, 60, 841, 390))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout.setSpacing(60)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.btnSchwa = QtGui.QPushButton(self.widget)
        self.btnSchwa.setText(_fromUtf8(""))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8("Logo/c4-90.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btnSchwa.setIcon(icon)
        self.btnSchwa.setIconSize(QtCore.QSize(400, 380))
        self.btnSchwa.setObjectName(_fromUtf8("btnSchwa"))
        self.horizontalLayout.addWidget(self.btnSchwa)
        self.btnCamera = QtGui.QPushButton(self.widget)
        self.btnCamera.setText(_fromUtf8(""))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8("Logo/c5-90.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btnCamera.setIcon(icon1)
        self.btnCamera.setIconSize(QtCore.QSize(400, 380))
        self.btnCamera.setObjectName(_fromUtf8("btnCamera"))
        self.horizontalLayout.addWidget(self.btnCamera)

        self.retranslateUi(Form)
        QtCore.QObject.connect(self.btnCamera, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.close)
        QtCore.QObject.connect(self.btnSchwa, QtCore.SIGNAL(_fromUtf8("clicked()")), Form.close)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "DAC position", None))
        self.label.setText(_translate("Form", "<html><head/><body><p><span style=\" font-size:18pt;\">Define current DAC position</span></p></body></html>", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Form = QtGui.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

