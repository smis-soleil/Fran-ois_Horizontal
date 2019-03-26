# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Micromode.ui'
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

class Ui_Microscopemode(object):
    def setupUi(self, Microscopemode):
        Microscopemode.setObjectName(_fromUtf8("Microscopemode"))
        Microscopemode.setWindowModality(QtCore.Qt.WindowModal)
        Microscopemode.resize(880, 480)
        self.label = QtGui.QLabel(Microscopemode)
        self.label.setGeometry(QtCore.QRect(10, 10, 880, 31))
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName(_fromUtf8("label"))
        self.layoutWidget = QtGui.QWidget(Microscopemode)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 60, 886, 390))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.layoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, -1)
        self.horizontalLayout.setSpacing(60)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.btnIRmode = QtGui.QPushButton(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.btnIRmode.setFont(font)
        self.btnIRmode.setMouseTracking(True)
        self.btnIRmode.setAutoFillBackground(False)
        self.btnIRmode.setText(_fromUtf8(""))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8("../../Desktop/Horizon 4.0/Script/Logo/b2-90.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btnIRmode.setIcon(icon)
        self.btnIRmode.setIconSize(QtCore.QSize(400, 380))
        self.btnIRmode.setObjectName(_fromUtf8("btnIRmode"))
        self.horizontalLayout.addWidget(self.btnIRmode)
        self.btnRamanmode = QtGui.QPushButton(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.btnRamanmode.setFont(font)
        self.btnRamanmode.setMouseTracking(True)
        self.btnRamanmode.setAutoFillBackground(False)
        self.btnRamanmode.setText(_fromUtf8(""))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8("../../Desktop/Horizon 4.0/Script/Logo/b3-90.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btnRamanmode.setIcon(icon1)
        self.btnRamanmode.setIconSize(QtCore.QSize(400, 380))
        self.btnRamanmode.setObjectName(_fromUtf8("btnRamanmode"))
        self.horizontalLayout.addWidget(self.btnRamanmode)
        self.label_2 = QtGui.QLabel(Microscopemode)
        self.label_2.setGeometry(QtCore.QRect(160, 90, 141, 31))
        self.label_2.setAlignment(QtCore.Qt.AlignCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(Microscopemode)
        self.label_3.setGeometry(QtCore.QRect(602, 90, 201, 31))
        self.label_3.setAlignment(QtCore.Qt.AlignCenter)
        self.label_3.setObjectName(_fromUtf8("label_3"))

        self.retranslateUi(Microscopemode)
        QtCore.QObject.connect(self.btnIRmode, QtCore.SIGNAL(_fromUtf8("clicked()")), Microscopemode.close)
        QtCore.QObject.connect(self.btnRamanmode, QtCore.SIGNAL(_fromUtf8("clicked()")), Microscopemode.close)
        QtCore.QMetaObject.connectSlotsByName(Microscopemode)

    def retranslateUi(self, Microscopemode):
        Microscopemode.setWindowTitle(_translate("Microscopemode", "Microscope mode", None))
        self.label.setText(_translate("Microscopemode", "<html><head/><body><p><span style=\" font-size:18pt;\">Select current configuration</span></p></body></html>", None))
        self.label_2.setText(_translate("Microscopemode", "<html><head/><body><p><span style=\" font-size:18pt;\">Infrared</span></p></body></html>", None))
        self.label_3.setText(_translate("Microscopemode", "<html><head/><body><p><span style=\" font-size:18pt;\">Raman / FL</span></p></body></html>", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Microscopemode = QtGui.QWidget()
    ui = Ui_Microscopemode()
    ui.setupUi(Microscopemode)
    Microscopemode.show()
    sys.exit(app.exec_())

