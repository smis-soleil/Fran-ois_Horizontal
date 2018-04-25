# --------------------------------------------------------------------

# Editeur: Matthias Glachant (glachant.matthias@gmail.com  )
# Date: 7/6/2016

# --------------------------------------------------------------------

import logging
import sys  # We need sys so that we can pass argv to QApplication
import time

import numpy as np
import serial
from PyQt4 import QtGui  # Import the PyQt4 module we'll need
import seabreeze.spectrometers as sb    # on importe le module seabreeze Ocean Optics

import Advancedparameters_ui
import DACposition_ui
import DACposition_simple
import Horizontal_ui  # This file holds our MainWindow and all design related things
import Micromode_ui

# create logger
logger = logging.getLogger('HM App:')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('HMApp_full.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s\t%(name)s\t(%(levelname)s)\t:\t%(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.info('--- Starting Application ---')


###----------DEF_FUNCTIONS--------------------------------------

def execution(ser, commande):
    logger.debug('Executing serial command: ' + commande)
    # Code to send a command to the controller
    # send the character to the device
    ser.write(bytes(commande + '\r\n', 'UTF-8'))
    out = ''
    # let's wait one second before reading output (let's give device time to answer)
    time.sleep(0.1)  ### is this necessary? maybe that is why the program is slow
    # Formation of the answer to something understable
    liste = []
    reponse = []
    while (ser.inWaiting() > 0):
        out = str(ser.read(1))
        liste.append(out)
    if (len(liste) > 0):
        del liste[len(liste) - 1]
        for car in liste:
            reponse.append(car[2])  # Answer is in a list of caracters --> a simple word (rep)
    rep = str()
    for i in range(len(reponse)):
        rep += reponse[i]
    return (rep)


def initialisation(motorNb): #Initilisation des moteurs
    execution(ser, 'INIT'+ str(motorNb))
    return()

def signe_str(motorNB):  # Le moteur 5 est monté à l'envers (str)
    if motorNB == 5:
        return '-'
    return '+'

def signe_int(motorNB):   # Le moteur 5 est monté à l'envers (entier)
    if motorNB == 5:
        return -1
    return 1

def position(motorNb, position):                # On définit la position désiré pour le moteur N
    execution(ser, 'PSET' + str(motorNb) + '=' + str(position))  # Set target position respectively relative travel distance for an axis (Relat)
    return()

def positionvalue(motorNb):
    position_ini = execution(ser, '?CNT' + str(motorNb))     # ?CNT > Read out current position counter for an axis > Mesure de l'Offset
    x = ''
    for i in range(len(position_ini)):
        y = position_ini[i]
        x = x + str(y)
    return x

def convert_str_int(nb, n):       #On convertit une chaîne de caractère en entier
    if nb == '': return ''
    nb = float(nb)
    nb = int(nb*n)
    return nb

def velocity(motorNb, velocity):                 #On applique la vélocitée désirée
    execution(ser, 'VVEL' + str(motorNb) + '=' + str(velocity))
    return()

def acceleration(motorNb, acceleration):        #On applique l'accélération désirée
    execution(ser, 'ACC' + str(motorNb) + '=' + str(acceleration))
    return()

def deceleration(motorNb, deceleration):       #On applique la décélération désirée
    execution(ser, 'DACC' + str(motorNb) + '=' + str(deceleration))
    return()

#def move(motorNB):
    #execution(ser, 'PGO' + str (motorNB)) #PGO > Start positioning for an axis
    #return()

def stop(motorNb):                     # On arrête le moteur N
    execution(ser, 'STOP' + str(motorNb))
    return()

def mouve(motorNb, value, mode):
    execution(ser, str(mode) + str(motorNb))
    position(motorNb, value)
    execution(ser, 'PGO' + str (motorNb))
    return()

# Ruby pressure scales and calculation-----------
x0 = 694.22

def dewaele(x):
    a = 1920
    b = 9.61
    P = (a / b) * ((x / x0) ** b - 1)
    return P


def mao_hydro(x):
    a = 1876
    b = 10.71
    P = (a / b) * ((x / x0) ** b - 1)
    return P


def mao_no_hydro(x):
    a = 1904
    b = 5
    P = (a / b) * ((x / x0) ** b - 1)
    return P


## TODO: simple table and plotting in a separate tab?

dacpos = ''
mode = ''
XYZmotors = ''
allmotors = ''
Ycam = 0
Yschwa = 1500000
Speed = 1
c = 1  # X,Y conversion parameter
c_z = 1 # Z conversion parameter

#----------------Initialisation scan------------------------

devices = sb.list_devices()  # On cherche un spectrometre
spec = sb.Spectrometer(devices[0])  # On va utiliser le premier specromètre trouvé qui correspond au Hr4000

# ------------------------CLASSES------------------------------------------------------------
class AdvancedparametersWindow(QtGui.QDialog, Advancedparameters_ui.Ui_Advanced_parameters_window):
    def __init__(self, parent=None):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)

        self.btnSeeparam.clicked.connect(self.seeparam)
        self.btnSaveparam.clicked.connect(self.saveparam)
        self.btnposition.clicked.connect(self.positionval)
        self.btnplus.clicked.connect(self.mouvementplus)
        self.btnmoins.clicked.connect(self.mouvementmoins)

    def positionval(self):
        motorNb = self.comboBox.currentIndex() + 1
        sgn = signe_int(motorNb)
        self.position.setText(str(sgn * int(positionvalue(motorNb))))

    def mouvementplus(self):
        motorNb = self.comboBox.currentIndex() + 1
        sgn = signe_int(motorNb)
        value = sgn * int(self.step.text())

        mouve(motorNb, value, 'RELAT')
        while (execution(ser, '?VACT' + str(motorNb)) != '0'):
            time.sleep(0.02)
        self.position.setText(str(sgn * int(positionvalue(motorNb))))

    def mouvementmoins(self):
        motorNb = self.comboBox.currentIndex() + 1
        sgn = signe_int(motorNb)
        value = sgn * -1 * int(self.step.text())

        mouve(motorNb, value, 'RELAT')
        while (execution(ser, '?VACT' + str(motorNb)) != '0'):
            time.sleep(0.02)
        self.position.setText(str(sgn * int(positionvalue(motorNb))))


    def saveparam(self):
        xyzacc = self.xyzAccNew.text()
        xyzdec = self.xyzDecNew.text()
        xyzspeed = self.xyzSpeedNew.text()
        largeacc = self.largeAccNew.text()
        largedec = self.largeDecNew.text()
        largespeed = self.largeSpeedNew.text()

        if (xyzacc != ''):
            acceleration(4, xyzacc)
            acceleration(5, xyzacc)
            acceleration(6, xyzacc)
        if (xyzdec != ''):
            deceleration(4, xyzdec)
            deceleration(5, xyzdec)
            deceleration(6, xyzdec)
        if (xyzspeed != ''):
            velocity(4, xyzspeed)
            velocity(5, xyzspeed)
            velocity(6, xyzspeed)

        if (largeacc != ''):
            acceleration(1, largeacc)
            acceleration(2, largeacc)
            acceleration(3, largeacc)
        if (largedec != ''):
            deceleration(1, largedec)
            deceleration(2, largedec)
            deceleration(3, largedec)
        if (largespeed != ''):
            velocity(1, largespeed)
            velocity(2, largespeed)
            velocity(3, largespeed)

        for i in range(1, 7):
            execution(ser, 'SAVEAXPA' + str(i))

    def seeparam(self):
        self.xyzAccAct.setText(execution(ser, '?ACC4'))
        self.xyzDecAct.setText(execution(ser, '?DACC4'))
        self.xyzSpeedAct.setText(execution(ser, '?VVEL4'))
        self.largeAccAct.setText(execution(ser, '?ACC1'))
        self.largeDecAct.setText(execution(ser, '?DACC1'))
        self.largeSpeedAct.setText(execution(ser, '?VVEL1'))

# Ask for the current DAC position
class DACpositionWindow(QtGui.QDialog, DACposition_ui.Ui_Form):
    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)

        self.btnCamera.clicked.connect(self.cameramode)
        self.btnSchwa.clicked.connect(self.schwamode)

    def cameramode(self):
        global dacpos
        dacpos = 'Cam'

    def schwamode(self):
        global dacpos
        dacpos = 'Schwa'

# Ask for the current microscope mode
class MicromodeWindow(QtGui.QDialog, Micromode_ui.Ui_Microscopemode):
    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        self.setupUi(self)

        self.btnIRmode.clicked.connect(self.IRmode)
        self.btnRamanmode.clicked.connect(self.Ramanmode)

    def IRmode(self):
        global mode
        mode = 'IR'

    def Ramanmode(self):
        global mode
        mode = 'Raman'

class MainHorizontalWindow(QtGui.QMainWindow, Horizontal_ui.Ui_MainWindow):
    def __init__(self):
        # Explaining super is out of the scope of this article
        # So please google it if you're not familar with it
        # Simple reason why we use it here is that it allows us to
        # access variables, methods etc in the design.py file
        super(self.__class__, self).__init__()
        self.setupUi(self)

        self.btnStop.clicked.connect(self.stopall)
        self.radioButton_1.clicked.connect(self.statemotor1)
        self.radioButton_2.clicked.connect(self.statemotor2)
        self.radioButton_3.clicked.connect(self.statemotor3)
        self.radioButton_4.clicked.connect(self.statemotor4)
        self.radioButton_5.clicked.connect(self.statemotor5)
        self.radioButton_6.clicked.connect(self.statemotor6)
        self.btn_apply.clicked.connect(self.convert_unite)
        if c == 1:
            self.btn_step.setChecked(True)
            self.btn_micrometer.setChecked(False)
        if c == 40:
            self.btn_micrometer.setChecked(True)
            self.btn_step.setChecked(False)
        self.btn_apply.clicked.connect(self.syz_speed)
        if Speed == 1:
            self.btn_x1.setChecked(True)
            self.btn_x5.setChecked(False)
        if Speed == 5:
            self.btn_x5.setChecked(True)
            self.btn_x1.setChecked(False)
        self.btnxh.clicked.connect(self.mouvementXh)
        self.btnxb.clicked.connect(self.mouvementXb)
        self.btnxh.clicked.connect(self.Xvalue)
        self.btnxb.clicked.connect(self.Xvalue)
        self.btnyh.clicked.connect(self.mouvementYh)
        self.btnyb.clicked.connect(self.mouvementYb)
        self.btnyh.clicked.connect(self.Yvalue)
        self.btnyb.clicked.connect(self.Yvalue)
        self.btnzh.clicked.connect(self.mouvementZh)
        self.btnzb.clicked.connect(self.mouvementZb)
        self.btnzh.clicked.connect(self.Zvalue)
        self.btnzb.clicked.connect(self.Zvalue)
        self.btnMicromode.clicked.connect(self.micromode)
        self.btnMicromode.clicked.connect(self.DACposition)
        self.btnReadX.clicked.connect(self.Xvalue)
        self.btnReadX.clicked.connect(self.Yvalue)
        self.btnReadX.clicked.connect(self.Zvalue)
        self.btnSetPosition.clicked.connect(self.Xsetvalue)
        self.btnSetPosition.clicked.connect(self.Xvalue)
        self.btnSetPosition.clicked.connect(self.Ysetvalue)
        self.btnSetPosition.clicked.connect(self.Yvalue)
        self.btnSetPosition.clicked.connect(self.Zsetvalue)
        self.btnSetPosition.clicked.connect(self.Zvalue)
        self.btnSaveP1.clicked.connect(self.SaveP1)
        self.btnSaveP2.clicked.connect(self.SaveP2)
        self.btnSaveP3.clicked.connect(self.SaveP3)
        self.btnGoP1.clicked.connect(self.GoP1)
        self.btnGoP2.clicked.connect(self.GoP2)
        self.btnGoP3.clicked.connect(self.GoP3)
        self.btnGoP1.clicked.connect(self.Xvalue)
        self.btnGoP1.clicked.connect(self.Yvalue)
        self.btnGoP1.clicked.connect(self.Zvalue)
        self.btnGoP2.clicked.connect(self.Xvalue)
        self.btnGoP2.clicked.connect(self.Yvalue)
        self.btnGoP2.clicked.connect(self.Zvalue)
        self.btnGoP3.clicked.connect(self.Xvalue)
        self.btnGoP3.clicked.connect(self.Yvalue)
        self.btnGoP3.clicked.connect(self.Zvalue)
        self.calcul.clicked.connect(self.pressurecalcul)
        self.btnsettings.clicked.connect(self.openparameters)
        if dacpos == 'Schwa':
            self.radioButton_Schwa.setChecked(True)
        if dacpos == 'Cam':
            self.radioButton_Camera.setChecked(True)
        if mode == 'IR':
            self.radioButton_IR.setChecked(True)
        if mode == 'Raman':
            self.radioButton_Raman.setChecked(True)
        self.btnSYZ.clicked.connect(self.SYZstatus)
        self.btnAll.clicked.connect(self.statusAll)
        self.btnAll.setStyleSheet("QPushButton {background : green}")
        self.btnAll.setText('Turn all OFF')
        self.labelSYZ.setText('Turn SYZ OFF')
        self.labelSYZ.setStyleSheet("QLabel {color : green}")
        self.microstep.setText('Step')
        self.microstep.setStyleSheet("QLabel {color : red}")
        self.btnSetzero.clicked.connect(self.SetZero)
        #self.btnSetzero2.clicked.connect(self.SetZero)
        self.btnScan.clicked.connect(self.SetZero)
        self.btnScan.clicked.connect(self.Scan)



        self.btnPosText.clicked.connect(self.PosTextChange)

    def PosTextChange(self):
        print(self.btnSaveP1.text())
        if self.btnSaveP1.text() == 'BKG':
            self.btnSaveP1.setText('Save P1')
            self.btnSaveP2.setText('Save P2')
            self.btnSaveP3.setText('Save P3')

            self.btnGoP1.setText('Go to P1')
            self.btnGoP2.setText('Go to P2')
            self.btnGoP3.setText('Go to P3')

        else:
            self.btnSaveP1.setText('BKG')
            self.btnSaveP2.setText('Sample')
            self.btnSaveP3.setText('Ruby')

            self.btnGoP1.setText('BKG')
            self.btnGoP2.setText('Sample')
            self.btnGoP3.setText('Ruby')

    def openparameters(self):
        AdvancedparametersWindow().setModal(True)
        AdvancedparametersWindow().show()
        AdvancedparametersWindow().exec_()

    def GoP1(self):
        xval = convert_str_int(self.Position1X.text(), c)
        yval = convert_str_int(self.Position1Y.text(), c)
        zval = convert_str_int(self.Position1Z.text(), c_z)

        Position = [xval, -int(yval), zval]
        for i in range(3):
            j = i + 4
            mouve(j, Position[i], 'ABSOL')

    def GoP2(self):
        xval = convert_str_int(self.Position2X.text(), c)
        yval = convert_str_int(self.Position2Y.text(), c)
        zval = convert_str_int(self.Position2Z.text(), c_z)

        Position = [xval, -int(yval), zval]
        for i in range(3):
            j = i + 4
            mouve(j, Position[i], 'ABSOL')

    def GoP3(self):
        xval = convert_str_int(self.Position3X.text(), c)
        yval = convert_str_int(self.Position3Y.text(), c)
        zval = convert_str_int(self.Position3Z.text(), c_z)

        Position = [xval, -int(yval), zval]
        for i in range(3):
            j = i + 4
            mouve(j, Position[i], 'ABSOL')

    def SaveP1(self):
        self.Position1X.setText(str(int(execution(ser, '?CNT4')) / c))  # Motor 4
        self.Position1Y.setText(str((-1 / c) * int(execution(ser, '?CNT5'))))
        self.Position1Z.setText(str(int(execution(ser, '?CNT6')) / c_z))

    def SaveP2(self):
        self.Position2X.setText(str(int(execution(ser, '?CNT4')) / c))  # Motor 4
        self.Position2Y.setText(str((-1 / c) * int(execution(ser, '?CNT5'))))
        self.Position2Z.setText(str(int(execution(ser, '?CNT6')) / c_z))

    def SaveP3(self):
        self.Position3X.setText(str(int(execution(ser, '?CNT4')) / c))  # Motor 4
        self.Position3Y.setText(str((-1 / c) * int(execution(ser, '?CNT5'))))
        self.Position3Z.setText(str(int(execution(ser, '?CNT6')) / c_z))

# Set the X,Y,S to the specified values
    def Xsetvalue(self):
        xval = convert_str_int(self.SetpositionX.text(), c)
        mouve(4, xval, 'ABSOL')

    def Ysetvalue(self):
        yval = convert_str_int(self.SetpositionY.text(), c)
        mouve(5, -yval, 'ABSOL')

    def Zsetvalue(self):
        zval =  convert_str_int(self.SetpositionZ.text(), c_z)
        mouve(6, zval, 'ABSOL')

# Read the X,Y,S values
    def Xvalue(self):
        while (execution(ser, '?VACT4') != '0'):
            time.sleep(0.02)
        self.positionX.setText(str(int(execution(ser, '?CNT4')) / c))

    def Yvalue(self):
        while (execution(ser, '?VACT5') != '0'):
            time.sleep(0.02)
        self.positionY.setText(str((-1 / c) * int(execution(ser, '?CNT5'))))

    def Zvalue(self):
        while (execution(ser, '?VACT6') != '0'):
            time.sleep(0.02)
        self.positionZ.setText(str(int(execution(ser, '?CNT6')) / c_z))

# Cassegrain and DAC long translations
    def micromode(self):
        if self.radioButton_IR.isChecked() and mode == 'Raman':
            mouve(1, -1500000, 'RELAT')
            mouve(2, 1500000, 'RELAT')
            global mode
            mode = 'IR'
        if self.radioButton_Raman.isChecked() and mode == 'IR':
            mouve(1, 1500000, 'RELAT')
            mouve(2, -1500000, 'RELAT')
            global mode
            mode = 'Raman'

    def DACposition(self):
        if self.radioButton_Camera.isChecked() and dacpos == 'Schwa':
            mouve(3, 0, 'ABSOL')
            global dacpos
            dacpos = 'Cam'
        if self.radioButton_Schwa.isChecked() and dacpos == 'Cam':
            mouve(3, 1500000, 'ABSOL')
            global dacpos
            dacpos = 'Schwa'

# Plus/minus X,Y,S movements relative to the current position
    def mouvementXh(self):
        xval = convert_str_int(self.stepX.text(), c)
        mouve(4, xval, 'RELAT')
    def mouvementXb(self):
        xval = convert_str_int(self.stepX.text(), c)
        mouve(4, -xval, 'RELAT')

    def mouvementYh(self):
        yval = convert_str_int(self.stepY.text(), c)
        mouve(5, -yval, 'RELAT')
    def mouvementYb(self):
        yval = convert_str_int(self.stepY.text(), c)
        mouve(5, yval, 'RELAT')

    def mouvementZh(self):
        zval = convert_str_int(self.stepZ.text(), c_z)
        mouve(6, zval, 'RELAT')
    def mouvementZb(self):
        zval = convert_str_int(self.stepZ.text(), c_z)
        mouve(6, -zval, 'RELAT')

    def convert_unite(self):
        if self.btn_step.isChecked() and c == 20 and c_z == 40:
            self.Position1X.setText(str(convert_str_int(self.Position1X.text(), 20)))
            self.Position2X.setText(str(convert_str_int(self.Position2X.text(), 20)))
            self.Position3X.setText(str(convert_str_int(self.Position3X.text(), 20)))
            self.Position1Y.setText(str(convert_str_int(self.Position1Y.text(), 20)))
            self.Position2Y.setText(str(convert_str_int(self.Position2Y.text(), 20)))
            self.Position3Y.setText(str(convert_str_int(self.Position3Y.text(), 20)))
            self.Position1Z.setText(str(convert_str_int(self.Position1Z.text(), 40)))
            self.Position2Z.setText(str(convert_str_int(self.Position2Z.text(), 40)))
            self.Position3Z.setText(str(convert_str_int(self.Position3Z.text(), 40)))
            self.positionX.setText(str(convert_str_int(self.positionX.text(), 20)))
            self.positionY.setText(str(convert_str_int(self.positionY.text(), 20)))
            self.positionZ.setText(str(convert_str_int(self.positionZ.text(), 40)))
            self.SetpositionX.setText(str(convert_str_int(self.SetpositionX.text(), 20)))
            self.SetpositionY.setText(str(convert_str_int(self.SetpositionY.text(), 20)))
            self.SetpositionZ.setText(str(convert_str_int(self.SetpositionZ.text(), 40)))
            self.stepX.setText(str(convert_str_int(self.stepX.text(), 20)))
            self.stepY.setText(str(convert_str_int(self.stepY.text(), 20)))
            self.stepZ.setText(str(convert_str_int(self.stepZ.text(), 40)))
            global c, c_z
            c = 1
            c_z = 1
            self.microstep.setText('Step')
            self.microstep.setStyleSheet("QLabel {color : red}")
        if self.btn_micrometer.isChecked() and c == 1 and c_z == 1:
            self.Position1X.setText(str(convert_str_int(self.Position1X.text(), 0.05)))
            self.Position2X.setText(str(convert_str_int(self.Position2X.text(), 0.05)))
            self.Position3X.setText(str(convert_str_int(self.Position3X.text(), 0.05)))
            self.Position1Y.setText(str(convert_str_int(self.Position1Y.text(), 0.05)))
            self.Position2Y.setText(str(convert_str_int(self.Position2Y.text(), 0.05)))
            self.Position3Y.setText(str(convert_str_int(self.Position3Y.text(), 0.05)))
            self.Position1Z.setText(str(convert_str_int(self.Position1Z.text(), 0.025)))
            self.Position2Z.setText(str(convert_str_int(self.Position2Z.text(), 0.025)))
            self.Position3Z.setText(str(convert_str_int(self.Position3Z.text(), 0.025)))
            self.positionX.setText(str(convert_str_int(self.positionX.text(), 0.05)))
            self.positionY.setText(str(convert_str_int(self.positionY.text(), 0.05)))
            self.positionZ.setText(str(convert_str_int(self.positionZ.text(), 0.025)))
            self.SetpositionX.setText(str(convert_str_int(self.SetpositionX.text(), 0.05)))
            self.SetpositionY.setText(str(convert_str_int(self.SetpositionY.text(), 0.05)))
            self.SetpositionZ.setText(str(convert_str_int(self.SetpositionZ.text(), 0.025)))
            self.stepX.setText(str(convert_str_int(self.stepX.text(), 0.05)))
            self.stepY.setText(str(convert_str_int(self.stepY.text(), 0.05)))
            self.stepZ.setText(str(convert_str_int(self.stepZ.text(), 0.025)))
            global c, c_z
            c = 20
            c_z = 40
            self.microstep.setText('Micro')
            self.microstep.setStyleSheet("QLabel {color : red}")

    def syz_speed(self):
        if self.btn_x1.isChecked() and Speed == 5:
            for i in range(3, 6):
                acceleration(i + 1, 5)
                deceleration(i + 1, 5)
                velocity(i + 1, 2)
            global Speed
            Speed = 1
        if self.btn_x5.isChecked() and Speed == 1:
            for i in range(3, 6):
                acceleration(i + 1, 20)
                deceleration(i + 1, 20)
                velocity(i + 1, 8)
            global Speed
            Speed = 5

    def pressurecalcul(self):
        # temp = float(self.t_value.text())
        method = self.pressure_method.currentText()
        x = float(self.wl_value.text())
        pressure = 0
        if method == "Dewaele":
            pressure = dewaele(x)
        if method == "Mao hydro":
            pressure = mao_hydro(x)
        if method == "Mao no-hydro":
            pressure = mao_no_hydro(x)
        self.p_value.setText((str(round(pressure, 2))))

    def color(self):
        if not (self.radioButton_4.isChecked()):
            if not (self.radioButton_5.isChecked()):
                if not (self.radioButton_6.isChecked()):
                    execution(ser, 'JACC4=2')
                    self.labelSYZ.setText('Turn SYZ ON')
                    self.labelSYZ.setStyleSheet("QLabel {color : red}")
                    if not (self.radioButton_1.isChecked()):
                        if not (self.radioButton_2.isChecked()):
                            if not (self.radioButton_3.isChecked()):
                                execution(ser, 'JACC5=2')
                                self.btnAll.setStyleSheet("QPushButton {background : red}")
                                self.btnAll.setText('Turn all ON')
        if (self.radioButton_4.isChecked()):
            if (self.radioButton_5.isChecked()):
                if (self.radioButton_6.isChecked()):
                    execution(ser, 'JACC4=1')
                    self.labelSYZ.setText('Turn SYZ OFF')
                    self.labelSYZ.setStyleSheet("QLabel {color : green}")
                    if (self.radioButton_1.isChecked()):
                        if (self.radioButton_2.isChecked()):
                            if (self.radioButton_3.isChecked()):
                                execution(ser, 'JACC5=1')
                                self.btnAll.setStyleSheet("QPushButton {background : green}")
                                self.btnAll.setText('Turn all OFF')

    def SYZstatus(self):
        if (execution(ser, '?JACC4') == '1'):  # Turn them off
            execution(ser, 'MOFF4')
            execution(ser, 'MOFF5')
            execution(ser, 'MOFF6')
            self.radioButton_4.setChecked(False)
            self.radioButton_5.setChecked(False)
            self.radioButton_6.setChecked(False)
            status = 'off'
        if (execution(ser, '?JACC4') == '2'):  # Turn them on
            execution(ser, 'MON4')
            execution(ser, 'MON5')
            execution(ser, 'MON6')
            self.radioButton_4.setChecked(True)
            self.radioButton_5.setChecked(True)
            self.radioButton_6.setChecked(True)
            status = 'on'
        if (status == 'on'):
            execution(ser, 'JACC4=1')
            self.labelSYZ.setText('Turn SYZ OFF')
            self.labelSYZ.setStyleSheet("QLabel {color : green}")
        if (status == 'off'):
            execution(ser, 'JACC4=2')
            self.labelSYZ.setText('Turn SYZ ON')
            self.labelSYZ.setStyleSheet("QLabel {color : red}")

    def statusAll(self):
        if (execution(ser, '?JACC5') == '1'):  # Turn them all off
            execution(ser, 'MOFF1')
            execution(ser, 'MOFF2')
            execution(ser, 'STOP3')
            execution(ser, 'MOFF4')
            execution(ser, 'MOFF5')
            execution(ser, 'MOFF6')
            self.radioButton_1.setChecked(False)
            self.radioButton_2.setChecked(False)
            self.radioButton_3.setChecked(False)
            self.radioButton_4.setChecked(False)
            self.radioButton_5.setChecked(False)
            self.radioButton_6.setChecked(False)
            status = 'all_off'
        if (execution(ser, '?JACC5') == '2'):  # Turn them all on
            execution(ser, 'MON1')
            execution(ser, 'MON2')
            execution(ser, 'MON3')
            execution(ser, 'MON4')
            execution(ser, 'MON5')
            execution(ser, 'MON6')
            self.radioButton_1.setChecked(True)
            self.radioButton_2.setChecked(True)
            self.radioButton_3.setChecked(True)
            self.radioButton_4.setChecked(True)
            self.radioButton_5.setChecked(True)
            self.radioButton_6.setChecked(True)
            status = 'all_on'
        if (status == 'all_on'):
            execution(ser, 'JACC5=1')
            self.btnAll.setStyleSheet("QPushButton {background : green}")
            self.btnAll.setText('Turn all OFF')
            execution(ser, 'JACC4=1')
            self.labelSYZ.setText('Turn SYZ OFF')
            self.labelSYZ.setStyleSheet("QLabel {color : green}")
        if (status == 'all_off'):
            execution(ser, 'JACC5=2')
            self.btnAll.setStyleSheet("QPushButton {background : red}")
            self.btnAll.setText('Turn all ON')
            execution(ser, 'JACC4=2')
            self.labelSYZ.setText('Turn SYZ ON')
            self.labelSYZ.setStyleSheet("QLabel {color : red}")

    def statemotor1(self):
        if self.radioButton_1.isChecked():
            execution(ser, 'MON1')
        else: execution(ser, 'MOFF1')
        self.color()

    def statemotor2(self):
        if self.radioButton_2.isChecked():
            execution(ser, 'MON2')
        else: execution(ser, 'MOFF2')
        self.color()

    def statemotor3(self):
        if self.radioButton_3.isChecked():
            execution(ser, 'MON3')
        else: execution(ser, 'MOFF3')
        self.color()

    def statemotor4(self):
        if self.radioButton_4.isChecked():
            execution(ser, 'MON4')
        else: execution(ser, 'MOFF4')
        self.color()

    def statemotor5(self):
        if self.radioButton_5.isChecked():
            execution(ser, 'MON5')
        else: execution(ser, 'MOFF5')
        self.color()

    def statemotor6(self):
        if self.radioButton_6.isChecked():
            execution(ser, 'MON6')
        else: execution(ser, 'MOFF6')
        self.color()

    def SetZero(self):
        execution(ser, 'CRES4')    # CRESn -> reset current position for an axis
        execution(ser, 'CRES5')
        execution(ser, 'CRES6')
        self.btnReadX.click()

    def stopall(self):
        execution(ser, 'JACC5=1')
        self.statusAll()



#--------------Scan-------------------------------------------------------------
##pas mal!  j'ai ajouté un peu des commentaires avec ## ...

    def Scan(self):

        intTval = convert_str_int(self.intT.text(), 1000)   # On lit la valeur du temps d'ingration choisi par l'utilisateur
        spec.integration_time_micros(intTval)  # Temps d'intégration en micro secondes * 1000 -> ms
        diameterZ = convert_str_int(self.Di.text(), 40)
        diameterY = convert_str_int(self.Di.text(), 18)
        pasY = convert_str_int(self.StepC.text(), 20)



        # LOM = []
        LIM = []    
        LY = []     ## lst_y ou y_values ou autre?! minuscules pour les variables et noms un peu plus descriptifs...
                    ##si on respecte les conventions, le code devien plus lisible !
        LZ = []


        for w in range(0, diameterZ, 600):                     # Boucle pour le mouvement latéral
            for i in range(0, diameterY, pasY):                    #Boucle mouvement vertical haut
                # L = spec.wavelengths()  # Liste des Longueurs d'ondes
                 #l = spec.intensities()  # Liste des Intensités
                y = pasY                     # Valeur du pas en step    ## pourquoi on a besoin de cette def, une variable est suffisant il me semble..
                mouve(5, y, 'RELAT')                # On se déplace par rapport à la valeur précédente
                LIM.append(max(spec.intensities()))       # On stocke l'intensité max à chaque point
                Y = execution(ser, '?CNT' + str(5))        # On lit la coordonée en Y
                Z = execution(ser, '?CNT' + str(6))        # On lit la coordonée en Z
                LY.append(Y)           # On stocke la coordonée en Y
                LZ.append(Z)           # On stocke la coordonée en Z    ## dans cette boucle Z sera toujours le meme...on a besoin d'une liste ?

                #print('+z', w, i)
                print(max(spec.intensities()))     # Indication de l'intensité max à chaque tour

            z = 300              # Valeur du pas de Z en step
            mouve(6, z, 'RELAT')       # on se déplace par rapport à la valeur précédente
            for i in range(0, diameterY, pasY):             # Boucle mouvement veritcal bas ## on peut eviter cette deuxieme boucle si on crée une boucle qui part d'une cordonnée negative
                # L = spec.wavelengths()  # Liste des Longeurs d'ondes
                #l = spec.intensities()  # Liste des Intensités
                y = pasY
                mouve(5, -y, 'RELAT')
                LIM.append(max(spec.intensities()))
                Y = execution(ser, '?CNT' + str(5))
                Z = execution(ser, '?CNT' + str(6))
                LY.append(Y)
                LZ.append(Z)
                #print('-z', w, i)
                print(max(spec.intensities()))
            z = 300
            mouve(6, z, 'RELAT')

        cLZ = np.asarray(LZ)  ##on peut pas utiliser les array numpy depuis le debut ?
        cLY = np.asarray(LY)
        cLim = np.asarray(LIM)

        def zoom_scan():
            s = 0
            choice = QtGui.QMessageBox.question(self, 'Precision Scan', "Do you want to get a more accurate point?",
                                                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
            if choice == QtGui.QMessageBox.Yes and s != 1:
                mouve(4, -300, 'RELAT')
                mouve(5, -200, 'RELAT')
                mouve(6, -400, 'RELAT')
                self.Di.setText(str(20))
                self.intT.setText(str(100))
                self.StepC.setText(str(2))
                s = 1

            else:
                pass






        print("Max Intensity", max(LIM))      # On affiche l'intensité maximale parmi toutes les intensitées mesurées
        Z_m = LZ[np.array(LIM).argmax()]      # On cherche la coordonée z correspondante
        Y_m = LY[np.array(LIM).argmax()]      # On cherche la coordonnée y correpsondante
        Y_max = convert_str_int(Y_m,1)        # On convertit les valeurs en entier
        Z_max = convert_str_int(Z_m,1)
        print("Z Max :", Z_max)               # On affiche les coordonnées du ruby
        print("Y Max : -", Y_max)
        self.Yscan.setText('-' + str(Y_max))
        self.Zscan.setText(str(Z_max))
        self.Imax.setText(str(max(LIM)))
        mouve(5, Y_max,'ABSOL')               # On se positione sur le ruby à la fin du scan
        mouve(6, Z_max, 'ABSOL')
        zoom_scan()




        #mouve(4, 4000, 'RELAT')
        #mouve()





#--------------------------Scan UI---------------------------------------



# ------Open the windows-----------------------------------------------
def setup_old():
    logger.info('STEP: initial setup')
    app = QtGui.QApplication(sys.argv)  # A new instance of QApplication
    logger.info('\tDAC Position')
    form2 = DACpositionWindow()  # We set the form to be our ExampleApp (design)
    form2.show()
    logger.info('\tOptical Arrangement')
    form3 = MicromodeWindow()  # We set the form to be our ExampleApp (design)
    form3.show()
    app.exec_()  # and execute the app

def setup():
    logger.info('STEP: setup')
    app = QtGui.QApplication(sys.argv)  # A new instance of QApplication
    logger.info('\tDAC Position')
    form2 = DACpositionWindow()  # We set the form to be our ExampleApp (design)
    form2.show()
    logger.info('\tOptical Arrangement')
    form3 = MicromodeWindow()  # We set the form to be our ExampleApp (design)
    form3.show()
    app.exec_()  # and execute the app


def main():
    logger.info('STEP: main')
    logger.info('\tlaunching main window...')
    app = QtGui.QApplication(sys.argv)
    form = MainHorizontalWindow()  # We set the form to be our ExampleApp (design)
    form.show()
    app.exec_()  # and execute the app


def discover_and_connect():
    logger.info('STEP: discover_and_connect')
    logger.info('\tfind serial port and connect')
    # Code to choose the serial port used
    ## this could be done much more elegantyl. what is the way of finding out the # of com ports in a machine?
    # consider this
    # import serial.tools.list_ports
    #
    # ports = list(serial.tools.list_ports.comports())
    # for p in ports:
    #     print p
    # then just loop over these

    # choixport = 1
    # port = 0
    ser = serial.Serial  # why is this here?
    for port in range(50):  # (choixport == 1):
        logger.debug('\tChecking COM' + str(port + 1))
        try:
            # Open the port
            ser = serial.Serial(port='COM' + str(port + 1),
                                baudrate=9600,
                                stopbits=serial.STOPBITS_ONE,
                                bytesize=serial.EIGHTBITS)

            # try to get an answer from the PS90 controller. The answer would identify the system
            #  see ?ASTAT documentation
            longueur = len(execution(ser, '?ASTAT'))
            if (longueur == 9 and str(ser.isOpen()) == 'True'):  # this could be done without the "str()"
                choixport = 2
                logger.debug('\t...attempting to connect to COM' + str(port + 1) + ': success.')
                logger.info('\tConnected to COM' + str(port + 1) + '.')
                break
            if not (longueur == 9):
                logger.debug('\t...attempting to connect to COM' + str(port + 1) + ': failed')
                choixport = 1
        except:
            if (port < 50):
                port += 1
            else:
                logger.error('*** Connection error. Make sure USB cable is connected.')
                # input()
                # exit()
    if (port >= 49):
        logger.error(
            '*** Connection error, did not find PS90 controller below port<50. Make sure USB cable is connected.')
        ### would be nice to add a window here or at least a pause
        logger.error('--- Shutting down ---')
        exit()

    return (ser)
    # Emptying the buffer
    # ser.reset_input_buffer()
    # ser.reset_output_buffer()


def init_motors(ser):
    logger.info('STEP: init_motors')
    relat = ('RELAT1', 'RELAT2', 'RELAT3', 'RELAT4', 'RELAT5', 'RELAT6')
    mon = ('MON1', 'MON2', 'MON3', 'MON4', 'MON5', 'MON6')
    for i in range(6):
        execution(ser, relat[i])
        execution(ser, mon[i])
        initialisation(i+1)

    # ---Variables---

if __name__ == '__main__':  # if we're running file directly and not importing it

    ser = discover_and_connect()
    # init_motors(ser)
    init_motors(ser)
    setup()  # run the two configuration functions
    main()  # run the main function
    ser.close()
    logger.info('--- Clean Exit ---')
