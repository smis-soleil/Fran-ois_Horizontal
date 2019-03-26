
#------- DESCRIPTION -------#
" This program is used to control PACE 5000 device."


#-------IMPORTS------------------------------

import serial.tools.list_ports
import serial
import time
from PyQt4 import QtGui
from PyQt4 import QtCore
import sys
import pyqtgraph as pg
import numpy as np

import G_window                    # Import the script of our graphical interface

#-------DEF FUNCTIONS------------------------

def execution(ser, commande):
    # Code to send a command to the controller
    # send the character to the device
    try:
        ser.write(bytes(commande + '\r\n', 'UTF-8'))
        print("command sent: ", (commande + '\r\n'))
        time.sleep(0.15)                                     # waiting for a response
        print(ser.outWaiting())
        print('? =', ser.inWaiting())                       # Look into the input buffer, if it's not empty it's mean that the device received our command
        out = ''
        liste = []
        for i in range(2):
            if ser.inWaiting()!=0:
                while (ser.inWaiting() > 0):                        # Receive and read the device's answer
                    out = str(ser.read(ser.inWaiting()))
                    print('Out = ', out)
                    liste.append(out)
                if (len(liste) > 0):                               # If we get an answer, we take only the usefull information before retruning it
                    rep = out.split(" ")[1]
                    rep = rep[0:-5]
                    print('reponse=', rep, "|")
                    print(type(rep))
                    if type(rep)==str:
                        return rep
            else:
                time.sleep(0.2)
        else:
            return '0'
    except AttributeError :
        info_pace = QtGui.QMessageBox.question(None, 'PACE 5000 Information', "Can't find any device, make sure RS232 cable is connected.")
        return ''

#------- CLASSES -------#

class MainWindow(QtGui.QMainWindow, G_window.Ui_MainWindow):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)


        self.btn_read.clicked.connect(self.readparam)        # We define the purpose of the buttons
        self.btn_stpm.clicked.connect(self.Stop_m)
        self.btn_strtm.clicked.connect(self.Start_m)


#------- INITIALISATION ---------

        self.line_setp1.setText(str(round(float(execution(ser, ":SOUR:PRES?")[0:4]), 3)))
        self.line_slewrate1.setText(str(round((float(execution(ser, ":SOUR:SLEW?")[0:-1]) * 60), 3)))


#-------FUNCTIONS-------

    def readparam(self):                                               # This function allow us to read the current parameters
        self.line_setp2.setText(str(round(float(execution(ser, ":SOUR:PRES?")[0:4]), 3)))
        self.line_slewrate2.setText(str(round((float(execution(ser, ":SOUR:SLEW?")[0:-1])*60), 3)))
        self.line_slewrate_act.setText(str(round((float(execution(ser, ":SENS:SLEW?")[0:5])*60), 3)))
        self.line_pressure2.setText(str(round(float(execution(ser, ":SENS:PRES?")[0:4]), 3)))

        #for i in range(100):                                                                                     # Fonction qui m'a permis de débugger
        #    try:
        #        print("i=",i)
        #        r=self.line_slewrate_act.setText(str(round(float(execution(ser, ":SENS:SLEW?")[0:8]) * 60, 3)))
        #        s= self.line_pressure2.setText(str(round(float(execution(ser, ":SENS:PRES?")[0:4]), 3)))
        #    except TypeError:
        #        print("error : ", r, s)






    def applyparam(self):                                            # This function allow us to apply new parameters
        set_p = self.line_setp1.text()                          # Those parameters are chosen by the user
        slew_rate = str(float((self.line_slewrate1.text()))/60)

        max_setpoint = 80
        if float(set_p)>max_setpoint:
            choice = QtGui.QMessageBox.question(self, 'Control Mode', "Are you sure to apply a set point higher than 80?",
                                                # We put a window to warn the user
                                                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

            if choice == QtGui.QMessageBox.Yes:
                execution(ser, ":SOUR:PRES " + set_p)
        if float(slew_rate)>2/60:
            choice = QtGui.QMessageBox.question(self, 'Control Mode',
                                                "Are you sure to apply a slew rate higher than 2?",
                                                # We put a window to warn the user
                                                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

            if choice == QtGui.QMessageBox.Yes:
                execution(ser, ":SOUR:PRES:SLEW " + slew_rate)
        else:
            execution(ser, ":SOUR:PRES " + set_p)
            execution(ser, ":SOUR:PRES:SLEW " + slew_rate)
            #execution(ser, ":SOUR:SLEW:MODE "+rate_mode)
            QtCore.QCoreApplication.processEvents()

    def Stop_m(self):                                         # Turns controller off
        execution(ser, ":OUTP:STAT OFF")                      # Print the state of the device
        self.line_state.setText('Measure Mode')
        self.line_state.setStyleSheet("QLabel {color : green}")


    def Start_m(self):                                       # Turns controller on
        self.applyparam()
        self.readparam()
        self.plot_p.clear()
        execution(ser, ":OUTP:STAT ON")
        while execution(ser, ":OUTP:STAT?")[0] == '1':
            self.line_state.setText('Control Mode')
            self.line_state.setStyleSheet("QLabel {color : blue}")
            x = 1
            pressure_data = []
            time_p =[]
            start = time.time()
            while x == 1 :

                self.line_slewrate_act.setText(str(round(float(execution(ser, ":SENS:SLEW?")[0:7])*60, 3)))
                self.line_pressure2.setText(str(round(float(execution(ser, ":SENS:PRES?")[0:4]), 3)))
                y_value = time.time()
                time_p.append(y_value-start)
                pressure_data.append(float(execution(ser, ":SENS:PRES?")[0:5]))
                QtCore.QCoreApplication.processEvents()  # Dirty way, I should try with a thread
                time.sleep(0.45)
                QtCore.QCoreApplication.processEvents()  # Dirty way, I should try with a thread
                self.plot_p.plot(y=pressure_data, x=time_p, pen=(19, 234, 201), symbolBrush=(19, 234, 201), symbol='h')
                if execution(ser, ":OUTP:STAT?")[0] == '0':
                    break


                #if round(float(execution(ser, ":SENS:PRES?")[0:4]), 4) == round(float(execution(ser, ":SOUR:PRES?")[0:4]), 4):
                #    execution(ser, ":OUTP:STAT OFF")
                #    self.line_state.setText('Measure Mode')
                #    self.line_state.setStyleSheet("QLabel {color : green}")
                #    print("List data = ", time_p)
                #    self.plot_p.plot(y=pressure_data, x=time_p, pen=(19,234,201), symbolBrush=(19,234,201), symbol='h')
                #    break#



#-------WINDOW-----------

def Main():
    self.btn_read.clicked.connect(self.readparam)  # We define the purpose of the buttons
    self.btn_stpm.clicked.connect(self.Stop_m)
    self.btn_strtm.clicked.connect(self.Start_m)

    # ------- INITIALISATION ---------

    self.line_setp1.setText(str(round(float(execution(ser, ":SOUR:PRES?")[0:4]), 3)))
    self.line_slewrate1.setText(str(round((float(execution(ser, ":SOUR:SLEW?")[0:-1]) * 60), 3)))


    # -------FUNCTIONS-------


def readparam(self):  # This function allow us to read the current parameters
    self.line_setp2.setText(str(round(float(execution(ser, ":SOUR:PRES?")[0:4]), 3)))
    self.line_slewrate2.setText(str(round((float(execution(ser, ":SOUR:SLEW?")[0:-1]) * 60), 3)))
    self.line_slewrate_act.setText(str(round((float(execution(ser, ":SENS:SLEW?")[0:5]) * 60), 3)))
    self.line_pressure2.setText(str(round(float(execution(ser, ":SENS:PRES?")[0:4]), 3)))

    # for i in range(100):                                                                                     # Fonction qui m'a permis de débugger
    #    try:
    #        print("i=",i)
    #        r=self.line_slewrate_act.setText(str(round(float(execution(ser, ":SENS:SLEW?")[0:8]) * 60, 3)))
    #        s= self.line_pressure2.setText(str(round(float(execution(ser, ":SENS:PRES?")[0:4]), 3)))
    #    except TypeError:
    #        print("error : ", r, s)


def applyparam(self):  # This function allow us to apply new parameters
    set_p = self.line_setp1.text()  # Those parameters are chosen by the user
    slew_rate = str(float((self.line_slewrate1.text())) / 60)

    max_setpoint = 80
    if float(set_p) > max_setpoint:
        choice = QtGui.QMessageBox.question(self, 'Control Mode', "Are you sure to apply a set point higher than 80?",
                                            # We put a window to warn the user
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

        if choice == QtGui.QMessageBox.Yes:
            execution(ser, ":SOUR:PRES " + set_p)
    if float(slew_rate) > 2 / 60:
        choice = QtGui.QMessageBox.question(self, 'Control Mode',
                                            "Are you sure to apply a slew rate higher than 2?",
                                            # We put a window to warn the user
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

        if choice == QtGui.QMessageBox.Yes:
            execution(ser, ":SOUR:PRES:SLEW " + slew_rate)
    else:
        execution(ser, ":SOUR:PRES " + set_p)
        execution(ser, ":SOUR:PRES:SLEW " + slew_rate)
        # execution(ser, ":SOUR:SLEW:MODE "+rate_mode)
        QtCore.QCoreApplication.processEvents()


def Stop_m(self):  # Turns controller off
    execution(ser, ":OUTP:STAT OFF")  # Print the state of the device
    self.line_state.setText('Measure Mode')
    self.line_state.setStyleSheet("QLabel {color : green}")


def Start_m(self):  # Turns controller on
    self.applyparam()
    self.readparam()
    self.plot_p.clear()
    execution(ser, ":OUTP:STAT ON")
    while execution(ser, ":OUTP:STAT?")[0] == '1':
        self.line_state.setText('Control Mode')
        self.line_state.setStyleSheet("QLabel {color : blue}")
        x = 1
        pressure_data = []
        time_p = []
        start = time.time()
        while x == 1:

            self.line_slewrate_act.setText(str(round(float(execution(ser, ":SENS:SLEW?")[0:7]) * 60, 3)))
            self.line_pressure2.setText(str(round(float(execution(ser, ":SENS:PRES?")[0:4]), 3)))
            y_value = time.time()
            time_p.append(y_value - start)
            pressure_data.append(float(execution(ser, ":SENS:PRES?")[0:5]))
            QtCore.QCoreApplication.processEvents()  # Dirty way, I should try with a thread
            time.sleep(0.45)
            QtCore.QCoreApplication.processEvents()  # Dirty way, I should try with a thread
            self.plot_p.plot(y=pressure_data, x=time_p, pen=(19, 234, 201), symbolBrush=(19, 234, 201), symbol='h')
            if execution(ser, ":OUTP:STAT?")[0] == '0':
                break





# -------WINDOW-----------

def Main():
    app = QtGui.QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()


def discover_connect():  # Allow us to find the COM port

    ports = list(serial.tools.list_ports.comports())  # List all the COM ports available

    for p in ports:  # Loop over available ports
        print("p = ", p)
        # Open the port
        ser = serial.Serial(port=p[0],  # Set the right parameters in order to communicate with the PACE 5000
                            baudrate=9600,
                            stopbits=serial.STOPBITS_ONE,
                            bytesize=serial.EIGHTBITS,
                            xonxoff=True)

        Longueur = len(execution(ser, ':SOUR:PRES?'))  # Query something to the device
        print('Exec : longueur = ', Longueur)
        if Longueur > 1 and str(ser.isOpen()) == 'True':  # If we get an answer it's the right device
            return (ser)  # ser is now the right COM port
        else:

            print('Connection error. Make sure RS232 cable is connected.')







if __name__ == '__main__':                  # If we start from this file launch the folowing
    ser = discover_connect()
    Main()
    ser.close()




