#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ZetCode PySide tutorial 

This program centers a window 
on the screen. 

author: Jan Bodnar
website: zetcode.com 
last edited: August 2011
"""

import sys
import numpy as np
from PySide import QtGui

class Example(QtGui.QWidget):
    
    def __init__(self):
        super(Example, self).__init__()
        
        self.initUI()
        
    def initUI(self):               
        
        QtGui.QToolTip.setFont(QtGui.QFont('SansSerif', 10))
        
        self.setToolTip('This is a <b>QWidget</b> widget')
        
        btn = QtGui.QPushButton('Button', self)
        btn.setToolTip('This is a <b>QPushButton</b> widget')
        btn.resize(btn.sizeHint())
        btn.move(50, 50)       

	# dynamically create button array
	bbGroup={}
	bb={}
	bbID={0:'Schwarzshield',
	    1:'Raman',
	    2:'DAC',
	    3:'S',
	    4:'Y',
	    5:'Z'
	}
	bbStackHPos = 0 # for storing the positions
	for i in range(len(bbID)):
		bbGroup[i] = QtGui.QButtonGroup(self) # each button has to belong to its own QButtonGroup to be independently clickable
		bb[i] = QtGui.QRadioButton(bbID[i], self)
		bbGroup[i].addButton(bb[i])
		bb[i].resize(bb[i].sizeHint())
		bb[i].move(10+bbStackHPos, 10)	
		bbStackHPos = bbStackHPos + bb[i].sizeHint().width()

	qst = QtGui.QStatusBar(self)
	qst.showMessage('Ready')
	#qst.resize(qst.sizeHint())
	qst.resize(240, 20)
	qst.move(5,120)
	qst.show()

        #menubar = self.QMenuBar()
        #fileMenu = menubar.addMenu('&File')
        #fileMenu.addAction(exitAction)
        
	#moods = [QtGui.QRadioButton("Happy"), QtGui.QRadioButton("Sad"), QtGui.QRadioButton("Angry")]
	#for i in range(len(moods)):
		#moods[i].resize(moods.sizeHint())
		#moods[i].move(20,50+i*10)
 
        self.resize(250, 150)
        self.center()
        
        self.setWindowTitle('Center')    
        self.show()
        
    def center(self):
        
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
        
def main():
    
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
