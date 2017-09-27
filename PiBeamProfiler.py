#!/usr/bin/env python
# Copyright (C) 2015 Anthony Ransford
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from picamera.array import PiRGBArray
from picamera import PiCamera
from PIL.ImageQt import ImageQt
from PyQt4 import QtGui, QtCore
import numpy as np
from scipy.misc.pilutil import toimage
from scipy.optimize import curve_fit
import time, sys
import matplotlib.animation as animation
import cv2

class proflayout(QtGui.QWidget):

    def __init__(self):
        super(proflayout, self).__init__()
	self.imageres = [640,480]
	self.zoom = 1
	self.getzoomgaps()
	self.fitting = True
	self.breakloop = False
	desktop = QtGui.QDesktopWidget()
	screensize = desktop.availableGeometry()
	self.screenres = [screensize.width(),screensize.height()]
	self.initCamera()
	self.initializeGUI()

    def initCamera(self):

	# initialize the camera
	self.camera = PiCamera()

	#set camera resolution, gain , sutter speed and framerate
	self.camera.resolution = (self.imageres[0], self.imageres[1])
	self.camera.framerate = 33 # in Hz
	self.camera.shutter_speed = 100 # in us
	self.camera.exposure_mode = 'off'
	self.camera.iso = 300

	#grab a reference to the raw camera capture
	self.rawCapture = PiRGBArray(self.camera, size=(self.imageres[0], self.imageres[1]))

	# allow the camera to warmup
	time.sleep(0.1)

    def initializeGUI(self):

        self.setWindowTitle('Beam Profiler')
	self.setGeometry(0, 0, self.screenres[0], self.screenres[1])
        layout = QtGui.QGridLayout()
	
	self.createPlots()
	self.setupPlots()

	self.expslider = QtGui.QSlider(QtCore.Qt.Vertical)
	self.expslider.setSingleStep(1)
	self.explabel = QtGui.QLabel('Exposure')
	self.expbar = QtGui.QProgressBar()
	self.expbar.setOrientation(QtCore.Qt.Vertical)
	self.expbar.setValue(65)
	self.videowindow = QtGui.QLabel(self)

	self.xwaist = QtGui.QLabel()
	self.ywaist = QtGui.QLabel()
	self.xwaist.setStyleSheet('color: #FF6600; font-weight: bold; font-family: Copperplate / Copperplate Gothic Light, sans-serif') 
	self.ywaist.setStyleSheet('color: #FF6600; font-weight: bold; font-family: Copperplate / Copperplate Gothic Light, sans-serif')
	self.zoominbutton = QtGui.QPushButton('Zoom In')
	self.zoomoutbutton = QtGui.QPushButton('Zoom Out')
	buttonsize = [int(self.screenres[1]/8 ), int(self.screenres[1]/4)]
	self.highresbutton = QtGui.QPushButton('1296x972')
	self.lowresbutton = QtGui.QPushButton('640x480')
	self.highresbutton.setCheckable(True)
	self.lowresbutton.setCheckable(True)
	self.lowresbutton.setChecked(True)
	self.highresbutton.setFixedSize(buttonsize[0],buttonsize[1])
	self.lowresbutton.setFixedSize(buttonsize[0],buttonsize[1])
	self.zoominbutton.setFixedSize(buttonsize[0],buttonsize[1])
	self.zoomoutbutton.setFixedSize(buttonsize[0],buttonsize[1])
        self.zoominbutton.toggled.connect(self.zoomin) 
	self.setupPlots()
        self.canvasrow = FigureCanvas(self.figurerow)
	self.canvascolumn = FigureCanvas(self.figurecolumn)

	self.expslider.valueChanged[int].connect(self.changeExposure)
        self.zoominbutton.clicked.connect(self.zoomin) 
        self.zoomoutbutton.clicked.connect(self.zoomout) 
	self.lowresbutton.clicked.connect(self.lowres)
	self.highresbutton.clicked.connect(self.highres)

	layout.addWidget(self.videowindow,   0,0,2,1)
	layout.addWidget(self.canvasrow,     2,0,2,1)
	layout.addWidget(self.canvascolumn,  0,1,2,1)
	layout.addWidget(self.expbar,        0,4,2,1)
	#withholds these widgets for tiny screens
	print self.screenres
	if not ((self.screenres[0] or self.screenres[1]) <= 400):
		layout.addWidget(self.lowresbutton,  1,2)
		layout.addWidget(self.highresbutton, 1,3)
		layout.addWidget(self.zoominbutton,  0,3)
		layout.addWidget(self.zoomoutbutton, 0,2)
		layout.addWidget(self.expslider,     0,5,2,1)
	layout.addWidget(self.xwaist,        2,1)
	layout.addWidget(self.ywaist,        3,1)

        self.setLayout(layout)

    def startCamera(self):
	# capture frames from the camera

	for frame in self.camera.capture_continuous(self.rawCapture, format="bgr", use_video_port=True):

		#start = time.time()

		# grab the raw NumPy array representing the imagef
		image = frame.array
		np.nan_to_num(image)


		#take the green part of the image
		greenimage = image[:,:,1]
		globmax = np.max(greenimage)

		#cv2 thingy
 		key = cv2.waitKey(1) & 0xFF
		
		#row and colum sum for live plots
		columnsum = greenimage.sum(axis=1)/40.0
		columnsum = columnsum[::-1]
		rowsum = greenimage.sum(axis=0)/40.0

		#subtract minumum value (background subtraction)
		columnsum = columnsum - np.min(columnsum)
		rowsum = rowsum - np.min(rowsum)
		columnampguess = columnsum.max()
		columncenterguess = np.argmax(columnsum)

		rowampguess = rowsum.max()
		rowcenterguess = np.argmax(rowsum)
		percexp = 100 * globmax/255.0
		self.expbar.setValue(percexp)

		rowampguess = rowsum.max()
		rowcenterguess = np.argmax(columnsum)
		coursecolumny, coursecolumnx = self.coursen(self.ypixels, columnsum, 3)
		courserowx, courserowy = self.coursen(self.xpixels, rowsum, 3)
		coursecolumny = np.nan_to_num(coursecolumny)
		coursecolumnx = np.nan_to_num(coursecolumnx)
		courserowy = np.nan_to_num(courserowy)
		courserowx = np.nan_to_num(courserowx)
		columnampguess = coursecolumnx.max()
		columncenterguess = np.argmax(coursecolumnx)
		#print 'init matrix time = ', time.time() - start
		#start = time.time()
		if self.fitting == True:
			try:
				popt1, pcov1 = curve_fit(self.func, courserowx,courserowy, p0=[rowampguess,rowcenterguess,200])
			except:
				popt1 = [0,0,1]

			try:
				popt2, pcov2 = curve_fit(self.func, coursecolumny, coursecolumnx, p0=[columnampguess,columncenterguess,200])
			except:
				popt2 = [0,0,1]
		else:
			popt1, popt2 = [[0,0,1],[0,0,1]]

		#print 'fitting time = ', time.time() - start

		#start = time.time()
		#print (popt1[0] - rowampguess), popt1[1] - rowcenterguess
		#print (popt2[0] - columnampguess), popt2[1] - columncenterguess
		#updates data for row and column plots, also mirrors column data
        	self.linesrow.set_xdata(courserowx)
        	self.linesrow.set_ydata(courserowy)

        	self.linescolumn.set_xdata(coursecolumnx)
        	self.linescolumn.set_ydata(coursecolumny)

		#updates data for fit row and column plots



        	self.linesrowfit.set_xdata(courserowx)
        	self.linesrowfit.set_ydata(self.func(courserowx, popt1[0],popt1[1],popt1[2]))

        	self.linescolumnfit.set_xdata(self.func(coursecolumny, popt2[0],popt2[1],popt2[2]))
        	self.linescolumnfit.set_ydata(coursecolumny)

        	#self.linescolumnfit.set_xdata(coursecolumnx)
        	#self.linescolumnfit.set_ydata(coursecolumny)

        	#self.linesrowfit.set_xdata(courserowx)
        	#self.linesrowfit.set_ydata(courserowy)


        	#draw data and flush
        	self.figurerow.canvas.draw()
        	self.figurerow.canvas.flush_events()

        	self.figurecolumn.canvas.draw()
        	self.figurecolumn.canvas.flush_events()
		

        	#update X and Y waist labels with scaled waists
		if self.imageres[0] == 640:
			f = 2
		else:
			f=1
		self.xwaist.setText('X = ' + str(np.abs(popt1[2]*2*5.875*f))[0:5] + 'um')
		self.ywaist.setText('Y = ' +str(np.abs(popt2[2]*2*5.875*f))[0:5]  + 'um')

		#print 'updating plots = ', time.time() - start

		#start = time.time()

		# convert RGB image np array to qPixmap and update canvas widget
		image = image[int(self.gaprow):self.imageres[0] - int(self.gaprow),int(self.gapcolumn):self.imageres[1] - int(self.gapcolumn)]
		qPixmap = self.nparrayToQPixmap(image)
		videoy = int(self.screenres[0]/2.1)
		videox = int(1.333 * videoy)
		self.videowindow.setPixmap(qPixmap.scaled(videox,videoy))
 
		# clear the stream in preparation for the next frame
		self.rawCapture.truncate(0)
		#print 'image updating = ', time.time() - start

    def createPlots(self):

        #Set up plot axes and figure positions
        self.figurerow, self.axrow = plt.subplots()
	#self.figurerow.gca().set_position([0,0,1,1])

        self.figurecolumn, self.axcolumn = plt.subplots()
	#self.figurecolumn.gca().set_position([0,0,1,1])

	#Create line objects for fast plot redrawing
        self.linesrow, = self.axrow.plot([],[],linewidth=2,color='purple')
        self.linesrowfit, = self.axrow.plot([],[],linestyle='--',linewidth=2,color='yellow')

        self.linescolumn, = self.axcolumn.plot([],[],linewidth=2,color='purple')
	self.linescolumnfit, = self.axcolumn.plot([],[],linestyle='--',linewidth=2,color='yellow')

    def setupPlots(self):

	self.xpixels = np.linspace(0,self.imageres[0],self.imageres[0])
	self.ypixels = np.linspace(0,self.imageres[1],self.imageres[1])

        self.axrow.set_xlim(0, self.imageres[0])
	self.axrow.set_ylim(0,300)

        self.axcolumn.set_xlim(0, 300)
	self.axcolumn.set_ylim(0,self.imageres[1])

	self.axrow.xaxis.set_ticks_position('none')
	self.axrow.yaxis.set_ticks_position('none')
	self.axrow.get_xaxis().set_visible(False)
	self.axrow.get_yaxis().set_visible(False)
	self.axrow.patch.set_visible(False)

	self.axcolumn.xaxis.set_ticks_position('none')
	self.axcolumn.yaxis.set_ticks_position('none')
	self.axcolumn.get_xaxis().set_visible(False)
	self.axcolumn.get_yaxis().set_visible(False)
	self.axcolumn.patch.set_visible(False)

    def changeExposure(self, value):
	scaledvalue = 0.5 * value**2 + 1
	self.camera.shutter_speed = int(scaledvalue)

    #gaussian function used in fitting routine
    def func(self, x, a, x0, sigma):
   	return a*np.exp(-(x-x0)**2/(2*sigma**2))

    #converts nparray to qpixmap
    def nparrayToQPixmap(self, arrayImage):
    	pilImage = toimage(arrayImage)
    	qtImage = ImageQt(pilImage)
    	qImage = QtGui.QImage(qtImage)
    	qPixmap = QtGui.QPixmap(qImage)
    	return qPixmap

    #to be added
    def zoomin(self):
	if self.zoom >= 10:
		self.zoom = 10
	else:
		self.zoom += 1
		self.resizePlots()

    def zoomout(self):
	if self.zoom <= 1:
		self.zoom = 1
	else:
		self.zoom -= 1
		self.resizePlots()

    def lowres(self):
	self.highresbutton.setChecked(False)
	self.breakloop = True
	self.imageres = [640, 480]
	time.sleep(1)
	self.camera.close()
	self.setupPlots()
	self.initCamera()
	self.startCamera()

    def highres(self):
	self.lowresbutton.setChecked(False)
	self.breakloop = True
	self.imageres = [1296,972]
	time.sleep(1)
	self.camera.close()
	self.setupPlots()
	self.initCamera()
	self.startCamera()

    def getzoomgaps(self):
	self.gaprow = self.imageres[0]*(self.zoom * 0.04)
	self.gapcolumn = self.imageres[1]*(self.zoom * 0.04)

    def resizePlots(self):
	self.getzoomgaps()
        self.axrow.set_xlim(self.gaprow, self.imageres[0] - self.gaprow)
	self.axrow.set_ylim(0,300)

        self.axcolumn.set_xlim(0, 300)
	self.axcolumn.set_ylim(self.gapcolumn,self.imageres[1] - self.gapcolumn)


    def coursen(self, xdata,ydata,points):
        newlength = int(len(xdata)/points)
        newxdata = []
        newydata = []
        j = 0
        for i in range(newlength):
        	i = points*(i)
        	newydata.append(np.mean(ydata[int(j):int(i)]))
        	newxdata.append(xdata[int((i + j)/2)])
        	j = i
    	return np.array(newxdata), np.array(newydata)
	
    def closeEvent(self, x):
	self.camera.close()

if __name__ == "__main__":

    a = QtGui.QApplication([])
    proflayoutwidget = proflayout()
    proflayoutwidget.show()
    proflayoutwidget.startCamera()
    sys.exit(a.exec_())


