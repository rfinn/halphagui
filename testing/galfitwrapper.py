#!/usr/bin/env python

import os
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Tophat2DKernel
import numpy as np
#import argparse
#import pyds9
from scipy.stats import scoreatpercentile
from PyQt5 import QtCore, QtGui, QtWidgets
#from PyQt5 import  QtWidgets
#from ginga.qtw.QtHelp import QtGui
#from PyQt5 import QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
from ginga import colors
from ginga.canvas.CanvasObject import get_canvas_types

from ginga.misc import log
from ginga.util.loader import load_data

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


from galfitWidget import Ui_Form as Ui_galfitWindow
from maskwrapper import my_cutout_image

import rungalfit as rg

# add a frame on top to show
# cutout, mask, model, residual

class galfitwindow(Ui_galfitWindow, QtCore.QObject):
    mask_saved = QtCore.pyqtSignal(str)
    def __init__(self, MainWindow, logger, image=None, haimage=None, sepath=None, config=None, threshold=0.05,snr=2,cmap='gist_heat_r'):
        super(galfitwindow, self).__init__()
        
        self.ui = Ui_galfitWindow()
        self.ui.setupUi(MainWindow)
        MainWindow.setWindowTitle('makin a mask...')
        self.MainWindow = MainWindow

        #self.readout = QtWidgets.QLabel('this is a test')
        #self.ui.gridLayout_2.addWidget(self.readout, 0,0,1,2)
        #self.readout.setText('this is another test')

        self.logger = logger
    def add_cutout_frames(self):
        # r-band cutout
        a = QtWidgets.QLabel('r-band')
        self.ui.cutoutsLayout.addWidget(a, 0, 0, 1, 1)
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.rcutout = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 4, 1)
        self.hacutout = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 4, 1)
        self.maskcutout = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 4, 1)
        self.maskcutout.key_pressed.connect(self.key_press_func)
    def connect_buttons(self):
        self.ui.quitButton.clicked.connect(self.quit_program)
        self.ui.helpButton.clicked.connect(self.print_help_menu)
        self.ui.runGalfitButton.clicked.connect(self.run_galfit)
    def key_press_func(self,text):
        key, x, y = text.split(',')
        self.xcursor = float(x)
        self.ycursor = float(y)
        self.cursor_value = self.maskdat[int(self.ycursor),int(self.xcursor)]
        print('cursor value = ',self.cursor_value, key)
        if key == 'a':
            self.add_object()
        elif key == 'r': 
            print('removing object')
            self.remove_object(int(self.cursor_value))
        elif key == 'o': 
            self.off_center()
        #elif key == 's': 
        #    self.set_box_size()
        #elif key == 't': 
        #    self.set_threshold()
        #elif key == 'n': 
        #    self.set_sesnr()
        elif key == 'h': 
            self.print_help_menu()
        elif key == 's': 
            self.save_mask()
        elif key == 'q':
            self.quit_program()
        else:
            print('did not understand that.  \n Try again!')

    def print_help_menu(self):
        print('Click on mask image, then enter:\n \t r to remove object in mask at the cursor position;'
              '\n \t a to mask additional pixels at cursor position;'
              '\n \t o if target is off center (and program is removing the wrong object);'
              #'\n \t s to change the size of the mask box;'
              #'\n \t t to adjust SE threshold (0=lots, 1=no deblend );'
              #'\n \t n to adjust SE SNR; '
              '\n \t h to print this menu; '
              '\n \t s to save output;'
              #'\n \t q to quit \n \n'
              'Display shortcuts (click on image to adjust):'
              '\n \t scroll = zoom'
              '\n \t ` = zoom to fit'
              '\n \t ALT-right_click = adjust contrast \n \n'
              'Click Red X to close window')
    def quit_program(self):
        self.clean_links()
        self.close_window()

class galfit_galaxy():

    def __init__(self, image=None, sigma_image=None, psf=None):
        
        self.image = image
        self.noise_image = noise_image
        self.psf = psf
    def initialize_galfit(self,convflag=True):
        '''
        GOAL: Preparing file to be run in galfit. Initialize galfit image parameters 

        INPUT: nsaid 

        OUTPUT: A definition of everything from galname to convflag, necessary for running galfit

        '''
        print 'self.psfimage = ',self.psf_image
        
        self.galfit = rg.galfit(galname=self.image_rootname,image=self.image, mask_image = self.mask_image, sigma_image=self.sigma_image,psf_image=self.psf_image,psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,convolution_size=self.convolution_size,magzp=self.magzp,pscale=self.pscale,ncomp=self.ncomp,convflag=convflag)
        
   def run_galfit(self,fitBA=1,fitPA=1):
        '''
        GOAL: take values from NSA tables in /Virgo 

        INPUT: nsaid, ba/pa = 1

        OUTPUT: several output files

        '''
        self.galfit.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
        self.galfit.set_sky(0)
        self.galfit.run_galfit()
        #self.galfit.display_results()
        #self.galfit.close_input_file()
        #self.galfit.print_params()
        #self.galfit.print_galfit_results()
    def get_galfit_results(self,printflag = False):
        '''
        GOAL: Grab results from galfit (xc, yc, mag, re, nsersic, BA, PA, sky, error, chi2nu) and parse them into self.filename

        INPUT: nsaid

        OUTPUT: 1Comp file of output from galfit, header for the output file

        '''

        self.filename = 'NSA-'+str(self.nsaid)+'-unwise-'+'w'+str(self.band)+'-1Comp-galfit-out.fits'
        t = rg.parse_galfit_1comp(self.filename)
        if printflag:
            self.gal1.print_galfit_results(self.filename)
        
        header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        self.xc, self.xc_err = t[0]
        self.yc, self.yc_err = t[1]
        self.mag, self.mag_err = t[2]
        self.re, self.re_err = t[3]
        self.nsersic, self.nsersic_err = t[4]
        self.BA, self.BA_err = t[5]
        self.PA, self.PA_err = t[6]
        self.sky, self.sky_err = t[7]
        self.error = t[8]
        self.chi2nu = t[9]
    def display_results(self):
        print '%%%%%%%%%%%%%%%%%%'
        print 'inside display_results'
        print 'self.galfit_flag = ',self.galfit_flag
        ### THIS SHOULD BE MOVED TO PARENT PROGRAM
        ### BECAUSE HOW U DISPLAY RESULTS WILL VARY WITH USAGE

        
if __name__ == "__main__":
    #catalog = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    #from halphamain import cutout_image
    #from halphamain import cutout_image
    logger = log.get_logger("galfitlog", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    MainWindow = QtWidgets.QWidget()
    ui = galfitwindow(MainWindow, logger)
    #ui.setupUi(MainWindow)
    #ui.test()

    MainWindow.show()
    sys.exit(app.exec_())

    
