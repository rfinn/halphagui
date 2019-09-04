#!/usr/bin/env python

'''
GOAL:
- open gui window to run galfit
- show image, model, residual

PROCEDURE

'''


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
    def __init__(self, MainWindow, logger, image=None, sigma_image=None, mask_image=None, psf=None,psf_oversampling=None, xmaxfit=None, ymaxfit=None, xminfit=1, yminfit=1, ncomp=1, convflag = True, convolution_size=None, fitallflag=False,xc=None, yc=None,mag=None,rad=None,nsersic=None, BA=None,PA=None):
        super(galfitwindow, self).__init__()

        # boiler plate gui stuff
        # ... at least I think it is...
        self.ui = Ui_galfitWindow()
        self.ui.setupUi(MainWindow)
        MainWindow.setWindowTitle('GALFIT!')
        self.MainWindow = MainWindow

        # name of input file
        self.image_name = image
        # store image data and header
        self.image_data, self.image_header = fits.getdata(self.image_name, header=True)

        # define root name for galfit output file
        self.image_rootname = self.image_name.split('.fits')[0]
        
        # define fit area as the complete size of the input image
        # this assumes that the user is providing a cutout rather than a large mosaic image
        if xmaxfit == None:
            self.xmaxfit = self.image_data.shape[0]
        else:
            self.xmaxfit = xmaxfit

        if ymaxfit == None:
            self.ymaxfit = self.image_data.shape[1]
        else:
            self.ymaxfit = ymaxfit
        self.xminfit = xminfit
        self.yminfit = yminfit
        
        # try to extract the magnitude ZP from image header
        # if not available, set to 25.
        try:
            self.magzp = self.image_header['PHOTZP']
        except KeyError:
            self.magzp = 25.

        # get pixel scale from image header
        # convert from degrees/pix to arcsec/pix
        self.pscale = abs(float(self.image_header['CD1_1'])*3600)
        
        self.sigma_image = sigma_image
        self.psf_image = psf
        self.psf_oversampling = psf_oversampling
        
        self.mask_image = mask_image
        if self.mask_image != None:
            self.mask_data = fits.getdata(self.mask_image)
            # create a masked array to display to the user
            self.image_data = np.ma.masked_array(self.image_data, mask=np.array(self.mask_data,'bool'))

        # set to false for fitting central object only
        self.fitallflag = fitallflag
        # number of components to fit
        # set to 1 for single component sersic fit
        self.ncomp = ncomp
        
        ###########################################################3
        # enable psf convolution in fit
        ###########################################################3
        self.convflag = convflag
        if convolution_size == None:
            self.convolution_size = min(self.xmaxfit, self.ymaxfit)
        else:
            self.convolution_size = convolution_size
            
        ###########################################################3
        # set up initial guesses for sersic parameters
        ###########################################################3
        if xc == None:
            self.xc=self.image_data.shape[0]/2
        else:
            self.xc=xc
        if yc == None:
            self.yc=self.image_data.shape[1]/2
        else:
            self.yc=yc
        if mag == None:
            self.mag=15.
        else:
            self.mag=mag
        if rad == None:
            self.re=5
        else:
            self.re=rad
        if nsersic == None:
            self.nsersic=2
        else:
            self.nsersic=nsersic
        if BA == None:
            self.BA= 0.6
        else:
            self.BA= BA
        if PA == None:
            self.PA=0
        else:
            self.PA=PA

        self.fitBA = 1
        self.fitPA = 1

        # define galfit output image the same way that rungalfit does
        self.output_image=self.image_rootname+'-'+ str(self.ncomp) +'Comp-galfit-out.fits'
        # continue with other functions
        self.logger = logger
        self.add_cutout_frames()
        self.connect_buttons()
        self.display_image()

        self.initialize_galfit()        
        #self.galfit = rg.galfit(galname=self.image_rootname,image=self.image, mask_image = self.mask_image, sigma_image=self.sigma_image,psf_image=self.psf_image,psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,convolution_size=self.convolution_size,magzp=self.magzp,pscale=self.pscale,ncomp=self.ncomp,convflag=convflag)


    def add_cutout_frames(self):
        # r-band cutout
        a = QtWidgets.QLabel('Image')
        self.ui.cutoutsLayout.addWidget(a, 0, 0, 1, 1)
        a = QtWidgets.QLabel('Model')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Residual')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.cutout_frame = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 4, 1)
        self.model_frame = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 4, 1)
        self.residual_frame = my_cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 4, 1)
        self.residual_frame.key_pressed.connect(self.key_press_func)
    def display_image(self):
        if self.mask_image == None:
            self.cutout_frame.load_file(self.image)
        else:
            # display image with masked values
            self.cutout_frame.load_image(self.image_data)

    def display_galfit_results(self):
        pass
    def connect_buttons(self):
        self.ui.quitButton.clicked.connect(self.quit_program)
        self.ui.helpButton.clicked.connect(self.print_help_menu)
        self.ui.runGalfitButton.clicked.connect(lambda: self.run_galfit(fitBA=self.fitBA, fitPA=self.fitPA))
    def key_press_func(self,text):
        key, x, y = text.split(',')
        self.xcursor = float(x)
        self.ycursor = float(y)
        if key == 'n':
            self.galfit.set_n()
        elif key == 'r':
            self.galfit.set_r()
        elif key == 'o':
            self.galfit.toggle_fitall()
        elif key == 'c':
            if ((self.xcursor > self.image_data.shape[0]) | (self.ycursor > self.image_data.shape[0]) | (self.xcursor < 1.) | (self.ycursor < 1.)):
                print('cursor value out of range.  put mouse over desired center and type c')
                
            self.galfit.set_center(self.xcursor, self.ycursor)
        elif key == 'f':
            self.galfit.print_fix_menu()
        elif key == 'a':
            self.galfit.toggle_asymmetry()
        elif key == 'R':
            self.galfit.reset_sersic_params()
        elif key == 'g':
            self.run_galfit()
        elif key == 'h': 
            self.print_help_menu()
        elif key == 'x':
            self.quit_program()
        else:
            print('did not understand that.  \n Try again!')

    def print_help_menu(self):
        print('What is wrong?\n n = adjust sersic \
        \n r = reset Re \
        \n o = nearby object (toggle fitall)\
        \n b = B/A \n p = PA \n m = mag\
        \n c = recenter\
        \n f = hold values fixed\
        \n a = toggle asymmetry parameter\
        \n R = reset to original values\
        \n g = go (run galfit)\
        \n x=quit \n \
        \nDisplay shortcuts (click on image to adjust):'+
        '\n \t scroll = zoom'+
        '\n \t ` = zoom to fit'+
        '\n \t ALT-right_click = adjust contrast \n \n'+
        'Click Red X to close window')
    def quit_program(self):
        self.close_window()

#class galfit_galaxy():


    def initialize_galfit(self,convflag=True):
        '''
        GOAL: Preparing file to be run in galfit. Initialize galfit image parameters 

        INPUT: nsaid 

        OUTPUT: A definition of everything from galname to convflag, necessary for running galfit

        '''
        print('self.psfimage = ',self.psf_image)
        
        self.galfit = rg.galfit(galname=self.image_rootname,image=self.image_name,
                                mask_image = self.mask_image,
                                sigma_image=self.sigma_image,psf_image=self.psf_image,
                                psf_oversampling=self.psf_oversampling,xminfit=self.xminfit,
                                yminfit=self.yminfit,xmaxfit=self.xmaxfit,ymaxfit=self.ymaxfit,
                                convolution_size=self.convolution_size,magzp=self.magzp,
                                pscale=self.pscale,ncomp=self.ncomp,convflag=convflag,
                                fitallflag = self.fitallflag)
        
    def run_galfit(self,fitBA=1,fitPA=1):
        '''
        GOAL: take values from NSA tables in /Virgo 

        INPUT: nsaid, ba/pa = 1

        OUTPUT: several output files

        '''

        self.galfit.set_sersic_params(xobj=self.xc,yobj=self.yc,mag=self.mag,rad=self.re,nsersic=self.nsersic,BA=self.BA,PA=self.PA,fitmag=1,fitcenter=1,fitrad=1,fitBA=fitBA,fitPA=fitPA,fitn=1,first_time=0)
        print('in galfitwrapper, run_galfit: fitBA = ',fitBA)
        self.galfit.set_sky(0)
        self.galfit.run_galfit()
        self.display_results()
        self.get_galfit_results(printflag=True)
        #self.galfit.close_input_file()
        #self.galfit.print_params()
        #self.galfit.print_galfit_results()
    def get_galfit_results(self,printflag = False):
        '''
        GOAL: Grab results from galfit (xc, yc, mag, re, nsersic, BA, PA, sky, error, chi2nu) and parse them into self.filename

        INPUT: nsaid

        OUTPUT: 1Comp file of output from galfit, header for the output file

        '''


        t = rg.parse_galfit_1comp(self.output_image)
        if printflag:
            self.galfit.print_galfit_results(self.output_image)
        
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
        self.model_data = fits.getdata(self.output_image,2)
        self.residual_data = fits.getdata(self.output_image,3)
        self.model_frame.load_image(self.model_data)
        self.residual_frame.load_image(self.residual_data)
        pass
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
    ui = galfitwindow(MainWindow, logger, image='MKW8-18216-R.fits', mask_image = 'MKW8-18216-R-mask.fits',
                      psf='MKW8_R.coadd-psf.fits',psf_oversampling=2.)
    #ui.setupUi(MainWindow)
    #ui.test()

    MainWindow.show()
    sys.exit(app.exec_())

    
