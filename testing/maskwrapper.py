#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to create a mask for a galaxy image to mask
out other objects within the cutout area.

USAGE:

from within ipython type:

   %run ~/github/HalphaImaging/uat_mask.py --image 'A1367-140231-R.fits'

you just need to run this on R-band images.

Interacting with the display is finicky.  This works fine when running
within ipython - not so much when running from the command line.  When running
from the command line, I am not able to interact with the figure.  This may
have something to do with setting block=False in show().


PROCEDURE:


REQUIRED MODULES:
   os
   astropy
   numpy
   argsparse
   matplotlib
   scipy

NOTES:
- rewrote using a class

'''

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Tophat2DKernel, convolve
from astropy.convolution.kernels import CustomKernel
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

#from maskGui import Ui_maskWindow
from maskWidget import Ui_Form as Ui_maskWindow
from halphaCommon import cutout_image, circle_pixels

defaultcat='default.sex.HDI.mask'

class my_cutout_image(QtCore.QObject):#QtCore.QObject):
    #mouse_clicked = QtCore.pyqtSignal(str)
    key_pressed = QtCore.pyqtSignal(str)
    def __init__(self,panel_name,ui, logger, row, col, drow, dcol,autocut_params='stddev'):
        #super(image_panel, self).__init__(panel_name,ui, logger, row, col, drow, dcol)
        # enable some user interaction
        #fi.get_bindings.enable_all(True)
        #super(my_cutout_image, self).__init__(panel_name,ui, logger, row, col, drow, dcol)
        #super(my_cutout_image,self).__init__(self)
        QtCore.QObject.__init__(self)
        self.logger = logger
        self.ui = ui
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params(autocut_params)
        fi.enable_autozoom('on')
        #fi.set_callback('drag-drop', self.drop_file)
        #fi.set_callback('none-move',self.cursor_cb)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        self.fitsimage = fi
        fi.show_focus_indicator(True)
        
        bd = fi.get_bindings()
        bd.enable_all(True)


        self.cutout = fi.get_widget()
        self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.fitsimage = fi
        
        self.fitsimage.set_callback('none-move',self.cursor_cb)
        #self.fitsimage.set_callback('cursor-down',self.cursor_down)
        self.readout = QtWidgets.QLabel('')
        ui.readoutGridLayout.addWidget(self.readout, 1, col, 1, 1)
        #self.ui.readoutLabel.setText('this is another test')
        self.fitsimage.set_callback('key-press',self.key_press_cb)
    def load_image(self, imagearray):
        #self.fitsimage.set_image(imagearray)
        self.fitsimage.set_data(imagearray)

    def load_file(self, filepath):
        image = load_data(filepath, logger=self.logger)
        self.image_wcs = WCS(filepath)
        self.fitsimage.set_image(image)
        #self.setWindowTitle(filepath)

    # adding cursor readout so we can identify the pixel values to change

    def cursor_cb(self, viewer, button, data_x, data_y):
        """This gets called when the data position relative to the cursor
        changes.
        """
        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + viewer.data_off),
                                    int(data_y + viewer.data_off))

        except Exception:
            value = None

        fits_x, fits_y = data_x + 1, data_y + 1
        
        try:
            text = "X: %.1f  Y: %.1f  Value: %.2f" % (fits_x, fits_y, float(value))
            #self.readout.setText('this is another test')
            self.readout.setText(text)
        except:
            pass
    '''    
    def cursor_down(self, fitsimage, event, data_x, data_y):
        print('clicked at ',data_x, data_y)
        #self.mouse_clicked.emit(str(data_x)+','+str(data_y))
    '''    
    def key_press_cb(self, canvas, keyname):
        print('key pressed! ',keyname)
        data_x, data_y = self.fitsimage.get_last_data_xy()
        self.key_pressed.emit(keyname+','+str(data_x)+','+str(data_y))
        #return self.imexam_cmd(self.canvas, keyname, data_x, data_y, func)

        
class maskwindow(Ui_maskWindow, QtCore.QObject):
    mask_saved = QtCore.pyqtSignal(str)
    def __init__(self, MainWindow, logger, image=None, haimage=None, sepath=None, config=None, threshold=0.05,snr=2,cmap='gist_heat_r'):
        super(maskwindow, self).__init__()
        
        self.ui = Ui_maskWindow()
        self.ui.setupUi(MainWindow)
        MainWindow.setWindowTitle('makin a mask...')
        self.MainWindow = MainWindow

        #self.readout = QtWidgets.QLabel('this is a test')
        #self.ui.gridLayout_2.addWidget(self.readout, 0,0,1,2)
        #self.readout.setText('this is another test')

        self.logger = logger

        ###  The lines below are for testing purposes
        ###  and should be removed before release.
        if image == None:
            image='MKW8-18216-R.fits'
        if haimage == None:
            haimage='MKW8-18216-CS.fits'
        if sepath == None:
            sepath=os.getenv('HOME')+'/github/halphagui/astromatic/'
        if config == None:
            config='default.sex.HDI.mask'
        self.image_name = image
        self.haimage_name = haimage
        print(self.image_name)
        print(self.haimage_name)
        print(sepath)
        self.sepath = sepath
        self.config = config
        self.threshold = threshold
        self.snr = snr
        self.snr_analysis = snr
        self.cmap = cmap
        self.xcursor_old = -99
        self.xcursor = -99
        self.mask_size = 20.
        # create name for output mask file
        t = self.image_name.split('.fit')
        self.mask_image=t[0]+'-mask.fits'
        self.mask_inv_image=t[0]+'-inv-mask.fits'
        print('saving image as: ',self.mask_image)
        
        # read in image and define center coords
        self.image, self.imheader = fits.getdata(self.image_name,header = True)
        self.ymax,self.xmax = self.image.shape
        self.xc = self.xmax/2.
        self.yc = self.ymax/2.

        self.v1,self.v2=scoreatpercentile(self.image,[5.,99.5])
        self.adjust_mask = True
        self.figure_size = (10,5)
        self.mask_size = 20. # side of square to mask out when user clicks on a pixel

        # set up array to store the user-created object masks

        self.usr_mask = np.zeros_like(self.image)
        print(self.image.shape, self.usr_mask.shape)
        # set off center flag as false by default
        self.off_center_flag = False

        # keep track of extra objects that the user deletes from mask

        self.deleted_objects = []
        self.link_files()
        self.add_cutout_frames()
        self.runse()
        self.display_cutouts()
        self.connect_buttons()
 
    def connect_buttons(self):
        #self.ui.msaveButton.clicked.connect(self.write_mask)
        self.ui.mquitButton.clicked.connect(self.quit_program)
        self.ui.mhelpButton.clicked.connect(self.print_help_menu)
        self.ui.mrunSEButton.clicked.connect(self.runse)
        #self.ui.msaveButton.clicked.connect(self.save_mask)
        #self.ui.mremoveButton.clicked.connect(self.remove_object)
        self.ui.boxSizeLineEdit.textChanged.connect(self.set_box_size)
        self.ui.seThresholdLineEdit.textChanged.connect(self.set_threshold)
        self.ui.seSNRLineEdit.textChanged.connect(self.set_sesnr)
        self.ui.seSNRAnalysisLineEdit.textChanged.connect(self.set_sesnr_analysis)
    def close_window(self):
        print('click red x to close window')
        #sys.exit()
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
        #self.maskcutout.mouse_clicked.connect(self.add_object)

        # this allows the user to press editing keys in any of the 3 image panels
        # not just in the mask panel
        self.maskcutout.key_pressed.connect(self.key_press_func)
        self.rcutout.key_pressed.connect(self.key_press_func)
        self.hacutout.key_pressed.connect(self.key_press_func)

    def display_cutouts(self):
        self.rcutout.load_file(self.image_name)
        self.rcutout.fitsimage.set_autocut_params('stddev')
        self.hacutout.load_file(self.haimage_name)
        self.display_mask()
    def display_mask(self):
        self.maskcutout.load_file(self.mask_image)
    def link_files(self):
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory
        sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('ln -s '+self.sepath+'/'+file+' .')
    def clean_links(self):
        # clean up symbolic links to sextractor files
        # sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('unlink '+file)

    def read_se_cat(self):
        sexout=fits.getdata('test.cat')
        self.xsex=sexout['XWIN_IMAGE']
        self.ysex=sexout['YWIN_IMAGE']
        self.fwhm = sexout['FWHM_IMAGE']
        dist=np.sqrt((self.yc-self.ysex)**2+(self.xc-self.xsex)**2)
        #   find object ID
        objIndex=np.where(dist == min(dist))
        objNumber=sexout['NUMBER'][objIndex]
        return objNumber[0] # not sure why above line returns a list

    def runse(self,galaxy_id = None):
        print('using a deblending threshold = ',self.threshold)
        print('se %s -c %s -CATALOG_NAME test.cat -CATALOG_TYPE FITS_1.0 -DEBLEND_MINCONT %f -DETECT_THRESH %f -ANALYSIS_THRESH %f'%(self.image_name,self.config,float(self.threshold),float(self.snr),float(self.snr_analysis)))
        os.system('se %s -c %s -CATALOG_NAME test.cat -CATALOG_TYPE FITS_1.0 -DEBLEND_MINCONT %f -DETECT_THRESH %f -ANALYSIS_THRESH %f'%(self.image_name,self.config,float(self.threshold),float(self.snr),float(self.snr)))
        self.maskdat = fits.getdata('segmentation.fits')
        # grow masked areas
        bool_array = np.array(self.maskdat.shape,'bool')
        #for i in range(len(self.xsex)):
            
            
        if self.off_center_flag:
            print('setting center object to objid ',self.galaxy_id)
            self.center_object = self.galaxy_id
        else:
            self.center_object = self.read_se_cat()
        self.maskdat[self.maskdat == self.center_object] = 0
        # add back the square masked areas that the user added
        self.maskdat = self.maskdat + self.usr_mask
        # remove objects that have already been deleted by user
        if len(self.deleted_objects) > 0:
            for objID in self.deleted_objects:
                self.maskdat[self.maskdat == objID] = 0.
        # write out mask image
        fits.writeto(self.mask_image,self.maskdat,header = self.imheader,overwrite=True)
        invmask = self.maskdat > 0.
        invmask = np.array(~invmask,'i')
        fits.writeto(self.mask_inv_image,invmask,header = self.imheader,overwrite=True)
        self.mask_saved.emit(self.mask_image)
        self.display_mask()


    def show_mask(self):
        if self.nods9:
            plt.close('all')
            self.fig = plt.figure(1,figsize=self.figure_size)
            plt.clf()
            plt.subplots_adjust(hspace=0,wspace=0)
            plt.subplot(1,2,1)
            plt.imshow(self.image,cmap='gray_r',vmin=self.v1,vmax=self.v2,origin='lower')
            plt.title('image')
            plt.subplot(1,2,2)
            #plt.imshow(maskdat,cmap='gray_r',origin='lower')
            plt.imshow(self.maskdat,cmap=self.cmap,origin='lower')
            plt.title('mask')
            plt.gca().set_yticks(())
            #plt.draw()
            plt.show(block=False)

    def key_press_func(self,text):
        key, x, y = text.split(',')
        self.xcursor = float(x)
        self.ycursor = float(y)
        try:
            self.cursor_value = self.maskdat[int(self.ycursor),int(self.xcursor)]
        except IndexError:
            print('out of bounds, try again')
        #print('cursor value = ',self.cursor_value, key)
        if key == 'a':
            self.add_circ_object()
        elif key == 'b':
            self.add_box_object()
        elif key == 'r': 
            print('removing object')
            self.remove_object(int(self.cursor_value))
        elif key == 'o': 
            self.off_center()
        elif key == 'g': 
            self.grow_mask()
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
        print('Click on mask or r/ha image, then enter:\n \t r to remove object in mask at the cursor position;'
              '\n \t a to add CIRCULAR mask at cursor position;'
              '\n \t b to add BOX mask at cursor position;'
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


    def add_box_object(self):
        '''
        this adds a square region
        '''
        print('adding pixels to the mask')
        # mask out a rectangle around click
        # size is given by mask_size
        xmin = int(self.xcursor) - int(0.5*self.mask_size)
        ymin = int(self.ycursor) - int(0.5*self.mask_size)
        xmax = int(self.xcursor) + int(0.5*self.mask_size)
        ymax = int(self.ycursor) + int(0.5*self.mask_size)
        
        # make sure cursor click is not outside of the image
        if ((self.xcursor >= self.xmax) or (self.xcursor <= 0) or (self.ycursor >= self.ymax) or (self.ycursor <= 0)):
            print('you clicked outside the image area')
            return
        
        # make sure mask dimensions are not outside of the image
        xmin = max(0,xmin)
        xmax = min(self.xmax,xmax)
        ymin = max(0,ymin)
        ymax = min(self.ymax,ymax)

        #print('xcursor, ycursor = ',self.xcursor, self.ycursor)
        mask_value = np.max(self.maskdat) + 1
        #print(xmin,xmax,ymin,ymax,self.mask_size)
        self.usr_mask[ymin:ymax,xmin:xmax] = mask_value*np.ones([ymax-ymin,xmax-xmin])
        self.maskdat = self.maskdat + self.usr_mask
        self.save_mask()
        print('added mask object '+str(mask_value))

    def add_circ_object(self):
        print('adding circular obj to the mask, with radius = ',self.mask_size)
        # mask out a rectangle around click
        # size is given by mask_size
        pixel_mask = circle_pixels(float(self.xcursor),float(self.ycursor),float(self.mask_size/2.),self.xmax,self.ymax)

        #print('xcursor, ycursor = ',self.xcursor, self.ycursor)
        mask_value = np.max(self.maskdat) + 1
        #print(xmin,xmax,ymin,ymax,self.mask_size)
        self.usr_mask[pixel_mask] = mask_value*np.ones_like(self.usr_mask)[pixel_mask]
        self.maskdat = self.maskdat + self.usr_mask
        self.save_mask()
        print('added mask object '+str(mask_value))
        
    def remove_object(self, objID):
        '''
        this will remove masked pixels near the cursor

        '''
        #objID = int(input('enter pixel value to remove object in mask'))
        xmin = int(self.xcursor) - int(0.5*self.mask_size)
        ymin = int(self.ycursor) - int(0.5*self.mask_size)
        xmax = int(self.xcursor) + int(0.5*self.mask_size)
        ymax = int(self.ycursor) + int(0.5*self.mask_size)

        if objID == 0:
            return
        else:
            self.maskdat[self.maskdat == objID] = 0.
            self.deleted_objects.append(objID)
        self.save_mask()
        self.display_mask()
    def set_threshold(self,t):
        '''
        adjust threshold used in SE deblending
         (0=lots, 1=no deblend )
        '''
        print('Adjust threshold for SE deblending')
        print('0=lots, 1=no deblend')
        #t = raw_input('enter new threshold')
        try:
            self.threshold = float(t)
            self.runse()

        except ValueError:
            pass
        
    def set_sesnr(self,t):
        #t = raw_input('enter new SNR')
        try:
            self.snr = float(t)
        except ValueError:
            pass
        #self.runse()
                
    def set_sesnr_analysis(self,t):
        #t = raw_input('enter new SNR')
        try:
            self.snr_analysis = float(t)
        except ValueError:
            pass

                
    def set_box_size(self,t):
        # change box size used for adding pixels to mask
        #print('current box size = '+str(self.mask_size))
        #t = input('enter new size for square area to be masked (in pixels)\n')
        try:
            self.mask_size = float(t)
        except:
            print('error reading input')

    def off_center(self):
        t = input('enter object number for target galaxy\n')
        self.off_center_flag = True
        self.galaxy_id = int(t)
        self.runse()

    def quit_program(self):
        self.clean_links()
        self.close_window()

    def save_mask(self):
        #super(maskwindow,self).mask_saved(event)
        print('saving mask: ',self.mask_image)
        #newfile = fits.PrimaryHDU()
        #newfile.data = self.maskdat
        #newfile.header = self.imheader
        #fits.writeto(self.mask_image, newfile.data, header = newfile.header, overwrite=True)

        # send message to main program that the mask is updated
        fits.writeto(self.mask_image, self.maskdat, header = self.imheader, overwrite=True)
        self.mask_saved.emit(self.mask_image)
        self.display_mask()
        
        #print(self.mask_image)
        self.mask_saved.emit(self.mask_image)


    def edit_mask(self):
        self.runse()
        while self.adjust_mask:    
            self.show_mask()
            self.print_menu()
            fits.writeto(self.mask_image,self.maskdat,header = self.imheader,overwrite=True)
            self.mask_saved.emit(self.mask_image)
    def grow_mask(self, size=7):
        # convolve mask with top hat kernel
        #kernel = Tophat2DKernel(5)
        '''
        Convolution: one way to grow the mask is to convolve the image with a kernel

        however, this does not preserve the pixels values of the original
        object, which come from the sextractor segmentation image.

        if the user wants to remove an object, it's much easier to do this
        by segmentation number rather than by pixels (as in the reverse of how we add objects
        to mask).

        Alternative: is to loop over just the masked pixels, and replace all pixels
        within a square region with the masked value at the central pixel.
        This will preserve the numbering from the segmentation image.

        Disadvantage: everything assumes a square shape after repeated calls.

        Alternative is currently implemented.
        
        '''
        #mykernel = np.ones([5,5])
        #kernel = CustomKernel(mykernel)
        #self.maskdat = convolve(self.maskdat, kernel)
        #self.maskdat = np.ceil(self.maskdat)
        nx,ny = self.maskdat.shape
        masked_pixels = np.where(self.maskdat > 0.)
        for i,j in zip(masked_pixels[0], masked_pixels[1]):
            rowmin = max(0,i-int(size/2))
            rowmax = min(nx,i+int(size/2))
            colmin = max(0,j-int(size/2))
            colmax = min(ny,j+int(size/2))
            if rowmax <= rowmin:
                # something is wrong, return without editing mask
                continue
            if colmax <= colmin:
                # something is wrong, return without editing mask
                continue
            #print(i,j,rowmin, rowmax, colmin, colmax)
            self.maskdat[rowmin:rowmax,colmin:colmax] = self.maskdat[i,j]*np.ones([rowmax-rowmin,colmax-colmin])
        self.display_mask()
        # save convolved mask as new mask
        self.save_mask()

        # plot mpl figure
        # this was for debugging purposes
        #plt.figure()
        #plt.imshow(self.maskdat)
        #plt.show()
        
# run sextractor on input image
# return segmentation image with central object removed



    


if __name__ == "__main__":
    #catalog = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    #from halphamain import cutout_image
    #from halphamain import cutout_image
    logger = log.get_logger("masklog", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    MainWindow = QtWidgets.QWidget()
    ui = maskwindow(MainWindow, logger)
    #ui.setupUi(MainWindow)
    #ui.test()

    MainWindow.show()
    sys.exit(app.exec_())

    
