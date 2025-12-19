#!/usr/bin/env python 
"""
USAGE:

Virgo 2019 INT data, running on laptop

cd github/halphagui

source venv/bin/activate

cd /data-pool/Halpha/halphagui-output-20230626

%run ~/github/halphagui/testing/halphamain.py --virgo --rimage /home/rfinn/data/reduced/virgo-coadds-feb2019-int/VF-118.1817+20.9822-INT-20190205-p001-r-shifted.fits --haimage /home/rfinn/data/reduced/virgo-coadds-feb2019-int/VF-118.1817+20.9822-INT-20190205-p001-Halpha.fits --filter inthalpha --psfdir /home/rfinn/data/reduced/psf-images/ --tabledir /home/rfinn/research/Virgo/tables-north/v1/ --auto


Testing after MVC semi-implementation.  I created a testing directory on the linux laptop, and this is command I used to run the gui:

(venv) (base) rfinn@virgof:~/research/Virgo-dev/halphagui-test$ 

python ~/github/halphagui/halphamain.py --virgo --tabledir ~/research/Virgo/tables-north/v2/ --rimage VF-145.781+31.887-HDI-20180313-p019-R.fits --haimage VF-145.781+31.887-HDI-20180313-p019-ha4.fits --csimage VF-145.781+31.887-HDI-20180313-p019-ha4-CS-ZP.fits --psfdir ~/research/Virgo-dev/halphagui-test/ --filter ha4 --prefix VF-145.781+31.887-HDI-20180313-p019

"""


import sys, os
sys.path.append(os.getcwd())
#sys.path.append(os.getenv('HOME')+'github/HalphaImaging/')
sys.path.append(os.getenv('HOME')+'/github/HalphaImaging/python3/')

import numpy as np
import platform

#from PyQt5 import  QtWidgets
#from PyQt5 import QtCore
#from PyQt5.Qtcore import  Qt
from PyQt5 import QtCore,QtWidgets, QtGui
#from ginga.qtw.QtHelp import QtGui #, QtCore
from halphav5 import Ui_MainWindow
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
#from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
from ginga import colors
from ginga.canvas.CanvasObject import get_canvas_types
from ginga.misc import log
from ginga.util.loader import load_data

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, FK5
import astropy.units as u
from astropy import nddata
from astropy.table import Table, Column
from astropy.visualization import simple_norm
from astropy.cosmology import WMAP9 as cosmo
from astropy.time import Time

# packages for ellipse fitting routine
# https://photutils.readthedocs.io/en/stable/isophote.html
#from photutils.isophote import EllipseGeometry
#from photutils.isophote import Ellipse


#from photutils import EllipticalAperture



import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from datetime import date
import time
# routines for measuring elliptical photometry
from photwrapper import ellipse

from maskwrapper import maskwindow

from galfitwrapper import galfitwindow

from buildpsf import psf_parent_image
from halphaCommon import cutout_image

from fit_profile import profile, dualprofile, rprofile, haprofile, ratio_error
# code from HalphaImaging repository
import filterratio as runse

# code for calculating redshift cutoffs of filter
# and for calculating the transmission correction for each galaxy
# based on where it falls ini the filter bandpass
from filter_transmission import filter_trace

from join_catalogs import join_cats, make_new_cats

import imutils

#from uat_mask import mask_image
# filter information

# INT filters from https://astro.ing.iac.es/filter/list.php?instrument=WFC
# lmin is center - 0.5 FWHM
# lmin is center + 0.5 FWHM
#
# Filter # 197 is INT197 = (6568, 95)
# Filter # 227 is INT227 = (6657, 80)
#
# BOK using the NOAO filter, so Halpha 4
# 
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.,'INT197':6520.5, 'INT227':6617}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.,'INT197':6615.5, 'INT227':6697}

# Force a specific toolkit on mac
macos_ver = platform.mac_ver()[0]
try:
    matplotlib.use('Qt5Agg')
except ImportError:
    print("WARNING! could not load Qt5Agg")
    
import matplotlib.pyplot as plt

### code to measure galaxy sizes from photutils segmentation image
from get_galaxy_size import getobjectsize

# default size for cutouts, multiple of NSA PETROTH90
cutout_scale = 14

# now in terms of R25 for
cutout_scale = 2.5

mask_scalefactor = 1 # number to multiple R24 by
######################################################
## FUNCTIONS
######################################################

def get_params_from_name(image_name):
    t = os.path.basename(image_name).split('-')
    #print(t)
    if len(t) == 5:
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]
    elif len(t) == 6: # meant to catch negative declinations
        telescope = t[3]
        dateobs = t[4]
        pointing = t[5]
    #else:
    #    print("ruh roh - trouble with get_params_from_name for image ",image_name, len(t))
        #print(image_name)
        #print(t)
    return telescope,dateobs,pointing

def get_params_from_name_uat(image_name):
    t = os.path.basename(image_name).split('_')
    #print(t)
    if len(t) == 7: # meant to catch negative declinations
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]+'_'+t[5]
    elif len(t) == 8:
        telescope = t[3]
        dateobs = t[4]
        pointing = t[4]+'_'+t[5]
        
    else:
        print("ruh roh - trouble getting info from ",image_name, len(t))
        print(image_name)
        print(t)
        return
    return telescope,dateobs,pointing

class psfimage():
    def __init__(self):
        fwhm = 5.6
        fwhm_arcsec = 2.0


class image_panel(QtCore.QObject):#(QtGui.QMainWindow,
    key_pressed = QtCore.pyqtSignal(str)
    def __init__(self,panel_name,ui,logger):
        super(image_panel, self).__init__()
        QtCore.QObject.__init__(self)
        self.ui = ui
        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()
        #self.figure = plt.figure()
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        #fi.set_autocut_params('histogram')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('once')
        fi.set_callback('drag-drop', self.drop_file)
        fi.set_callback('none-move',self.cursor_cb)
        # not sure how to add tab for multiple images.
        # going to try using keystroke instead
        #fi.add_callback('add-channel',self.add_channel)
        #fi.add_callback('channel-change', self.focus_cb)

        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        #fi.set_figure(self.figure)
        self.fitsimage = fi
        self.fitsimage.set_callback('key-press',self.key_press_cb)
        # enable some user interaction
        #fi.get_bindings.enable_all(True)
        bd = fi.get_bindings()
        bd.enable_all(True)
        
        w = fi.get_widget()
        #w.resize(512, 512)

        # add scrollbar interface around this viewer
        si = ScrolledView(fi)


        panel_name.addWidget(w,2,0,7,1)
        #panel_name.setMinimumSize(QtCore.QSize(512, 512))     
        # canvas that we will draw on
        canvas = self.dc.DrawingCanvas()
        canvas.enable_draw(True)
        canvas.enable_edit(True)
        canvas.set_drawtype('rectangle', color='lightblue')
        canvas.set_surface(fi)
        #canvas.rectangle(.5,.5,10,0)
        self.canvas = canvas
        # add canvas to view
        #fi.add(canvas)
        private_canvas = fi.get_canvas()
        private_canvas.add(canvas)
        canvas.register_for_cursor_drawing(fi)
        canvas.add_callback('draw-event', self.draw_cb)
        canvas.set_draw_mode('draw')
        canvas.ui_set_active(True)
        self.canvas = canvas

        self.drawtypes = canvas.get_drawtypes()
        self.drawtypes.sort()

        # add a color bar
        #fi.show_color_bar(True)
        fi.show_focus_indicator(True)

        self.readout = QtWidgets.QLabel('this is a test')
        panel_name.addWidget(self.readout)
        self.readout.setText('this is another test')
        
        #wdrawcolor = QtGui.QComboBox()
        wdrawcolor = QtWidgets.QComboBox()
        for name in self.drawcolors:
            wdrawcolor.addItem(name)
        index = self.drawcolors.index('lightblue')
        wdrawcolor.setCurrentIndex(index)
        wdrawcolor.activated.connect(self.set_drawparams)
        self.wdrawcolor = wdrawcolor
        
        wdrawtype = QtWidgets.QComboBox()
        for name in self.drawtypes:
            wdrawtype.addItem(name)
        index = self.drawtypes.index('rectangle')
        wdrawtype.setCurrentIndex(index)
        wdrawtype.activated.connect(self.set_drawparams)
        self.wdrawtype = wdrawtype
        ui.wclear.clicked.connect(self.clear_canvas)

        #ui.wopen.clicked.connect(self.open_file)

        '''
        # add little mode indicator that shows keyboard modal states
        fi.show_mode_indicator(True, corner='ur')
        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))
        wclear = QtGui.QPushButton("Clear Canvas")
        wclear.clicked.connect(self.clear_canvas)
        wopen = QtGui.QPushButton("Open File")
        

        #wgalaxies = QtGui.QPushButton("Mark Galaxies")
        #wgalaxies.clicked.connect(self.mark_galaxies)
        
        #gridLayout = QtWidgets.QGridLayout()
        panel_name.addWidget(wopen)#,4,4,4,1)
        panel_name.addWidget(wgalaxies)#,4,4,4,2)
        panel_name.addWidget(wdrawcolor)#,4,4,4,3)
        panel_name.addWidget(wdrawtype)#,4,4,4,4)
        #panel_name.addWidget(wclear)
        '''

        #self.add_cutouts()

 
    def add_channel(self):
        print('adding a channel')
    def key_press_cb(self, canvas, keyname):
        #print('key pressed! ',keyname)
        self.key_pressed.emit(keyname)
        
        
    def set_drawparams(self, kind):
        index = self.wdrawtype.currentIndex()
        kind = self.drawtypes[index]
        index = self.wdrawcolor.currentIndex()
        fill = (self.wfill.checkState() != 0)
        alpha = self.walpha.value()

        params = {'color': self.drawcolors[index],
                  'alpha': alpha,
                  }
        if kind in ('circle', 'rectangle', 'polygon', 'triangle',
                    'righttriangle', 'ellipse', 'square', 'box'):
            params['fill'] = fill
            params['fillalpha'] = alpha

        self.canvas.set_drawtype(kind, **params)
    def clear_canvas(self):
        self.canvas.delete_all_objects()
        
        
    def load_file(self, filepath):
        image = load_data(filepath, logger=self.logger)
        # RAF - updating for new WCS call
        header = fits.getheader(filepath)
        self.image_wcs = WCS(header)
        self.fitsimage.set_image(image)
        #t = fits.getdata(image)
        #v1 = np.scoreatpercentile(image,1)
        #v2 = np.scoreatpercentile(image,99.5)
        #self.fitsimage.auto_levels(autocuts=(v1,v2))
        #self.setWindowTitle(filepath)
        self.coadd_filename = filepath

    def open_file(self):
        res = QtWidgets.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                ".", "FITS files (*.fits)")
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    def drop_file(self, fitsimage, paths):
        fileName = paths[0]
        self.load_file(fileName)

        
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
        
        # Calculate WCS RA
        #print(fits_x, fits_y)
        #ra_txt, dec_txt = self.image_wcs.wcs_pix2world(fits_x, fits_y,1)
        try:
            # NOTE: image function operates on DATA space coords            
            ra_txt, dec_txt = self.image_wcs.wcs_pix2world(fits_x, fits_y,1)
        except Exception as e:
            self.logger.warning("Bad coordinate conversion: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'
        try:
            text = "RA: %.6f  DEC: %.4f  X: %.2f  Y: %.2f  Value: %.3f" % ((1.*ra_txt), (1.*dec_txt), fits_x, fits_y, float(value))
            self.readout.setText(text)
        except:
            pass
        # WCS stuff is not working so deleting for now...
        #text = "X: %.2f  Y: %.2f  Value: %s" % (fits_x, fits_y, value)


    def set_mode_cb(self, mode, tf):
        self.logger.info("canvas mode changed (%s) %s" % (mode, tf))
        if not (tf is False):
            self.canvas.set_draw_mode(mode)
        return True

    def draw_cb(self, canvas, tag):
        obj = canvas.get_object_by_tag(tag)
        obj.add_callback('pick-down', self.pick_cb, 'down')
        obj.add_callback('pick-up', self.pick_cb, 'up')
        obj.add_callback('pick-move', self.pick_cb, 'move')
        obj.add_callback('pick-hover', self.pick_cb, 'hover')
        obj.add_callback('pick-enter', self.pick_cb, 'enter')
        obj.add_callback('pick-leave', self.pick_cb, 'leave')
        obj.add_callback('pick-key', self.pick_cb, 'key')
        obj.pickable = True
        obj.add_callback('edited', self.edit_cb)

    def pick_cb(self, obj, canvas, event, pt, ptype):
        self.logger.info("pick event '%s' with obj %s at (%.2f, %.2f)" % (
            ptype, obj.kind, pt[0], pt[1]))
        return True

    def edit_cb(self, obj):
        self.logger.info("object %s has been edited" % (obj.kind))
        return True

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.deleteLater()
        
    def add_galaxies(self):
        ax = self.figure.gca()
        r = patches.Rectangle((wd * 0.10, ht * 0.10), wd * 0.6, ht * 0.5, ec='b',
                      fill=False)
        ax.add_patch(r)

class output_table_view():
    """ methods that interact with the gui  """
    
    def update_gui_table(self):

        self.ui.tableWidget.setColumnCount(len(self.table.columns))
        self.ui.tableWidget.setRowCount(len(self.table))
        
        self.ui.tableWidget.setHorizontalHeaderLabels(self.table.colnames)

        '''
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(4, item)
        '''
        if self.igal is not None:
            #print(self.ui.commentLineEdit.text())
            self.table['COMMENT'][self.igal] = str(self.ui.commentLineEdit.text())
            
        for col, c in enumerate(self.table.columns):
            #item = self.ui.tableWidget.horizontalHeaderItem(col)
            #item.setText(_translate("MainWindow", self.table.columns[col].name))
            for row in range(len(self.table[c])):
                item = self.table[row][col]
                self.ui.tableWidget.setItem(row,col,QtWidgets.QTableWidgetItem(str(item)))
        #item = self.tableWidget.horizontalHeaderItem(0)
        #item.setText(_translate("MainWindow", "ID"))
        self.write_fits_table()

    def update_gui_table_cell(self,row,col,item):
        # row will be igal, so that's easy
        # how to make it easy to get the right column number?
        #
        # easiest is to pass the column name, and
        # then match the column name to find the

        colmatch = False
        for i,c in enumerate(self.table.colnames):
            if c == col:
                ncol = i
                colmatch = True
                break
        if colmatch:    
            self.ui.tableWidget.setItem(row,ncol,QtWidgets.QTableWidgetItem(str(item)))
            self.table[row][col]=item
        else:
            print('could not match column name ',col)
        self.write_fits_table()

        
class create_output_table(output_table_view):
    """
    output table that stores all of the measured values for each galaxy in FOV
    """
    def initialize_results_table(self, prefix=None,virgo=False,uat=False,nogui=False):
        #print("in initialize_results_table, value of uat = ", uat)
        #print("in initialize_results_table, value of virgo = ", virgo)        
        self.nogui = nogui
        
        '''
        Data to store:
        - NSAID
        * AGC number
        - RA
        - DEC
        - filter_ratio
        - cutout_size
        - xmin:xmax,ymin:ymax from parent image
        - ha_flag -- boolean
        - ha_class -- category of Halpha emission
        - psf_fwhm
        - galfit re
        - galfit n
        - galfit BA
        - galfit PA
        - galfit (xc,yc)
        - galfit (RA,DEC) - translate the pixel coords to RA and Dec of galaxy center
        - galfit mag
        - galfit sky
        - ellipse PA
        - ellipse BA
        - ellipse Gini
        - ellipse skynoise
        - ellipse mag R
        - ellipse mag Ha
        - ellipse SFR Ha
        - profiles Re r
        - profiles Re Ha
        - becky inner ssfr
        - becky outer ssfr
        - becky C30
        - becky C70
        '''
        ## define fits table output name
        # get directory after Users - this should be username for 
        user = os.getenv('USER')
        today = date.today()
        str_date_today = today.strftime('%Y-%b-%d')
        if prefix is None:
            self.output_table = 'halpha-data-'+user+'-'+str_date_today+'.fits'
        else:
            self.output_table = prefix+'-data-'+user+'-'+str_date_today+'.fits'
        ## check for existing table
        ##
        ## load if it exists

        # why is this block at the beginning???
        #if virgo:
        #    self.create_table_virgo()
        #elif uat:
        #    print("just to check, I am running create_table_uat")
        #    self.create_table_uat()
        #else:
        #    self.create_table()


        
        if os.path.exists(self.output_table):
            if virgo:
                self.read_table_virgo()
            elif uat:
                self.read_table_uat()
            else:
                self.read_table()
                self.agc2 = fits.getdata(self.prefix+'-agc-matched.fits')
        ## if not, create table                
        else:
            if virgo:
                self.create_table_virgo()
            elif uat:
                self.create_table_uat()
            else:
                self.create_table()
            # call other methods to add columns to the table
            self.add_part1()
            # skipping for now b/c this will have to be different for virgo
            #self.add_nsa()
            self.add_cutout_info()
            self.add_galfit_r()
            #self.add_galfit_ha()            
            self.add_ellipse()
            self.add_profile_fit()
            self.add_photutils()
            self.add_statmorph()            
        if not self.nogui:
            self.update_gui_table()
                
    def read_table():
        ''' read in output from previous run, if it exists'''
        self.table = Table(fits.getdata(self.output_table))
        self.gredshift = self.table['REDSHIFT']
        self.ngalaxies = len(self.gredshift)
        self.ra = self.table['NSA_RA']*self.table['NSA_FLAG']+ self.table['AGC_RA']*(~self.table['NSA_FLAG'])
        self.dec = self.table['NSA_DEC']*self.table['NSA_FLAG']+ self.table['AGC_DEC']*(~self.table['NSA_FLAG'])
        self.gradius = self.table['SERSIC_TH50']*self.table['NSA_FLAG']/self.pixelscale + 100.*np.ones(self.ngalaxies)*(~self.table['NSA_FLAG'])
        self.gzdist = self.table['ZDIST']
        charar1 = np.chararray(self.ngalaxies)
        charar1[:] = 'N'
        charar2 = np.chararray(self.ngalaxies)
        charar2[:] = '-A'
        self.galid=np.zeros(self.ngalaxies, dtype='U15')
        for i in np.arange(self.ngalaxies):
            self.galid[i] = 'N'+str(self.table['NSAID'][i])+'-A'+str(self.table['AGCNUMBER'][i])
        # read in nsa2
        self.nsa2 = fits.getdata(self.prefix+'-nsa-matched.fits')
        # read in agc2
        self.agc2 = fits.getdata(self.prefix+'-agc-matched.fits')                                                                              
        ## if not, create table
    def read_table_virgo(self):
        self.table = Table(fits.getdata(self.output_table))
        self.gredshift = self.table['REDSHIFT']
        self.ngalaxies = len(self.gredshift)
        self.ra = self.table['RA']
        self.dec = self.table['DEC']
        self.gradius = self.table['radius']
        self.gzdist = self.table['ZDIST']
        self.galid=self.table['VFID']
        self.NEDname=self.table['NEDname']
        self.gprefix=self.table['prefix']                

    def read_table_uat(self):
        self.table = Table(fits.getdata(self.output_table))
        self.gredshift = self.table['vopt']/3.e5
        self.ngalaxies = len(self.gredshift)
        self.ra = self.table['RA']
        self.dec = self.table['DEC']
        self.gradius = self.table['a']
        self.gzdist = self.table['vopt']/3.e5
        self.galid= self.table['AGCnr']
        self.NEDname= self.table['AGCnr']
        self.haflag= self.table['AGCnr']        
        try:
            self.gprefix= self.table['prefix']
        except:
            print("WARNING: no prefix in gal table")
            self.gprefix = None
        
    def create_table_virgo(self):
        # updating this part for virgo filament survey 

        self.table = self.defcat.cat['VFID','RA','DEC','vr','radius','NEDname','prefix']
        self.table['VFID'].description = 'ID from Virgo Filament catalog'                
        self.table['RA'].unit = u.deg
        self.table['RA'].description = 'RA from VF catalog'        
        self.table['DEC'].unit = u.deg
        self.table['DEC'].description = 'DEC from VF catalog'                
        self.table['vr'].unit = u.km/u.s
        self.table['vr'].description = 'recession velocity from VF catalog'
        self.table['radius'].unit = u.arcsec
        self.table['radius'].description = 'radius from VF catalog'        
        self.ngalaxies = len(self.table)
        #print('number of galaxies = ',self.ngalaxies)
        self.haflag = np.zeros(self.ngalaxies,'bool')
        self.galid = self.table['VFID']
        self.NEDname = self.table['NEDname']                
        self.gredshift = self.defcat.cat['vr']/3.e5
        self.gzdist = self.defcat.cat['vr']/3.e5

        ##
        # update this to use the SMA_SB24
        ##
        #self.gradius = self.defcat.cat['radius']/self.pixelscale
        self.gradius = self.radius_arcsec/self.pixelscale

        
        self.ra = self.defcat.cat['RA']
        self.dec = self.defcat.cat['DEC']        
        c1 = Column(self.haflag, name='HAflag', description='Halpha flag')
        c2 = Column(self.gredshift, name='REDSHIFT', description='redshift')
        c3 = Column(self.gredshift, name='ZDIST', description='redshift')        
        self.table.add_columns([c1,c2,c3])

    def create_table_uat(self):
        # updating this part for uat Halpha Groups

        self.table = self.defcat.cat['AGCnr','RA','DEC','vopt','a','b','v21','hiflux',]
        self.table['AGCnr'].description = 'Number in AGC catalog'
        self.table['RA'].unit = u.deg
        self.table['RA'].description = 'RA from AGC catalog'        
        self.table['DEC'].unit = u.deg
        self.table['DEC'].description = 'DEC from AGC catalog'                
        self.table['vopt'].unit = u.km/u.s
        self.table['vopt'].description = 'optical recession velocity from AGC catalog'
        self.table['v21'].unit = u.km/u.s
        self.table['v21'].description = 'HI recession velocity from AGC catalog'
        self.table['a'].unit = u.arcmin
        self.table['a'].description = 'AGC semi-major axis'
        self.table['b'].unit = u.arcmin
        self.table['b'].description = 'AGC semi-minor axis'
        self.ngalaxies = len(self.table)
        #print('number of galaxies = ',self.ngalaxies)
        self.haflag = np.zeros(self.ngalaxies,'bool')
        self.galid = self.table['AGCnr']
        self.NEDname = None #self.table['NEDname']                
        self.gredshift = self.defcat.cat['vopt']/3.e5
        self.gzdist = self.defcat.cat['vopt']/3.e5

        ##
        # update this to use the SMA_SB24
        ##
        #self.gradius = self.defcat.cat['radius']/self.pixelscale
        self.gradius = self.radius_arcsec/self.pixelscale

        
        self.ra = self.defcat.cat['RA']
        self.dec = self.defcat.cat['DEC']        
        c1 = Column(self.haflag, name='HAflag', description='Halpha flag')
        c2 = Column(self.gredshift, name='REDSHIFT', description='redshift')
        c3 = Column(self.gredshift, name='ZDIST', description='redshift')        
        self.table.add_columns([c1,c2,c3])
        
    def create_table(self):
        # updating this part to make use of NSA and AGC catalogs
        # not going to make this backward compatible, meaning you need to enter both catalogs
        # probably just being lazy, but it's giving me a headache...

        # much better approach would be to match entire NSA and AGC catalogs
        # and then just use that,
        # but forging ahead for now.



        # this returns new nsa and agc catalogs, that are row matched to the joined table
        if self.agcflag:
            self.nsa2, self.nsa_matchflag, self.agc2, self.agc_matchflag = make_new_cats(self.nsa.cat, self.agc.cat)
            # write matched tables
            self.nsa2.write(self.prefix+'-nsa-matched.fits', format='fits', overwrite=True)
            self.agc2.write(self.prefix+'-agc-matched.fits', format='fits', overwrite=True)

            self.ngalaxies = len(self.nsa_matchflag)        
            # create arrays that we need for other parts of the programs, like ra, dec, size
            self.ra = self.nsa2['RA']*self.nsa_matchflag + self.agc2['RA']*(~self.nsa_matchflag)

            self.dec = self.nsa2['DEC']*self.nsa_matchflag + self.agc2['DEC']*(~self.nsa_matchflag)
 
            self.gradius = self.nsa2['SERSIC_TH50']*self.nsa_matchflag/self.pixelscale + 100.*np.ones(self.ngalaxies)*(~self.nsa_matchflag)
            charar1 = np.chararray(self.ngalaxies)
            charar1[:] = 'N'
            charar2 = np.chararray(self.ngalaxies)
            charar2[:] = '-A'
            self.galid=np.zeros(self.ngalaxies, dtype='U15')
            for i in np.arange(self.ngalaxies):
                self.galid[i] = 'N'+str(self.nsa2['NSAID'][i])+'-A'+str(self.agc2['AGCnr'][i])
            voptflag = self.agc2['vopt'] > 0.
            agcredshift = self.agc2['vopt']/3.e5*voptflag + self.agc2['v21']/3.e5*(~voptflag)
            self.gredshift = self.nsa2['Z']*self.nsa_matchflag + agcredshift*(~self.nsa_matchflag)
        else:
            self.nsa2 = self.nsa.cat
            self.nsa_matchflag = np.ones(len(self.nsa2))
            self.ra = self.nsa2['RA']

            self.dec = self.nsa2['DEC']

            self.gradius = self.nsa2['SERSIC_TH50']*self.nsa_matchflag/self.pixelscale 
            self.ngalaxies = len(self.nsa_matchflag)
            self.galid = self.nsa2['NSAID']*self.nsa_matchflag
            self.gredshift = self.nsa2['Z']*self.nsa_matchflag
            
        # number of galaxies in the joined table
        self.haflag = np.zeros(self.ngalaxies,'bool')
        c0 = Column(self.galid,name='ID')
        c1 = Column(self.gredshift,name='REDSHIFT')
        c2 = Column(self.nsa2['NSAID'], name='NSAID',dtype=np.int32, description='NSAID')
        c3 = Column(self.nsa_matchflag, name='NSA_FLAG',dtype='bool', description='NSA_FLAG')
        c4 = Column(self.nsa2['RA'], name='NSA_RA',dtype='f', unit=u.deg)
        c5 = Column(self.nsa2['DEC'], name='NSA_DEC',dtype='f', unit=u.deg)
        self.table = Table([c0,c1,c2,c3,c4,c5])
        if self.agcflag:
            c5 = Column(self.agc2['AGCnr'], name='AGCNUMBER',dtype=np.int32, description='AGC ID NUMBER')
            c6 = Column(self.agc_matchflag, name='AGC_FLAG',dtype='bool', description='AGC_FLAG')
            c7 = Column(self.agc2['RA'], name='AGC_RA',dtype='f', unit=u.deg)
            c8 = Column(self.agc2['DEC'], name='AGC_DEC',dtype='f', unit=u.deg)
            self.table.add_columns([c5,c6,c7,c8])
        
    def add_part1(self):
        g1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_RA', unit=u.deg,description='R-band center RA from galfit')
        g2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_DEC', unit=u.deg,description='R-band center DEC from galfit')
        g3 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HRA', unit=u.deg,description='HA center RA from galfit')
        g4 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HDEC', unit=u.deg,description='HA center DEC from galfit')
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_RA', unit=u.deg,description='R-band center RA from photutil centroid')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_DEC', unit=u.deg,description='R-band center DEC from photutil centroid')
        e1a = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HRA', unit=u.deg,description='Halpha center RA from photutil centroid')
        e2a = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HDEC', unit=u.deg,description='Halpha center DEC from photutil centroid')
        c9 = Column(self.haflag, name='HA_FLAG',description='shows HA emission')
        c10 = Column(np.ones(self.ngalaxies,'f'),name='FILT_COR',unit='', description='max filt trans/trans at gal z')
        c11 = Column(np.zeros(self.ngalaxies,'f'),name='R_FWHM',unit=u.arcsec, description='R FWHM in arcsec')
        c12 = Column(np.zeros(self.ngalaxies,'f'),name='H_FWHM',unit=u.arcsec, description='HA FWHM in arcsec')
        c13 = Column(np.zeros(self.ngalaxies,dtype='|S40'),name='POINTING', description='string specifying year and pointing')
        c14 = Column(np.zeros(self.ngalaxies,dtype='|S3'),name='TEL', description='telescope/instrument')
        c15 = Column(np.zeros(self.ngalaxies,dtype='i'),name='DATE-OBS', description='string specifying date of observation')                        

        self.table.add_columns([g1,g2,g3,g4,e1,e2,e1a,e2a,c9,c10,c11,c12,c13,c14,c15])
    def add_nsa(self):
        # add some useful info from NSA catalog (although matching to NSA could be done down the line)
        r = 22.5 - 2.5*np.log10(self.nsa2['NMGY'][:,4])
        c11 = Column(r,name='NSA_RMAG',unit=u.mag,description='NSA r mag')
        c12 = Column(self.nsa2['SERSIC_TH50'],name='SERSIC_TH50', unit=u.arcsec,description='NSA SERSIC_TH50')
        c13 = Column(self.nsa2['SERSIC_N'],name='SERSIC_N',description='NSA SERSIC index')
        c14 = Column(self.nsa2['SERSIC_BA'],name='SERSIC_BA',description='NSA SERSIC B/A')
        c15 = Column(self.nsa2['SERSIC_PHI'],name='SERSIC_PHI', unit=u.deg,description='NSA SERSIC PHI')
        self.gzdist = self.nsa2['ZDIST']*self.nsa_matchflag + self.gredshift*~self.nsa_matchflag
        c16 = Column(self.gzdist,name='ZDIST',description='NSA ZDIST')
        self.table.add_columns([c11,c12,c13,c14,c15,c16,])
    def add_cutout_info(self):
        # cutout region in coadded images
        c1 = Column(np.zeros(len(self.table),dtype='U22'), name='BBOX',description='location of galaxy cutout in mosaic')
        # R-band scale factor for making continuum-subtracted image
        c2 = Column(np.zeros(len(self.table),'f'), name='FILTER_RATIO',description='R/Ha ratio used in cont subtraction')
        # r-band ZP
        c3 = Column(np.zeros(len(self.table),'f'), name='RZP',description='R-band ZP')
        # Halpha ZP        
        c4 = Column(np.zeros(len(self.table),'f'), name='HZP',description='Halpha ZP')
        c5 = Column(np.zeros(len(self.table),'f'), name='PIXSCALE',description='Pixel scale')                

        self.table.add_columns([c1,c2,c3,c4,c5])
    def add_galfit_r(self):
        ##############################################3
        ### GALFIT R-BAND FITS
        ##############################################3

        fields = ['XC','YC','MAG','RE','N','BA','PA']
        units = ['pixel','pixel','mag','arcsec',None,'deg',None]
        descriptions = ['R-band center from galfit (pix)',\
                        'R-band center from galfit (pix)',\
                        'R-band mag from galfit',\
                        # this is currently written in pixels - need to write out in arcsec
                        #'R-band effective radius from galfit (pix)',\
                        'R-band effective radius from galfit (arcsec)',\
                        'R-band sersic index from galfit',\
                        'R-band axis ratio from galfit',\
                        'R-band position angle from galfit']
        i=0
        for f,unit in zip(fields,units):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR',description='err in '+descriptions[i])
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f, unit=unit,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR', unit=unit,description='err in '+descriptions[i])
            #print(c1)
            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_SKY',unit=u.adu,description='sky from galfit')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_CHISQ',description='chisq of galfit sersic model')
        #c3 = Column(np.zeros(self.ngalaxies,'f'), name='GAL_GINI')
        #c4 = Column(np.zeros(self.ngalaxies), name='GAL_GINI2')
        #c5 = Column(np.zeros(self.ngalaxies,'f'), name='GAL_ASYM')
        #c6 = Column(np.zeros(self.ngalaxies,'f'), name='GAL_ASYM2')
        self.table.add_columns([c1,c2])#,c3,c4,c5,c6])

    def add_galfit_2comp_r(self):
        # galfit sersic parameters from 2 comp fit
        c16 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_2SERSIC',description='galfit R-band 2comp fit')
        c17 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_2SERSIC_ERR',description='galfit R-band 2comp fit errors')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_2SERSIC_ERROR',description='galfit R-band 2comp fit num err flag')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_2SERSIC_CHISQ',description='galfit R-band 2comp chi sq')
        self.table.add_columns([c16,c17,c18,c19])

    def add_galfit_1comp_with_asymmetry_r(self):

        # galfit 1 comp with asymmetry
        c16 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_SERSASYM',description='galfit R-band 1comp sersic w/asymmetry')
        c17 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_SERSASYM_ERR')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_ERROR',description='galfit R-band 1comp sersic w/asymmetry num err flag')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_CHISQ',description='galfit R-band 1comp sersic w/asymmetry chi sq')
        c20 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_RA',unit='deg',description='RA from galfit R-band 1comp sersic w/asymmetry')
        c21 = Column(np.zeros(self.ngalaxies), name='GAL_SERSASYM_DEC',unit='deg',description='DEC from galfit R-band 1comp sersic w/asymmetry')
        self.table.add_columns([c16,c17,c18,c19,c20,c21])
    def add_galfit_ha(self):
        ##############################################
        ### GALFIT Halpha FITS
        ##############################################

        fields = ['XC','YC','MAG','RE','N','BA','PA']
        units = ['pixel','pixel','mag','arcsec',None,'deg',None]
        descriptions = ['HA center from galfit (pix)',\
                        'HA center from galfit (pix)',\
                        'HA mag from galfit',\
                        # currently written in pixels
                        # need to convert to arcsec
                        'HA effective radius from galfit (arcsec)',\
                        'HA sersic index from galfit',\
                        'HA axis ratio from galfit',\
                        'HA position angle from galfit']
        i=0
        for f,unit in zip(fields,units):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f+'_ERR',description='err in '+descriptions[i])
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f, unit=unit,description=descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_H'+f+'_ERR', unit=unit,description='err in '+descriptions[i])

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HSKY',description='galfit HA sky')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_HCHISQ',description='galfit chisq of HA model')
        self.table.add_columns([c1,c2])#,c3,c4,c5,c6])

    def add_galfit_2comp_ha(self):
        # galfit sersic parameters from 2 comp fit
        c16 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_H2SERSIC',description='galfit HA 2-comp fit')
        c17 = Column(np.zeros((self.ngalaxies,15),'f'), name='GAL_H2SERSIC_ERR')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_H2SERSIC_ERROR',description='galfit HA 2-comp num error code')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_H2SERSIC_CHISQ',description='galfit HA 2-comp chisq')
        self.table.add_columns([c16,c17,c18,c19])

    def add_galfit_1comp_with_asymmetry_ha(self):
        # galfit 1 comp with asymmetry
        c16 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_HSERSASYM',description='galfit HA model w/asym')
        c17 = Column(np.zeros((self.ngalaxies,10),'f'), name='GAL_HSERSASYM_ERR')
        c18 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_ERROR',description='galfit HA asym num error code')
        c19 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_CHISQ',description='galfit HA asym chisq')
        c20 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_RA',unit='deg')
        c21 = Column(np.zeros(self.ngalaxies), name='GAL_HSERSASYM_DEC',unit='deg')
        self.table.add_columns([c16,c17,c18,c19,c20,c21])

    def add_ellipse(self):
        #####################################################################
        # ellipse output
        # xcentroid, ycentroid, eps, theta, gini, sky_centroid, area, background_mean, source_sum, source_sum_err
        #####################################################################
        e0 = Column(np.zeros(self.ngalaxies,'bool'), name='BADGAL',description='bad galaxy flag - maybe partial coverage')
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_EPS',description='axis ratio from ellipse')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_THETA', unit=u.degree,description='position angle from ellipse')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_GINI',description='gini coeff from ellipse')
        e6 = Column(np.zeros(self.ngalaxies), name='ELLIP_HGINI',description='gini coeff method 2')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_M20',description='M20 for r image')
        e8 = Column(np.zeros(self.ngalaxies), name='ELLIP_HM20',description='M20 for Halpha image ')
        e9 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_UNMASKED_AREA',description='unmasked source area from photutils')
        e9b = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_TOTAL_AREA',description='total source area from photutils')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM', unit = u.erg/u.s/u.cm**2,description='total flux from ellipse')
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_SUM_MAG', unit = u.mag,description='mag from ellipse')
        e12 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM',description='asym from ellipse')
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_ASYM_ERR')
        e14 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM', unit=u.erg/u.s/u.cm**2,description='HA flux from ellipse')
        e15 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HSUM_MAG', unit=u.mag,description='HA mag from ellipse')
        e16 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM',description='HA asymmetry from ellipse')
        e17 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HASYM_ERR')
        e18 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKYNOISE',description='R skynoise in 1E-17 erg/s/cm^2/arcsec^2')
        e19 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKYNOISE',description='HA skynoise in 1E-17  erg/s/cm^2/arcsec^2')
        e20 = Column(np.zeros(self.ngalaxies,'e'), name='R_SKY',description='R sky level in ADU')
        e21 = Column(np.zeros(self.ngalaxies,'e'), name='H_SKY',description='HA sky level in ADU')

        # photutils radii
        e22 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R30',description='photutils R flux frac 30')
        e23 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R50',description='photutils R flux frac 50')
        e24 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_R90',description='photutils R flux frac 90')
        
        e25 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR30',description='photutils Halpha flux frac 30')
        e26 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR50',description='photutils Halpha flux frac 50')
        e27 = Column(np.zeros(self.ngalaxies,'f'), name='ELLIP_HR90',description='photutils Halpha flux frac 90')        

        
        self.table.add_columns([e0,e1,e2,e3,e4,e5,e6,e7,e8, e9, e9b,e10, e11, e12, e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24,e25,e26,e27])
    def add_profile_fit(self):
        #####################################################################
        # profile fitting using galfit geometry
        #####################################################################
        #
        # r-band parameters
        # 
        self.fields_r = ['R24','R25','R26','R_F25','R24V','R25V','R_F50','R_F75','M24','M25','M26', 'F_30R24','F_R24','C30',\
                    'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        self.units_r = [u.arcsec,u.arcsec,u.arcsec,u.arcsec,u.arcsec,\
                   u.arcsec,u.arcsec,u.arcsec,\
                   u.mag, u.mag, u.mag, \
                   u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,'',\
                   u.arcsec,u.erg/u.s/u.cm**2,u.arcsec, u.arcsec,'',u.mag
                   ]
        self.descriptions= ['isophotal radius at 24mag/sqarc AB',\
                            'isophotal radius at 25mag/sqarc AB',\
                            'isophotal radius at 26mag/sqarc AB',\
                            'radius that encloses 25% of total flux',\
                            'isophotal radius at 24mag/sqarc Vega',\
                            'isophotal radius at 24mag/sqarc Vega',\
                            'radius that encloses 50% of total flux',\
                            'radius that encloses 75% of total flux',\
                            'isophotal mag within R24',\
                            'isophotal mag within R25',\
                            'isophotal mag within R26',\
                            'flux within 30% of R24',\
                            'flux within R24',\
                            'C30 = flux w/in 0.3 r24 / flux w/in r24',\
                            'petrosian radius: where sb is 0.2 times mean sb',\
                            'flux enclosed within 2xpetro radius',\
                            'radius enclosing 50% of petrosian flux',\
                            'radius enclosing 90% of petrosian flux',\
                            '90% petro radius / 50% petro radius',\
                            'magnitude of petrosian flux']
        i=0
        for f,unit in zip(self.fields_r,self.units_r):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f,description='galfit '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f, unit=unit,description='galfit '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        #
        # Halpha parameters
        # 
        self.fields_ha = ['R16','R17',\
                  'R_F25','R_F50','R_F75',\
                  'M16','M17', \
                  'F_30R24','F_R24','C30',\
                  'R_F95R24','F_TOT',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG'
                  ]
        self.units_ha = [u.arcsec,u.arcsec,\
                 u.arcsec,u.arcsec, u.arcsec, \
                 u.mag, u.mag, \
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2, '',\
                 u.arcsec,u.erg/u.s/u.cm**2,\
                 u.arcsec,u.erg/u.s/u.cm**2,u.arcsec, u.arcsec,'',u.mag]
        self.descriptions_ha= ['HA isophotal radius at 16erg/s/cm^2',\
                            'HA isophotal radius at 17erg/s/cm^2',\
                            'HA radius that encloses 25% of total flux',\
                            'HA radius that encloses 50% of total flux',\
                            'HA radius that encloses 75% of total flux',\
                            'HA isophotal radius at 16erg/s/cm^s',\
                            'HA isophotal radius at 17erg/s/cm^2',\
                            'HA flux within 30% of R-band R24',\
                            'HA flux within R-band R24',\
                            'HA C30 = flux w/in 0.3 R-band r24 / flux w/in R-band r24',\
                            'HA flux within 30% of R-band R24',\
                            'HA total flux',\
                            'petrosian radius: where sb is 0.2 times mean sb',\
                            'flux enclosed within 2xpetro radius',\
                            'radius enclosing 50% of petrosian flux',\
                            'radius enclosing 90% of petrosian flux',\
                            '90% petro radius / 50% petro radius',\
                            'magnitude of petrosian flux']
        i=0
        for f,unit in zip(self.fields_ha,self.units_ha):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f,description='galfit '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f, unit=unit,description='galfit '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='GAL_'+'H'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1           
        f='GAL_'+'LOG_SFR_HA'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=u.M_sun/u.yr,description='log10 of HA SFR in Msun/yr')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR',unit=u.M_sun/u.yr,description='error in log10 of HA SFR in Msun/yr')
        c3 = Column(np.zeros(self.ngalaxies,'bool'),name=f+'_FLAG')        
        self.table.add_column(c1)
        self.table.add_column(c2)
        self.table.add_column(c3)        
        
        f='GAL_'+'SSFR_IN'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_column(c1)
        self.table.add_column(c2)
        f='GAL_'+'SSFR_OUT'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_column(c1)
        self.table.add_column(c2)
    def add_photutils(self):
        #####################################################################
        # profile fitting using photutils geometry
        #####################################################################
        #
        # r-band parameters
        #
        i=0
        for f,unit in zip(self.fields_r,self.units_r):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='ellipse '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=unit,description='ellipse '+self.descriptions[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        #
        # Halpha parameters
        #
        i=0
        for f,unit in zip(self.fields_ha,self.units_ha):
            if unit == None:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f,description='ellipse '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f+'_ERR')
            else:
                c1 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f, unit=unit,description='ellipse '+self.descriptions_ha[i])
                c2 = Column(np.zeros(self.ngalaxies,'f'),name='H'+f+'_ERR', unit=unit)

            self.table.add_column(c1)
            self.table.add_column(c2)
            i += 1
        f='LOG_SFR_HA'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f, unit=u.M_sun/u.yr,description='log10 of HA SFR in Msun/yr')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR',unit=u.M_sun/u.yr)
        c3 = Column(np.zeros(self.ngalaxies,'bool'),name=f+'_FLAG')
        print('testing: colname = ',f+'_FLAG')
        self.table.add_columns([c1,c2,c3])

        ######################################################################
        ### LAST TWO QUANTITIES, I SWEAR!
        ######################################################################        
        
        f='SSFR_IN'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(r) within 0.3 R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_columns([c1,c2])

        f='SSFR_OUT'
        c1 = Column(np.zeros(self.ngalaxies,'f'),name=f,description='F(HA)/F(R) outside 0.3R24')
        c2 = Column(np.zeros(self.ngalaxies,'f'),name=f+'_ERR')
        self.table.add_columns([c1,c2])


        self.add_flags()
        
        self.table.add_column(Column(np.zeros(self.ngalaxies,dtype='U50'), name='COMMENT'))
        #print(self.table)
        
    def add_statmorph(self):
        #####################################################################
        # statmorph output
        #####################################################################

        # rband area
        e1 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_XCENTROID', unit='pixel',description='xcentroid from ellipse')
        e2 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_YCENTROID', unit='pixel',description='ycentroid from ellipse')
        e3 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RPETRO_CIRC', unit='arcsec',description='rpetro circ')
        e4 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RPETRO_ELLIP', unit='arcsec',description='rpetro ellip')
        e5 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_RHALF_ELLIP', unit='arcsec',description='rhalf ellip')
        e6 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R20', unit='arcsec',description='R20')
        e7 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_R80', unit='arcsec',description='R80')
        e8 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_GINI',description='statmorph gini')
        e8b = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_M20',description='statmorph M20')        
        e9 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_F_GM20',description='statmorph F(G,M20)')
        e10 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_S_GM20',description='statmorph S(G,M20)')
        e11 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_C',description='statmorph concentration')
        e12 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_A',description='statmorph asymmetry')
        e13 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_S',description='statmorph smoothness')
        e14 = Column(np.zeros(self.ngalaxies,'bool'), name='SMORPH_FLAG',description='statmorph flag')                 

        ## Halpha parameters
        h1 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HXCENTROID', unit='pixel',description='xcentroid from ellipse')
        h2 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HYCENTROID', unit='pixel',description='ycentroid from ellipse')
        h3 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRPETRO_CIRC', unit='arcsec',description='rpetro circ')
        h4 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRPETRO_ELLIP', unit='arcsec',description='rpetro ellip')
        h5 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HRHALF_ELLIP', unit='arcsec',description='rhalf ellip')
        h6 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR20', unit='arcsec',description='R20')
        h7 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HR80', unit='arcsec',description='R80')
        h8 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HGINI',description='statmorph gini')
        h8b = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HM20',description='statmorph M20')        
        h9 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HF_GM20',description='statmorph F(G,M20)')
        h10 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HS_GM20',description='statmorph S(G,M20)')
        h11 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HC',description='statmorph concentration')
        h12 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HA',description='statmorph asymmetry')
        h13 = Column(np.zeros(self.ngalaxies,'f'), name='SMORPH_HS',description='statmorph smoothness')
        h14 = Column(np.zeros(self.ngalaxies,'bool'), name='SMORPH_HFLAG',description='statmorph flag')                 
        
        
        self.table.add_columns([e1,e2,e3,e4,e5,e6,e7,e8,e8b, e9, e10, e11, e12, e13,e14,\
                                h1,h2,h3,h4,h5,h6,h7,h8,h8b, h9, h10, h11, h12, h13,h14])

    def add_flags(self):
        '''
        these are common comments that the user will be able to select
        '''
        names = ['CONTSUB_FLAG','MERGER_FLAG','SCATLIGHT_FLAG','ASYMR_FLAG','ASYMHA_FLAG','OVERSTAR_FLAG','OVERGAL_FLAG','PARTIAL_FLAG','EDGEON_FLAG','NUC_HA']
        descriptions =  ['Halpha Emission','Cont Sub Prob','merger/tidal','scattered light','asymmetric R-band', 'asymmetric Ha','foreground star', 'foreground gal','galaxy is edge-on','galaxy is only partially covered by mosaic','nuclear ha emission'] 
        for i,n in enumerate(names):
            #print(n)
            c = Column(np.zeros(self.ngalaxies,'bool'),name=n,description=descriptions[i])
            self.table.add_column(c)
            
    def append_column(self, variable, var_name, var_dtype=None, var_unit=None):
        if (var_dtype != None) & (var_unit != None):
            z = Column(variable, name=var_name, dtype = var_dype, unit  = var_unit)
        elif (var_dtype != None) & (var_unit == None):
            z = Column(variable, name=var_name, dtype = var_dype)
        elif (var_dtype == None) & (var_unit != None):
            z = Column(variable, name=var_name, unit  = var_unit)
        else:
            z = Column(variable, name=var_name)
        self.table.add_column(z)
    def write_fits_table(self):
        if (self.igal is not None) & (not self.auto):
            #print(self.ui.commentLineEdit.text())
            t = str(self.ui.commentLineEdit.text())
            if len(t) > 1:
                self.table['COMMENT'][self.igal] = t
                self.update_gui_table_cell(self.igal, 'COMMENT',t)
        #fits.writeto('halpha-data-'+user+'-'+str_date_today+'.fits',self.table, overwrite=True)
        if self.prefix is not None:
            # this is not working when running gui - need to feed in the r-band image name
            try:
                telescope,dateobs,p = get_params_from_name(self.prefix)
            except UnboundLocalError:
                if self.uat:
                    telescope,dateobs,p = get_params_from_name_uat(self.rcoadd_fname)
                else:
                    telescope,dateobs,p = get_params_from_name(self.rcoadd_fname)
                #print(f"telescope={telescope}, dateobs={dateobs}, p={p}")
            for i in range(len(self.table)):
                self.table['POINTING'][i] = self.prefix
                self.table['TEL'][i] = telescope
                self.table['DATE-OBS'] = dateobs
        self.table.write(self.output_table, format='fits', overwrite=True)

    
class uco_table():
    '''
    table for collecting positions of objects that are not in NSA or AGC catalogs
    '''
    def initialize_uco_arrays(self):
        # columns: id, x, y, ra, dec
        user = os.getenv('USER')
        today = date.today()
        str_date_today = today.strftime('%Y-%b-%d')
        self.uco_output_table = 'halpha-uco-data-'+user+'-'+str_date_today+'.fits'
        if os.path.exists(self.uco_output_table):
            self.uco_table = Table(fits.getdata(self.uco_output_table))
            self.uco_id = self.uco_table['ID'].tolist()
            self.uco_ra = self.uco_table['RA'].tolist()
            self.uco_dec = self.uco_table['DEC'].tolist()
            self.uco_x = self.uco_table['X'].tolist()
            self.uco_y = self.uco_table['Y'].tolist()

        ## if not, create table
        else:
            
            self.uco_id = []
            self.uco_ra = []
            self.uco_dec = []
            self.uco_x = []
            self.uco_y = []
    def create_uco_table(self):
        c1 = Column(np.array(self.uco_id), name='ID',dtype=np.int32, description='ID')
        c2 = Column(np.array(self.uco_ra), name='RA',dtype='f', unit=u.deg)
        c3 = Column(np.array(self.uco_dec), name='DEC',dtype='f', unit=u.deg)
        c4 = Column(np.array(self.uco_x), name='X',dtype='f', unit=u.pixel)
        c5 = Column(np.array(self.uco_y), name='Y',dtype='f', unit=u.pixel)
        self.uco_table = Table([c1,c2,c3,c4,c5])
    def write_uco_table(self):
        self.create_uco_table()

        self.uco_table.write(self.uco_output_table, format='fits',overwrite=True)

class hagui_methods():
    """ class to handle Halpha model, in Model-View-Controller design"""

    def read_rcoadd(self): # model
        #print('reading rband image ',self.rcoadd_fname)
        self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)
        try:
            self.pixelscale = np.abs(float(self.r_header['PIXSCAL1'])) # convert deg/pix to arcsec/pixel                        
        except KeyError:
            try:
                self.pixelscale = np.abs(float(self.r_header['CD1_1']))*3600. # convert deg/pix to arcsec/pixel
            except KeyError:
                self.pixelscale = np.abs(float(self.r_header['PC1_1']))*3600. # Siena pipeline from astronometry.net
        #self.pixelscale = np.abs(float(coadd_header['CD1_1']))*3600          
        #self.psf.psf_image_name = 'MKW8_R.coadd-psf.fits'

        # check for weight image
        # assuminging naming conventions from swarp
        # image = NRGb161_R.coadd.fits
        # weight = NRGb161_R.coadd.weight.fits
        
        weight_image = self.rcoadd_fname.split('.fits')[0]+'.weight.fits'
        if weight_image.find('-shifted') > -1:
            weight_image = self.rcoadd_fname.replace('shifted','weight-shifted')

        if os.path.exists(weight_image):
            self.rweight = weight_image
            self.rweight_flag = True
        else:
            self.rweight_flag = False
        # RAF - updating how WCS is constructed
        header = fits.getheader(self.rcoadd_fname)
        self.coadd_wcs= WCS(header)#OF R IMAGE, SO THAT HA MATCHES WCS OF R, SO THEY'RE THE SAME
    def read_hacoadd(self): # model
        #print(self.hacoadd_fname)
        self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)
        self.hacoadd_cs_fname = self.hacoadd_fname.split('.fits')[0]+'-CS.fits'
        #self.psf.psf_image_name = 'MKW8_R.coadd-psf.fits'

        # check for weight image
        # assuminging naming conventions from swarp
        # image = NRGb161_R.coadd.fits
        # weight = NRGb161_R.coadd.weight.fits
        
        weight_image = self.hacoadd_fname.split('.fits')[0]+'.weight.fits'
        if os.path.exists(weight_image):
            self.haweight = weight_image
            self.haweight_flag = True
        else:
            self.haweight_flag = False
        
    def load_rcoadd(self): # model
        self.coadd.load_file(self.rcoadd_fname)
        self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)
        try:
            self.pixelscale = abs(float(self.r_header['PIXSCAL1'])) # in deg per pixel
        except KeyError:
            try:
                self.pixelscale = abs(float(self.r_header['CD1_1']))*3600. # in deg per pixel
            except KeyError:
                self.pixelscale = abs(float(self.r_header['PC1_1']))*3600. # in deg per pixel
        #self.psf.psf_image_name = 'MKW8_R.coadd-psf.fits'

        # check for weight image
        # assuminging naming conventions from swarp
        # image = NRGb161_R.coadd.fits
        # weight = NRGb161_R.coadd.weight.fits
        
        weight_image = self.rcoadd_fname.split('.fits')[0]+'.weight.fits'
        if os.path.exists(weight_image):
            self.rweight = weight_image
            self.rweight_flag = True
        else:
            self.rweight_flag = False

        self.coadd_wcs= WCS(self.rcoadd_fname)#OF R IMAGE, SO THAT HA MATCHES WCS OF R, SO THEY'RE THE SAME
        self.check_previous_r_psf()
        #not sure what the following line is suppose to do
        # ephot['SMMA_MOMENT'][~noradius_flag]
    def load_hacoadd(self): # model
        self.coadd.load_file(self.hacoadd_fname)
        self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)

        #self.psf.psf_image_name = 'MKW8_R.coadd-psf.fits'
        weight_image = self.hacoadd_fname.split('.fits')[0]+'.weight.fits'
        if os.path.exists(weight_image):
            self.haweight = weight_image
            self.haweight_flag = True
        else:
            self.haweight_flag = False
        self.check_previous_ha_psf()
    def get_zcut(self): # model
        self.zmax=(((lmax[self.hafilter])/6563.)-1)
        self.zmin=(((lmin[self.hafilter])/6563.)-1)
    def check_previous_ha_psf(self): # MVC - model
        """
        check for data from previous runs, including
        - psf file
        - data table with results for all/subset of galaxies
        - filter ratio
        """

        # check for halpha psf
        basename = os.path.basename(self.hacoadd_fname)
        psf_image_name = basename.split('.fits')[0]+'-psf.fits'
        if os.path.exists(psf_image_name):
            self.psf_haimage_name = psf_image_name
        else:
            self.psf_haimage_name = None
    def check_previous_r_psf(self): # MVC - model
        """
        check for data from previous runs, including
        - psf file
        - data table with results for all/subset of galaxies
        - filter ratio
        """
        # check for rband psf
        basename = os.path.basename(self.rcoadd_fname)
        psf_image_name = basename.split('.fits')[0]+'-psf.fits'
        if os.path.exists(psf_image_name):
            self.psf_image_name = psf_image_name
        else:
            self.psf_image_name = None

    def find_galaxies(self): # MVC - model
        print('getting galaxies in FOV')
        flag = self.get_gal_list()
        if flag == False:
            #print("WARNING: could not find galaxies!")
            return


        ##
        # rewrite to get the size from JM's ephot, +/- 3*SMA_SB24
        #
        # 3 is too small for most galaxies - forgot to convert radius to diameter
        # ran with 2*2.5, and this was tight for a lot of galaxies.
        # increasing size by 50% (2.5*1.5 = 15/4 = 3.75 so rounding to 4)
        ##
        scale = 2*2.5
        # ROSE - DON'T CHANGE THIS ANYMORE!!!!!!!
        scale = 2*3 
        #scale = 2*4        

        # NOW USING SB25, so going to decrease scale factor to 2
        scale = 2*1.9

        if self.uat:
            scale = 3.5



        
        # this is the total length in one dimension
        # the scale factor includes an extra factor of 2 to convert radius to diameter
        self.cutout_sizes = self.radius_arcsec*scale/self.pixelscale # in pixels


        
        
        #start_time = time.perf_counter()        
        #print('getting object sizes from segmentation image')
        #if 'HDI' is in self.rcoadd_fname:
        #    scale = 3.5
        #else:
        #    scale = 1.75
        #self.cutout_sizes = getobjectsize(self.rcoadd_fname,np.array(self.gximage,'i'),np.array(self.gyimage,'i'),scale=scale)
        #print('...done with segmentation image in {:2f} sec'.format(time.perf_counter()-start_time))
        #print('...sorry for the wait')
        #print('\t cutout sizes = ',self.cutout_sizes)
              
        # convert to arcsec
        self.cutout_sizes_arcsec = self.cutout_sizes*self.pixelscale*u.arcsec


        #######################################################
        # OK, so here is the thing...
        # we really should decrease the cutout size as the
        # image gets bigger.  We don't need arcminutes around
        # big galaxies to get a proper sky sample.
        #######################################################    

        smallflag = self.cutout_sizes_arcsec < self.global_min_cutout_size# make the smallest cutout size equal to 60

        self.cutout_sizes_arcsec[smallflag] = self.global_min_cutout_size
        self.cutout_sizes[smallflag] = self.global_min_cutout_size.value/self.pixelscale
        

        # this step is described in the notebook cutout-sizing in havirgo/notebooks
        bigflag = self.radius_arcsec > 240.
        self.cutout_sizes_arcsec[bigflag] = self.cutout_sizes_arcsec[bigflag] - (2.5/60)*(self.radius_arcsec[bigflag] - 240.)*u.arcsec
        self.cutout_sizes[bigflag] = self.cutout_sizes_arcsec[bigflag].value/self.pixelscale
        
        #print('\t cutout sizes arcsec = ',self.cutout_sizes_arcsec)
        # check to see if any are outside the allowed range


        
        # ##
        # # SKIPPING MIN/MAX RESET - would rather correct values in source catalogs
        # ##
        # for i,s in enumerate(self.cutout_sizes_arcsec):
        #     # the cutouts for the biggest galaxies are getting cut short,
        #     # so commenting this part out
            
        #     #if s > self.global_max_cutout_size:
        #     #    self.cutout_sizes_arcsec[i] = self.global_max_cutout_size
        #     #    self.cutout_sizes[i] = self.global_max_cutout_size.value/self.pixelscale
            
        #     if s < self.global_min_cutout_size:
        #         self.cutout_sizes_arcsec[i] = self.global_min_cutout_size
        #         self.cutout_sizes[i] = self.global_min_cutout_size.value/self.pixelscale
                
        # #print('after comparing with max/min size limits:')
        # #print('\t cutout sizes arcsec = ',self.cutout_sizes_arcsec)        
        # # set up the output table that will store results from various fits

        #print("Value of self.uat = ",self.uat)
        self.initialize_results_table(prefix=self.prefix,virgo=self.virgo,uat=self.uat,nogui=self.auto)

        # testing
        
        # TODO - this needs to be moved into the view class
        if not self.auto:
            # plot location of galaxies in the coadd image
            self.mark_galaxies()

            # populate a button that contains list
            # of galaxies in the field of view,
            # user can select from list to set the active galaxy
            for name in self.galid:
                self.ui.wgalid.addItem(str(name))
            print(len(self.galid),' galaxies in FOV')
            self.ui.wgalid.activated.connect(self.select_galaxy)

        # get transmission correction for each galaxy
        # add prefix to figure that is created by filter_trace.get_trans_correction        
        figfile = f"{self.prefix}-galaxies_in_filter.png"        
        self.filter_correction= self.filter_trace.get_trans_correction(self.table['REDSHIFT'],outfile=figfile)
        #print(self.filter_correction)
        self.table['FILT_COR'] = self.filter_correction

        #TODO - add this to the NB image header!!!

        
    def get_gal_list(self): # MVC - model
        
        """
        # get list of NSA galaxies on image viewer
        #
        # for reference:
        # self.r, header_r = fits.getdata(self.rcoadd_fname,header=True)
        # self.ha, header_ha = fits.getdata(self.hacoadd_fname, header=True)
        #
        """
        
        n2,n1 = self.r.data.shape #should be same for Ha too, maybe? IDK

        #keepflag = self.defcat.galaxies_in_fov(self.coadd_wcs, nrow=n2,ncol=n1,zmin=self.zmin,zmax=self.zmax,virgoflag=self.virgo)
        #keepflag = self.defcat.galaxies_in_fov(self.coadd_wcs, nrow=n2,ncol=n1,zmin=self.zmin,zmax=self.zmax)
        #try:
        #    keepflag = self.defcat.galaxies_in_fov(self.coadd_wcs, nrow=n2,ncol=n1,image_name=self.rcoadd_fname,zmin=self.zmin,zmax=self.zmax)
        #except AttributeError:
        try:
            self.defcat
        except AttributeError:
            print(f'Warning: no catalog defined. self.defcat probably not defined')
            print('Make sure you set the path to the NSA/Virgo parent catalog\n')
            return False


        #keepflag = self.defcat.galaxies_in_fov(self.coadd_wcs, nrow=n2,ncol=n1,image_name=self.rcoadd_fname,zmin=self.zmin,zmax=self.zmax)
        #print("finished running galaxies_in_fov")
        #print()
        try:
            print(f'min and max redshift = {self.zmin:.3f}, {self.zmax:.3f}')
            keepflag = self.defcat.galaxies_in_fov(self.coadd_wcs, nrow=n2,ncol=n1,image_name=self.rcoadd_fname,zmin=self.zmin,zmax=self.zmax)
            try:
                print('number of galaxies in FOV = ',np.sum(keepflag))
            except TypeError:
                print("WARNING: trouble counting galaxies in FOV - make sure you have the filter correct.")
            
        except AttributeError:
            print('problem finding galaxies.')
            print('make sure you selected a filter!')
            return False
        
        #print('keepflag = ',keepflag)
            
        # check weight image to make sure the galaxy actually has data
        # reject galaxies who have zero in the weight image
        px,py = self.coadd_wcs.wcs_world2pix(self.defcat.cat['RA'],self.defcat.cat['DEC'],0)
        px = np.nan_to_num(px, nan=-99)
        py = np.nan_to_num(py, nan=-99)

        self.xpixel = px
        self.ypixel = py
        
        if self.rweight_flag and self.haweight_flag:
            rweight = fits.getdata(self.rweight)
            haweight = fits.getdata(self.haweight)
            
            # multiply weights
            # result will equal zero if exposure is zero in either image
            weight = rweight * haweight
            
            # check location of pixels to make sure weight is not zero
            # this will have the length = # of galaxies that have keepflag True
            offimage = (weight[np.array(py[keepflag],'i'),np.array(px[keepflag],'i')] == 0)
            
            # store another array that has the indices in original keepflag array
            # where keepflag = True
            # need to take [0] element because np.where is weird
            keepindex = np.where(keepflag)[0]
            # change value of keepflag for galaxies that are off the image
            keepflag[keepindex[offimage]] = np.zeros(len(keepindex[offimage]),'bool')
            #print(offimage)
        # cut down NSA catalog to keep information only for galaxies within FOV
        #print('number of galaxies in FOV = ',sum(keepflag))
        if (np.sum(keepflag) == 0) | (keepflag is None):
            print('WARNING: no NSA galaxies in FOV')
            print('\t make sure you have selected the right filter!')
            return
        if self.prefix is None:
            print('you need to set the prefix')
            print('try again!')
            return
        self.defcat.cull_catalog(keepflag,self.prefix)
    
        if self.virgo: # why don't we do this for uat?
            self.nsa.cull_catalog(keepflag,self.prefix)
        if self.virgo | self.uat:
            self.radius_arcsec = self.radius_arcsec[keepflag]
            self.BA = self.BA[keepflag]
            self.PA = self.PA[keepflag]
        
        #self.gra=self.nsa.cat.RA
        #self.gdec=self.nsa.cat.DEC
        #self.gradius=self.nsa.cat.SERSIC_TH50
        #self.galid=self.nsa.cat.NSAID
        #self.gredshift = self.nsa.cat.Z
        #self.gzdist= self.nsa.cat.ZDIST
        self.gximage,self.gyimage =self.coadd_wcs.wcs_world2pix(self.defcat.cat['RA'],self.defcat.cat['DEC'],0)
        # set up a boolean array to track whether Halpha emission is present

        if self.agcflag:
            #try:
            #    px,py = self.coadd_wcs.wcs_world2pix(self.agc.cat.RA,self.agc.cat.DEC,0)
            #except AttributeError:
            #    px,py = self.coadd_wcs.wcs_world2pix(self.agc.cat.radeg,self.agc.cat.decdeg,0)

            keepagc = self.agc.galaxies_in_fov(self.coadd_wcs, nrow=n2,ncol=n1, agcflag=True,image_name=self.rcoadd_fname,zmin=self.zmin,zmax=self.zmax)

            print('number of AGC galaxies in FOV = ',sum(keepagc))

            #if self.uat: # testing to see if this fixes the large radius problem
            print(f"length of radius_arcsec = {len(self.radius_arcsec)}, len(keepagc)={len(keepagc)}")
            #    self.radius_arcsec = self.radius_arcsec[keepagc]
            #    self.BA = self.BA[keepagc]
            #    self.PA = self.PA[keepagc]
            
            self.agc.cull_catalog(keepagc, self.prefix)
            if sum(keepagc) == 0:
                self.agcflag = False
            else:
                try:
                    self.agcximage,self.agcyimage =self.coadd_wcs.wcs_world2pix(self.agc.cat['RA'],self.agc.cat['DEC'],0)
                except AttributeError:
                    self.agcximage,self.agcyimage =self.coadd_wcs.wcs_world2pix(self.agc.cat['radeg'],self.agc.cat['decdeg'],0)        
            #print(self.agcximage)
            #print(self.agcyimage)


        return True
    def initialize_output_arrays(self): # MVC - model
        ngal = len(self.galid)

        # galfit output
        self.gal_xc = np.zeros((ngal,2),'f')
        self.gal_xc = np.zeros((ngal,2),'f')
        self.gal_mag = np.zeros((ngal,2),'f')
        self.gal_n = np.zeros((ngal,2),'f')
        self.gal_re = np.zeros((ngal,2),'f')
        self.gal_PA = np.zeros((ngal,2),'f')
        self.gal_BA = np.zeros((ngal,2),'f')
        self.gal_sky = np.zeros((ngal,2),'f')
        
    def link_files(self): # MVC- model
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory
        sextractor_files=['default.sex.HDI','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            if not os.path.exists(file):
                os.system('ln -s '+self.sepath+'/'+file+' .')
    def clean_links(self): # MVC - model
        # clean up symbolic links to sextractor files
        # sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        sextractor_files=['default.sex.HDI','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('unlink '+file)

        # remove catalog
    def get_filter_ratio(self): # MVC - model
        #
        # get ratio of Halpha to Rband filters
        # 
        # cannabalizing HalphaImaging/uat_find_filter_ratio.py
        #

        ##
        # Look for ZP info in headers
        ##

        try:
            header = fits.getheader(self.rcoadd_fname)
            #rZP = header['PHOTZP']
            #header = fits.getheader(self.hacoadd_fname)
            #hZP = header['PHOTZP']
            #dm = hZP-rZP
            #ave = 10**(dm/2.5) # f2/f1

            ave = header['FRATIOZP']
            run_se = False
            print("Getting Filter Ratio from r-band image header")
        except KeyError:
            print("Getting Filter Ratio: FRATIOZP not found in r-band image; will run SE")
            self.link_files()
            current_dir = os.getcwd()
            image_dir = os.path.dirname(self.rcoadd_fname)
            #os.chdir(image_dir)

            ZP1, zp1flag, ZP2, zp2flag = runse.run_sextractor(self.rcoadd_fname, self.hacoadd_fname)
            if zp1flag and zp2flag:
                # print("got ZP ratio")
                zpargs = (ZP1, ZP2)
            else:
                zpargs = None
            t = runse.make_plot(self.rcoadd_fname, self.hacoadd_fname, return_flag = True, plotdir = current_dir,zps = zpargs)
            if len(t) == 2:
                ave, std = t
                fzpratio = None
            elif len(t) == 3:
                ave, std, fzpratio = t
            print(f"filter ratio = {ave:.4f} +/- {std:.4f}")
            run_se = True
            #plt.show()
            #os.chdir(current_dir)
        self.filter_ratio = ave
        self.reset_ratio = ave
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio

        self.subtract_images()
        if not self.auto:
            self.ui.filterRatioLineEdit.setText(str(self.filter_ratio))
        #if not self.auto:
        #    self.setup_ratio_slider()
        if run_se:
            self.clean_links()

            # the following lines remove the SE files, but we shouldn't do this
            # it takes a while to run SE.  Should keep the catalogs
            images = [self.rcoadd_fname, self.hacoadd_fname]
            for im in images:
                catfile = os.path.basename(im).split('.fits')[0]+'.cat'
                os.remove(catfile)

            # add filter ratio to image headers
            # add ratio to r-band image headers
            r,header = fits.getdata(self.rcoadd_fname,header=True)
            header.set('FLTRATIO',ave)
            header.set('FLTR_ERR',std)
            header.set('HAIMAGE',os.path.basename(self.hacoadd_fname))
            if fzpratio is not None:
                header.set('FRATIOZP',fzpratio)
            
            fits.writeto(self.rcoadd_fname,r,header=header,overwrite=True)
                
    def subtract_images(self,overwrite=False): # MVC - model
        # use the ZP subtracted images
        self.hacoadd_cs_fname = self.hacoadd_fname.split('.fits')[0]+'-CS-ZP.fits'
        
        if os.path.exists(self.hacoadd_cs_fname) and not overwrite:
            # don't need to subtract images
            self.halpha_cs = fits.getdata(self.hacoadd_cs_fname)
        else:
            try:
                self.halpha_cs = self.ha - self.filter_ratio*self.r
                fits.writeto(self.hacoadd_cs_fname,self.halpha_cs,header=self.ha_header,overwrite=True)           
            except AttributeError:

                print('WARNING: no filter ratio')

        # TODO - needs to be moved to view and then called from controller
        if not self.auto:
            # display continuum subtracted Halpha image in the large frame
            self.coadd.fitsimage.set_autocut_params('zscale')
            self.coadd.fitsimage.set_data(self.halpha_cs)
        
    def get_galaxy_cutout(self): 

        # scale cutout size according to NSA Re
        #size = cutout_scale*self.gradius[self.igal]
        # updating to use the size from segmentation image from photutils
        size_arcsec = self.cutout_sizes_arcsec[self.igal]
        # set the min size to 100x100 pixels (43"x43")
        #print('size, global min, global max = ',size_arcsec,self.global_min_cutout_size,self.global_max_cutout_size)
        size = max(self.global_min_cutout_size,size_arcsec)

        ##
        # cutouts for biggest galaxies are too small, so commenting this out
        ##
        #if size > self.global_max_cutout_size:
        #    size = self.global_max_cutout_size

        
        #print('new cutout size = ',size, self.igal, self.gradius[self.igal])

        self.bad_galaxy = False # this will get set to True if the noise in the sky is zero, which happens if gal is in top right corner of INT
        
        self.reset_size = size.value
        self.cutout_size_arcsec = size
        self.cutout_size = size.value/self.pixelscale
        self.mincutout_size = 0.2*self.cutout_size
        self.maxcutout_size = 3.*self.cutout_size

        

        # updating on march 5, 2021 to make cutout directory name unique
        # to prevent accidentally writing over data
        # add date and pointing

        if self.verbose:
            print("\n in get_cutouts, accessing halpha header\n")
        header = self.ha_header

        # date is stored in epoch for HDI data

        if self.verbose:
            print("\n in get_cutouts, getting date\n")

        try:
            # store time 
            t = header['EPOCH']
            # convert to year, month,day
            t = Time(t,format='decimalyear')
        except KeyError:
            # try format expected for INT data
            try:
                t = Time(header['DATE-OBS'],format='isot')
            except KeyError:
                # putting in a place holder b/c older (2017)
                # HDI coadds don't have epoch in the header.
                t = Time('2017-05-20T00:00:00.0',format='isot')
        dateobs = t.iso.split()[0].replace('-','')

        # get instrument
        if self.verbose:
            print("\n in get_cutouts, getting instrument\n")
        
        try:
            instrument = self.ha_header['INSTRUME']
            if instrument.find('hdi') > -1:
                instrument='HDI'
            elif instrument.find('Mosaic') > -1:
                instrument='MOSAIC'
        except KeyError:
            instrument='INT'

        if self.verbose:
            print("\n in get_cutouts, getting object\n")
            
        # read in object
        o = header['OBJECT']
        if instrument == '90prime':
            # this should be the VFID of the primary target
            pointing = o.split('_')[0]
            instrument = 'BOK'
        else:
            try:

                if (o.find('197') > -1) | (o.find('227') > -1): #INT, may data
                    pnumber = o.split('-')[1]
                else:
                    # try to identify format
                    #print(o,len(o))
                    #print(o)
                    if len(o.split()) > 1:
                        split_string=' '
                        pnumber = int(o.split(split_string)[1])
                        #print('object names contain ',split_string)        
                    elif len(o.split('-')) > 1:
                        split_string='-'
                        pnumber = int(o.split(split_string)[1])
                        #print('object names contain ',split_string)
                    elif len(o.split('_')) > 1:
                        split_string='_'
                        #print('object names contain ',split_string)                     
                        pnumber = int(o.split(split_string)[1])
                        #print('testing')
                if (o.find('lm') > -1)| (o.find('LM') > -1):
                    # low-mass pointing
                    prefix = "lmp"
                else:
                    prefix = "p"
            
                pointing = "{}{:03d}".format(prefix,pnumber)

            except ValueError:
                pointing=self.rcoadd_fname.replace('R.fits','')
            except UnboundLocalError:
                pointing = o


            if self.rcoadd_fname.find('nNGC5846') > -1:
                pointing = os.path.basename(self.rcoadd_fname.replace('R.fits',''))
        # first pass of mask
        # radial profiles
        # save latest of mask
        if self.virgo:
            if self.verbose:
                print("\n in get_cutouts, getting NEDname\n")
            
            nedname = self.NEDname[self.igal].replace(" ","")
            nedname = nedname.replace("[","")
            nedname = nedname.replace("]","")
            nedname = nedname.replace("/","")

            instrument,dateobs,pointing = get_params_from_name(self.prefix)
            
            
            cprefix = "{}-{}-{}-{}-{}".format(self.galid[self.igal],nedname,instrument,dateobs,pointing)
            
        else:
            cprefix = "{}-{}-{}-{}".format(self.galid[self.igal],instrument,dateobs,pointing)
            
    

        
        ## create the output directory to store the cutouts and figures
        if self.verbose:
            print("\n in get_cutouts, creating output directory\n")
        
        self.cutoutdir = 'cutouts/'+cprefix+'/'
        if not os.path.exists('cutouts'):
            os.mkdir('cutouts')
        if not os.path.exists(self.cutoutdir):
            os.mkdir(self.cutoutdir)

        if self.uat:
            self.cutout_name_r = self.cutoutdir+cprefix+'-R.fits'
            self.cutout_name_ha =self.cutoutdir+cprefix+'-CS.fits'
            self.cutout_name_hawc= self.cutoutdir+cprefix+'-Ha.fits'
        else:
            self.cutout_name_r = self.cutoutdir+cprefix+'-R.fits'
            self.cutout_name_ha =self.cutoutdir+cprefix+'-CS.fits'
            self.cutout_name_hawc= self.cutoutdir+cprefix+'-Ha.fits'

        position = SkyCoord(ra=self.ra[self.igal],dec=self.dec[self.igal],unit='deg')

        if self.verbose:
            print("\n in get_cutouts, calling Cutout2D\n")
        
        try:
            ###########################################################
            # TODONE - this part should be moved to get_galaxy_cutout
            #print(self.cutout_size_arcsec)
            
            #get cutout, and require entire image to be on parent image
            self.cutoutR = Cutout2D(self.r.data, position, self.cutout_size_arcsec, wcs=self.coadd_wcs, mode='trim') 
            ############################################################

        except nddata.utils.PartialOverlapError:# PartialOverlapError:
            print('galaxy is only partially covered by mosaic - skipping ',self.galid[self.igal])
            return
        except nddata.utils.NoOverlapError:# PartialOverlapError:
            print('galaxy is not covered by mosaic - skipping ',self.galid[self.igal])
            return

        t = self.cutout_name_r.split('.fit')
        self.mask_image_name=t[0]+'-mask.fits'

        if self.verbose:
            print("\n in get_cutouts, writing cutouts\n")

        self.write_cutouts()

        #self.display_cutouts() - moved this to calling function in controller


        
    def write_cutouts(self): # MVC - model
        #print(ymin,ymax,xmin,xmax)
        #print('in write_cutouts')
        try:
            self.table['FILTER_RATIO'][self.igal] = self.filter_ratio
        except AttributeError:
            print('Warning - could not set filter ratio in table')
            print('I think something is wrong')


        w = WCS(self.rcoadd_fname)
        try:
            ((ymin,ymax),(xmin,xmax)) = self.cutoutR.bbox_original
        except AttributeError:
            print('make sure you have selected a galaxy and saved the cutout')
            return

        if self.verbose:
            print("\n in write_cutouts, calling Cutout2D\n")

        # add other image parameters to table        
        self.table['RZP'][self.igal] = self.r_header['PHOTZP']
        
        self.table['HZP'][self.igal] = self.ha_header['PHOTZP']        
        self.table['PIXSCALE'][self.igal] = self.pixelscale


        ###################################################
        # saving R-band Cutout as fits image
        ###################################################  
        newfile = fits.PrimaryHDU()
        # specify cutout from the region returned with Cutout2D,
        # which is saved with self.cutoutR
        if self.verbose:
            print("\nCutout region: ",ymin,ymax,xmin,xmax)
        newfile.data = self.r[ymin:ymax,xmin:xmax]
        # add sky subtraction here

        if self.verbose:
            print("\n in write_cutouts, writing header\n")

        newfile.header = self.r_header
        newfile.header.update(w[ymin:ymax,xmin:xmax].to_header())
        newfile.header.set('REDSHIFT',float('{:.6f}'.format(self.gredshift[self.igal])))
        newfile.header.set('ZDIST',float('{:.6f}'.format(self.gzdist[self.igal])))
        newfile.header.set('ID',str(self.galid[self.igal]))
        newfile.header.set('SERSTH50',float('{:.2f}'.format(self.gradius[self.igal])))
        #print('trying to add exptime now')
        # set the exposure time to 1 sec
        # for some reason, in the coadded images produced by swarp, the gain has been corrected
        # to account for image units in ADU/s, but the exptime was not adjusted.
        # this will impact galfit magnitudes if we don't correct it here.
        # alternatively, we could fix it right after running swarp
        newfile.header['EXPTIME']=1.0

        # subtract the sky, but we need the mask - no we don't!
        skysub_data,rmed,rstd = imutils.subtract_median_sky(newfile.data,getstd=True)

        if rmed is not np.nan:
            newfile.data -= rmed
        if self.verbose:
            print("\n in write_cutouts, add sky stats to header\n")

        try:
            newfile.header.set('SKYMED',rmed)
        except ValueError:
            print("error writing SKYMED ",rmed)
            print("setting header value to zero")            
            newfile.header.set('SKYMED',0)
        try:
            newfile.header.set('SKYSTD',rstd)
        except ValueError:
            print("error writing SKYSTD ",rstd)
            print("setting header value to zero")
            newfile.header.set('SKYSTD',0)
            

        # need to add a case to handle when sky std is zero.
        # This happens when object is in wacky region of INT.
        if rstd == 0:
            # return to caller and move to the next galaxy
            self.bad_galaxy = True
            return
        
        # add PSF image

        # add coadded image

        #print('saving r-band cutout')

        if self.verbose:
            print("\n in write_cutouts, saving cutout\n")

        print(f"subtracted {rmed:.3f} ADU from r-band cutout",rstd)            
        fits.writeto(self.cutout_name_r, newfile.data, header = newfile.header, overwrite=True)

        ###################################################
        # saving Ha CS Cutout as fits image
        ###################################################  
        try:
            t = self.halpha_cs
            newfile1 = fits.PrimaryHDU()
            newfile1.data = self.halpha_cs[ymin:ymax,xmin:xmax]
            newfile1.header = self.ha_header
            newfile1.header.update(w[ymin:ymax,xmin:xmax].to_header())
            newfile1.header.set('REDSHIFT',float('{:.6f}'.format(self.gredshift[self.igal])))
            newfile1.header.set('ZDIST',float('{:.6f}'.format(self.gzdist[self.igal])))
            newfile1.header.set('ID',str(self.galid[self.igal]))

            newfile1.header.set('SERSIC_TH50',float('{:.2f}'.format(self.gradius[self.igal])))
            newfile1.header['EXPTIME']=1.0

            # subtract sky from CS Halpha image
            skysub_hdata,hmed,hstd = imutils.subtract_median_sky(newfile1.data,getstd=True)

            if hmed is not np.nan:
                newfile1.data -= hmed
            

            try:
                newfile1.header.set('SKYMED',hmed)
            except ValueError:
                print("error writing SKYMED ",hmed)
                print("setting header value to zero")            
                newfile1.header.set('SKYMED',0)
                try:
                    newfile1.header.set('SKYSTD',hstd)
                except ValueError:
                    print("error writing SKYSTD ",hstd)
                    print("setting header value to zero")
                    newfile1.header.set('SKYSTD',0)
        

            # print sky stats
            print(f"subtracted {hmed:.3f} ADU from halpha cutout",hstd)        

        
            fits.writeto(self.cutout_name_ha, newfile1.data, header = newfile1.header, overwrite=True)
            #print('saving halpha cutout')

        except AttributeError:
            print("WARNING: No continuum subtracted image.  Be sure to run get_filter_ratio")

        ###################################################
        # saving Ha with continuum Cutout as fits image
        ###################################################        
        newfile2 = fits.PrimaryHDU()
        newfile2.data = self.ha[ymin:ymax,xmin:xmax]
        newfile2.header = self.ha_header
        newfile2.header.update(w[ymin:ymax,xmin:xmax].to_header())
        newfile2.header.set('REDSHIFT',float('{:.6f}'.format(self.gredshift[self.igal])))
        newfile2.header.set('ZDIST',float('{:.6f}'.format(self.gzdist[self.igal])))
        newfile2.header.set('ID',str(self.galid[self.igal]))

        newfile2.header.set('SERSIC_TH50',float('{:.2f}'.format(self.gradius[self.igal])))
        newfile2.header['EXPTIME']=1.0

        # subtract sky from Halpha image
        skysub_hdata,hmed,hstd = imutils.subtract_median_sky(newfile2.data,getstd=True)

        newfile2.data -= hmed
        newfile2.header.set('SKYMED',hmed)
        newfile2.header.set('SKYSTD',hstd)
        
        fits.writeto(self.cutout_name_hawc, newfile2.data, header = newfile2.header, overwrite=True)

        # save bounding box info
        ((ymin,ymax),(xmin,xmax)) = self.cutoutR.bbox_original
        bbox = '[{:d}:{:d},{:d}:{:d}]'.format(int(xmin),int(xmax),int(ymin),int(ymax))
        self.table['BBOX'][self.igal] = bbox

        # this should be the end of the write function
        
    def build_psf(self): # MVC - model
        # check to see if R-band PSF images exist
        coadd_header = fits.getheader(self.rcoadd_fname)


        basename = os.path.basename(self.rcoadd_fname)
        psf_image_name = basename.split('.fits')[0]+'-psf.fits'
        if psf_image_name.find('-shifted') > -1:
            psf_image_name = psf_image_name.replace('-shifted','')
        basename = os.path.basename(self.hacoadd_fname)
        psf_image_name_ha = basename.split('.fits')[0]+'-psf.fits'
        psf_image_name = os.path.join(self.psfdirectory,psf_image_name)
        psf_image_name_ha = os.path.join(self.psfdirectory,psf_image_name_ha).replace('-CS','')
        if self.verbose:
            print('\nPSF NAME = ',psf_image_name,'\n')
        if os.path.exists(psf_image_name):
            # get fwhm from the image header
            # and oversampling
            print("LOADING EXISTING PSF IMAGE")
            header = fits.getheader(psf_image_name)
            self.psf_data = fits.getdata(psf_image_name)

            # add attributes
            self.psf = psfimage()
            
            self.psf.fwhm = header['FWHM'] # in pixels
            self.psf.fwhm_arcsec = self.psf.fwhm*self.pixelscale
            self.oversampling = float(header['OVERSAMP'])
            # if psf is in another directory, create a link to the current directory
            # this will avoid having a long filename b/c galfit does not handle long filenames

            ##
            # GALFIT does not handle the long image name - it craps out
            #
            # that must be why I copied the psf image to r-psf.fits
            # need to figure out another way so that I can run them in parallel
            ##
            #self.psf_image_name = psf_image_name

            ##
            # changing to use the psf without copying to current directory
            ##
            outname = os.path.basename(psf_image_name)
            command = f'ln -s {psf_image_name} {outname}'
            #print('running: ',command)
            os.system(command)
            #self.psf_image_name = 'r-psf.fits'
            self.psf_image_name = outname


            
        else:
            if self.verbose:
                print('oversampling = ',self.oversampling)
            print('PSF RESULTS FOR R-BAND COADDED IMAGE')
            self.psf = psf_parent_image(image=self.rcoadd_fname, size=21, nstars=100, oversampling=self.oversampling)
            self.psf.run_all()
            self.psf_image_name = self.psf.psf_image_name
            
        print('\nPSF NAME = ',psf_image_name_ha,'\n')
        if os.path.exists(psf_image_name_ha):
            # get fwhm from the image header
            # and oversampling
            print("LOADING EXISTING HALPHA PSF IMAGE")            
            header = fits.getheader(psf_image_name_ha)
            self.hapsf = psfimage()
            self.hapsf.fwhm = header['FWHM'] # in pixels
            self.hapsf.fwhm_arcsec = self.hapsf.fwhm*self.pixelscale

            # if psf is in another directory, create a link to the current directory
            # this will avoid having a long filename b/c galfit does not handle long filenames

            #self.psf_haimage_name = psf_image_name
            #command = 'cp {} ha-psf.fits'.format(psf_image_name)
            #print('running: ',command)
            #os.system(command)
            #self.psf_haimage_name = 'ha-psf.fits'

            outname = os.path.basename(psf_image_name_ha)
            command = f'ln -s {psf_image_name_ha} {outname}'
            #print('running: ',command)
            os.system(command)
            self.psf_haimage_name = outname
            
        else:
            print('PSF RESULTS FOR HA COADDED IMAGE')
            self.hapsf = psf_parent_image(image=self.hacoadd_fname, size=21, nstars=100, oversampling=self.oversampling)
            self.hapsf.run_all()
            self.psf_haimage_name = self.hapsf.psf_image_name

    def add_psf_to_table(self): # MVC - model
        fields = ['R_FWHM','H_FWHM']
        values = [self.psf.fwhm_arcsec,self.hapsf.fwhm_arcsec]
        for i,f in enumerate(fields):
            for j in range(len(self.table)):
                self.table[f][j]=values[i]
        pass
    def run_galfit(self, ncomp=1, asym=0, ha=0): # MVC - model?
        if self.psf_image_name is None:
            print('WARNING: psf could not be found')
            print('Please run build_psf')
            return
        if self.psf is not None:
            self.add_psf_to_table()
        self.ncomp = ncomp
        self.asym=asym
        self.galha = ha
        print('running galfit with ',ncomp,' components')
        #self.gwindow = QtWidgets.QWidget()
        '''
        if self.testing:
            self.ncomp = ncomp
            self.galfit = galfitwindow(self.gwindow, self.logger, image = 'MKW8-18037-R.fits', mask_image = 'MKW8-18037-R-mask.fits', psf='MKW8_R.coadd-psf.fits', psf_oversampling=2, ncomp=ncomp)
        else:
            self.galfit = galfitwindow(self.gwindow, self.logger, image = self.cutout_name_r, mask_image = self.mask_image_name, psf=self.psf.psf_image_name, psf_oversampling = self.oversampling, ncomp=ncomp)
        self.galfit.model_saved.connect(self.galfit_save)        
        self.galfit.setupUi(self.gwindow)

        self.gwindow.show()
        '''
        try:
            if ha:
                self.galimage = self.cutout_name_ha
                # setup psfimage
                psf = self.psf_haimage_name
                psf_oversampling = self.oversampling
            else:
                self.galimage = self.cutout_name_r
                # setup psfimage
                psf = self.psf_image_name
                psf_oversampling = self.oversampling
        except AttributeError:
            print('make sure you selected a galaxy')
            return
        try:
            if not self.auto:
                self.gwindow = QtWidgets.QWidget()
            #self.gwindow.aboutToQuit.connect(self.galfit_closed)

            
            if (ncomp == 1) & (asym == 0):
                #self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=self.nsa.rmag[self.igal], BA = self.nsa.cat.SERSIC_BA[self.igal], PA=self.nsa.cat.SERSIC_PHI[self.igal],nsersic=self.nsa.cat.SERSIC_N[self.igal], convolution_size=80)
                print('GALFIT psf image = ',psf)
                if self.auto:
                    self.galfit = galfitwindow(None, None, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, rad=self.gradius[self.igal],mag=10, BA = .8, PA=0,nsersic=2, convolution_size=80,auto=self.auto)
                else:
                    self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=14, BA = .8, PA=0,nsersic=2, convolution_size=80,auto=self.auto)                
            elif (ncomp == 1) & (asym == 1):
                #self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=self.nsa.rmag[self.igal], BA = self.nsa.cat.SERSIC_BA[self.igal], PA=self.nsa.cat.SERSIC_PHI[self.igal],nsersic=self.nsa.cat.SERSIC_N[self.igal], convolution_size=80,asym=1)
                self.galfit = galfitwindow(self.gwindow, self.logger, image = self.galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, mag=14, BA = .8, PA=5,nsersic=2, convolution_size=80,asym=1)                
            elif ncomp == 2:
                # use results from 1 component fit as input
                try:
                    mag = self.table['GAL_MAG'][self.igal]
                    re = self.table['GAL_RE'][self.igal]
                    BA = self.table['GAL_BA'][self.igal]
                    PA = self.table['GAL_BA'][self.igal]
                
                except KeyError:
                    print('WARNING!!!!')
                    print('trouble reading galfit results from data table')
                    print('make sure you run 1 component fit first')
                    return
                ########################
                # assume bulge contains 20% of light for initial guess
                ########################
                mag_disk = mag+.25
                mag_bulge = mag + 1.75
        
                ########################
                # require n=1 for disk, n=4 for bulge (allow bulge to vary)
                # also start PA=0 and BA=1 for bulge
                ########################
                nsersic_disk=1
                nsersic_bulge=4
        
                ########################
                # set re=1.5*re_initial for disk
                # set re = 0.5*re_initial for bulge
                ########################
                re1=1.2*re
                re2=.5*re
                    
                self.galfit = galfitwindow(self.gwindow, self.logger, image = galimage, mask_image = self.mask_image_name, psf=psf, psf_oversampling = psf_oversampling, ncomp=ncomp, rad=re1, mag=mag_disk, BA=BA, PA=PA,nsersic=nsersic_disk,nsersic2=nsersic_bulge,mag2=mag_bulge, rad2=re2, fitn=False, fitn2=True, convolution_size=80)

            if not self.auto:
                self.galfit.model_saved.connect(self.galfit_save)        
                self.galfit.setupUi(self.gwindow)
                self.gwindow.show()
            else:
                self.galfit_save(None)
        except ValueError:
            print('WARNING - ERROR RUNNING GALFIT!!!')
            print('Make sure you have measured the PSF and made a mask!')
            print("error:", sys.exc_info()[0])
        #except:
        #    print("Unexpected error:", sys.exc_info()[0])
        #    raise
    def galfit_save(self,msg): # MVC - model?
        #print('galfit model saved!!!',msg)
        if self.testing:
            self.ncomp = int(msg)
            print('ncomp = ',self.ncomp)
        if self.galha:
            prefix = 'GAL_H'
        else:
            prefix = 'GAL_'
        if (self.ncomp == 1) & (self.asym == 0):
            if self.galha:
                self.galfit_hresults = self.galfit.galfit_results
            else:
                self.galfit_results = self.galfit.galfit_results
            fields = ['XC','YC','MAG','RE','N','BA','PA']
            values = np.array(self.galfit.galfit_results[:-2])[:,0].tolist()
            for i,f in enumerate(fields):
                colname = prefix+f
                if i == 3: # multiply radius by pixel scale
                    self.table[colname][self.igal]=values[i]*self.pixelscale
                else:
                    self.table[colname][self.igal]=values[i]
            fields = ['XC','YC','MAG','RE','N','BA','PA']
            print("testing, galfit_results = ", self.galfit.galfit_results[:-2])
            values = np.array(self.galfit.galfit_results[:-2])[:,1].tolist()
            for i,f in enumerate(fields):
                colname = prefix+f+'_ERR'
                if i == 3: # convert radius from pixels to arcsec
                    self.table[colname][self.igal]=values[i]*self.pixelscale
                else:
                    self.table[colname][self.igal]=values[i]
            fields = ['SKY','CHISQ']
            values = [self.galfit.galfit_results[-2],self.galfit.galfit_results[-1]]
            for i,f in enumerate(fields):
                colname = prefix+f
                self.table[colname][self.igal]=values[i]
            wcs = WCS(self.galimage)
            ra,dec = wcs.wcs_pix2world(self.galfit.galfit_results[0][0],self.galfit.galfit_results[1][0],0)
            self.table[prefix+'RA'][self.igal]=ra
            self.table[prefix+'DEC'][self.igal]=dec
            #self.update_gui_table()

            # save results to output table
            #self.table['GAL_SERSIC'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,0]
            #self.table['GAL_SERSIC_ERR'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,1]
            #self.table['GAL_SERSIC_SKY'][self.igal] = (self.galfit.galfit_results[-2])
            #self.table['GAL_SERSIC_CHISQ'][self.igal] = (self.galfit.galfit_results[-1])
            #print(self.table[self.igal])
        elif (self.ncomp == 1) & (self.asym == 1):
            if self.galha:
                self.galfit_haasym_results = self.galfit.galfit_results
            else:
                self.galfit_asym_results = self.galfit.galfit_results
            #print('writing results for galfit w/asymmetry')
            self.table[prefix+'SERSASYM'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,0]
            self.table[prefix+'SERSASYM_ERR'][self.igal] = np.array(self.galfit.galfit_results[:-2])[:,1]
            self.table[prefix+'SERSASYM_ERROR'][self.igal] = (self.galfit.galfit_results[-2])
            self.table[prefix+'SERSASYM_CHISQ'][self.igal] = (self.galfit.galfit_results[-1])
            wcs = WCS(self.galimage)
            ra,dec = wcs.wcs_pix2world(self.galfit.galfit_results[0][0],self.galfit.galfit_results[1][0],0)
            self.table[prefix+'SERSASYM_RA'][self.igal]=ra
            self.table[prefix+'SERSASYM_DEC'][self.igal]=dec
            self.update_gui_table()
        elif self.ncomp == 2:
            self.galfit_results2 = self.galfit.galfit_results
            #print(self.galfit_results2)
            self.table['GAL_2SERSIC'][self.igal] = np.array(self.galfit_results2[:-2])[:,0]
            self.table['GAL_2SERSIC_ERR'][self.igal] = np.array(self.galfit_results2[:-2])[:,1]
            self.table['GAL_2SERSIC_ERROR'][self.igal] = (self.galfit_results2[-2])
            self.table['GAL_2SERSIC_CHISQ'][self.igal] = (self.galfit_results2[-1])
            #print(self.table[self.igal])
        if not self.auto:
            self.update_gui_table()
        

    def galfit_ellip_phot(self): # MVC - model
        '''
        use galfit ellipse parameters as input for photutils elliptical photometry

        '''
        ### CLEAR R-BAND CUTOUT CANVAS
        #self.rcutout.canvas.delete_all_objects()

        ### MAKE SURE GALFIT 1 COMP MODEL WAS RUN

        try:
            xc,yc,mag,re,n,BA,pa = np.array(self.galfit_results[:-3])[:,0].tolist()
            #print('GALFIT PA = ',xc,yc,mag,re,n,BA,pa )
        except AttributeError:
            print('Warning - galfit 1 comp fit results not found!')
            print('Make sure you run galfit, then try again.')
            return
        
        ### FIT ELLIPSE
        #
        if self.auto:
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = None,image2_filter=self.hafilter, filter_ratio=self.filter_ratio,psf=self.psf_image_name,psf_ha=self.psf_haimage_name)
        else:
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = self.rcutout,image2_filter=self.hafilter, filter_ratio=self.filter_ratio,psf=self.psf_image_name,psf_ha=self.psf_haimage_name)
        #fields = ['XC','YC','MAG','RE','N','BA','PA']

        # TRANSFORM THETA
        # GALFIT DEFINES THETA RELATIVE TO Y AXIS
        # PHOT UTILS DEFINES THETA RELATIVE TO X AXIS
        THETA = pa + 90 # in degrees
        #print('THETA = ',THETA)
        self.e.run_with_galfit_ellipse(xc,yc,BA=BA,THETA=THETA)
        self.e.plot_profiles()
        #os.chdir(current_dir)

        '''
        fields = ['ASYM','ASYM_ERR','ASYM2','ASYM2_ERR']
        values = [self.e.asym, self.e.asym_err, self.e.asym2,self.e.asym2_err]
        for i,f in enumerate(fields):
            colname = 'GAL_'+f
            self.table[colname][self.igal]=values[i]
        self.update_gui_table()
        '''

        # fit profiles
        self.fit_profiles(prefix='GAL')
        # save results
        self.write_profile_fits(prefix='GAL_')

        # TODO - this is calling a view function - should not do that.  need to separate model and view, or call both from controller
        if not self.auto:
            self.draw_ellipse_results(color='cyan')

    def photutils_ellip_phot(self):
        #current_dir = os.getcwd()
        #image_dir = os.path.dirname(self.rcoadd_fname)
        #os.chdir(image_dir)

        ### CLEAR R-BAND CUTOUT CANVAS
        #self.rcutout.canvas.delete_all_objects()

        ### FIT ELLIPSE
        #
        if self.verbose:
            print()
            print("running photutils_ellip_phot")
            print()
        if self.auto:
            ra = self.objparams[0]
            dec = self.objparams[1]            
            
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = None,image2_filter='16', filter_ratio=self.filter_ratio, psf=self.psf_image_name,psf_ha=self.psf_haimage_name,objra=ra,objdec=dec)
        else:
            ra = self.defcat.cat['RA'][self.igal]
            dec = self.defcat.cat['DEC'][self.igal]            
            self.e = ellipse(self.cutout_name_r, image2=self.cutout_name_ha, mask = self.mask_image_name, image_frame = self.rcutout,image2_filter='16', filter_ratio=self.filter_ratio, psf=self.psf_image_name,psf_ha=self.psf_haimage_name,objra=ra,objdec=dec )
        self.e.run_for_gui(runStatmorphFlag=True)

        if self.verbose:
            print("calling plot_profiles\n")
        self.e.plot_profiles()
        if self.verbose:
            print("finished plot_profiles\n")
        
        #os.chdir(current_dir)
        if self.verbose:
            print("saving data\n")

        ### SAVE DATA TO TABLE
        fields = ['BADGAL','XCENTROID','YCENTROID','EPS','THETA','GINI','HGINI',\
                  'M20','HM20',
                  'UNMASKED_AREA','TOTAL_AREA',\
                  'SUM','SUM_MAG','ASYM','ASYM_ERR',\
                  'HSUM','HSUM_MAG','HASYM','HASYM_ERR']#,'SUM_ERR']
        values = [self.bad_galaxy,self.e.xcenter, self.e.ycenter,self.e.eps, np.degrees(self.e.theta), \
                  self.e.cat.gini[self.e.objectIndex],self.e.cat2.gini[self.e.objectIndex],\
                  self.e.M20_1,self.e.M20_2,\
                  self.e.cat[self.e.objectIndex].area.value*self.pixelscale*self.pixelscale,\
                  #self.e.cat[self.e.objectIndex].segment_area.value*self.pixelscale*self.pixelscale,\
                  self.e.masked_pixel_area*self.pixelscale*self.pixelscale,\
                  self.e.source_sum_erg, self.e.source_sum_mag,self.e.asym, self.e.asym_err, \
                  self.e.source_sum2_erg,self.e.source_sum2_mag,self.e.asym2,self.e.asym2_err]
        for i,f in enumerate(fields):
            if i == 0:
                colname = f
            else:
                colname = 'ELLIP_'+f
            #print(colname)
            #self.table[colname][self.igal]=float('%.2e'%(values[i]))            
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
        # update sky noise
        fields = ['R_SKYNOISE','H_SKYNOISE']
        values = [self.e.im1_skynoise/1.e-17,self.e.im2_skynoise/1.e-17]
        #print("before writing skynoise: ",values)
        for i,f in enumerate(fields):
            #print(values[i])
            self.table[f][self.igal] = values[i]
        fields = ['R_SKY','H_SKY']
        values = [self.e.sky,self.e.sky2]
        #print("before writing sky values: ",values)
        for i,f in enumerate(fields):
            #print(values[i])
            self.table[f][self.igal] = values[i]

        
        # what are we doing here?
        wcs = WCS(self.cutout_name_r)
        ra,dec = wcs.wcs_pix2world(self.e.xcenter,self.e.ycenter,0)
        
        ra2,dec2 = wcs.wcs_pix2world(self.e.cat2.xcentroid[self.e.objectIndex],self.e.cat2.ycentroid[self.e.objectIndex],0)        
        self.table['ELLIP_RA'][self.igal]=ra
        self.table['ELLIP_DEC'][self.igal]=dec
        self.table['ELLIP_HRA'][self.igal]=ra2
        self.table['ELLIP_HDEC'][self.igal]=dec2
        self.write_fits_table()
        # TODONE - write out phot table
        colnames = ['area',
                    'background_mean',
                    'bbox_xmax',
                    'bbox_xmin',
                    'bbox_ymax',
                    'bbox_ymin',
                    'cxx',
                    'cxy',
                    'cyy',
                    'eccentricity',
                    'ellipticity',
                    'elongation',
                    'equivalent_radius',
                    'fwhm',
                    'gini',
                    'inertia_tensor',
                    'kron_flux',
                    'kron_fluxerr',
                    'kron_radius',
                    'local_background',
                    'moments', 
                    'moments_central',
                    'orientation',
                    'perimeter',
                    'segment_flux',
                    'segment_fluxerr',
                    'semimajor_sigma',
                    'semiminor_sigma',
                    'xcentroid',
                    'ycentroid']


        # this is hanging when trying to calculate for cat2 - not sure why
        # skipping this for now.
        # if self.verbose:
        #     print("Calculating fluxfrac_radius 30 for cat2")
            
        # r30 = self.e.cat2.fluxfrac_radius(0.3)*self.pixelscale*u.arcsec/u.pixel
        # if self.verbose:
        #     print("Calculating fluxfrac_radius 50 for cat2")
        
        # r50 = self.e.cat2.fluxfrac_radius(0.5)*self.pixelscale*u.arcsec/u.pixel
        # if self.verbose:
        #     print("Calculating fluxfrac_radius 90 for cat2")
        
        # r90 = self.e.cat2.fluxfrac_radius(0.9)*self.pixelscale*u.arcsec/u.pixel

        # if self.verbose:
        #     print()
        #     print("adding r30 to e.cat2")
        #     print()
        

        # self.e.cat2.add_extra_property('PHOT_R30',r30)

        # if self.verbose:
        #     print()
        #     print("adding r50 to e.cat2")
        #     print()
        
        # self.e.cat2.add_extra_property('PHOT_R50',r50)

        # if self.verbose:
        #     print()
        #     print("adding r90 to e.cat2")
        #     print()
        
        # self.e.cat2.add_extra_property('PHOT_R90',r90)

        if self.verbose:
            print("Calculating fluxfrac_radius")
        # calculate fractional radii, but these are circular, and in pixels
        r30 = self.e.cat.fluxfrac_radius(0.3)*self.pixelscale*u.arcsec/u.pixel
        r50 = self.e.cat.fluxfrac_radius(0.5)*self.pixelscale*u.arcsec/u.pixel
        r90 = self.e.cat.fluxfrac_radius(0.9)*self.pixelscale*u.arcsec/u.pixel

        if self.verbose:
            print()
            print("done calculating fluxfrac_radius")
        # calculate fractional radii, but these are circular, and in pixels

        if self.verbose:
            print()
            print("adding extra properties to e.cat")
            print()
        self.e.cat.add_extra_property('PHOT_R30',r30)
        self.e.cat.add_extra_property('PHOT_R50',r50)
        self.e.cat.add_extra_property('PHOT_R90',r90)



        if self.verbose:
            print("writing radii to main table")
        # write these out to the main table
        fields = ['R30','R50','R90',\
                  'HR30','HR50','HR90']
        
        #values = [self.e.cat.PHOT_R30[self.e.objectIndex].value,\
        #          self.e.cat.PHOT_R50[self.e.objectIndex].value,\
        #          self.e.cat.PHOT_R90[self.e.objectIndex].value,\
        #          self.e.cat2.PHOT_R30[self.e.objectIndex].value,\
        #          self.e.cat2.PHOT_R50[self.e.objectIndex].value,\
        #          self.e.cat2.PHOT_R90[self.e.objectIndex].value]
        values = [self.e.cat.PHOT_R30[self.e.objectIndex].value,\
                  self.e.cat.PHOT_R50[self.e.objectIndex].value,\
                  self.e.cat.PHOT_R90[self.e.objectIndex].value,\
                  0,0,0]
        for i,f in enumerate(fields):
            colname = 'ELLIP_'+f
            #print(colname,values[i])
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i].value))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except AttributeError:
                print("AttributeError: ",colname, values[i])
                print("sorry for the shit show...")
                self.table[colname][self.igal]=float('%.4e'%(values[i]))                
            #except:
            #    print("problem writing table element",colname,values[i])

        if self.verbose:
            print("writing fits table\n")
        self.write_fits_table()
        
        if self.e.statmorph_flag:
            print()
            print("running statmorph")
            self.write_statmorph()
            if self.verbose:
                print("writing fits table after running statmorph")
            self.write_fits_table()
        #c1 = Column(data=np.array(r30[self.e.objectIndex]),name='PHOTR30',unit='arcsec',description='photutils fluxfrac_radius')
        #c2 = Column(data=np.array(r50[self.e.objectIndex]),name='PHOTR50',unit='arcsec',description='photutils fluxfrac_radius')
        #c3 = Column(data=r90[self.e.objectIndex],name='PHOTR90',unit='arcsec',description='photutils fluxfrac_radius')
        #qtable.add_columns([c1,c2,c3])

        if self.verbose:
            print("setting up photutil_tab.fits")
        qtable = self.e.cat[self.e.objectIndex].to_table(colnames)
        
        phot_table_name = self.cutout_name_r.replace('.fits','-photutil_tab.fits')
        qtable = Table(qtable)
        qtable.write(phot_table_name,format='fits',overwrite=True)
        
       
        if not self.auto:
            self.update_gui_table()

        # convert theta to degrees, and subtract 90 to get angle relative to y axis
        #self.e.theta = np.degrees(self.e.theta) - 90
        # fit profiles

        if args.verbose:
            print("in photutils_ellip_phot, fitting profiles")
        self.fit_profiles()
        # save results
        self.write_profile_fits()

        
        if not self.auto:
            self.draw_ellipse_results(color='magenta')
    
    def write_statmorph(self):
        #########################################################
        ## ADD STATMORPH PARAMETERS
        #########################################################

        
        # write these out to the main table
        fields = ['XCENTROID','YCENTROID',\
                  'RPETRO_CIRC','RPETRO_ELLIP','RHALF_ELLIP',\
                  'R20','R80',\
                  'GINI','M20','F_GM20','S_GM20',\
                  'C','A','S','FLAG']
        
        values = [self.e.morph.xc_centroid,\
                  self.e.morph.yc_centroid,\
                  self.e.morph.rpetro_circ*self.pixelscale,\
                  self.e.morph.rpetro_ellip*self.pixelscale,\
                  self.e.morph.rhalf_ellip*self.pixelscale,\
                  self.e.morph.r20*self.pixelscale,\
                  self.e.morph.r80*self.pixelscale,\
                  self.e.morph.gini,\
                  self.e.morph.m20,\
                  self.e.morph.gini_m20_bulge,\
                  self.e.morph.gini_m20_merger,\
                  self.e.morph.concentration,\
                  self.e.morph.asymmetry,\
                  self.e.morph.smoothness,\
                  self.e.morph.flag]
                  
        for i,f in enumerate(fields):
            colname = 'SMORPH_'+f
            #print(colname)
            try:
                self.table[colname][self.igal]=float('%.4e'%(values[i]))
            except KeyError:
                print("KeyError: ",colname)
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()
            except TypeError:
                print("TypeError: ",colname, values[i])
                print("sorry for the shit show...")
                print("\ntable column names: \n",self.table.colnames)
                sys.exit()

        if self.e.statmorph_flag2:
            ## Add Halpha values
            values = [self.e.morph2.xc_centroid,\
                      self.e.morph2.yc_centroid,\
                      self.e.morph2.rpetro_circ*self.pixelscale,\
                      self.e.morph2.rpetro_ellip*self.pixelscale,\
                      self.e.morph2.rhalf_ellip*self.pixelscale,\
                      self.e.morph2.r20*self.pixelscale,\
                      self.e.morph2.r80*self.pixelscale,\
                      self.e.morph2.gini,\
                      self.e.morph2.m20,\
                      self.e.morph2.gini_m20_bulge,\
                      self.e.morph2.gini_m20_merger,\
                      self.e.morph2.concentration,
                      self.e.morph2.asymmetry,
                      self.e.morph2.smoothness,\
                      self.e.morph2.flag]

            for i,f in enumerate(fields):
                colname = 'SMORPH_H'+f
                #print(colname)
                try:
                    self.table[colname][self.igal]=float('%.4e'%(values[i]))
                except KeyError:
                    print("KeyError: ",colname)
                    print("\ntable column names: \n",self.table.colnames)
                    sys.exit()
                except TypeError:
                    print("TypeError: ",colname, values[i])
                    print("sorry for the shit show...")
                    print("\ntable column names: \n",self.table.colnames)
                    sys.exit()

        
    def fit_profiles(self,prefix=None):
        #current_dir = os.getcwd()
        #image_dir = os.path.dirname(self.rcoadd_fname)
        #os.chdir(image_dir)
        if prefix is None:
            rphot_table = self.cutout_name_r.split('.fits')[0]+'_phot.fits'
            haphot_table = self.cutout_name_ha.split('.fits')[0]+'_phot.fits'
        else:
            rphot_table = self.cutout_name_r.split('.fits')[0]+'-'+prefix+'_phot.fits'
            haphot_table = self.cutout_name_ha.split('.fits')[0]+'-'+prefix+'_phot.fits'

        self.rfit = rprofile(self.cutout_name_r, rphot_table, label='R')
        self.rfit.becky_measurements()
        self.hafit = haprofile(self.cutout_name_ha, haphot_table, label=r"$H\alpha$")
        #if self.hafit.fit_flag:
        self.hafit.becky_measurements()
        self.hafit.get_r24_stuff(self.rfit.iso_radii[self.rfit.isophotes == 24.][0][0])

        # TODO - is this mixing model and view????  ARGGGHHHHHH
        both = dualprofile(self.rfit,self.hafit)
        try:
            both.make_3panel_plot()
        except:
            print()            
            print("problem making 3panel plot - weird values...")
            print()
    def write_profile_fits(self,prefix=None): # MVC - model
        fields = ['R24','R25','R26','R24V','R25V',\
                  'R_F25','R_F50','R_F75',\
                  'M24','M25','M26',\
                  'F_30R24','F_R24','C30',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        d = self.rfit
        values = [d.iso_radii[0],d.iso_radii[1],d.iso_radii[2],d.iso_radii[3],d.iso_radii[4],\
                  d.flux_radii[0],d.flux_radii[1],d.flux_radii[2],\
                  d.iso_mag[0],d.iso_mag[1],d.iso_mag[2],\
                  d.flux_30r24,d.flux_r24,d.c30,\
                  d.petrorad,d.petroflux_erg,d.petror50,d.petror90,d.petrocon,d.petromag
                  ]
        for i,f in enumerate(fields):
            if prefix is None:
                colname = f
            else:
                colname = prefix+f
            #print(colname, values[i])
            self.table[colname][self.igal]=float('%.4e'%(values[i][0]))
            self.table[colname+'_ERR'][self.igal]=float('%.4e'%(values[i][1]))
            
        fields = ['R16','R17','R_F25','R_F50','R_F75','M16','M17','F_30R24','F_R24','C30','R_F95R24','F_TOT',\
                  'PETRO_R','PETRO_FLUX','PETRO_R50','PETRO_R90','PETRO_CON','PETRO_MAG']
        d = self.hafit
        values = [d.iso_radii[0],d.iso_radii[1],\
                  d.flux_radii[0],d.flux_radii[1],d.flux_radii[2],\
                  d.iso_mag[0],d.iso_mag[1],\
                  d.flux_30r24,d.flux_r24,d.c30,d.flux_95r24, d.total_flux,\
                  d.petrorad,d.petroflux_erg,d.petror50,d.petror90,d.petrocon,d.petromag
                  ]
        for i,f in enumerate(fields):
            if prefix is None:
                colname = 'H'+f
            else:
                colname = prefix+'H'+f

            #print(colname,values[i])
            self.table[colname][self.igal]=float('{:.4e}'.format(values[i][0]))
            self.table[colname+'_ERR'][self.igal]=float('{:.4e}'.format(values[i][1]))

        # SFR conversion from Kennicutt and Evans (2012)
        # log (dM/dt/Msun/yr) = log(Lx) - logCx
        logCx = 41.27
        #print(len(self.hafit.total_flux),len(self.gzdist))
        L = self.hafit.total_flux*(4.*np.pi*cosmo.luminosity_distance(self.gzdist[self.igal]).cgs.value**2)
        #print(L)
        detect_flag = L > 0
        self.sfr = np.zeros(len(L),'d')
        self.sfr[detect_flag] = np.log10(L[detect_flag]) - logCx
        if prefix is None:
            colname='LOG_SFR_HA'
        else:
            colname=prefix+'LOG_SFR_HA'
        #print('sfr = ',self.sfr)
        #print(self.sfr[0], self.sfr[1])
        self.table[colname][self.igal]=float('%.4e'%(self.sfr[0]))
        self.table[colname+'_ERR'][self.igal]=float('%.4e'%(self.sfr[1]))
        self.table[colname+'_FLAG'][self.igal]=detect_flag[0]
        # inner ssfr
        a = self.hafit.flux_30r24
        b = self.rfit.flux_30r24
        self.inner_ssfr = a[0]/b[0]
        self.inner_ssfr_err = ratio_error(a[0],b[0],a[1],b[1])
        if prefix is None:
            colname='SSFR_IN'
        else:
            colname = prefix+'SSFR_IN'
        self.table[colname][self.igal]=float('%.4e'%(self.inner_ssfr))
        self.table[colname+'_ERR'][self.igal]=float('%.4e'%(self.inner_ssfr_err))
        # outer ssfr
        c = self.hafit.flux_r24
        d = self.rfit.flux_r24
        self.outer_ssfr = (c[0] - a[0])/(d[0] - b[0])
        self.outer_ssfr_err = ratio_error(c[0] - a[0],d[0] - b[0],np.sqrt(a[1]**2 + c[1]**2),np.sqrt(b[1]**2 + d[1]**2))
        if prefix is None:
            colname='SSFR_OUT'
        else:
            colname=prefix+'SSFR_OUT'
        self.table[colname][self.igal]=float('%.4e'%(self.outer_ssfr))
        self.table[colname+'_ERR'][self.igal]=float('%.4e'%(self.outer_ssfr_err))
        self.write_fits_table()        
        if not self.auto:
            self.update_gui_table()

class hagui_interactive():
    """ 
    class to handle all gui setup and visualizations 

    I think this needs to inherit the model class, because buttons are connected to model functions

    Or the functions that are connecting the buttons should be called from the controller class

    """

    def setup_gui(self): # view
        #print(MainWindow)

        self.ui = Ui_MainWindow()        
        self.ui.setupUi(MainWindow)
        #self.ui.setGeometry(0,0,400,300)
        #self.ui.setFont(QtGui.Qfont('Arial',10))
        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()
        self.add_coadd_frame(self.ui.leftLayout)
        self.add_cutout_frames()
        self.connect_setup_menu()
        self.connect_ha_menu()
        self.connect_halpha_type_menu()
        self.connect_comment_menu()
        #self.connect_buttons()
        #self.add_image(self.ui.gridLayout_2)
        #self.add_image(self.ui.gridLayout_2)
        self.connect_buttons()
    def connect_buttons(self): #view
        self.ui.wmark.clicked.connect(self.find_galaxies)
        #self.ui.editMaskButton.clicked.connect(self.edit_mask)
        self.ui.makeMaskButton.clicked.connect(self.make_mask)
        #self.ui.saveCutoutsButton.clicked.connect(self.write_cutouts)
        self.ui.fitEllipseGalfitButton.clicked.connect(self.galfit_ellip_phot)
        self.ui.fitEllipsePhotutilsButton.clicked.connect(self.photutils_ellip_phot)
        self.ui.wfratio.clicked.connect(self.get_filter_ratio)
        self.ui.resetRatioButton.clicked.connect(self.reset_cutout_ratio)
        self.ui.resetSizeButton.clicked.connect(self.reset_cutout_size)
        self.ui.prefixLineEdit.textChanged.connect(self.set_prefix)
        #self.ui.fitEllipseButton.clicked.connect(self.fit_ellipse_phot)
        self.ui.galfitButton.clicked.connect(lambda: self.run_galfit(ncomp=1))
        self.ui.galfitAsymButton.clicked.connect(lambda: self.run_galfit(ncomp=1,asym=1))
        self.ui.galfitHaButton.clicked.connect(lambda: self.run_galfit(ncomp=1,ha=1))
        self.ui.galfitHaAsymButton.clicked.connect(lambda: self.run_galfit(ncomp=1,asym=1,ha=1))
        self.ui.galfit2Button.clicked.connect(lambda: self.run_galfit(ncomp=2))
        self.ui.psfButton.clicked.connect(self.build_psf)
        self.ui.saveButton.clicked.connect(self.write_fits_table)
        self.ui.clearCutoutsButton.clicked.connect(self.clear_cutouts)
        self.ui.filterRatioLineEdit.returnPressed.connect(self.set_filter_ratio)
        self.ui.cutoutSizeLineEdit.returnPressed.connect(self.set_cutout_size)
    def add_coadd_frame(self,panel_name): # view
        logger = log.get_logger("example1", log_stderr=True, level=40)
        self.coadd = image_panel(panel_name, self.ui,logger)
        self.coadd.key_pressed.connect(self.key_press_func)
        #self.coadd.add_cutouts()

    def add_cutout_framesv4(self):# with halphav4 # view
        # r-band cutout
        self.rcutout_label = QtWidgets.QLabel('r-band')
        self.ui.cutoutsLayout.addWidget(self.rcutout_label, 0, 0, 1, 1)
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.rcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 8, 1,autocut_params='histogram')
        self.hacutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 8, 1,autocut_params='stddev')
        self.maskcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 8, 1,autocut_params='stddev')
    def add_cutout_frames(self): # view
        # r-band cutout
        self.rcutout_label = QtWidgets.QLabel('r-band')
        drow = 20
        self.ui.cutoutsLayout.addWidget(self.rcutout_label, 0, 0, 1, 2)
        temp_label = QtWidgets.QLabel('')
        self.ui.cutoutsLayout.addWidget(temp_label, 0, 2, 1, 2)
        self.nsa_label = QtWidgets.QLabel('NSA ID')
        self.ui.cutoutsLayout.addWidget(self.nsa_label, 0, 4, 1, 2)
        temp_label = QtWidgets.QLabel('')
        self.ui.cutoutsLayout.addWidget(temp_label, 0, 6, 1, 2)
        self.agc_label = QtWidgets.QLabel('AGC Number')
        self.ui.cutoutsLayout.addWidget(self.agc_label, 0, 8, 1, 2)
        self.rcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, drow, 10,autocut_params='stddev')
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, drow+1, 0, 1, 2)
        self.hacutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, drow+2, 0, drow, 10,autocut_params='stddev')

        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 2*drow+2, 0, 1, 2)
        self.maskcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,2*drow+3, 0, drow, 10)
        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)



    def clear_cutouts(self): # view - should this be in the controller class?
        self.rcutout.canvas.delete_all_objects()
        self.hacutout.canvas.delete_all_objects()
    def connect_setup_menu(self): # view or controller - these are calling model functions/quantities
        self.ui.actionR_coadd.triggered.connect(self.get_rcoadd_file)
        self.ui.actionHa_coadd_2.triggered.connect(self.get_hacoadd_file)
        self.ui.actionNSA_catalog_path.triggered.connect(self.getnsafile)
        self.ui.actionAGC_catalog_path.triggered.connect(self.getagcfile)
        
    def connect_ha_menu(self): # view
        #print('working on this')
        #extractAction.triggered.connect(self.close_application)

        self.ui.actionhalpha4.triggered.connect(lambda: self.set_hafilter('4'))
        self.ui.actionhalpha8.triggered.connect(lambda: self.set_hafilter('8'))
        self.ui.actionhalpha12.triggered.connect(lambda: self.set_hafilter('12'))
        self.ui.actionhalpha16.triggered.connect(lambda: self.set_hafilter('16'))
        self.ui.actioninthalpha.triggered.connect(lambda: self.set_hafilter('inthalpha'))
        self.ui.actionintha6657.triggered.connect(lambda: self.set_hafilter('intha6657'))
        self.ui.actionsienaha.triggered.connect(lambda: self.set_hafilter('sienaha'))
        
    def connect_halpha_type_menu(self): # view
        ha_types = ['Ha Emission','No Ha']
        for name in ha_types:
            self.ui.haTypeComboBox.addItem(str(name))
        self.ui.haTypeComboBox.activated.connect(self.set_halpha_type)
    def connect_comment_menu(self): # view
        comment_types = ['Cont Sub Prob','merger/tidal','scat light','asym R', 'asym Ha','fore. star', 'fore. gal','edge-on','part cov','nuc ha']
        for name in comment_types:
            self.ui.commentComboBox.addItem(str(name))
        self.ui.commentComboBox.activated.connect(self.set_comment)
    def set_prefix(self,prefix): # MVC - probably model? - no this is from the gui
        self.prefix = prefix
        #print('prefix for output files = ',self.prefix)
    def set_prefix_on_gui(self,prefix): # MVC - probably model? - no this is from the gui
        """  use this if the prefix is provided by the user - this will fill in the box with the provided prefix """
        self.ui.prefixLineEdit.setText(prefix)        
        
    def mark_galaxies(self): # MVC - view or controller, b/c this relies on model quantities
        #
        # using code in TVMark.py as a guide for adding shapes to canvas
        #
        #

        objlist = []
        markcolor='cyan'
        markwidth=1
        #size = cutout_scale*self.gradius
        size = self.cutout_sizes
        #size[size > self.global_max_cutout_size] = self.global_max_cutout_size
        #size[size < self.global_min_cutout_size] = self.global_min_cutout_size
        print(f"DEBUGGING: in mark_galaxies: len(gximage) = {len(self.gximage)}")
        for i,x in enumerate(self.gximage):
            print(f"{i}, {self.galid[i]}, cutout_size = {size[i]:.2f}, gximage[i]={self.gximage[i]:.1f}, agcximage={self.agcximage[i]:.1f}") 
            obj = self.coadd.dc.Box(
                x=x, y=self.gyimage[i], \
                xradius=size[i]/2,yradius=size[i]/2, \
                #xradius=100,yradius=100, \
                color=markcolor, linewidth=markwidth)
            glabel = self.coadd.dc.Text(x-50,self.gyimage[i]+60,\
                                        str(self.galid[i]), color=markcolor)
            objlist.append(obj)
            objlist.append(glabel)
        if self.agcflag:
            for i,x in enumerate(self.agcximage):
                #print(x,self.agcyimage[i],self.agc.cat.AGCNUMBER[i])
                obj = self.coadd.dc.Box(
                    x=x, y=self.agcyimage[i], xradius=75,\
                    yradius=75, color='purple', linewidth=markwidth)
                glabel = self.coadd.dc.Text(x-40,self.agcyimage[i]+40,\
                                        str(self.agc.cat['AGCnr'][i]), color='purple')
                objlist.append(obj)
                objlist.append(glabel)
            
        self.markhltag = self.coadd.canvas.add(self.coadd.dc.CompoundObject(*objlist))
        self.coadd.fitsimage.redraw()
    def setup_ratio_slider(self): # MVC - view
        self.ui.ratioSlider.setRange(0,100)
        self.ui.ratioSlider.setValue(50)
        self.ui.ratioSlider.setSingleStep(1)
        #self.ui.ratioSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self.ui.ratioSlider.setFocusPolicy(QtCore.StrongFocus)
        self.ui.ratioSlider.valueChanged.connect(self.ratio_slider_changed)
    def setup_cutout_slider(self): # MVC - view
        self.ui.cutoutSlider.setRange(0,100)
        self.ui.cutoutSlider.setValue(50)
        self.ui.cutoutSlider.setSingleStep(1)
        #self.ui.ratioSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self.ui.ratioSlider.setFocusPolicy(QtCore.StrongFocus)
        self.ui.cutoutSlider.valueChanged.connect(self.cutout_slider_changed)
            
    def clear_comment_field(self): # MVC - view?
        self.ui.commentLineEdit.clear()
        
    def display_cutouts(self): # MVC - view

        # TODONE  - split into separate model and view functions
        #################################################
        # this part needs to move to the view class
        #################################################
        # only do this when running gui (as opposed to automatically)
        self.ui.cutoutSizeLineEdit.setText(str(self.cutout_size))
        self.reset_cutout_size()
        self.reset_cutout_ratio()
        
        ############################################################
        # TODONE: this is the only part that should be in view
        #cutoutHa = Cutout2D(self.ha.data, position, self.size, wcs=self.coadd_wcs, mode = 'trim')
        ((ymin,ymax),(xmin,xmax)) = self.cutoutR.bbox_original
        bbox = '[{:d}:{:d},{:d}:{:d}]'.format(int(xmin),int(xmax),int(ymin),int(ymax))        
        #print(ymin,ymax,xmin,xmax)
        
        ############################################################                

        self.rcutout.load_image(self.r[ymin:ymax,xmin:xmax])
        self.hacutout.load_image(self.halpha_cs[ymin:ymax,xmin:xmax])
        #cutoutR.plot_on_original(color='white')
        
        # ###################################################################################
        # TODO - this should not be called in the model b/c model does not interact with view
        # ###################################################################################
        self.update_gui_table_cell(self.igal,'BBOX',str(bbox))
        
    def draw_ellipse_results(self, color='cyan'): # MVC - view
        # mark r24
        markcolor=color#, 'yellow', 'cyan']
        markwidth=1
        #print('inside draw_ellipse_results')
        image_frames = [self.rcutout, self.hacutout]
        radii = self.rfit.iso_radii[:,0][0:2]
        objlist = []
        for i,im in enumerate(image_frames):
            for r in radii:
                #print('r = ',r)
                r = r/self.pixelscale
                #print('r = ',r)
                obj =im.dc.Ellipse(self.e.xcenter,self.e.ycenter,r, r*(1-self.e.eps), rot_deg = np.degrees(self.e.theta), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            if i == 1: # add R17 for Halpha image
                r = self.hafit.iso_radii[:,0][1]/self.pixelscale
                obj =im.dc.Ellipse(self.e.xcenter,self.e.ycenter,r, r*(1-self.e.eps), rot_deg = np.degrees(self.e.theta), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            self.markhltag = im.canvas.add(im.dc.CompoundObject(*objlist))
            im.fitsimage.redraw()
        # mark R17 in halpha image
        
    def display_mask(self, mask_image_name): # MVC - controller or view??? putting in view for now
        t = self.cutout_name_r.split('.fit')
        self.mask_image_name=t[0]+'-mask.fits'
        
        if not self.auto:
            if os.path.exists(self.mask_image_name):
                self.maskcutout.load_file(self.mask_image_name)
                self.mask_image_exists = True
                
            else:
                # clear mask frame
                self.maskcutout.fitsimage.clear()

        
        #self.mask_image = mask_image_name


        
class hacontroller():
    """ class to handle interactions with the user and send requests to model and view  """
    def get_rcoadd_file(self): # view or controller? controller b/c interact w/user
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            self.rcoadd_fname = fname[0]
            #print(self.rcoadd_fname)
            self.load_rcoadd()
            #self.le.setPixmap(QPixmap(fname))
    def get_hacoadd_file(self): # view or controller? controller
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            
            self.hacoadd_fname = fname[0]
            #print(self.hacoadd_fname)
            self.load_hacoadd()
            #self.le.setPixmap(QPixmap(fname))
    def getnsafile(self): # view or controller?
        """ get NSA file location from the user """
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            self.nsa_fname = fname[0]
            self.nsa = galaxy_catalog(self.nsa_fname,nsa=True)
        print('Got NSA catalog with {} lines!'.format(len(self.nsa.cat)))
        self.defcat = self.nsa
        #self.le.setPixmap(QPixmap(fname))
    def getagcfile(self): # controller
        """ get AGC file location from the user """        
        fname = QtWidgets.QFileDialog.getOpenFileName()
        if len(fname[0]) < 1:
            print('invalid filename')
        else:
            self.agc_fname = fname[0]
            self.agc = galaxy_catalog(self.agc_fname,agc=True)
            self.agcflag = True
            self.defcat = self.agc
    def set_hafilter(self,filterid): # model or controller? I think this is a mix...
        #print('setting ha filter to ',filterid)

        # this is in controller
        self.hafilter = filterid
        print('halpha filter = ',self.hafilter)
        
        # this part should be in the model, so this should be a method in model
        # and then call model.set_filter() which would execute the following lines

        self.filter_trace = filter_trace(self.hafilter)
        self.zmin = self.filter_trace.minz_trans10
        self.zmax = self.filter_trace.maxz_trans10
        #self.get_zcut()
        
    def set_halpha_type(self,hatype): # controller
        self.halpha_type = hatype
        print(hatype)
        if int(hatype) == 0:
            self.haflag[self.igal]=True
            self.update_gui_table_cell(self.igal,'HA_FLAG',str(True))
        try:
            if int(hatype) == 0:
                self.haflag[self.igal]=True
                self.update_gui_table_cell(self.igal,'HA_FLAG',str(True))
        except AttributeError:
            print('make sure you selected a galaxy')
    def set_filter_ratio(self,ratio): # controller?
        try:
            self.filter_ratio = float(ratio)
            self.subtract_images()
            self.ui.filterRatioLineEdit.setText(str(ratio))
        except:
            print('did not understand.  try again')
    def set_cutout_size(self,size): # controller?
        try:
            self.cutout_size_arcsec = int(size)*u.arcsec
            self.cutout_size_arcsec = int(size)            
            self.display_cutouts()
            self.ui.cutoutSizeLineEdit.setText(str(size))            
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
            self.ui.cutoutSizeLineEdit.setText(str(self.cutout_size))            
    def key_press_func(self,key): # MVC - is this controller or view? controller
        """ define keys to that control behavior of image display """
        print(key)
        if key == 'r':
            z = self.coadd.fitsimage.settings.get_setting('zoomlevel')
            print('zoom = ',z)
            p = self.coadd.fitsimage.settings.get_setting('pan')
            print('pan = ',p)

            self.coadd.fitsimage.set_data(self.r)
            self.coadd.fitsimage.zoom_to(z.value)
            self.coadd.fitsimage.panset_xy(p.value[0],p.value[1])
            self.coadd.canvas.redraw()
        elif key == 'h':
            try:
                self.coadd.fitsimage.set_data(self.halpha_cs)
            except AttributeError:
                print('no continuum subtracted image yet - get filter ratio')
                self.coadd.fitsimage.set_data(self.ha)
        elif key == 'u': # unidentified object!
            x,y = self.coadd.fitsimage.get_last_data_xy()
            x = float(x)
            y = float(y)
            self.uco_x.append(x)
            self.uco_y.append(y)
            ra,dec = self.coadd_wcs.wcs_pix2world(x,y,0)
            self.uco_ra.append(ra)
            self.uco_dec.append(dec)
            if len(self.uco_id) == 0:
                self.uco_id.append(1)
            else:
                self.uco_id.append(np.max(self.uco_id)+1)
            self.write_uco_table()
        elif key == 'down':
            '''
            up arrow will go to previous galaxy in the list
            '''
            #print(self.igal)
            if self.igal == (len(self.ra)-1):
                self.igal = 0
            else:
                self.igal += 1
                self.ui.wgalid.setCurrentIndex(self.igal)
                self.select_galaxy(self.igal)
            #print(self.igal)
        elif key == 'up':
            '''
            down arrow will go to previous galaxy in the list
            '''
            if self.igal == 0: 
                self.igal = len(self.ra)-1
            else:
                self.igal -= 1
                self.ui.wgalid.setCurrentIndex(self.igal)
                self.select_galaxy(self.igal)
    def ratio_slider_changed(self, value): # MVC - controller
        #print(self.minfilter_ratio, self.maxfilter_ratio, self.filter_ratio)
        delta = self.maxfilter_ratio - self.minfilter_ratio
        self.filter_ratio = self.minfilter_ratio + (delta)/100.*self.ui.ratioSlider.value()
        #print(value,' ratio slider changed to', round(self.filter_ratio,4))
        try:
            self.subtract_images(overwrite=True)
            self.display_cutouts()
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
    def cutout_slider_changed(self, value): # MVC - controller
        #print(self.minfilter_ratio, self.maxfilter_ratio, self.filter_ratio)
        delta = self.maxcutout_size - self.mincutout_size
        #self.cutout_size = self.mincutout_size + (delta)/100.*value
        self.cutout_size = self.cutout_sizes[self.igal]
        self.cutout_size_arcsec = self.cutout_sizes_arcsec[self.igal]        
        #print(value,' ratio slider changed to', round(self.filter_ratio,4))
        try:
            self.display_cutouts()
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
    def select_galaxy(self,id): # MVC - view or controller?
        print()
        print()
        print('selecting a galaxy')
        self.igal = self.ui.wgalid.currentIndex()
        if self.virgo:
            self.rcutout_label.setText('r-band '+str(self.defcat.cat['VFID'][self.igal]))
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]
            #print("new galaxy params = ",self.objparams)
            #print("compare lengths of catalogs ",len(self.defcat.cat),len(self.BA))
            print()
        elif self.uat:
            # removed this an we resolved the issue where we got the wrong galaxy in the viewer.
            #self.igal = self.igal-1 # why do we need this???
            self.rcutout_label.setText('r-band '+str(self.defcat.cat['AGCnr'][self.igal]))
            self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.defcat.cat['b'][self.igal]/self.defcat.cat['a'][self.igal],0]
            #print("new galaxy params = ",self.objparams)
            #print("compare lengths of catalogs ",len(self.defcat.cat),len(self.BA))
            print()
            
        else:
            self.rcutout_label.setText('r-band '+str(self.nsa2['NSAID'][self.igal]))
            self.objparams = None
        self.rcutout_label.show()
                                   
                                   
        print('active galaxy = ',self.igal)
        # when galaxy is selected from list, trigger
        # cutout imaages
        self.get_galaxy_cutout()

        # will need to also call new view functions
        self.display_cutouts()
        #self.display_mask()

        # the following two lines are the same as self.clear_cutouts
        # seems odd that this is called here - why are we clearing the cutouts after we displayed the cutouts???
        # need to test this when actually running gui, as opposed to running in auto mode
        self.clear_cutouts()
        #self.rcutout.canvas.delete_all_objects()
        #self.hacutout.canvas.delete_all_objects()            
        self.clear_comment_field()
        #################################################



    def reset_cutout_size(self): # MVC - controller
        self.cutout_size = self.reset_size
        self.update_images()
        #self.ui.cutoutSlider.setValue(50)
    def reset_cutout_ratio(self): # MVC - controller
        self.filter_ratio = self.reset_ratio
        self.update_images()
        #self.ui.ratioSlider.setValue(50)
        
    def update_images(self): # MVC - controller b/c calls model and view
        self.subtract_images(overwrite=True)
        #self.display_cutouts()
    def make_mask(self,objparams=None): # MVC - is this controller, or view?
        # TODO - break off view functions into a method within the haview class
        #current_dir = os.getcwd()
        #image_dir = os.path.dirname(self.rcoadd_fname)
        #os.chdir(image_dir)
        try:
            print()
            print("using ellipe parameters to unmask central region - woo hoo!")
            print()
            objparams = self.objparams
            #print("\t ",self.objparams)
        except AttributeError:
            print()
            print("problem getting objparams for masking routing")
            print()
            pass
    
    
        try:
            self.write_cutouts()
        except AttributeError:
            print('are you rushing to make a mask w/out selecting galaxies?')
            print('try selecting filter, then selecting galaxies')
            return
        #try:
        self.mwindow = QtWidgets.QWidget()
        print()
        #print("initiating mask window")
        #print("\t object params = ",objparams)
        try:
            self.mui = maskwindow(self.mwindow, self.logger, image = self.cutout_name_r, haimage=self.cutout_name_ha, sepath='~/github/halphagui/astromatic/',objparams=objparams)
        
            self.mui.mask_saved.connect(self.display_mask)
            self.mui.setupUi(self.mwindow)
            self.mwindow.show()
        except AttributeError:
            print('Hey - make sure you selected a galaxy!')
        #os.chdir(current_dir)
        
    def set_comment(self,comment): # MVC - controller or model??
        """ what is this doing? """
        col_names = ['CONTSUB_FLAG','MERGER_FLAG','SCATLIGHT_FLAG','ASYMR_FLAG','ASYMHA_FLAG','OVERSTAR_FLAG','OVERGAL_FLAG','EDGEON_FLAG','PARTIAL_FLAG','NUC_HA']
        self.table[col_names[int(comment)]][self.igal] = not(self.table[col_names[int(comment)]][self.igal])
        if not self.auto:
            self.update_gui_table_cell(self.igal,col_names[int(comment)],str(True))

        
        

    def closeEvent(self, event):
        # send signal that window is closed

        # write the data table
        print('writing fits table')
        self.table.write_fits_table()
        #event.accept()

class hafunctions(Ui_MainWindow, create_output_table, uco_table, hagui_methods, hagui_interactive, hacontroller):
    """ Main class for the halpha image analysis  """
    def __init__(self,MainWindow, logger, sepath=None, args=None):
        #testing=False,nebula=False,virgo=False,laptop=False,pointing=None,prefix=None,auto=False,obsyear=None,psfdir=None,rimage=None,haimage=None,csimage=None,filter=None,tabledir=None):
        super(hafunctions, self).__init__()
        self.auto = args.auto
        self.cscoadd_fname = None
        #if self.auto:
        #    matplotlib.use('TkAgg')
        #self.obsyear = args.obsyear

        # initialize psf as None
        self.psf = None
        if not(self.auto):
            self.setup_gui()
        self.prefix = args.prefix
        if (self.prefix is not None) and not(self.auto):
            self.set_prefix_on_gui(self.prefix)
        self.testing = args.testing
        self.draco = args.draco
        self.nebula = args.nebula        
        self.laptop = args.laptop
        self.virgo = args.virgo
        self.uat = args.uat
        self.verbose = args.verbose
        # this is the oversampling that I use when creating the PSF images
        self.oversampling = 2        
        if sepath == None:
            self.sepath = os.getenv('HOME')+'/github/halphagui/astromatic/'
        else:
            self.sepath = sepath
        if args.psfdir is None:
            self.psfdirectory = os.getcwd()
        else:
            self.psfdirectory = args.psfdir
        self.igal = None
        ############################################################
        ### CHECK TO SEE IF IMAGE NAMES ARE SPECIFIED
        ############################################################
        if args.rimage is not None:
            self.rcoadd_fname = args.rimage
        if args.haimage is not None:
            self.hacoadd_fname = args.haimage
        if args.csimage is not None:
            self.cscoadd_fname = args.csimage
        if args.filter is not None:
            self.filter = args.filter
            self.set_hafilter(args.filter)
        if args.tabledir is not None:
            self.tabledir = args.tabledir

        if (args.rimage is None):
            ############################################################
            ### CONFIGURATION SETUP FOR RUNNING ON DIFFERENT COMPUTERS
            ############################################################
        
            if self.draco & self.virgo:
                self.setup_draco_virgo()
                self.setup_virgo(pointing=pointing)
            elif self.laptop & self.virgo:
                self.setup_laptop_virgo()
                self.setup_virgo(pointing=pointing)
            elif self.nebula & self.virgo: 
                self.setup_nebula_virgo()
                if pointing is not None:
                    print('GOT A POINTING NUMBER FOR VIRGO FIELD')
                    self.setup_virgo(pointing=pointing)
                else:
                    self.setup_virgo()
            elif self.nebula:
                self.setup_nebula()
            elif self.testing:
                self.setup_testing()
        else:

            if not self.auto:
                # load rband image
                self.load_rcoadd()
            
                if args.haimage is not None:
                    # load Halpha
                    self.load_hacoadd()
            
            #print('running build_psf')
            #self.build_psf()

            
        ############################################################
        ### SET PARAMETERS FOR VIRGO VS HALPHA GROUPS PROJECT
        ### (WHICH USES THE NSA)
        ############################################################

        if self.virgo:
            self.setup_virgo_catalogs()
            
            self.defcat = self.vf
            self.def_label = 'VF.v1.'
            self.def_label = 'VF.v2.'            
            self.radius_label = 'radius'
            self.global_min_cutout_size = 60*u.arcsec
            self.global_max_cutout_size = 480*u.arcsec # 9 arcmin

        # this is for halpha groups analysis
        # commenting this out for now
        '''
        else:
            self.defcat = self.nsa
            self.def_label = 'NSAID'
            self.radius_label = 'PETROTH90'
            self.global_min_cutout_size = 100
            self.global_max_cutout_size = 250
        '''

        if self.uat:
            self.setup_uat_catalogs()
            self.defcat = self.agc
            self.global_min_cutout_size = 60*u.arcsec
            self.global_max_cutout_size = 360*u.arcsec # 6 arcmin
            
        self.initialize_uco_arrays()
        if self.auto:
            self.auto_run()

    def auto_run(self):
        # run the analysis without starting the gui
        # read in coadded images
        self.read_rcoadd()
        self.read_hacoadd()
        # set filter
        self.set_hafilter(self.filter)
        
        # get filter ratio
        # this will look for filter ratio in r-band image header
        # will run SE if ratio is not found in header
        self.get_filter_ratio()
        
        # subtract images
        # this will look for CS-ZP.fits image
        # if not found, it will subtract images according to the ZP ratio of r and halpha images
        self.subtract_images(overwrite=True) # checks out ok
        
        # measure psf
        # function will check if psf image already exists
        self.build_psf()

        # get galaxies
        self.find_galaxies()

        self.write_fits_table()
        # analyze galaxies
        # the default catalog has been trimmed to only include
        # galaxies in FOV

        # skipping for now
        # just making cutouts and getting galaxies in FOV

        # need to get a list of RA, DEC, SMA, BA, PA to feed into masking routine
        if self.verbose:
            print("starting processing of each galaxy ",len(self.gximage))
        for i in range(len(self.gximage)):
        #for i in [2]: # for testing
            self.igal = i
            # get cutouts
            self.auto_gal()
            if self.verbose:
                print("writing fits table after running auto_gal")
            self.write_fits_table()
            if self.verbose:
                print(f"##########################\nFinished galaxy {i+1}/{len(self.gximage)}")
        
    def auto_gal(self):
        # run the analysis on an individual galaxy
        if self.verbose:
            print("in auto_gal \n")
        self.bad_galaxy = False
        # create cutout
        if self.verbose:
            print("\ngetting galaxy cutouts\n")
        self.get_galaxy_cutout()

        if self.verbose:
            print("\nfinished cutouts\n")
        
        if self.bad_galaxy:
            print("\nzero std in sky - this is not real, so skipping galaxy ",self.cutout_name_r)
            print()
            return
        # create mask
        self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]

        if self.verbose:
            print("initiating maskwindow\n")

        self.mui = maskwindow(None, None, image = self.cutout_name_r, haimage=self.cutout_name_ha, \
                              sepath='~/github/halphagui/astromatic/',auto=self.auto,\
                              objparams=self.objparams,unmaskellipse=True)
                                  
        # subtract the sky, using the mask image, and resave cutouts
        
        # run galfit

        # try skipping galfit for Messier 109, VFID1167...
        if self.igal == 1167:
            # skip galfit, galaxy is just too darn big and parameters are ridiculous
            print("skipping galfit b/c output is nonsense...")
        else:
            
            try:
                self.run_galfit(ncomp=1,ha=False)
            except:
                print('WARNING: problem running galfit ellip phot',self.cutout_name_r)

            # don't need this - the geometry is so similar
            self.galfit_ellip_phot()            
            try:
                if self.verbose:
                    print("running galfit_ellip_phot ")
                self.galfit_ellip_phot()
            except:
                print("##################################")
                print('WARNING: problem running galfit ellip phot',self.cutout_name_r)
                print("##################################")                

        # run galfit ellip phot
        # use try in case fit fails
        #self.galfit_ellip_phot()

        
        # run phot util ellip phot
        # use try in case fit fails
        #self.photutils_ellip_phot()

        # RF - debugging on 07/07/2024 - not getting size measurements in output file
        # taking command out of try/except
        self.photutils_ellip_phot()
        #try:
        #    self.photutils_ellip_phot()
        #except:
        #    print('\nWARNING: problem running photutils ellip phot\n',self.cutout_name_r)
        
            
      
    def setup_testing(self):
        #self.hacoadd_fname = os.getenv('HOME')+'/research/halphagui_test/MKW8_ha16.coadd.fits'
        #self.hacoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_ha16.coadd.fits'
        self.hacoadd_fname = os.getenv('HOME')+'/research/VirgoFilaments/Halpha/virgo-coadds-2017/pointing-1_ha4.coadd.fits'

        #self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)
        self.haweight = self.hacoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.haweight_flag = True
        #self.rcoadd_fname = os.getenv('HOME')+'/research/halphagui_test/MKW8_R.coadd.fits'
        #self.rcoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_R.coadd.fits'
        self.rcoadd_fname = os.getenv('HOME')+'/research/VirgoFilaments/Halpha/virgo-coadds-2017/pointing-1_R.coadd.fits'
        #self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)

        #self.rweight = self.rcoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.rweight_flag = True
        #self.pixelscale = abs(float(self.r_header['CD1_1']))*3600. # in deg per pixel
        self.nsa_fname = os.getenv('HOME')+'/research/NSA/nsa_v0_1_2.fits'
        self.nsa = galaxy_catalog(self.nsa_fname,nsa=True)
        self.agc_fname = os.getenv('HOME')+'/research/AGC/agcnorthminus1.2019Sep24.fits'
        self.agc = galaxy_catalog(self.agc_fname,agc=True)
        self.agcflag = True
        #self.coadd.load_file(self.rcoadd_fname)
        #self.filter_ratio = 0.0416
        self.filter_ratio = 0.0406
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        self.load_hacoadd()
        self.load_rcoadd()        
        self.subtract_images()
        #self.setup_ratio_slider()
        #self.setup_cutout_slider()
    def setup_nebula(self):
        self.hacoadd_fname = '/mnt/astrophysics/reduced/20150418/MKW8_ha16.coadd.fits'
        #self.hacoadd_fname = '/mnt/qnap_home/rfinn/Halpha/reduced/virgo-coadds-2017/pointing-1_ha4.coadd.fits'
        #self.hacoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_ha16.coadd.fits'
        self.load_hacoadd()
        #self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)
        #self.haweight = self.hacoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.haweight_flag = True
        self.rcoadd_fname ='/mnt/astrophysics/reduced/20150418/MKW8_R.coadd.fits'
        #self.rcoadd_fname = '/mnt/qnap_home/rfinn/Halpha/reduced/virgo-coadds-2017/pointing-1_R.coadd.fits'
        #self.rcoadd_fname = os.getenv('HOME')+'/research/HalphaGroups/reduced_data/HDI/20150418/NRGs27_R.coadd.fits'
        #self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)
        self.load_rcoadd()
        #self.rweight = self.rcoadd_fname.split('.fits')[0]+'.weight.fits'
        #self.rweight_flag = True
        #self.pixelscale = abs(float(self.r_header['CD1_1']))*3600. # in deg per pixel
        #self.nsa_fname = '/mnt/qnap_home/share/catalogs/nsa_v0_1_2.fits'
        self.nsa_fname = '/mnt/astrophysics/catalogs/nsa_v0_1_2.fits'
        self.nsa = galaxy_catalog(self.nsa_fname,nsa=True)
        #self.agc_fname = '/mnt/qnap_home/share/catalogs/agcnorthminus1.2019Sep24.fits'
        self.agc_fname = '/mnt/astrophysics/catalogs/agcnorthminus1.2019Sep24.fits'
        self.agc = galaxy_catalog(self.agc_fname,agc=True)
        self.agcflag = True
        #self.coadd.load_file(self.rcoadd_fname)
        self.filter_ratio = 0.0422 #MKW8
        #self.filter_ratio = 0.0406
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        self.subtract_images()
        #self.setup_ratio_slider()
        #self.setup_cutout_slider()
    def setup_nebula_virgo(self):
        self.imagedir =  '/mnt/astrophysics/reduced/virgo-coadds-2017/'
        self.tabledir= '/mnt/astrophysics/catalogs/Virgo/tables-north/v1/'
    def setup_laptop_virgo(self):
        if self.obsyear == '2018':
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-2018/'
        elif self.obsyear == '2020':
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-feb2020/'
        elif self.obsyear == '2019':
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-feb2019-int/'
        else:
            self.imagedir =  '/home/rfinn/data/reduced/virgo-coadds-2017/'
        self.tabledir= '/home/rfinn/research/Virgo/tables-north/v1/'
    def setup_draco_virgo(self):
        self.imagedir =  '/data-pool/Halpha/coadds/all-virgo-coadds/'
        self.tabledir= '/home/rfinn/research/Virgo/tables-north/v2/'
    def setup_virgo(self,pointing=None):
        ''' construct image names from input.  only works for data reduced before 2021 '''
        if pointing is None:
            self.hacoadd_fname = self.imagedir+'pointing-3_ha4.coadd.fits'
            self.rcoadd_fname = self.imagedir+'pointing-3_R.coadd.fits'
        else:
            #print('got a pointing')
            if self.obsyear == '2020':
                self.hacoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_ha4.coadd.fits'
                self.rcoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_r.coadd.fits'
            else:
                self.hacoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_ha4.coadd.fits'
                self.rcoadd_fname = self.imagedir+'pointing-'+str(pointing)+'_R.coadd.fits'
            

        if not self.auto:
            self.load_hacoadd()
            self.load_rcoadd()
            try:
                self.subtract_images()
                #self.setup_ratio_slider()
                #self.setup_cutout_slider()
            except:
                pass
        self.setup_virgo_catalogs()

    def setup_virgo_catalogs(self):
        ## UPDATES TO USE VIRGO FILAMENT MASTER TABLE
        #self.vf_fname = self.tabledir+'vf_north_v1_main.fits'
        #self.nsa_fname = self.tabledir+'vf_north_v1_nsa_v0.fits'
        # RF 2023-03-24: updating to use the v2 catalogs
        self.vf_fname = os.path.join(self.tabledir,'vf_v2_main.fits')
        self.nsa_fname = os.path.join(self.tabledir,'vf_v2_nsa_v0.fits')
        ephot_fname = os.path.join(self.tabledir,'vf_v2_legacy_ephot.fits')        
        self.vf = galaxy_catalog(self.vf_fname,virgo=True)

        self.nsa = galaxy_catalog(self.nsa_fname,virgo=True)

        ##
        # get sizes for galaxies - will use this to unmask central region
        # need to cut this catalog based on keepflag
        ##

        ephot = Table.read(ephot_fname)

        
        #self.radius_arcsec = ephot['SMA_SB24']

        bad_sb25 = ephot['SMA_SB25'] == 0

        self.radius_arcsec = ephot['SMA_SB25']*(~bad_sb25) + 1.35*ephot['SMA_SB24']*bad_sb25
        # OK, I know what you are thinking, I can't possibly be changing this again...

        # use SMA_SB25 instead of SB24 - this should work better for both high and low SB galaxies
        # if SMA_SB25 is not available use 1.35*SMA_SB24

        # for galaxies with SMA_SB24=0, set radius to value in main table 
        noradius_flag = self.radius_arcsec == 0
        self.radius_arcsec[noradius_flag] = self.vf.cat['radius'][noradius_flag]

        # also save BA and PA from John's catalog
        # use the self.radius_arcsec for the sma
        self.BA = np.ones(len(self.radius_arcsec))
        self.PA = np.zeros(len(self.radius_arcsec))
        
        self.BA[~noradius_flag] = ephot['BA_MOMENT'][~noradius_flag]
        self.PA[~noradius_flag] = ephot['PA_MOMENT'][~noradius_flag]
        
        self.RA = self.vf.cat['RA']
        self.DEC = self.vf.cat['DEC']        
        
        self.agcflag = False
        self.nsaflag = False

        
        #self.coadd.load_file(self.rcoadd_fname)
        #self.filter_ratio = 0.0422 #MKW8
        # filter ratio for ha4
        self.filter_ratio = 0.0426
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        

    def setup_uat_catalogs(self):
        """
        use the latest AGC as the source catalog
        """

        # on RF laptop, tabledir is /Users/rfinn/research/
        self.agc_fname = os.path.join(self.tabledir,'AGC/agc.allsky.210720.fits')
            

        self.agc = galaxy_catalog(self.agc_fname,virgo=False,agc=True)

        self.nsa_fname = os.path.join(self.tabledir,'NSA/nsa_v1_0_1.fits')
        if not os.path.exists(self.nsa_fname):
            self.nsa_fname = os.path.join(self.tabledir,'NSA/nsa_v0_1_2.fits')            
            
        self.nsa = galaxy_catalog(self.nsa_fname,virgo=False, agc=False, nsa=True)        
        #self.agc.check_ra_colname() # should check automatically with agc=True
        ##
        # get sizes for galaxies - will use this to unmask central region
        # need to cut this catalog based on keepflag
        ##

        #agc['a'] is the blue semi-major diameter in arcmin
        self.radius_arcsec = self.agc.cat['a']*60
        
        noradius_flag = self.radius_arcsec == 0
        self.radius_arcsec[noradius_flag] = 60 # set size of galaxies with no A value to 60 arcsec

        # also save BA and PA from John's catalog
        # use the self.radius_arcsec for the sma
        self.BA = np.ones(len(self.radius_arcsec))
        self.PA = np.zeros(len(self.radius_arcsec))
        
        self.BA[~noradius_flag] = self.agc.cat['b'][~noradius_flag]/self.agc.cat['a'][~noradius_flag]

        
        self.RA = self.agc.cat['RA']
        self.DEC = self.agc.cat['DEC']        
        
        self.agcflag = True
        self.nsaflag = True

        
        # filter ratio for ha4
        self.filter_ratio = 0.0426
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
		

        
class galaxy_catalog():
    agcflag: bool
    

    def __init__(self,catalog,nsa=False,agc=False,virgo=False,sizecat=None):
        self.cat = Table.read(catalog)
        #self.cat = Table(self.cat)
        self.catalog_name = catalog
        self.agcflag = agc
        self.nsaflag = nsa
        self.virgoflag = virgo
        if self.agcflag:
            self.check_ra_colname()
        if sizecat is not None:
            self.sizecat = sizecat
        else:
            self.sizecat = None
    def check_ra_colname(self):
        """
        GOAL:
        make sure the catalog has RA and DEC
        columns that are named RA and DEC
        
        this is set up to rename the AGC fields radeg/decdeg to the more standard RA/DEC

        PARAMS:
        * self

        METHOD:
        * will edit the column names of self.cat (the galaxy catalog)

        """
        try:
            t = self.cat['RA']
        except AttributeError:
            print('defining new catalog columns for RA/DEC')
            self.cat.rename_column('radeg','RA')
            self.cat.rename_column('decdeg','DEC')            

        except KeyError:
            print('defining new catalogs columns for RA/DEC')
            #print(self.cat.colnames)
            self.cat.rename_column('radeg','RA')
            self.cat.rename_column('decdeg','DEC')            
            
    def galaxies_in_fov(self,wcs,nrow=None,ncol=None,zmin=None,zmax=None,image_name = None,weight_image=None, agcflag=None,virgoflag=None):
        """
        GOAL: get galaxies in FOV

        PROCEDURE:
        * transforms catalog coords to image coords using wcs

        PARAMS:
        * wcs - of image
        
        OPTIONAL PARAMS:
        * nrow
        * ncol
        * zmin - apply redshift cut to galaxies, e.g. that fall within halpha filter window
        * zmax - apply redshift cut to galaxies, e.g. that fall within halpha filter window
        
        """
        #print('in galaxies in fov, nrow,ncol = ',nrow,ncol) # debug
        #print(f"self.nsa flag is {self.nsa}")

        if agcflag is None:
            agcflag = self.agcflag
        if virgoflag is None:
            virgoflag = self.virgoflag
        if (nrow is None) | (ncol is None):
            print('need image dimensions')
            return None
        else:
            print("")
            print(f"size of input image = ({nrow},{ncol})")

        ###########################################################################
        # use astropy.WCS.wcs.footprint_contains to get galaxies w/in FOV of image
        ###########################################################################        
        # 
        # this can replace the method below, where I transform all the coordinates
        # to pixels.  However, astropy now returns nans for objects that are far
        # from the field center, and this is causing errors downstream.
        # So footprint_contains should be more robust.
        #
        coords = SkyCoord(ra=self.cat['RA'],dec=self.cat['DEC'],unit='deg') 
        self.keepflag = wcs.footprint_contains(coords)
        print(f"number of galaxies based on keepflag  = {np.sum(self.keepflag)}")       


        ###########################################################################
        # old approach using wcs_world2pix and comparing with image dimensions
        ###########################################################################
        #px,py = wcs.wcs_world2pix(self.cat['RA'],self.cat['DEC'],0)

        # nanflag = np.isnan(px) | np.isnan(py)
        # print("number of nans in transformed coordinates = ",np.sum(nanflag))

        # print(f"min/max px of catalog galaxies = {np.min(px)} - {np.max(px)}")
        # print(f"min/max py of catalog galaxies = {np.min(py)} - {np.max(py)}")
        # print()
        # print("resulting coords from world2pix = ",px[0:10],py[0:10])
        # print("")
        # #print('in galaxies_in_fov: px={},py={}'.format(px,py))
        # colflag = (px < ncol) & (px >0)

        # rowflag = (py < nrow) & (py > 0)
        # print(f"number of galaxies within range of columns = {np.sum(colflag)}")
        # print(f"number of galaxies within range of rows = {np.sum(rowflag)}")        

        # print(f"number of galaxies within range of rows/columns  = {np.sum(rowflag & colflag)}")
        
        #keepflag=(px < ncol) & (px >0) & (py < nrow) & (py > 0)

        #keepflag = rowflag & colflag 
        

        # WCS returns nans for objects that are too far from the central coordinate
        # so these are also objects that will NOT be within the image FOV

        
        # replace keepflag of 
        #nanflag = np.isnan(self.cat['RA']) | np.isnan(self.cat['DEC'])
        #print("number of nans in RA/DEC coordinates = ",np.sum(nanflag))

        # check number of galaxies in fov
        if self.keepflag is None:
            print("WARNING: found no galaxies in FOV")
            return
        else:
            print(f"found {np.sum(self.keepflag)} after RA/DEC cuts")
            print()


                     

        ###########################################################################
        # check weight image to make sure galaxy is in good part of image
        ###########################################################################
            
        # should also check the weight image and remove galaxies with weight=0
        # this won't take care of images with partial exposures, but we can deal with that later...
        # TODO - how to handle images with partial exposures, meaning only part of galaxy is in FOV?
        
        
        imagename = image_name
        if imagename is not None:
            if imagename.find('shifted.fits') > -1:
                weightimage = imagename.replace('-r-shifted.fits','-r.weight-shifted.fits')
            else:
                weightimage = imagename.replace('.fits','.weight.fits')

            #if os.path.exists(weightimage):
            if 'MOS' not in imagename: # TODO: not sure why I am skipping MOS.  should just check to see if weightimage exists?
                if os.path.exists(weightimage):
                    px,py = wcs.wcs_world2pix(self.cat['RA'][self.keepflag],self.cat['DEC'][self.keepflag],0)
                    print()
                    print("cross checking object locations with weight image")
                    print()
                    whdu = fits.open(weightimage)
                    # just check center position?
                    int_px = np.array(px,'i')
                    int_py = np.array(py,'i')        
                    centerpixvals = whdu[0].data[int_py,int_px]
                    # weight image will have value > 0 if there is data there
                    weightflag = centerpixvals > 0
                    self.keepflag[self.keepflag] = self.keepflag[self.keepflag] & weightflag

        self.keepflag = self.apply_redshift_cut(zmin=zmin,zmax=zmax, agcflag=agcflag,virgoflag=virgoflag)
        
        return self.keepflag
    
    def apply_redshift_cut(self,zmin=None,zmax=None,image_name = None,weight_image=None, agcflag=False,virgoflag=False):
        ###########################################################################
        # get redshift cut
        ###########################################################################
        print(f"\nApplying redshift cut: zmin={zmin:.4f}, zmax={zmax:.4f}\n")
        # initialize value of zFlag
        zFlag = np.zeros(len(self.cat), 'bool')
        print(f"redshift of objects in FOV = ",self.cat['vopt'][self.keepflag].data/3.e5)
        #try: # should really edit the catalogs to have the same redshift/vel column name

        print(f"value of agcflag = {agcflag}")
        if agcflag:
            print("\t using the AGC velocities")
            zFlag1 = (self.cat['vopt']/3.e5 > zmin) & (self.cat['vopt']/3.e5 < zmax)
            zFlag2 = (self.cat['v21']/3.e5 > zmin) & (self.cat['v21']/3.e5 < zmax)
            zFlag = zFlag1 | zFlag2
            return (zFlag & self.keepflag)
        else:
            try:
                if self.virgoflag:
                    #print('virgo, right?')
                    print("\t using the Virgo velocities")                
                    zFlag = (self.cat['vr']/3.e5 > zmin) & (self.cat['vr']/3.e5 < zmax)
                elif self.nsaflag:
                    print("\t using the NSA velocities")                
                    zFlag = (self.cat.Z > zmin) & (self.cat.Z < zmax)
                print('number of galaxies on image, after z cut = ',np.sum(zFlag & self.keepflag))
                return (zFlag & self.keepflag)

            except AttributeError:
                print('AttributeError')
                print('make sure you selected the halpha filter')
                return self.keepflag


    def cull_catalog(self, keepflag,prefix):
        self.cat = self.cat[keepflag]
        if self.nsaflag:
            self.rmag = 22.5 - 2.5*np.log10(self.cat.NMGY[:,4])
        
        if self.nsaflag:
            outfile = prefix+'_nsa.fits'
            fits.writeto(outfile,self.cat, overwrite=True)
        elif self.agcflag:
            outfile = prefix+'_agc.fits'
            try:
                fits.writeto(outfile,self.cat, overwrite=True)
            except:
                self.cat.write(outfile, overwrite=True)
        elif self.virgoflag:
            #print('virgo, right???')
            outfile = prefix+'_virgo_cat.fits'
            #print('culled catalog = ',outfile)
            self.cat.write(outfile,format='fits',overwrite=True)
            # cull ephot
            
        
if __name__ == "__main__":
    ## RUNNING AS THE MAIN PROGRAM
    


    #####################################
    ## SETUP COMMAND-LINE PARAMETERS
    #####################################
    import argparse    
    parser = argparse.ArgumentParser(description ='Run gui for analyzing Halpha images')

    parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    
    parser.add_argument('--rimage',dest = 'rimage', default=None,help='r-band image')
    parser.add_argument('--haimage',dest = 'haimage', default=None,help='Halpha image')
    parser.add_argument('--csimage',dest = 'csimage', default=None,help='Continuum-subtracted Halpha image')    
    parser.add_argument('--filter',dest = 'filter', default=None,help='filter.  should be ha4, inthalpha, or intha6657.')
    parser.add_argument('--tabledir',dest = 'tabledir', default=None,help='table directory. something like /home/rfinn/research/Virgo/tables-north/v1/')
    parser.add_argument('--psfdir',dest = 'psfdir', default=None,help='set this to the directory containing PSF images')        
    parser.add_argument('--prefix',dest = 'prefix', default='v17p03',help='prefix associated with the coadded image.  Default is v17p03. required when running auto.')
    parser.add_argument('--auto',dest = 'auto', action='store_true',default=False,help='set this to process the images automatically, without the gui')
    
    parser.add_argument('--virgo',dest = 'virgo', action='store_true',default=False,help='set this if running on virgo data.  The virgo filaments catalog will be used as input.')
    parser.add_argument('--uat',dest = 'uat', action='store_true',default=False,help='set this if running on uat halpha groups.  The AGC (210720) will be used as the parent catalog.')     
    parser.add_argument('--draco',dest = 'draco', action='store_true',default=False,help='set this if running on draco.')   
    parser.add_argument('--nebula',dest = 'nebula', action='store_true',default=False,help='set this if running on open nebula virtual machine.  catalog paths will be set accordingly.')
    parser.add_argument('--laptop',dest = 'laptop', action='store_true',default=False,help="custom setting for running on Rose's laptop. catalog paths will be set accordingly.")
    
    parser.add_argument('--obsyear',dest = 'obsyear', default=None,help='year that data were taken.  this finds the right image directory if you are building the image name in pieces..  ')
    parser.add_argument('--pointing',dest = 'pointing', default=None,help='Pointing number that you want to load.  ONLY FOR VIRGO DATA, and only if you are buildling the image name in pieces.')
    
    parser.add_argument('--testing',dest = 'testing', action='store_true',default=False,help='set this if running on open nebula virtual machine')
    parser.add_argument('--verbose',dest = 'verbose', action='store_true',default=False,help='set this for extra print statements')    
        
    args = parser.parse_args()
    
    #catalog = os.getenv('HOME')+'/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    #print(sys.argv)
    logger = log.get_logger("example1", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    #custom_font = QtGui.QFont()
    #custom_font.setWeight(12)
    #QtWidgets.QApplication.setFont(custom_font,"QLabel")
    #QtWidgets.QApplication.setFont(custom_font,"QLabel")    
    MainWindow = QtWidgets.QMainWindow()
    #MainWindow.setStyleSheet("fusion")
    #font = MainWindow.font()
    #font.setPointSize(7)
    sepath = os.getenv('HOME')+'/github/halphagui/astromatic/'

    # OLD CODE TO MAKE USE OF ARGV INPUTS
    #if int(sys.argv[1]) == 0:
    #    ui = hafunctions(MainWindow, logger, sepath = sepath, testing=False)
    #elif int(sys.argv[1]) == 1:
    #    ui = hafunctions(MainWindow, logger, sepath = sepath, testing=True)
    #elif int(sys.argv[1]) == 2:
        # load default directories for virgo machine on open nebula
    #    ui = hafunctions(MainWindow, logger, sepath = sepath, testing=False,nebula=True)
    #ui.setupUi(MainWindow)
    #ui.test()

    #################################
    ## UPDATED TO USE ARGPARSE
    #################################
    if not args.auto:
        ui = hafunctions(MainWindow, logger, sepath = sepath, args=args)        
        #ui = hafunctions(MainWindow, logger, sepath = sepath, testing=args.testing,nebula=args.nebula,virgo=args.virgo,laptop=args.laptop,pointing=args.pointing,rimage=args.rimage,haimage=args.haimage,csimage=args.csimage,filter=args.filter,tabledir=args.tabledir,psfdir=args.psfdir)
        #if os.getenv("HOME").find('/home/') > -1: # special case for when running on linux laptop
        #    x = MainWindow.width()*.7
        #    y = MainWindow.height()*.7
        #    MainWindow.resize(x,y)
        MainWindow.show()
        sys.exit(app.exec_())
    else:
        # run functions non-interactively
        #ui = hafunctions(MainWindow, logger, sepath = sepath, testing=args.testing,nebula=args.nebula,virgo=args.virgo,laptop=args.laptop,pointing=args.pointing,auto=args.auto,prefix=args.prefix,obsyear=args.obsyear,rimage=args.rimage,haimage=args.haimage,csimage=args.csimage,filter=args.filter,tabledir=args.tabledir,psfdir=args.psfdir)
        ui = hafunctions(MainWindow, logger, sepath = sepath, args=args)
        pass
