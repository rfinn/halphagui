import sys, os
sys.path.append(os.getcwd())
#sys.path.append('/Users/rfinn/github/HalphaImaging/')
sys.path.append('/Users/rfinn/github/HalphaImaging/python3/')

import numpy as np

from PyQt5 import  QtWidgets
from PyQt5 import QtCore
#from PyQt5.Qtcore import  Qt
from ginga.qtw.QtHelp import QtGui #, QtCore
from halphav3 import Ui_MainWindow
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.mplw.ImageViewCanvasMpl import ImageViewCanvas
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

# packages for ellipse fitting routine
# https://photutils.readthedocs.io/en/stable/isophote.html
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
from photutils import EllipticalAperture

from matplotlib import pyplot as plt
import matplotlib.patches as patches

from maskwrapper import maskwindow
from halphaCommon import cutout_image
# code from HalphaImaging repository
import uat_sextractor_2image as runse
#from uat_mask import mask_image
# filter information
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.,'INT197':6540.5}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.,'INT197':6615.5}



class image_panel(QtGui.QMainWindow):
    def __init__(self,panel_name,ui,logger):
        super(image_panel, self).__init__()
        self.ui = ui
        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()
        #self.figure = plt.figure()
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        fi.set_callback('drag-drop', self.drop_file)
        fi.set_callback('none-move',self.cursor_cb)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        #fi.set_figure(self.figure)
        self.fitsimage = fi

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
        
        wdrawcolor = QtGui.QComboBox()
        for name in self.drawcolors:
            wdrawcolor.addItem(name)
        index = self.drawcolors.index('lightblue')
        wdrawcolor.setCurrentIndex(index)
        wdrawcolor.activated.connect(self.set_drawparams)
        self.wdrawcolor = wdrawcolor
        
        wdrawtype = QtGui.QComboBox()
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
        self.image_wcs = WCS(filepath)
        self.fitsimage.set_image(image)
        #self.setWindowTitle(filepath)
        self.coadd_filename = filepath

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
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
        except:
            print('invalid value')
        # WCS stuff is not working so deleting for now...
        #text = "X: %.2f  Y: %.2f  Value: %s" % (fits_x, fits_y, value)
        self.readout.setText(text)

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


class hafunctions(Ui_MainWindow):

    def __init__(self,MainWindow, logger):
        super(hafunctions, self).__init__()
        #print(MainWindow)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(MainWindow)
        self.prefix=''

        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()
        self.add_coadd_frame(self.ui.leftLayout)
        self.add_cutout_frames()
        self.connect_setup_menu()
        self.connect_ha_menu()
        self.connect_halpha_type_menu()
        #self.connect_buttons()
        #self.add_image(self.ui.gridLayout_2)
        #self.add_image(self.ui.gridLayout_2)
        self.connect_buttons()
    def connect_buttons(self):
        self.ui.wmark.clicked.connect(self.find_galaxies)
        #self.ui.editMaskButton.clicked.connect(self.edit_mask)
        self.ui.makeMaskButton.clicked.connect(self.make_mask)
        self.ui.saveCutoutsButton.clicked.connect(self.write_cutouts)
        self.ui.profileButton.clicked.connect(self.plot_profiles)
        self.ui.wfratio.clicked.connect(self.get_filter_ratio)
        self.ui.resetRatioButton.clicked.connect(self.reset_cutout_ratio)
        self.ui.resetSizeButton.clicked.connect(self.reset_cutout_size)
        self.ui.prefixLineEdit.textChanged.connect(self.set_prefix)
        self.ui.fitEllipseButton.clicked.connect(self.fit_ellipse)
        self.setup_testing()
    def setup_testing(self):
        self.hacoadd_fname = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_ha16.coadd.fits'
        self.ha, self.ha_header = fits.getdata(self.hacoadd_fname, header=True)
        self.rcoadd_fname = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_R.coadd.fits'
        self.r, self.r_header = fits.getdata(self.rcoadd_fname, header=True)
        self.nsa_fname = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
        self.nsa = galaxy_catalog(self.nsa_fname)
        self.coadd.load_file(self.rcoadd_fname)
        self.filter_ratio = 0.0422
        self.reset_ratio = self.filter_ratio
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio
        self.subtract_images()
        self.setup_ratio_slider()
        self.cutout_size = 100
        self.setup_cutout_slider()

    def add_coadd_frame(self,panel_name):
        logger = log.get_logger("example1", log_stderr=True, level=40)
        self.coadd = image_panel(panel_name, self.ui,logger)
        #self.coadd.add_cutouts()

    def add_cutout_frames(self):
        # r-band cutout
        a = QtWidgets.QLabel('r-band')
        self.ui.cutoutsLayout.addWidget(a, 0, 0, 1, 1)
        a = QtWidgets.QLabel('CS Halpha')
        self.ui.cutoutsLayout.addWidget(a, 0, 1, 1, 1)
        a = QtWidgets.QLabel('Mask')
        self.ui.cutoutsLayout.addWidget(a, 0, 2, 1, 1)

        #self.ui.cutoutsLayout.addWidget(self.cutout, row, col, drow, dcol)
        self.rcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 0, 4, 1)
        self.hacutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1, 1, 4, 1)
        self.maskcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger,1, 2, 4, 1)

    def connect_setup_menu(self):
        self.ui.actionR_coadd.triggered.connect(self.get_rcoadd_file)
        self.ui.actionHa_coadd_2.triggered.connect(self.get_hacoadd_file)
        self.ui.actionNSA_catalog_path.triggered.connect(self.getnsafile)
    def get_rcoadd_file(self):
        fname = QtGui.QFileDialog.getOpenFileName()
        self.rcoadd_fname = fname[0]
        print(self.rcoadd_fname)
        #self.le.setPixmap(QPixmap(fname))
    def get_hacoadd_file(self):
        fname = QtGui.QFileDialog.getOpenFileName()
        self.hacoadd_fname = fname[0]
        print(self.hacoadd_fname)
        #self.le.setPixmap(QPixmap(fname))
    def getnsafile(self):
        fname = QtGui.QFileDialog.getOpenFileName()
        self.nsa_fname = fname[0]
        self.nsa = galaxy_catalog(self.nsa_fname)
        #self.le.setPixmap(QPixmap(fname))
		
    def connect_ha_menu(self):
        #print('working on this')
        #extractAction.triggered.connect(self.close_application)

        self.ui.actionhalpha4.triggered.connect(lambda: self.set_hafilter('4'))
        self.ui.actionhalpha8.triggered.connect(lambda: self.set_hafilter('8'))
        self.ui.actionhalpha12.triggered.connect(lambda: self.set_hafilter('12'))
        self.ui.actionhalpha16.triggered.connect(lambda: self.set_hafilter('16'))
    def set_hafilter(self,filterid):
        self.hafilter = filterid
        print('halpha filter = ',self.hafilter)
        self.get_zcut()
    def get_zcut(self):
        self.zmax=(((lmax[self.hafilter])/6563.)-1)
        self.zmin=(((lmin[self.hafilter])/6563.)-1)
    def connect_halpha_type_menu(self):
        ha_types = ['Ha Emission','No Ha','Cont Sub Problem']
        for name in ha_types:
            self.ui.haTypeComboBox.addItem(str(name))
        self.ui.haTypeComboBox.activated.connect(self.select_galaxy)
    def set_halpha_type(self,hatype):
        self.halpha_type = hatype

    def set_prefix(self,prefix):
        self.prefix = prefix
        #print('prefix for output files = ',self.prefix)
        
    def find_galaxies(self):
        #
        # get list of NSA galaxies on image viewer
        #
        # for reference:
        # self.r, header_r = fits.getdata(self.rcoadd_fname,header=True)
        # self.ha, header_ha = fits.getdata(self.hacoadd_fname, header=True)
        #
        n2,n1 = self.r.data.shape #should be same for Ha too, maybe? IDK
        n4,n3 = self.ha.data.shape 
    
        self.coadd_wcs= WCS(self.rcoadd_fname)#OF R IMAGE, SO THAT HA MATCHES WCS OF R, SO THEY'RE THE SAME
        px,py = self.coadd_wcs.wcs_world2pix(self.nsa.cat.RA,self.nsa.cat.DEC,1)
        onimageflag=(px < n1) & (px >0) & (py < n2) & (py > 0)
        try:
            zFlag = (self.nsa.cat.Z > self.zmin) & (self.nsa.cat.Z < self.zmax)
        except AttributeError:
            print('AttributeError')
            print('make sure you selected the halpha filter')
            return
        keepflag=zFlag & onimageflag
        self.gra=self.nsa.cat.RA[keepflag]
        self.gdec=self.nsa.cat.DEC[keepflag]
        self.gradius=self.nsa.cat.SERSIC_TH50[keepflag]
        self.galid=self.nsa.cat.NSAID[keepflag]
        self.gximage = px[keepflag]
        self.gyimage = py[keepflag]

        self.gredshift = self.nsa.cat.Z[keepflag]
        self.gzdist = self.nsa.cat.ZDIST[keepflag]

        # populate a button that contains list
        #print('nsa galaxies on fov:')
        #print(self.galid)
        for name in self.galid:
            self.ui.wgalid.addItem(str(name))
        print(len(self.galid),' galaxies in FOV')
        self.ui.wgalid.activated.connect(self.select_galaxy)

        # plot location of galaxies in the coadd image
        self.mark_galaxies()
        
    def mark_galaxies(self):
        #
        # using code in TVMark.py as a guide for adding shapes to canvas
        #
        #
        print('testing')
        objlist = []
        markcolor='cyan'
        markwidth=1
        for i,x in enumerate(self.gximage):
            obj = self.coadd.dc.Box(
                x=x, y=self.gyimage[i], xradius=8*self.gradius[i], yradius=8*self.gradius[i], color=markcolor,
                linewidth=markwidth)
            glabel = self.coadd.dc.Text(x-4*self.gradius[i],self.gyimage[i]+8.5*self.gradius[i],str(self.galid[i]), color=markcolor)
            objlist.append(obj)
            objlist.append(glabel)
        self.markhltag = self.coadd.canvas.add(self.coadd.dc.CompoundObject(*objlist))
        self.coadd.fitsimage.redraw()
    def get_filter_ratio(self):
        #
        # get ratio of Halpha to Rband filters
        # 
        # cannabalizing HalphaImaging/uat_find_filter_ratio.py
        #
        current_dir = os.getcwd()
        image_dir = os.path.dirname(self.rcoadd_fname)
        os.chdir(image_dir)
        runse.run_sextractor(self.rcoadd_fname, self.hacoadd_fname)
        ave, std = runse.make_plot(self.rcoadd_fname, self.hacoadd_fname, return_flag = True, image_dir = current_dir)
        print(ave,std)
        os.chdir(current_dir)
        self.filter_ratio = ave
        self.reset_ratio = ave
        self.minfilter_ratio = self.filter_ratio - 0.12*self.filter_ratio
        self.maxfilter_ratio = self.filter_ratio + 0.12*self.filter_ratio

        self.subtract_images()
        self.setup_ratio_slider()
        
    def subtract_images(self):
        self.halpha_cs = self.ha - self.filter_ratio*self.r
        # display continuum subtracted Halpha image in the large frame        
        self.coadd.fitsimage.set_data(self.halpha_cs)
    def setup_ratio_slider(self):
        self.ui.ratioSlider.setRange(0,100)
        self.ui.ratioSlider.setValue(50)
        self.ui.ratioSlider.setSingleStep(1)
        #self.ui.ratioSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self.ui.ratioSlider.setFocusPolicy(QtCore.StrongFocus)
        self.ui.ratioSlider.valueChanged.connect(self.ratio_slider_changed)
    def ratio_slider_changed(self, value):
        #print(self.minfilter_ratio, self.maxfilter_ratio, self.filter_ratio)
        delta = self.maxfilter_ratio - self.minfilter_ratio
        self.filter_ratio = self.minfilter_ratio + (delta)/100.*self.ui.ratioSlider.value()
        #print(value,' ratio slider changed to', round(self.filter_ratio,4))
        try:
            self.subtract_images()
            self.display_cutouts()
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
    def setup_cutout_slider(self):
        self.ui.cutoutSlider.setRange(0,100)
        self.ui.cutoutSlider.setValue(50)
        self.ui.cutoutSlider.setSingleStep(1)
        #self.ui.ratioSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self.ui.ratioSlider.setFocusPolicy(QtCore.StrongFocus)
        self.ui.cutoutSlider.valueChanged.connect(self.cutout_slider_changed)
    def cutout_slider_changed(self, value):
        #print(self.minfilter_ratio, self.maxfilter_ratio, self.filter_ratio)
        delta = self.maxcutout_size - self.mincutout_size
        self.cutout_size = self.mincutout_size + (delta)/100.*value
        #print(value,' ratio slider changed to', round(self.filter_ratio,4))
        try:
            self.display_cutouts()
        except:
            print('Trouble plotting cutouts')
            print('make sure galaxy is selected')
            
    def select_galaxy(self,id):
        print('selecting a galaxy')
        self.igal = self.ui.wgalid.currentIndex()
        print('active galaxy = ',self.igal)
        # when galaxy is selected from list, trigger
        # cutout imaages

        size = 16*self.gradius[self.igal]
        print('new cutout size = ',size, self.igal, self.gradius[self.igal])
        self.reset_size = size
        self.cutout_size = u.Quantity((size, size), u.arcsec)
        self.mincutout_size = 0.2*self.cutout_size
        self.maxcutout_size = 3.*self.cutout_size
        
        self.reset_cutout_size()
        self.reset_cutout_ratio()
        self.display_cutouts()
        # first pass of mask
        # radial profiles
        # save latest of mask
        try:
            self.cutout_name_r = self.prefix+'-'+str(self.galid[self.igal])+'-R.fits'
            self.cutout_name_ha = self.prefix+'-'+str(self.galid[self.igal])+'-CS.fits'
        except:
            self.cutout_name_r = str(self.galid[self.igal])+'-R.fits'
            self.cutout_name_ha =str(self.galid[self.igal])+'-CS.fits'
        self.rcutout.canvas.delete_all_objects()
    def display_cutouts(self):
        position = SkyCoord(ra=self.gra[self.igal],dec=self.gdec[self.igal],unit='deg')
        
        try:
            self.cutoutR = Cutout2D(self.r.data, position, self.cutout_size, wcs=self.coadd_wcs, mode='trim') #require entire image to be on parent image
            #cutoutHa = Cutout2D(self.ha.data, position, self.size, wcs=self.coadd_wcs, mode = 'trim')
            ((ymin,ymax),(xmin,xmax)) = self.cutoutR.bbox_original
            #print(ymin,ymax,xmin,xmax)
            self.rcutout.load_image(self.r[ymin:ymax,xmin:xmax])
            self.hacutout.load_image(self.halpha_cs[ymin:ymax,xmin:xmax])
            #cutoutR.plot_on_original(color='white')
        except nddata.utils.PartialOverlapError:# PartialOverlapError:
            print('galaxy is only partially covered by mosaic - skipping ',self.galid[self.igal])
            return
        except nddata.utils.NoOverlapError:# PartialOverlapError:
            print('galaxy is not covered by mosaic - skipping ',self.galid[self.igal])
            return

    def reset_cutout_size(self):
        self.cutout_size = self.reset_size
        self.update_images()
        self.ui.cutoutSlider.setValue(50)
    def reset_cutout_ratio(self):
        self.filter_ratio = self.reset_ratio
        self.update_images()
        self.ui.ratioSlider.setValue(50)
    def update_images(self):
        self.subtract_images()
        self.display_cutouts()

        
    def write_cutouts(self):
        #print(ymin,ymax,xmin,xmax)
        w = WCS(self.rcoadd_fname)
        try:
            ((ymin,ymax),(xmin,xmax)) = self.cutoutR.bbox_original
        except AttributeError:
            print('make sure you have selected a galaxy and saved the cutout')
            return
        newfile = fits.PrimaryHDU()
        newfile.data = self.r[ymin:ymax,xmin:xmax] 
        newfile.header = self.r_header
        newfile.header.update(w[ymin:ymax,xmin:xmax].to_header())
        newfile.header.set('REDSHIFT',float('{:.6f}'.format(self.gredshift[self.igal])))
        newfile.header.set('ZDIST',float('{:.6f}'.format(self.gzdist[self.igal])))
        newfile.header.set('NSAID',float('{:d}'.format(self.galid[self.igal])))
        newfile.header.set('SERSIC_TH50',float('{:.2f}'.format(self.gradius[self.igal])))
        fits.writeto(self.cutout_name_r, newfile.data, header = newfile.header, overwrite=True)

        # saving Ha Cutout as fits image
        newfile1 = fits.PrimaryHDU()
        newfile1.data = self.halpha_cs[ymin:ymax,xmin:xmax]
        newfile1.header = self.ha_header
        newfile1.header.update(w[ymin:ymax,xmin:xmax].to_header())
        newfile.header.set('REDSHIFT',float('{:.6f}'.format(self.gredshift[self.igal])))
        newfile.header.set('ZDIST',float('{:.6f}'.format(self.gzdist[self.igal])))
        newfile.header.set('NSAID',float('{:d}'.format(self.galid[self.igal])))
        newfile.header.set('SERSIC_TH50',float('{:.2f}'.format(self.gradius[self.igal])))
        fits.writeto(self.cutout_name_ha, newfile1.data, header = newfile1.header, overwrite=True)
        
    def make_mask(self):
        current_dir = os.getcwd()
        image_dir = os.path.dirname(self.rcoadd_fname)
        os.chdir(image_dir)
        try:
            self.write_cutouts()
        except AttributeError:
            print('are you rushing to make a mask w/out selecting galaxies?')
            print('try selecting filter, then selecting galaxies')
            return

        self.mwindow = QtWidgets.QWidget()
        self.mui = maskwindow(self.mwindow, self.logger, image = self.cutout_name_r, haimage=self.cutout_name_ha, sepath='~/github/HalphaImaging/astromatic/')
        self.mui.mask_saved.connect(self.display_mask)
        self.mui.setupUi(self.mwindow)
        self.mwindow.show()
        
        #print('make mask')
        #m.edit_mask()
        #m.clean_links()
        #self.display_mask()
        os.chdir(current_dir)
        

    def display_mask(self, mask_image_name):
        t = self.cutout_name_r.split('.fit')
        self.mask_image_name=t[0]+'-mask.fits'
        #self.mask_image = mask_image_name
        self.maskcutout.load_file(self.mask_image_name)

    def fit_ellipse(self):
        # https://github.com/astropy/photutils-datasets/blob/master/notebooks/isophote/isophote_example4.ipynb
        print(self.cutout_name_r)
        try:
            rdata, rheader = fits.getdata(self.cutout_name_r, header=True)
        except AttributeError:
            print('make sure you have selected a galaxy')
            return
        ### CLEAR R-BAND CUTOUT CANVAS

        self.rcutout.canvas.delete_all_objects()

        ### CREATE MASKED ARRAY IF MASK IS AVAILABLE
        #
        # following https://docs.scipy.org/doc/numpy/reference/maskedarray.generic.html
        # read in mask, convert to boolean
        # 
        # boolmask = np.array(mask,'bool')
        #
        # create new array:
        # newarray = np.ma.array(rdata, mask = np.array(self.mask_data, 'bool'))
        #
        # mask should have bad values = True
        #
        # NOTE! According to photutils documentation, using masked arrays
        # slows down the process considerably
        # https://github.com/astropy/photutils-datasets/blob/master/notebooks/isophote/isophote_example4.ipynb
        #
        # for details on ellipse fitting, including
        # holding EPS and PA fixed
        # step size in sma
        # https://photutils.readthedocs.io/en/stable/_modules/photutils/isophote/ellipse.html
        
        ### GENERATE GUESS ELLIPSE  ###
        # use center of cutout as first get for xcenter
        xcenter = int(rdata.shape[0]/2.)
        ycenter = int(rdata.shape[1]/2.)
        sma = 4.*self.gradius[self.igal]
        eps = 0.5
        pa = np.pi/2.
        guess = EllipseGeometry(x0=xcenter,y0=ycenter,sma=sma,eps = eps, pa = pa)
        aper = EllipticalAperture((guess.x0, guess.y0),guess.sma, guess.sma*(1 - guess.eps), guess.pa)

        ### DRAW ELLIPSE ON R-BAND CUTOUT
        #
        markcolor='cyan'
        markwidth=1

        ### FIT ELLIPSE
        #
        ellipse = Ellipse(rdata, guess)
        isolist = ellipse.fit_image(sma0 = 5, step=2, fix_pa = True, fix_eps = True)

        ### DRAW RESULTING FIT ON R-BAND CUTOUT
        iso = isolist.get_closest(5*self.gradius[self.igal])
        
        obj = self.coadd.dc.Ellipse(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps), rotdeg = np.degrees(iso.pa), color=markcolor,linewidth=markwidth)
        self.markhltag = self.rcutout.canvas.add(obj)
        self.rcutout.fitsimage.redraw()

        
        smas = np.linspace(10, np.max(isolist.sma), 5)
        objlist = []
        for sma in smas:
            iso = isolist.get_closest(sma)
            obj = self.coadd.dc.Ellipse(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps), rotdeg = np.degrees(iso.pa), color=markcolor,linewidth=markwidth)
            objlist.append(obj)
        self.markhltag = self.rcutout.canvas.add(self.coadd.dc.CompoundObject(*objlist))
        self.rcutout.fitsimage.redraw()
        
    def plot_profiles(self):
        print('edit mask')

class galaxy_catalog():

    def __init__(self,catalog):
        self.cat = fits.getdata(catalog)

if __name__ == "__main__":
    catalog = '/Users/rfinn/research/NSA/nsa_v0_1_2.fits'
    #gcat = galaxy_catalog(catalog)
    logger = log.get_logger("example1", log_stderr=True, level=40)
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = hafunctions(MainWindow, logger)
    #ui.setupUi(MainWindow)
    #ui.test()

    MainWindow.show()
    sys.exit(app.exec_())
