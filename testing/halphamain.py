import sys, os
sys.path.append(os.getcwd())
sys.path.append('/Users/rfinn/github/HalphaImaging/')

from PyQt5 import  QtWidgets
#from PyQt5.Qtcore import  Qt
from ginga.qtw.QtHelp import QtGui, QtCore
from halphav2 import Ui_MainWindow
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
from matplotlib import pyplot as plt
import matplotlib.patches as patches

# code from HalphaImaging repository
import uat_sextractor_2image as runse
# filter information
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.,'INT197':6540.5}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.,'INT197':6615.5}



class image_panel(QtGui.QMainWindow):
    def __init__(self,panel_name,ui, logger):
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
        '''
        # Calculate WCS RA
        try:
            # NOTE: image function operates on DATA space coords
            image = viewer.get_image()
            if image is None:
                # No image loaded
                return
            ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,format='str', coords='fits')
        except Exception as e:
            self.logger.warning("Bad coordinate conversion: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'

        text = "RA: %s  DEC: %s  X: %.2f  Y: %.2f  Value: %s" % (ra_txt, dec_txt, fits_x, fits_y, value)
        '''
        # WCS stuff is not working so deleting for now...
        text = "X: %.2f  Y: %.2f  Value: %s" % (fits_x, fits_y, value)
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
class cutout_image():
    def __init__(self,panel_name,ui, logger, index):
        #super(image_panel, self).__init__()
        self.logger = logger
        self.ui = ui
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        #fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        #fi.set_callback('drag-drop', self.drop_file)
        fi.set_bg(0.2, 0.2, 0.2)
        #fi.ui_set_active(True)
        self.fitsimage = fi
        bd = fi.get_bindings()
        bd.enable_all(True)
        self.cutout = fi.get_widget()
        self.ui.cutoutsLayout.addWidget(self.cutout, 0, index, 1, 1)
        self.fitsimage = fi
    def load_image(self, imagearray):
        #self.fitsimage.set_image(imagearray)
        self.fitsimage.set_data(imagearray)


class hafunctions(Ui_MainWindow):

    def __init__(self,MainWindow, logger):
        super(hafunctions, self).__init__()
        print(MainWindow)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(MainWindow)


        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()
        self.add_coadd_frame(self.ui.leftLayout)
        self.add_cutout_frames()
        self.connect_setup_menu()
        self.connect_ha_menu()
        #self.connect_buttons()
        #self.add_image(self.ui.gridLayout_2)
        #self.add_image(self.ui.gridLayout_2)
        self.ui.wmark.clicked.connect(self.mark_galaxies)
        self.ui.maskButton.clicked.connect(self.edit_mask)
        self.ui.wfratio.clicked.connect(self.get_filter_ratio)
        self.ui.resetButton.clicked.connect(self.reset_cutout_values)

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
        self.rcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 0)
        self.hacutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 1)
        self.maskcutout = cutout_image(self.ui.cutoutsLayout,self.ui, self.logger, 2)

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
        print('working on this')
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

    def mark_galaxies(self):
            
        print('working on this...')
        # get list of NSA galaxies on image viewer
        #self.r, header_r = fits.getdata(self.rcoadd_fname,header=True)
        #self.ha, header_ha = fits.getdata(self.hacoadd_fname, header=True)
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

        self.gredshift = self.nsa.cat.Z[keepflag]
        gzdist = self.nsa.cat.ZDIST[keepflag]

        # populate a button that contains list
        print('nsa galaxies on fov:')
        print(self.galid)
        for name in self.galid:
            self.ui.wgalid.addItem(str(name))
        print(len(self.galid),' galaxies in FOV')
        self.ui.wgalid.activated.connect(self.select_galaxy)

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
        self.reset_cutout_values()
        self.display_cutouts()
        # first pass of mask
        # radial profiles
        # save latest of mask

    def display_cutouts(self):
        position = SkyCoord(ra=self.gra[self.igal],dec=self.gdec[self.igal],unit='deg')
        
        try:
            cutoutR = Cutout2D(self.r.data, position, self.cutout_size, wcs=self.coadd_wcs, mode='trim') #require entire image to be on parent image
            #cutoutHa = Cutout2D(self.ha.data, position, self.size, wcs=self.coadd_wcs, mode = 'trim')
            ((ymin,ymax),(xmin,xmax)) = cutoutR.bbox_original
            #print(ymin,ymax,xmin,xmax)
            self.rcutout.load_image(self.r[ymin:ymax,xmin:xmax])
            self.hacutout.load_image(self.halpha_cs[ymin:ymax,xmin:xmax])
            #cutoutR.plot_on_original(color='white')
        except nddata.utils.PartialOverlapError:# PartialOverlapError:
            print('galaxy is only partially covered by mosaic - skipping ',IDNUMBER[i])
            return
        except nddata.utils.NoOverlapError:# PartialOverlapError:
            print('galaxy is not covered by mosaic - skipping ',IDNUMBER[i])
            return

    def reset_cutout_values(self):
        self.cutout_size = self.reset_size
        self.filter_ratio = self.reset_ratio
        self.subtract_images()
        self.display_cutouts()
        self.ui.cutoutSlider.setValue(50)
        self.ui.ratioSlider.setValue(50)

    def edit_mask(self):
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
