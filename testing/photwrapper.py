#!/usr/bin/env python

'''
This is going to be the wrapper to do photometry on the detected galaxies.

Would like to build this totally on photutils.

Useful references:


https://photutils.readthedocs.io/en/stable/segmentation.html
- detecting sources
- creating a segmentation image
- getting source properties (including total flux, Gini coefficient!!!)
- defining elliptical apertures for sources

'''

from photutils import detect_threshold, detect_sources
from photutils import source_properties
from photutils import Background2D, MedianBackground
from photutils import EllipticalAperture
from photutils.utils import calc_total_error
from photutils.isophote import EllipseGeometry, Ellipse

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import astropy.units as u
from astropy.io import fits
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize


from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

import numpy as np

import time
start_time = time.time()

# read in image and mask

# identify source for photometry

# run detect to detect source

# estimate ellipse parameters from source properties

# run 

class ellipse():

    def __init__(self, image, image2 = None, mask = None, image_frame=None, use_mpl=False):
        # use mpl for testing purposes, before integrating with pyqt gui
        # image2 is intended to be the Halpha image - use apertures from r-band but measure on Halpha
        self.image, self.header = fits.getdata(image, header=True)
        if mask != None:
            self.mimage, self.mheader = fits.getdata(mask,header=True)
            self.masked_image = np.ma.array(self.image, mask = np.array(self.mimage,'bool'))
        else:
            # for testing, skipping masking b/c it slows it down A LOT!!!
            self.masked_image = self.image
        self.image_frame = image_frame
        self.use_mpl = use_mpl
    def detect_objects(self, snrcut=2):
        self.threshold = detect_threshold(self.masked_image, nsigma=snrcut)
        self.segmentation = detect_sources(self.masked_image, self.threshold, npixels=10)
        self.cat = source_properties(self.masked_image, self.segmentation)
        #self.tbl = self.cat.to_table()
    def find_central_object(self):
        xdim,ydim = self.masked_image.shape
        distance = np.sqrt((self.cat.xcentroid.value - xdim/2.)**2 + (self.cat.ycentroid.value - ydim/2.)**2)
        # save object ID as the row in table with source that is closest to center
        self.objectID = np.arange(len(distance))[(distance == min(distance))][0]
        #print(self.objectID)
    def get_ellipse_guess(self, r=2.5):
        obj = self.cat[self.objectID]
        self.xcenter = obj.xcentroid.value
        self.ycenter = obj.ycentroid.value
        self.position = (obj.xcentroid.value, obj.ycentroid.value)

        self.sma = obj.semimajor_axis_sigma.value * r
        self.start_size = self.sma
        self.b = obj.semiminor_axis_sigma.value * r
        self.eps = 1 - self.b/self.sma
        t = obj.orientation.value
        if t < 0:
            self.theta = np.radians(180. - t)
        else:
            self.theta = obj.orientation.to(u.rad).value # orientation in radians
        # EllipticalAperture gives rotation angle in radians from +x axis, CCW
        self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)
    def draw_guess_ellipse(self):
        ### DRAW INITIAL ELLIPSE ON R-BAND CUTOUT
        #
        markcolor='magenta'
        markwidth=1
        obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,self.sma, self.sma*(1-self.eps), rotdeg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
        self.markhltag = self.image_frame.canvas.add(obj)
        self.image_frame.fitsimage.redraw()

    def draw_guess_ellipse_mpl(self):
        ### DRAW INITIAL ELLIPSE ON R-BAND CUTOUT
        #
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm , origin='lower')
        self.aperture.plot(color='white', lw=1.)
        plt.show()

    def fit_ellipse(self):
        ### FIT ELLIPSE
        #
        self.ellipse = Ellipse(self.masked_image, self.guess)
        self.isolist = self.ellipse.fit_image()#sfix_pa = True, step=.5)#, fix_eps=True, fix_center=True)
        self.table = self.isolist.to_table()
    def draw_fit_results(self):
        ### DRAW RESULTING FIT ON R-BAND CUTOUT
        markcolor='cyan'
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma), np.max(self.isolist.sma), 8)
            objlist = []
            for sma in smas:
                iso = self.isolist.get_closest(sma)
                obj = self.coadd.dc.Ellipse(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps), rotdeg = np.degrees(iso.pa), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            self.markhltag = self.image_panel.canvas.add(self.coadd.dc.CompoundObject(*objlist))
            self.image_panel.fitsimage.redraw()
        else:
            print('problem fitting ellipse')
    def draw_fit_results_mpl(self):
        
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm, origin='lower')
        apertures = []
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma)+2, np.max(self.isolist.sma), 12)
            objlist = []
            for sma in smas:
                iso = self.isolist.get_closest(sma)
                #print(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps),  np.degrees(iso.pa))
                apertures.append(EllipticalAperture((iso.x0,iso.y0),iso.sma, iso.sma*(1-iso.eps), theta = np.degrees(iso.pa)))
            for aperture in apertures:
                aperture.plot(color='white',lw=1.5)
        plt.show()



if __name__ == '__main__':
    image = 'MKW8-18216-R.fits'
    mask = 'MKW8-18216-R-mask.fits'
    image = 'MKW8-18037-R.fits'
    mask = 'MKW8-18037-R-mask.fits'
    image = 'r-18045-R.fits'
    mask = 'r-18045-R-mask.fits'
    e = ellipse(image,mask=mask, use_mpl=True)
    print('detect objects')
    e.detect_objects()
    print('find central')
    e.find_central_object()
    print('get guess')
    e.get_ellipse_guess()
    print('draw guess')
    e.draw_guess_ellipse_mpl()
    print('fit ellipse')
    e.fit_ellipse()
    print('plot results')
    e.draw_fit_results_mpl()

    print("--- %s seconds ---" % (time.time() - start_time))
