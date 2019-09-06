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
from photutils import aperture_photometry

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

    def __init__(self, image, image2 = None, mask = None, image_frame=None, use_mpl=False, napertures=20):
        # use mpl for testing purposes, before integrating with pyqt gui
        # image2 is intended to be the Halpha image - use apertures from r-band but measure on Halpha
        self.image, self.header = fits.getdata(image, header=True)
        self.image_name = image
        # image 2 is designed to be the Halpha image, but it can be any second
        # image whereby you define the ellipse geometry using image 1, and
        # measure the photometry on image 1 and image 2
        #
        # self.image2_flag is True is image2 is provided
        if image2 != None:
            self.image2_name = image2
            self.image2,self.header2 = fits.getdata(image2, header=True)
            self.image2_flag = True
        else:
            self.image2_flag = False

        # the mask should identify all pixels in the cutout image that are not associated with the target galaxy
        # these will be ignored when defining the shape of the ellipse and when measuring the photometry
        #
        # self.mask_flag is True if a mask is provided
        if mask != None:
            self.mask_image, self.mask_header = fits.getdata(mask,header=True)
            self.mask_flag = True
            self.mask_image_bool = np.array(self.mask_image,'bool')
            self.masked_image = np.ma.array(self.image, mask = self.mask_image_bool) 
        else:
            self.mask_flag = False
            self.masked_image = self.image
        
        # image frame for plotting inside a gui
        # like if this is called from halphamain.py
        self.image_frame = image_frame

        # alternatively, for plotting with matplotlib
        # use this if running this code as the main program
        self.use_mpl = use_mpl
        self.napertures = napertures

    def run_for_gui(self):
        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_surface_brightness()
        self.write_phot_tables()
        if self.use_mpl:
            self.draw_phot_results_mpl()
        else:
            self.draw_phot_results()
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
        self.gini = obj.gini
        self.sky_centroid = obj.sky_centroid
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
        obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,self.sma, self.sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
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
                obj = self.image_frame.dc.Ellipse(iso.x0,iso.y0,iso.sma, iso.sma*(1-iso.eps), rot_deg = np.degrees(iso.pa), color=markcolor,linewidth=markwidth)
                objlist.append(obj)
            self.markhltag = self.image_frame.canvas.add(self.coadd.dc.CompoundObject(*objlist))
            self.image_frame.fitsimage.redraw()
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

    def measure_phot(self):
        # alternative is to use ellipse from detect
        # then create apertures and measure flux

        # rmax is max radius to measure ellipse
        # could cut this off based on SNR
        # or could cut this off based on enclosed flux?
        rmax = 1.5*self.sma
        self.apertures_a = np.linspace(2,rmax,20)
        self.apertures_b = (1.-self.eps)*self.apertures_a
        self.flux1 = np.zeros(len(self.apertures_a),'f')
        if self.image2_flag:
            self.flux2 = np.zeros(len(self.apertures_a),'f')
        for i in range(len(self.apertures_a)):
            ap = EllipticalAperture((self.xcenter, self.ycenter),self.apertures_a[i],self.apertures_b[i],self.theta)#,ai,bi,theta) for ai,bi in zip(a,b)]

            if self.mask_flag:
                self.phot_table1 = aperture_photometry(self.image, ap, mask=self.mask_image_bool)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, mask=self.mask_image_bool)
            else:
                # subpixel is the method used by Source Extractor
                self.phot_table1 = aperture_photometry(self.image, ap, method = 'subpixel', subpixels=5)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, method = 'subpixel', subpixels=5)
            self.flux1[i] = self.phot_table1['aperture_sum'][0]
            if self.image2_flag:
                self.flux2[i] = self.phot_table2['aperture_sum'][0]

    def calc_surface_brightness(self):
        # calculate surface brightness in each aperture

        area = np.pi*self.apertures_a*self.apertures_b # area of each ellipse

        # first aperture is calculated differently
        self.surface_brightness1 = np.zeros(len(self.apertures_a),'f')
        self.surface_brightness1[0] = self.flux1[0]/area[0]
        # outer apertures need flux from inner aperture subtracted
        for i in range(1,len(area)):
            self.surface_brightness1[i] = (self.flux1[i] - self.flux1[i-1])/(area[i]-area[i-1])

        # repeat for image 2 if it is provided
        if self.image2_flag:
            self.surface_brightness2 = np.zeros(len(self.apertures_a),'f')
            self.surface_brightness2[0] = self.flux2[0]/area[0]
            for i in range(1,len(area)):
                self.surface_brightness2[i] = (self.flux2[i] - self.flux2[i-1])/(area[i]-area[i-1])

    def write_phot_tables(self):
        # write out photometry for r-band
        # radius enclosed flux
        outfile = open(self.image_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

        outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
        outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
        outfile.write('# radius enclosed_flux \n')
        for i in range(len(self.apertures_a)):
            outfile.write('%.2f %.3e %.3e \n'%(self.apertures_a[i],self.flux1[i],self.surface_brightness1[i]))
        outfile.close()

        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            outfile = open(self.image2_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles
    
            outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
            outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
            outfile.write('# radius enclosed_flux \n')
            for i in range(len(self.apertures_a)):
                outfile.write('%.2f %.3e %.3e \n'%(self.apertures_a[i],self.flux2[i],self.surface_brightness2[i]))
            outfile.close()
    def draw_phot_results(self):
        ### DRAW RESULTING FIT ON R-BAND CUTOUT
        markcolor='cyan'
        objlist=[]
        markwidth=1.5
        for sma in self.apertures_a:
            obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,sma, sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
            objlist.append(obj)
            #print(self.xcenter,self.ycenter,sma, sma*(1-self.eps), self.theta, np.degrees(self.theta))
        self.markhltag = self.image_frame.canvas.add(self.image_frame.dc.CompoundObject(*objlist))
        self.image_frame.fitsimage.redraw()
        
    def draw_phot_results_mpl(self):
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm, origin='lower')

        apertures = []
        for sma in self.apertures_a:
            apertures.append(EllipticalAperture((self.xcenter,self.ycenter),sma, sma*(1-self.eps), theta = self.theta))
            
        for aperture in apertures:
            aperture.plot(color='white',lw=1.5)
        plt.show()
    def plot_profiles(self):
        plt.figure(figsize=(10,4))
        plt.subplots_adjust(wspace=.3)
        plt.subplot(1,2,1)
        plt.plot(self.apertures_a,self.flux1,'bo')
        plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Enclosed flux')
        if self.image2_flag:
            plt.subplot(1,2,2)
            plt.plot(self.apertures_a,self.flux2,'bo')
            plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Enclosed flux')
        plt.show()
        plt.savefig(self.image_name.split('.fits')[0]+'-enclosed-flux.png')
if __name__ == '__main__':
    image = 'MKW8-18216-R.fits'
    mask = 'MKW8-18216-R-mask.fits'
    #image = 'MKW8-18037-R.fits'
    #mask = 'MKW8-18037-R-mask.fits'
    #image = 'r-18045-R.fits'
    #mask = 'r-18045-R-mask.fits'
    e = ellipse(image,mask=mask, use_mpl=True)
    ## print('detect objects')
    ## e.detect_objects()
    ## print('find central')
    ## e.find_central_object()
    ## print('get guess')
    ## e.get_ellipse_guess()
    ## print('draw guess')
    ## e.draw_guess_ellipse_mpl()
    ## print('fit ellipse')
    ## e.fit_ellipse()
    ## print('plot results')
    ## e.draw_fit_results_mpl()

    print("--- %s seconds ---" % (time.time() - start_time))
