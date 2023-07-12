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
from photutils.morphology import gini

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import astropy.units as u
from astropy.io import fits
from astropy.table import Table, Column
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import scipy.ndimage as ndi
import statmorph
from statmorph.utils.image_diagnostics import make_figure


from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

import numpy as np

import time
start_time = time.time()

import matplotlib
matplotlib.use('Qt5Agg')


# modules in halphagui
import imutils

## filter information
## from https://www.noao.edu/kpno/mosaic/filters/
central_wavelength = {'4':6620.52,'8':6654.19,'12':6698.53,'16':6730.72,'R':6513.5,'r':6292.28,'inthalpha':6568.,'intha6657':6657,'intr':6240} # angstrom
dwavelength = {'4':80.48,'8':81.33,'12':82.95,'16':81.1,'R':1511.3,'r':1475.17,'inthalpha':95.,'intha6657':80,'intr':1347} # angstrom

# read in image and mask

# identify source for photometry

# run detect to detect source

# estimate ellipse parameters from source properties

# run 

class ellipse():
    '''
    class to run photometry routines on image

    INPUT
    * image         - primary image (this is usually rband)
    * image2        - image2 is designed to be the Halpha image, 
                      but it can be any second image whereby you define 
                      the ellipse geometry using image 1, and
                      measure the photometry on image 1 and image 2
    * mask          - mask to apply when measuring photometry.  
                      this is usually created from the r-band image
    * image_frame   - image frame for plotting inside a gui; like if this is called from halphamain
    * use_mpl       - use mpl for testing purposes, before integrating with pyqt gui
    * napertures    - number of apertures for measuring photmetry in. default is 20.
    * image2_filter - this is used to calculate the flux levels in the second filter.  
                      this should be one of standard halpha filters in our dictionary 
                      ('4','inthalpha','intha6657').
    * filter_ratio  - ratio of flux in image2/image1
    * psf           - used by statmorph (not sure if this is the image or the image name)
    * psf_ha        - used by statmorph

    '''
    def __init__(self, image, image2 = None, mask = None, image_frame=None, use_mpl=False, napertures=20,image2_filter=None, filter_ratio=None,psf=None,psf_ha=None):
        '''  inputs described above '''

        self.image, self.header = fits.getdata(image, header=True)
        self.image_name = image

        # get image dimensions - will use this to determine the max sma to measure
        self.ximage_max, self.yimage_max = self.image.shape
        
        # image 2 is designed to be the Halpha image, but it can be any second
        # image whereby you define the ellipse geometry using image 1, and
        # measure the photometry on image 1 and image 2
        #
        # self.image2_flag is True is image2 is provided
        if image2 is not None:
            self.image2_name = image2
            self.image2,self.header2 = fits.getdata(image2, header=True)
            self.image2_flag = True
        else:
            self.image2_flag = False
        self.image2_filter = image2_filter
        self.filter_ratio = filter_ratio
        # will use the gain to calculate the noise in the image
        try:
            self.gain = self.header['GAIN']
        except KeyError:
            print("WARNING: no GAIN keyword in header. Setting gain=1")
            self.gain = 1.
        self.psf = psf
        self.psf_ha = psf_ha
        # the mask should identify all pixels in the cutout image that are not
        # associated with the target galaxy
        # these will be ignored when defining the shape of the ellipse and when measuring the photometry
        #
        # self.mask_flag is True if a mask is provided
        if mask is not None:
            self.mask_image, self.mask_header = fits.getdata(mask,header=True)
            self.mask_flag = True
            self.boolmask = np.array(self.mask_image,'bool')
            self.masked_image = np.ma.array(self.image, mask = self.boolmask)
            if self.image2_flag:
                self.masked_image2 = np.ma.array(self.image2, mask = self.boolmask)
        else:
            print('not using a mask')
            self.mask_flag = False
            self.masked_image = self.image
            if self.image2_flag:
                self.masked_image2 = self.image2
        # image frame for plotting inside a gui
        # like if this is called from halphamain.py
        self.image_frame = image_frame

        # alternatively, for plotting with matplotlib
        # use this if running this code as the main program
        self.use_mpl = use_mpl
        self.napertures = napertures
        # assuming a typical fwhm 
        self.fwhm = 3.5
    def get_noise_in_aper(self, flux, area):
        ''' calculate the noise in an area '''
        noise_e = np.sqrt(flux*self.gain + area*self.sky_noise*self.gain)
        noise_adu = noise_e/self.gain
        return noise_adu

    def run_for_gui(self):
        ''' 
        batch all of the functions that we run for the gui, including:

        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_sb()
        self.convert_units()
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()
        '''
        
        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()
        self.measure_phot()
        self.calc_sb()
        self.convert_units()
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_tables()
        self.write_phot_fits_tables()
        self.get_sky_noise()

        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()
    def run_with_galfit_ellipse(self, xc,yc,BA=1,THETA=0):
        '''
        replicating run_for_gui(), but taking input ellipse geometry from galfit

        '''
        self.detect_objects()
        self.find_central_object()
        self.get_ellipse_guess()

        ### RESET ELLIPSE GEOMETRY USING GALFIT VALUES
        self.xcenter = float(xc)
        self.ycenter = float(yc)
        self.position = (self.xcenter, self.ycenter)
        # leave sma as it is defined in get_ellipse_guess
        # so that we measure photometry over the same approximate region
        # in practice this could be different areas if the ellipticity is very different
        # self.sma = r

        # need to reset b to be consistent with galfit ellipticity
        BA = float(BA)
        THETA = float(THETA)
        #print('THETA inside phot wrapper',THETA, BA)
        self.b = BA*self.sma
        self.eps = 1 - BA
        #print(self.b,self.eps,self.sma,BA)
        t = THETA
        if t < 0:
            self.theta = np.radians(180. + t)
        else:
            self.theta = np.radians(t) # orientation in radians
        # EllipticalAperture gives rotation angle in radians from +x axis, CCW
        self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)

        ### AND NOW BACK TO OUR REGULAR PROGRAMMING
        #print('measuring phot')
        self.measure_phot()
        #print('measuring phot')
        self.calc_sb()
        #print('measuring converting units')
        self.convert_units()
        #print('writing table')
        self.get_image2_gini()
        self.get_asymmetry()
        self.write_phot_fits_tables(prefix='GAL')
        #if self.use_mpl:
        #    self.draw_phot_results_mpl()
        #else:
        #    self.draw_phot_results()
    def detect_objects(self, snrcut=1.5,npixels=10):
        ''' 
        run photutils detect_sources to find objects in fov.  
        you can specify the snrcut, and only pixels above this value will be counted.
        
        this also measures the sky noise as the mean of the threshold image
        '''
        if self.mask_flag:
            self.threshold = detect_threshold(self.image, nsigma=snrcut,mask=self.boolmask)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels, mask=self.boolmask)
            self.cat = source_properties(self.image, self.segmentation, mask=self.boolmask)
        else:
            self.threshold = detect_threshold(self.image, nsigma=snrcut)
            self.segmentation = detect_sources(self.image, self.threshold, npixels=npixels)
            self.cat = source_properties(self.image, self.segmentation)
        # get average sky noise per pixel
        # threshold is the sky noise at the snrcut level, so need to divide by this
        self.sky_noise = np.mean(self.threshold)/snrcut
        #self.tbl = self.cat.to_table()
    def get_sky_noise(self):
        '''
        * get the noise in image1 and image2 
        * noise is stored as SKYERR in image header
          - units of sky noise are erg/s/cm^2/arcsec^2
        '''
        # get sky noise for image 1
        if self.mask_flag:
            threshold = detect_threshold(self.image, nsigma=1,mask=self.boolmask)
        else:
            threshold = detect_threshold(self.image, nsigma=snrcut)

        # add sky noise to image 1 header
        sky_noise_erg = np.mean(threshold)*self.uconversion1/self.pixel_scale**2
        print('r sky noise = ',sky_noise_erg)
        self.header.set('SKYERR',float('{:.2f}'.format(sky_noise_erg)),'sky noise in erg/s/cm^2/arcsec^2')
        # save files
        fits.writeto(self.image_name,self.image,header=self.header,overwrite=True)
        self.im1_skynoise = sky_noise_erg
        # get sky noise for image 2
        if self.image2 is not None:
            if self.mask_flag:
                threshold = detect_threshold(self.image2, nsigma=1,mask=self.boolmask)
            else:
                threshold = detect_threshold(self.image2, nsigma=snrcut)
            # add sky noise to image 2 header
            sky_noise_erg = np.mean(threshold)*self.uconversion2/self.pixel_scale**2        
            self.header2.set('SKYERR',float('{:.2f}'.format(sky_noise_erg)),'sky noise in erg/s/cm^2/arcsec^2')

            fits.writeto(self.image2_name,self.image2,header=self.header2,overwrite=True)
            self.im2_skynoise = sky_noise_erg
            print('ha sky noise = ',sky_noise_erg)

    def find_central_object(self):
        ''' 
        find the central object in the image and get its objid in segmentation image.
        object is stored as self.objectIndex
        '''
        xdim,ydim = self.image.shape
        distance = np.sqrt((self.cat.xcentroid.value - xdim/2.)**2 + (self.cat.ycentroid.value - ydim/2.)**2)
        # save object ID as the row in table with source that is closest to center
        self.objectIndex = np.arange(len(distance))[(distance == min(distance))][0]
        #print(self.objectIndex)
    def run_statmorph(self):
        '''
        run statmorph on image1 and image2 (if provided).

        results are stored as self.morph and self.morph2

        summary figures are save as XX statmorph-r.pdf and statmore-ha.pdf
        '''
        # show original
        #plt.figure()
        #plt.imshow(self.segmentation.data)
        # need segmentation map of object only
        segmap = self.segmentation.data == self.cat.id[self.objectIndex]
        segmap_float = ndi.uniform_filter(np.float64(segmap), size=10)
        segmap = segmap_float > 0.5
        #plt.figure()
        #plt.imshow(segmap, origin='lower', cmap='gray')
        if self.psf is None:
            source_morphs = statmorph.source_morphology(self.image, segmap, gain=self.gain)
        else:
            source_morphs = statmorph.source_morphology(self.image, segmap, gain=self.gain, psf=self.psf)
        self.morph = source_morphs[0]
        fig = make_figure(self.morph)
        figname = self.image_name.split('.fits')[0]
        fig.savefig(figname+'statmorph-r.pdf')
        if self.psf_ha is None:
            source_morphs2 = statmorph.source_morphology(self.image2, segmap, gain=self.gain)
        else:
            source_morphs2 = statmorph.source_morphology(self.image2, segmap, gain=self.gain, psf=self.psf_ha)
        self.morph2 = source_morphs2[0]
        fig2 = make_figure(self.morph2)
        fig2.savefig(figname+'statmorph-ha.pdf')
    def get_image2_gini(self, snrcut=1.5):
        ''' 
        calculate gini coefficient for image2 using pixels that are associated with r-band object ID

        this also calculates the sum and mag of the pixels associated with the central galaxy 
        (not sure why this is done together...)
        
        '''
        if self.mask_flag:
            self.threshold2 = detect_threshold(self.image2, nsigma=snrcut, mask=self.boolmask)
            self.segmentation2 = detect_sources(self.image2, self.threshold2, npixels=10,mask=self.boolmask)
            self.cat2 = source_properties(self.image2, self.segmentation2, mask=self.boolmask)
        else:
            self.threshold2 = detect_threshold(self.image2, nsigma=snrcut)
            self.segmentation2 = detect_sources(self.image2, self.threshold2, npixels=10)
            self.cat2 = source_properties(self.image2, self.segmentation2)

        '''
        select pixels associated with rband image in the segmentation
        AND
        pixels that are above the SNR cut in the Halpha image (image2)
        '''
        self.gini_pixels = (self.segmentation.data == self.cat.id[self.objectIndex]) & (self.segmentation2.data > 0.)

        #self.tbl = self.cat.to_table()
        self.gini2 = gini(self.image2[self.gini_pixels])
        self.source_sum2 = np.sum(self.image2[self.gini_pixels])
        self.source_sum2_erg = self.uconversion1*self.source_sum2
        self.source_sum2_mag = self.magzp2 - 2.5*np.log10(self.source_sum2)
    def get_asymmetry(self):
        '''
        * goal is to measure the assymetry of the galaxy about its center
        * going to measure asymmetry from pixels in the segmentation image only, so

        '''

        # for pixels in segmentation image of central object
        # (can't figure out a way to do this without looping
        # calculate delta_x and delta_y from centroid

        self.object_pixels = self.segmentation.data == self.cat.id[self.objectIndex]

        xc = self.cat.xcentroid[self.objectIndex].value
        yc = self.cat.ycentroid[self.objectIndex].value
        row,col = np.where(self.object_pixels)

        grid_size = 3
        sum_diff = np.zeros((grid_size,grid_size),'f')
        source_sum = np.zeros((grid_size,grid_size),'f')
        for dxc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
            for dyc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                drow = np.array((row-(yc+dyc)),'i')
                dcol = np.array((col-(xc+dxc)),'i')
                row2 = np.array(((yc+dyc) -1*drow),'i')
                col2 = np.array(((xc+dxc) -1*dcol),'i')
                sum_diff[dyc,dxc] = np.sum(np.abs(self.masked_image[row,col] - self.masked_image[row2,col2]))
                # divide by the sum of the original pixel values for object
                source_sum[dyc,dxc] = np.sum(self.image[self.object_pixels])
        asym = sum_diff/source_sum
        #print(asym)
        self.asym = np.min(asym)
        self.asym_err = np.std(asym)
        r,c = np.where(asym == np.min(asym))
        self.asym_center = np.array([r+yc,c+xc])

        print('asymmetry = {:.3f}+/-{:.3f}'.format(self.asym,self.asym_err))
        
        if self.image2_flag:
            self.object_pixels2 = (self.segmentation.data == self.cat.id[self.objectIndex]) & (self.segmentation2.data > 0.)

            xc = self.cat.xcentroid[self.objectIndex].value
            yc = self.cat.ycentroid[self.objectIndex].value
            row,col = np.where(self.object_pixels2)
            sum_diff = np.zeros((grid_size,grid_size),'f')
            source_sum = np.zeros((grid_size,grid_size),'f')


            # looks like I am measuring asymmetry myself here?
            # doesn't statmorph do this? - looks like I am not using it
            for dxc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                for dyc in np.arange(int(-1*(grid_size/2)),int(grid_size/2)+1):
                    drow = np.array((row-(yc+dyc)),'i')
                    dcol = np.array((col-(xc+dxc)),'i')
                    row2 = np.array(((yc+dyc) -1*drow),'i')
                    col2 = np.array(((xc+dxc) -1*dcol),'i')
                    sum_diff[dyc,dxc] = np.sum(np.abs(self.masked_image2[row,col] - self.masked_image2[row2,col2]))
                    # divide by the sum of the original pixel values for object
                    source_sum[dyc,dxc] = np.sum(self.image2[self.object_pixels2])
            asym2 = sum_diff/source_sum
            
            #print('asym2 = ',asym2,r,c)
            # measure halpha asymmetry at pixel where R-band asymmetry is minimum
            try:
                self.asym2 = asym2[r,c][0]
                #print(self.asym2)
            except IndexError:
                try:
                    r,c = np.where(asym == np.min(asym2))
                
                    self.asym2 = asym2[r,c][0]
                except IndexError:
                    self.asym2 = np.nan
                    self.asym2_err = np.nan
                    self.asym2_center = np.nan
                    return
            self.asym2_err = np.std(asym2)
            r,c = np.where(asym == np.min(asym2))
            self.asym2_center = np.array([r+yc,c+xc])
            #print('asymmetry2 = ',self.asym2)
            print('asymmetry = {:.3f}+/-{:.3f}'.format(self.asym2,self.asym2_err))
            '''
            # use all the same images as for r-band measurement
            self.object_pixels2 = (self.segmentation.data == self.cat.id[self.objectIndex])# & (self.segmentation2.data > 0.)

            xc = self.cat.xcentroid[self.objectIndex].value
            yc = self.cat.ycentroid[self.objectIndex].value
            row,col = np.where(self.object_pixels2)

            drow = np.array((row-yc),'i')
            dcol = np.array((col-xc),'i')
            row2 = np.array((yc -1*drow),'i')
            col2 = np.array((xc -1*dcol),'i')
            sum_diff = np.sum(np.abs(self.masked_image2[row,col] - self.masked_image2[row2,col2]))
            # divide by the sum of the original pixel values for object
            source_sum = np.sum(self.image2[self.object_pixels2])
        
        
            self.asym2b = sum_diff/source_sum
            print('asymmetry2 = ',self.asym2b)
            '''
        
    def get_ellipse_guess(self, r=2.5):
        '''
        this gets the guess for the ellipse geometry from the detection catalog 
        '''
        obj = self.cat[self.objectIndex]
        self.xcenter = obj.xcentroid.value
        self.ycenter = obj.ycentroid.value
        self.position = (obj.xcentroid.value, obj.ycentroid.value)

        self.sma = obj.semimajor_axis_sigma.value * r
        self.start_size = self.sma
        self.b = obj.semiminor_axis_sigma.value * r
        self.eps = 1 - self.b/self.sma
        self.gini = obj.gini
        self.source_sum = self.cat[self.objectIndex].source_sum
        self.sky_centroid = obj.sky_centroid
        # orientation is angle in radians, CCW relative to +x axis
        t = obj.orientation.value
        #print('inside get_ellipse_guess, orientation = ',obj.orientation)
        if t < 0: # convert to positive angle wrt +x axis
            self.theta = np.pi+obj.orientation.to(u.rad).value
        else:
            self.theta = obj.orientation.to(u.rad).value # orientation in radians
        # EllipticalAperture gives rotation angle in radians from +x axis, CCW
        self.aperture = EllipticalAperture(self.position, self.sma, self.b, theta=self.theta)
        # EllipseGeometry using angle in radians, CCW from +x axis
        self.guess = EllipseGeometry(x0=self.xcenter,y0=self.ycenter,sma=self.sma,eps = self.eps, pa = self.theta)
    def draw_guess_ellipse(self):
        ''' DRAW INITIAL ELLIPSE ON R-BAND CUTOUT '''
        #
        markcolor='magenta'
        markwidth=1
        obj = self.image_frame.dc.Ellipse(self.xcenter,self.ycenter,self.sma, self.sma*(1-self.eps), rot_deg = np.degrees(self.theta), color=markcolor,linewidth=markwidth)
        self.markhltag = self.image_frame.canvas.add(obj)
        self.image_frame.fitsimage.redraw()

    def draw_guess_ellipse_mpl(self):
        ''' DRAW INITIAL ELLIPSE ON R-BAND CUTOUT '''
        #
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        plt.imshow(self.masked_image, cmap='Greys_r', norm=norm , origin='lower')
        self.aperture.plot(color='white', lw=1.)
        plt.show()

    def fit_ellipse(self):
        ''' FIT ELLIPSE '''
        #
        # create instance of photutils.Ellipse
        # https://photutils.readthedocs.io/en/stable/isophote.html
        self.ellipse = Ellipse(self.masked_image, self.guess)
        self.isolist = self.ellipse.fit_image()#sfix_pa = True, step=.5)#, fix_eps=True, fix_center=True)
        self.table = self.isolist.to_table()
        
    def draw_fit_results(self):
        ''' DRAW RESULTING FIT ON R-BAND CUTOUT '''
        markcolor='cyan'
        if len(self.isolist) > 5:
            smas = np.linspace(np.min(self.isolist.sma), np.max(self.isolist.sma), 3)
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
        ''' draw fit results in matplotlib figure '''
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
    def show_seg_aperture(self):
        ''' matplotlib plotting to show apertures   '''
        tbl1 = self.cat.to_table()
        cat = self.cat
        r=3.
        apertures = []
        for obj in cat:
            position = np.transpose((obj.xcentroid.value, obj.ycentroid.value))
            a = obj.semimajor_axis_sigma.value * r
            b = obj.semiminor_axis_sigma.value * r
            theta = obj.orientation.to(u.rad).value
            #print(theta)
            apertures.append(EllipticalAperture(position, a, b, theta=theta))
    
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
        ax1.imshow(self.image, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        #cmap = segm_deblend.make_cmap(random_state=12345)
        ax2.imshow(self.segmentation.data, origin='lower')
        ax2.set_title('Segmentation Image')
        for aperture in apertures:
            aperture.plot(axes=ax1, color='white', lw=1.5)
            aperture.plot(axes=ax2, color='white', lw=1.5)
        plt.show()
    def measure_phot(self):
        '''
        # alternative is to use ellipse from detect
        # then create apertures and measure flux

        # rmax is max radius to measure ellipse
        # could cut this off based on SNR
        # or could cut this off based on enclosed flux?
        # or could cut off based on image dimension, and do the cutting afterward
        
        #rmax = 2.5*self.sma
        '''
        
        '''
        this is how becky set the apertures
        a = [0]
        for i in range(1,500):
        a.append(a[i-1] + hwhm + (hwhm*i*.1))
        
        '''
        # rmax is set according to the image dimensions
        # look for where the semi-major axis hits the edge of the image
        # could by on side (limited by x range) or on top/bottom (limited by y range)
        # 
        #print('xcenter, ycenter, theta = ',self.xcenter, self.ycenter,self.theta)
        rmax = np.min([(self.ximage_max - self.xcenter)/abs(np.cos(self.theta)),\
                       (self.yimage_max - self.ycenter)/abs(np.sin(self.theta))])
        #print('print rmax, ximage_max, image.shape = ',rmax,self.ximage_max,self.image.shape)
        '''
        this is how becky set the apertures
        a = [0]
        for i in range(1,500):
        a.append(a[i-1] + hwhm + (hwhm*i*.1))
        
        '''

        index = np.arange(80)
        apertures = (index+1)*.5*self.fwhm*(1+(index+1)*.1)
        # cut off apertures at edge of image
        self.apertures_a = apertures[apertures < rmax]
        #print('number of apertures = ',len(self.apertures_a))
        #self.apertures_a = np.linspace(3,rmax,40)
        self.apertures_b = (1.-self.eps)*self.apertures_a
        self.area = np.pi*self.apertures_a*self.apertures_b # area of each ellipse


        self.flux1 = np.zeros(len(self.apertures_a),'f')
        self.flux1_err = np.zeros(len(self.apertures_a),'f')
        if self.image2_flag:
            self.flux2 = np.zeros(len(self.apertures_a),'f')
            self.flux2_err = np.zeros(len(self.apertures_a),'f')
        self.allellipses = []
        for i in range(len(self.apertures_a)):
            #print('measure phot: ',self.xcenter, self.ycenter,self.apertures_a[i],self.apertures_b[i],self.theta)
            #,ai,bi,theta) for ai,bi in zip(a,b)]
            # EllipticalAperture takes rotation angle in radians, CCW from +x axis
            ap = EllipticalAperture((self.xcenter, self.ycenter),self.apertures_a[i],self.apertures_b[i],self.theta)#,ai,bi,theta) for ai,bi in zip(a,b)]
            self.allellipses.append(ap)

            if self.mask_flag:
                self.phot_table1 = aperture_photometry(self.image, ap, mask=self.boolmask)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, mask=self.boolmask)
            else:
                # subpixel is the method used by Source Extractor
                self.phot_table1 = aperture_photometry(self.image, ap, method = 'subpixel', subpixels=5)
                if self.image2_flag:
                    self.phot_table2 = aperture_photometry(self.image2, ap, method = 'subpixel', subpixels=5)
            self.flux1[i] = self.phot_table1['aperture_sum'][0]
            
            # calculate noise
            self.flux1_err[i] = self.get_noise_in_aper(self.flux1[i], self.area[i])
            if self.image2_flag:
                self.flux2[i] = self.phot_table2['aperture_sum'][0]
                self.flux2_err[i] = self.get_noise_in_aper(self.flux2[i], self.area[i])
    def draw_phot_apertures(self):
        ''' matplotlib plotting to show apertures   '''
        tbl1 = self.cat.to_table()
        cat = self.cat
        r=3.
        apertures = []
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
        ax1.imshow(self.image, origin='lower', cmap='Greys_r', norm=norm)
        ax1.set_title('Data')
        #cmap = segm_deblend.make_cmap(random_state=12345)
        ax2.imshow(self.segmentation.data, origin='lower')
        ax2.set_title('Segmentation Image')
        for aperture in self.allellipses:
            aperture.plot(axes=ax1, color='white', lw=1.5)
            aperture.plot(axes=ax2, color='white', lw=1.5)
        plt.show()

    def calc_sb(self):
        # calculate surface brightness in each aperture

        # first aperture is calculated differently
        self.sb1 = np.zeros(len(self.apertures_a),'f')
        self.sb1_err = np.zeros(len(self.apertures_a),'f')

        self.sb1[0] = self.flux1[0]/self.area[0]
        self.sb1_err[0] = self.get_noise_in_aper(self.flux1[0], self.area[0])/self.area[0]
        # outer apertures need flux from inner aperture subtracted
        for i in range(1,len(self.area)):
            self.sb1[i] = (self.flux1[i] - self.flux1[i-1])/(self.area[i]-self.area[i-1])
            self.sb1_err[i] = self.get_noise_in_aper((self.flux1[i] - self.flux1[i-1]),(self.area[i]-self.area[i-1]))/(self.area[i]-self.area[i-1])

        # calculate SNR to follow Becky's method of cutting off analysis where SNR = 2
        self.sb1_snr = np.abs(self.sb1/self.sb1_err)
        # repeat for image 2 if it is provided
        if self.image2_flag:
            self.sb2 = np.zeros(len(self.apertures_a),'f')
            self.sb2_err = np.zeros(len(self.apertures_a),'f')
            self.sb2[0] = self.flux2[0]/self.area[0]
            self.sb2_err[0] = self.get_noise_in_aper(self.flux2[0], self.area[0])/self.area[0]
            for i in range(1,len(self.area)):
                self.sb2[i] = (self.flux2[i] - self.flux2[i-1])/(self.area[i]-self.area[i-1])
                self.sb2_err[i] = self.get_noise_in_aper((self.flux2[i] - self.flux2[i-1]),(self.area[i]-self.area[i-1]))/(self.area[i]-self.area[i-1])
            self.sb2_snr = np.abs(self.sb2/self.sb2_err)

    def convert_units(self):
        '''
        ###########################################################
        ### SET UP INITIAL PARAMETERS TO CALCULATE CONVERSION
        ### FROM ADU/S TO PHYSICAL UNITS
        ###########################################################
        '''

        self.pixel_scale = imutils.get_pixel_scale(self.header)
        try:
            self.magzp = float(self.header['PHOTZP'])

        except:
            print("WARNING: no PHOTZP keyword in image header. \nAssuming ZP=22.5")
            self.magzp = 22.5
        #print('mag zp = ',self.magzp)
        filter = self.header
        # multiply by bandwidth of filter to convert from Jy to erg/s/cm^2
        bandwidth1 = 3.e8*dwavelength['R']*1.e-10/(central_wavelength['R']*1.e-10)**2
        # need to figure out how to adjust automatically
        bandwidth1 = 3.e8*dwavelength['r']*1.e-10/(central_wavelength['r']*1.e-10)**2        
        self.uconversion1 = 3631.*10**(self.magzp/-2.5)*1.e-23*bandwidth1
        if self.image2_filter:
            bandwidth2 = 3.e8*dwavelength[self.image2_filter]*1.e-10/(central_wavelength[self.image2_filter]*1.e-10)**2
            try:
                self.magzp2 = float(self.header2['PHOTZP'])
                self.uconversion2 = 3631.*10**(self.magzp2/-2.5)*1.e-23*bandwidth2
            except:
                # use 25 as default ZP if none is provided in header
                self.uconversion2 = 3631.*10**(25/-2.5)*1.e-23*bandwidth2
        if self.filter_ratio is not None:
            if self.image2_flag:
                self.uconversion2b = self.filter_ratio*self.uconversion1
        else:
            self.uconversion2b = None
            
        ###########################################################
        ### CONVERT UNITS TO
        ### FLUX -> ERG/S/CM^2
        ### FLUX -> MAG
        ### SURFACE BRIGHTNESS -> ERG/S/CM^2/ARCSEC^2
        ### SURFACE BRIGHTNESS -> MAG/ARCSEC^2
        ###########################################################
        self.sky_noise_erg = self.sky_noise*self.uconversion1/self.pixel_scale**2
        self.flux1_erg = self.uconversion1*self.flux1
        self.flux1_err_erg = self.uconversion1*self.flux1_err
        self.source_sum_erg = self.uconversion1*self.source_sum
        self.source_sum_mag = self.magzp - 2.5*np.log10(self.source_sum)
        self.mag1 = self.magzp - 2.5*np.log10(self.flux1)
        self.mag1_err = self.mag1 - (self.magzp - 2.5*np.log10(self.flux1 + self.flux1_err))
        self.sb1_erg_sqarcsec = self.uconversion1*self.sb1/self.pixel_scale**2
        self.sb1_erg_sqarcsec_err = self.uconversion1*self.sb1_err/self.pixel_scale**2
        self.sb1_mag_sqarcsec = self.magzp - 2.5*np.log10(self.sb1/self.pixel_scale**2)
        self.sb1_mag_sqarcsec_err = self.sb1_mag_sqarcsec - (self.magzp - 2.5*np.log10((self.sb1 + self.sb1_err)/self.pixel_scale**2))
        if self.image2_flag:
            self.flux2_erg = self.uconversion2*self.flux2
            self.flux2_err_erg = self.uconversion2*self.flux2_err

            self.mag2 = self.magzp2 - 2.5*np.log10(self.flux2)
            self.mag2_err = self.mag2 - (self.magzp2 - 2.5*np.log10(self.flux2 + self.flux2_err))
            self.sb2_erg_sqarcsec = self.uconversion2*self.sb2/self.pixel_scale**2
            self.sb2_erg_sqarcsec_err = self.uconversion2*self.sb2_err/self.pixel_scale**2
            self.sb2_mag_sqarcsec = self.magzp2 - 2.5*np.log10(self.sb2/self.pixel_scale**2)
            # add error to flux and calculate magnitude, then take difference with original mag
            self.sb2_mag_sqarcsec_err = self.sb2_mag_sqarcsec - (self.magzp2 - 2.5*np.log10((self.sb2+self.sb2_err)/self.pixel_scale**2))
            # this next set uses the filter ratio and the filter 1 flux conversion to
            # convert narrow-band flux (filter 2) to physical units.
            if self.uconversion2b:
                conversion = self.uconversion2b
                self.flux2b_erg = conversion*self.flux2
                self.flux2b_err_erg = conversion*self.flux2_err
                self.sb2b_erg_sqarcsec = conversion*self.sb2/self.pixel_scale**2
                self.sb2b_erg_sqarcsec_err = conversion*self.sb2_err/self.pixel_scale**2
                self.sb2b_mag_sqarcsec = self.magzp2 - 2.5*np.log10(conversion*self.sb2/self.pixel_scale**2)
                self.sb2b_mag_sqarcsec_err = self.sb2b_mag_sqarcsec - (self.magzp2 - 2.5*np.log10(conversion*(self.sb2+self.sb2_err)/self.pixel_scale**2))
                
    def write_phot_tables(self):
        '''
        write out photometry for image and image2 in ascii format
        '''
        
        # radius enclosed flux
        outfile = open(self.image_name.split('.fits')[0]+'-phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

        #outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
        #outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
        outfile.write('# radius flux flux_err sb sb_err sb_snr flux_erg flux_erg_err mag mag_err sb_ergsqarc sb_err_ergsqarc sb_magsqarc sb_err_magsqarc \n')
        for i in range(len(self.apertures_a)):
            s='%.2f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e '% \
                          (self.apertures_a[i],self.flux1[i],self.flux1_err[i],\
                           self.sb1[i], self.sb1_err[i], \
                          self.sb1_snr[i], \
                          self.flux1_erg[i], self.flux1_err_erg[i],\
                          self.mag1[i], self.mag1_err[i], \
                          self.sb1_erg_sqarcsec[i],self.sb1_erg_sqarcsec[i], \
                          self.sb1_mag_sqarcsec[i],self.sb1_mag_sqarcsec[i])
            s=s+'\n'
            outfile.write(s)
        outfile.close()

        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            outfile = open(self.image2_name.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles
    
            #outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
            #outfile.write('# %.2f %.2f %.2f %.2f \n'%(self.xcenter,self.ycenter,self.eps,self.theta))
            s = '# radius flux flux_err sb sb_err sb_snr flux_erg flux_erg_err mag mag_err sb_ergsqarc sb_err_ergsqarc sb_magsqarc sb_err_magsqarc'
            if self.uconversion2b:
                s = s +' fluxb_erg fluxb_erg_err sbb_ergsqarc sbb_err_ergsqarc sbb_magsqarc sbb_err_magsqarc'
            s = s + '\n'
            outfile.write(s)
            for i in range(len(self.apertures_a)):
                s = '%.2f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e '% \
                              (self.apertures_a[i],self.flux2[i],self.flux2_err[i],\
                               self.sb2[i], self.sb2_err[i],self.sb2_snr[i],\
                              self.flux2_erg[i], self.flux2_err_erg[i],\
                              self.mag2[i], self.mag2_err[i], \
                              self.sb2_erg_sqarcsec[i],self.sb2_erg_sqarcsec[i], \
                              self.sb2_mag_sqarcsec[i],self.sb2_mag_sqarcsec[i])
                if self.uconversion2b:
                    s=s+' %.3e %.3e %.3e %.3e %.3e %.3e'% \
                      (self.flux2b_erg[i], self.flux2b_err_erg[i],\
                      self.sb2b_erg_sqarcsec[i],self.sb2b_erg_sqarcsec[i], \
                      self.sb2b_mag_sqarcsec[i],self.sb2b_mag_sqarcsec[i])
                s = s+'\n'
                outfile.write(s)

            outfile.close()
    def write_phot_fits_tables(self, prefix=None):
        ''' write out photometry for image and image2 in fits format '''

        if prefix is None:
             outfile = self.image_name.split('.fits')[0]+'-phot.fits'
        else:
             outfile = self.image_name.split('.fits')[0]+'-'+prefix+'-phot.fits'
        print('photometry outfile = ',outfile)

        data = [self.apertures_a*self.pixel_scale,self.apertures_a, \
             self.flux1,self.flux1_err,\
             self.sb1, self.sb1_err, \
             self.sb1_snr, \
             self.flux1_erg, self.flux1_err_erg,\
             self.mag1, self.mag1_err, \
             self.sb1_erg_sqarcsec,self.sb1_erg_sqarcsec_err, \
             self.sb1_mag_sqarcsec,self.sb1_mag_sqarcsec_err]

        names = ['sma_arcsec','sma_pix','flux','flux_err',\
                 'sb', 'sb_err', \
                 'sb_snr', \
                 'flux_erg', 'flux_erg_err',\
                 'mag', 'mag_err', \
                 'sb_erg_sqarcsec','sb_erg_sqarcsec_err', \
                 'sb_mag_sqarcsec','sb_mag_sqarcsec_err']

        units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '',\
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,\
                 u.mag,u.mag,\
                 u.erg/u.s/u.cm**2/u.arcsec**2,u.erg/u.s/u.cm**2/u.arcsec**2,\
                 u.mag/u.arcsec**2,u.mag/u.arcsec**2]


        #self.sky_noise,self.sky_noise_erg]
        #'sky_noise_ADU_sqpix','sky_noise_erg_sqarcsec']
        #u.adu/u.s/u.pixel**2,u.erg/u.s/u.cm**2/u.arcsec**2]        
        columns = []
        for i in range(len(data)):
            columns.append(Column(data[i],name=names[i],unit=units[i]))
        
        t = Table(columns)
        t.write(outfile, format='fits', overwrite=True)

        if self.image2_flag:
            # write out photometry for h-alpha
            # radius enclosed flux
            if prefix is None:
                outfile = self.image2_name.split('.fits')[0]+'-phot.fits'
            else:
                outfile = self.image2_name.split('.fits')[0]+'-'+prefix+'-phot.fits'
    
            data = [self.apertures_a*self.pixel_scale,self.apertures_a, \
                self.flux2,self.flux2_err,\
                self.sb2, self.sb2_err, \
                self.sb2_snr, \
                self.flux2_erg, self.flux2_err_erg,\
                self.mag2, self.mag2_err, \
                self.sb2_erg_sqarcsec,self.sb2_erg_sqarcsec_err, \
                self.sb2_mag_sqarcsec,self.sb2_mag_sqarcsec_err]
            names = ['sma_arcsec','sma_pix','flux','flux_err',\
                'sb', 'sb_err', \
                'sb_snr', \
                'flux_erg', 'flux_erg_err',\
                'mag', 'mag_err', \
                'sb_erg_sqarcsec','sb_erg_sqarcsec_err', \
                'sb_mag_sqarcsec','sb_mag_sqarcsec_err']
            units = [u.arcsec,u.pixel,u.adu/u.s,u.adu/u.s, \
                 u.adu/u.s/u.pixel**2, u.adu/u.s/u.pixel**2, '',\
                 u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,\
                 u.mag,u.mag,\
                 u.erg/u.s/u.cm**2/u.arcsec**2,u.erg/u.s/u.cm**2/u.arcsec**2,\
                 u.mag/u.arcsec**2,u.mag/u.arcsec**2]
            columns = []
            for i in range(len(data)):
                columns.append(Column(data[i],name=names[i],unit=units[i]))
            if self.uconversion2b:
                data = [self.flux2_erg, self.flux2_err_erg,\
                      self.sb2_erg_sqarcsec,self.sb2_erg_sqarcsec_err, \
                      self.sb2_mag_sqarcsec,self.sb2_mag_sqarcsec_err]
                names = ['flux2_erg', 'flux2_err_erg',\
                         'sb2_erg_sqarcsec','sb2_erg_sqarcsec_err', \
                         'sb2_mag_sqarcsec','sb2_mag_sqarcsec_err']
                units = [u.erg/u.s/u.cm**2,u.erg/u.s/u.cm**2,\
                         u.erg/u.s/u.cm**2/u.arcsec**2,u.erg/u.s/u.cm**2/u.arcsec**2,\
                         u.mag/u.arcsec**2,u.mag/u.arcsec**2]
                for i in range(len(data)):
                    columns.append(Column(data[i],name=names[i],unit=units[i]))
            t = Table(columns)
            t.write(outfile, format='fits', overwrite=True)

            
    def draw_phot_results(self):
        ''' DRAW RESULTING FIT ON R-BAND CUTOUT, for gui '''
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
        ''' draw results in matplotlib figure '''
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
        ''' enclosed flux and surface brightness profiles, save figure '''
        plt.figure(figsize=(10,4))
        plt.subplots_adjust(wspace=.3)
        plt.subplot(2,2,1)
        #plt.plot(self.apertures_a,self.flux1,'bo')
        plt.errorbar(self.apertures_a,self.flux1,self.flux1_err,fmt='b.')
        plt.title('R-band')
        #plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Enclosed flux')
        plt.gca().set_yscale('log')
        if self.image2_flag:
            plt.subplot(2,2,2)
            plt.errorbar(self.apertures_a,self.flux2,self.flux2_err,fmt='b.')
            #plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Enclosed flux')
            plt.title('H-alpha')
            plt.gca().set_yscale('log')
        # plot surface brightness vs radius
        plt.subplot(2,2,3)
        #plt.plot(self.apertures_a,self.flux1,'bo')
        plt.errorbar(self.apertures_a,self.sb1,self.sb1_err,fmt='b.')
        plt.xlabel('semi-major axis (pixels)')
        plt.ylabel('Surface Brightess')
        plt.gca().set_yscale('log')
        if self.image2_flag:
            plt.subplot(2,2,4)
            plt.errorbar(self.apertures_a,self.sb2,self.sb2_err,fmt='b.')
            plt.xlabel('semi-major axis (pixels)')
            plt.ylabel('Surface Brightness')
            plt.gca().set_yscale('log')
        #plt.show()
        plt.savefig(self.image_name.split('.fits')[0]+'-enclosed-flux.png')

        
if __name__ == '__main__':
    image = 'MKW8-18216-R.fits'
    mask = 'MKW8-18216-R-mask.fits'
    image2 = 'MKW8-18216-CS.fits'
    nsaid='18045'
    prefix = 'MKW8-'
    nsaid='110430'
    nsaid='157146'
    prefix = 'NRGs27-'
    image = prefix+nsaid+'-R.fits'
    mask = prefix+nsaid+'-R-mask.fits'
    image2 = prefix+nsaid+'-CS.fits'
    # testing on 2017 pointing 1
    # second galaxy has clear halpha but profile is not fit
    # want to make sure we record some size
    image = 'v17p01-N119230-A742747-R.fits'
    rphot_table = 'v17p01-N119230-A742747-R_phot.fits'
    image2 = 'v17p01-N119230-A742747-CS.fits'
    haphot_table = 'v17p01-N119230-A742747-CS_phot.fits'
    mask = 'v17p01-N119230-A742747-R-mask.fits'
    image = 'v17p01-N118647-A8219-R.fits'
    rphot_table = 'v17p01-N118647-A8219-R_phot.fits'
    image2 = 'v17p01-N118647-A8219-CS.fits'
    haphot_table = 'v17p01-N118647-A8219-CS_phot.fits'
    mask = 'v17p01-N118647-A8219-R-mask.fits'
    myfilter = '4'
    myratio = .0406
    
    # testing on 2019 pointing 1
    # second galaxy has clear halpha but profile is not fit
    # want to make sure we record some size
    image = 'VFID3623-CGCG118-019-v19p001-R.fits'
    
    rphot_table = 'VFID3623-CGCG118-019-v19p001-R-phot.fits'
    image2 = 'VFID3623-CGCG118-019-v19p001-CS.fits'
    haphot_table = 'VFID3623-CGCG118-019-v19p001-CS-phot.fits'
    mask = 'VFID3623-CGCG118-019-v19p001-R-mask.fits'
    
    myfilter = 'inthalpha'
    myratio = .0356
    #image = 'MKW8-18037-R.fits'
    #mask = 'MKW8-18037-R-mask.fits'
    #image = 'r-18045-R.fits'
    #mask = 'r-18045-R-mask.fits'
    try:
        e = ellipse(image,mask=mask, image2=image2, use_mpl=True,image2_filter=myfilter, filter_ratio=myratio)
    except FileNotFoundError:
        print("so sorry, but no images were loaded")
        print("try e = ellipse(imagename) to start")
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

