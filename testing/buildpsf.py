'''
GOAL:
- create a psf image from an R-band mosaic


PROCEDURE
- use sextractor to identify objects
- use CLASS_STAR and flags to identify stars that are not saturated
- create a table of x,y coordinates (astropy.table.Table)
- extract stars - photutils.psf.extract_stars
- build psf - photutils.EPSFBuilder
- save psf image (e.g. for input into galfit)

sounds easy enough, right?


INPUT:
- image
- saturation level

OUTPUT:
- psf image

REFERENCES:
- following instructions found here:
https://photutils.readthedocs.io/en/stable/epsf.html#build-epsf


'''
import os
import numpy as np

from photutils import EPSFBuilder
from photutils.psf import extract_stars
from photutils import centroid_com, centroid_1dg, centroid_2dg

from astropy.nddata import NDData
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
from astropy.stats import gaussian_sigma_to_fwhm

from matplotlib import pyplot as plt

class psf_parent_image():
    def __init__(self, image=None, max_good=None, size=25, se_config = 'default.sex.HDI', sepath=None,nstars=100, pixelscale=0.43,oversampling=None):
        self.image_name = image
        self.data, self.header = fits.getdata(self.image_name, header=True)
        self.sat_level = max_good
        self.size = size
        self.config = se_config
        self.sextractor_files =['default.sex.HDI','default.param','default.conv','default.nnw']
        # path to source extractor files
        # default is probably fine if you have github in ~/github
        if sepath == None:
            self.sepath=os.getenv('HOME')+'/github/halphagui/astromatic/'
        else:
            self.sepath = sepath
        # number of stars to use to determine the psf
        self.nstars = nstars
        # default pixelscale is set for HDI camera
        self.pixelscale = pixelscale
        if oversampling == None:
            self.oversampling = 2
        else:
            self.oversampling = oversampling
    def link_files(self):
        # these are the sextractor files that we need
        # set up symbolic links from sextractor directory to the current working directory        
        for file in self.sextractor_files:
            os.system('ln -s '+self.sepath+'/'+file+' .')
            
    def clean_links(self):
        # clean up symbolic links to sextractor files
        # sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        for file in self.sextractor_files:
            os.system('unlink '+file)
    def run_all(self):
        self.runse()
        self.read_se_table()
        self.find_stars()
        self.extract_stars()
        self.show_stars()
        self.build_psf()
        self.show_psf()
        self.measure_fwhm()
        self.save_psf_image()

    def runse(self):
        self.link_files()
        s = 'sex %s -c %s  -SATUR_LEVEL 40000.0'%(self.image_name,self.config)
        #print(s)
        os.system(s)
        ###################################
        # Read in Source Extractor catalog
        ###################################

        secat = fits.getdata('test_cat.fits',2)
        ###################################
        # get median fwhm of image
        ###################################
        fwhm = np.median(secat['FWHM_IMAGE'])*self.pixelscale

        #############################################################
        # rerun Source Extractor catalog with updated SEEING_FWHM
        #############################################################
        s = 'sex %s -c %s  -SATUR_LEVEL 40000.0 -SEEING_FWHM %s'%(self.image_name,self.config,str(fwhm))

        os.system(s)


    def read_se_table(self):
        self.secat = fits.getdata('test_cat.fits',2)
        pass
    def find_stars(self):
        # select objects with CLASS_STAR > 0.98
        # and FLUX_MAX 
        x = self.secat['X_IMAGE']
        y = self.secat['Y_IMAGE']

        # select based on star/class separator and no flags
        flag1 = (self.secat['CLASS_STAR'] > 0.95) & (self.secat['FLAGS'] == 0)

        # remove stars that are near the edge
        # being conservative by using the 10x full image size as the buffer rather than half image size
        # in part to compensate from zeros that sometimes are seen around perimeter after swarp
        flag2 = ((x > 5.*self.size) & (x < (self.data.shape[1] -1 - 5.*self.size)) &
                (y > 5.*self.size) & (y < (self.data.shape[0] -1 - 5.*self.size)))

        # remove stars with close neighbors
        c = SkyCoord(self.secat['ALPHA_J2000']*u.deg,self.secat['DELTA_J2000']*u.deg, frame='icrs')
        # returns index of closest match (not itself), separation in deg,
        # and 3d sep (not sure what this means if I don't provide redshift
        idx, sep2d, dist3d = c.match_to_catalog_sky(c,nthneighbor=2)
        # make sure stars don't have a neighbor within 15" (picked that fairly randomly...)
        flag3 = sep2d > 15./3600.*u.deg
        star_flag = flag1 & flag2 & flag3
        # select 25 objects with FLUX_MAX closest to median value
        # in other words, select 25 most central objects
        fm = self.secat['FLUX_MAX'][star_flag]
        x = x[star_flag]
        y = y[star_flag]
        sorted_indices = fm.argsort()
        print('number of potential stars', sum(flag1), sum(flag2))
        # select stars within +/- nstar/2 from a threshold
        # where threshold is the percentile ranking according to peak flux
        threshold = .65
        lower = int(threshold*len(sorted_indices)) - int((self.nstars)/2)
        upper = int(threshold*len(sorted_indices)) + int((self.nstars)/2)
        #lower = int(.25*len(sorted_indices)) 
        #upper = int(.85*len(sorted_indices)) 
        print('number of psf stars = ',upper-lower+1)
        self.xstar = x[lower:upper]
        self.ystar = y[lower:upper]
        print(self.xstar.shape)

        self.stars_tbl = Table()
        self.stars_tbl['x'] = self.xstar
        self.stars_tbl['y'] = self.ystar

    def extract_stars(self):
        # just in case sky was not properly subtracted
        mean_val, median_val, std_val = sigma_clipped_stats(self.data, sigma=2.)  
        self.data -= median_val
        nddata = NDData(data=self.data)  
        self.stars = extract_stars(nddata, self.stars_tbl, size=self.size)  
        # check to make sure stars don't have zeros
        # Virgo 2017 pointing-4_R.coadd.fits is giving me trouble b/c a lot of stars have zeros
        keepflag = np.ones(len(self.stars),'bool')
        for i,s in enumerate(self.stars):
            if len(np.where(s.data == 0)[0]) > 1:
                keepflag[i] = False
        self.keepflag2 = keepflag
        self.stars = extract_stars(nddata, self.stars_tbl[keepflag], size=self.size)  
        #self.stars = self.stars[keepflag]
        if len(self.stars) < self.nstars:
            self.nstars = len(self.stars)
    def show_stars(self):
        # round up to integer
        nrows = int(np.ceil(np.sqrt(len(self.stars))))
        ncols = int(np.ceil(np.sqrt(len(self.stars))))
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10), squeeze=True)
        ax = ax.ravel()
        for i in range(self.nstars):
            norm = simple_norm(self.stars[i], 'log', percent=99.)
            ax[i].imshow(self.stars[i], norm=norm, origin='lower', cmap='viridis')
        plt.show()
        plt.savefig('allstars.png')
    def build_psf(self):
        if self.oversampling == None:
            epsf_builder = EPSFBuilder(maxiters=12, progress_bar=False, smoothing_kernel=None, recentering_func = centroid_com)
        else:
            epsf_builder = EPSFBuilder(oversampling=self.oversampling, maxiters=13, progress_bar=False,  recentering_func = centroid_com, smoothing_kernel=None)  
        self.epsf, self.fitted_stars = epsf_builder(self.stars)
    def show_psf(self):
        norm = simple_norm(self.epsf.data, 'log', percent=99.)
        plt.figure()
        plt.imshow(self.epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        plt.show()
        plt.savefig('psf.png')
    def measure_fwhm(self):
        nx,ny = self.epsf.data.shape
        x,y = np.mgrid[:nx,:ny]
        # gaussian model with initial guess for center and std
        # assume the PSF is in the center of the image
        # assume a std of 4
        p_init = models.Gaussian2D(x_mean=nx/2.,y_mean=ny/2., x_stddev=4.,y_stddev=4.)
        fit_p = fitting.LevMarLSQFitter()
        t = fit_p(p_init,x,y,self.epsf.data)
        self.x_stddev = t.x_stddev.value
        self.y_stddev = t.y_stddev.value
        # compute total radial std
        # correct for oversampling of PSF
        self.std = np.sqrt(self.x_stddev**2 + self.y_stddev**2)/self.oversampling
        # convert gaussian std to fwhm
        self.fwhm = self.std*gaussian_sigma_to_fwhm
        print('image fwhm = %.2f pix (%.2f arcsec)'%(self.fwhm, self.fwhm*self.pixelscale))
    def save_psf_image(self):
        # save the psf file
        # outfile = append '-psf' to the input image name
        # just use the basename (filename), and drop any path information
        # this will have the effect of saving the psf image in the current directory,
        # rather than the directory where the image is
        basename = os.path.basename(self.image_name)
        self.psf_image_name = basename.split('.fits')[0]+'-psf.fits'
        print('PSF image name = ',self.psf_image_name)
        fits.writeto(self.psf_image_name,self.epsf.data, overwrite=True)

        # update image header
        # I know this is a kludge, but not sure how to access the default header that it writes
        # BEFORE it gets written

        data, header = fits.getdata(self.psf_image_name, header=True)

        header.append(card=('FWHM', float('{:.2f}'.format(self.fwhm)), 'PSF fwhm in pixels'))
        header.append(card=('STD', float('{:.2f}'.format(self.std)), 'PSF STD in pixels'))
        header.append(card=('OVERSAMP', self.oversampling, 'PSF oversampling'))
        fits.writeto(self.psf_image_name, data, header=header, overwrite=True)
        
if __name__ == '__main__':
    #image = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_R.coadd.fits'
    image = '/Users/rfinn/research/VirgoFilaments/Halpha/virgo-coadds-2017/pointing-4_R.coadd.fits'
    p = psf_parent_image(image=image, size=21, nstars=100, oversampling=2)
    p.runse()
    p.read_se_table()
    p.find_stars()
    p.extract_stars()
    p.show_stars()
    p.build_psf()
    p.show_psf()
    p.measure_fwhm()
    p.save_psf_image()
    '''
    getting some weird noise at the bottom of the psf image

    need to read more about how it handles the edge

    not sure that this will impact results significantly, so moving on...

    '''
    
