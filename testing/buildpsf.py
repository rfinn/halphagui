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
https://photutils.readthedocs.io/en/stable/epsf.html#build-epsf


'''
import os
from photutils import EPSFBuilder
from photutils.psf import extract_stars
from astropy.nddata import NDData
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.visualization import simple_norm

class parent_image():
    def __init__(self, image=None, max_good=None, size=25, se_config = 'default.sex.HDI', sepath=None,nstars=25):
        self.image_name = image
        self.data, self.header = fits.getdata(self.image_name, header=True)
        self.sat_level = max_good
        self.size = size
        self.config = se_config
        self.sextractor_files =['default.sex.HDI','default.param','default.conv','default.nnw']
        if sepath == None:
            self.sepath=os.getenv('HOME')+'/github/halphagui/astromatic/'
        else:
            self.sepath = sepath
        self.nstars = nstars
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
    def runse(self):
        self.link_files()
        s = 'sex %s -c %s  -SATUR_LEVEL 40000.0'%(self.image_name,self.config)
        print(s)
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
        hsize = (self.size - 1) / 2
        flag2 = ((x > hsize) & (x < (self.data.shape[1] -1 - hsize)) &
                (y > hsize) & (y < (self.data.shape[0] -1 - hsize)))

        # remove stars with close neighbors
        c = SkyCoord(self.secat['
        star_flag = flag1 & flag2
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
        lower = int(threshold*len(sorted_indices)) - int((self.nstars-1)/2)
        upper = int(threshold*len(sorted_indices)) + int((self.nstars-1)/2)
        print('number of psf stars = ',upper-lower)
        self.xstar = x[lower:upper+1]
        self.ystar = y[lower:upper+1]

        self.stars_tbl = Table()
        self.stars_tbl['x'] = self.xstar
        self.stars_tbl['y'] = self.ystar

    def extract_stars(self):
        nddata = NDData(data=self.data)  
        self.stars = extract_stars(nddata, self.stars_tbl, size=25)  
        pass
    def show_stars(self):
        nrows = 5
        ncols = 5
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10), squeeze=True)
        ax = ax.ravel()
        for i in range(nrows*ncols):
            norm = simple_norm(self.stars[i], 'log', percent=99.)
            ax[i].imshow(self.stars[i], norm=norm, origin='lower', cmap='viridis')
        plt.show()
    def build_psf(self):
        epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=False)  
        self.epsf, self.fitted_stars = epsf_builder(self.stars)
    def show_psf(self):
        norm = simple_norm(self.epsf.data, 'log', percent=99.)
        plt.figure()
        plt.imshow(self.epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar()
        plt.show()
    def save_psf_image(self):
        pass

if __name__ == '__main__':
    image = '/Users/rfinn/research/HalphaGroups/reduced_data/HDI/20150418/MKW8_R.coadd.fits'
    p = parent_image(image=image)
    p.runse()
    p.read_se_table()
    p.find_stars()
    p.extract_stars()
    p.show_stars()
    p.build_psf()
    p.show_psf()
    
