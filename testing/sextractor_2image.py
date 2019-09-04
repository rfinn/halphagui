#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to run SExtractor in two image mode, to detect objects on the R-band image and
  apply the same apertures to the Halpha image.

PROCEDURE:
  - copy default setup files to current directory
  - Run sextractor on each image

EXAMPLE:
   In the directory containing all flattened objects with fixed headers to run sextractor type in the command line:
      '/Users/alfalfa/Github/HalphaImaging/uat_sextractor_2image.py --s'(or whatever the path is to where this program is stored)

    from within ipython:
    
    %run ~/github/HalphaImaging/uat_sextractor_2image.py --image1 A1367_R
   ...: .coadd.fits --image2 A1367_ha12.coadd.fits --plot

    -You want to run rimage against the rimage to create .cat file for R, then use R image with the Haimage (image 2) to create .cat file for Ha.



WHAT THIS CODE DOES:
INPUT/OUPUT:
REQUIRED MODULES:
EXTRA NOTES:
- updated 8/14/19 to encase code in two functions.  this makes it easier to import into other programs, like the Halpha gui.

WRITTEN BY:
Rose Finn, 04 Jan 2017

'''
import glob
import os
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import subprocess
import numpy as np



def run_sextractor(image1,image2, default_se_dir = '/Users/rfinn/github/HalphaImaging/astromatic'):
    # get magnitude zeropoint for image 1 and image 2
    header1 = fits.getheader(image1)
    header2 = fits.getheader(image2)
    try:
        ZP1 = header1['PHOTZP']
        zp1flag = True
    except KeyError:
        print('no PHOTZP found in image 1 header.  Too bad :(')
        print('did you run getzp.py?')
        zp1flag = False
    try:
        ZP2 = header2['PHOTZP']
        zp2flag = True
    except KeyError:
        print('no PHOTZP found in image 2 header.  Too bad :(')
        print('did you run getzp.py?')
        zp2flag = False

    print('RUNNING SEXTRACTOR')
    t = image1.split('.fits')
    froot1 = t[0]
    if zp1flag:
        os.system('sex ' + image1+','+image1 + ' -c default.sex.hdi -CATALOG_NAME ' + froot1 + '.cat -MAG_ZEROPOINT '+str(ZP1))
    else:
        os.system('sex ' + image1+','+image1 + ' -c default.sex.hdi -CATALOG_NAME ' + froot1 + '.cat')
    os.rename('check.fits', froot1 + 'check.fits')
    # run on second image
    t = image2.split('.fits')
    froot2 = t[0]
    if zp2flag:
        os.system('sex ' + image1+','+image2 + ' -c default.sex.hdi -CATALOG_NAME ' + froot2 + '.cat -MAG_ZEROPOINT '+str(ZP2))
    else:
        os.system('sex ' + image1+','+image2 + ' -c default.sex.hdi -CATALOG_NAME ' + froot2 + '.cat')
    os.rename('check.fits', froot2 + 'check.fits')

def make_plot(image1, image2, return_flag = False, image_dir = './'):
    from matplotlib import pyplot as plt
    t = image1.split('.fits')
    froot1 = t[0]
    cat1 = fits.getdata(froot1+'.cat',2)
    t = image2.split('.fits')
    froot2 = t[0]
    cat2 = fits.getdata(froot2+'.cat',2)
    plt.figure(figsize=(6,4))
    plt.subplots_adjust(bottom=.2,left=.15,right=.95)
    x = cat2.FLUX_AUTO
    y = cat2.FLUX_AUTO/cat1.FLUX_AUTO
    flag = cat1.FLAGS == 0
    plt.plot(x[flag],y[flag],'k.')
    #plt.axis([0,500,.02,.08])
    x1,x2 = plt.xlim()
    plt.xlim(.2,x2)
    ave = np.median(y[(x> 50) & flag])
    std = np.std(y[(x > 50) & flag])
    plt.axhline(y=ave)
    print('%.4f (%.4f)'%(ave,std))
    plt.ylabel('$Flux (Halpha)/Flux(R) $',fontsize=20)
    plt.xlabel('$Flux(R) \ (ADU)$',fontsize=20)
    plt.text(20,.07,'$ ratio = %.4f (%.4f)$'%(ave,std),fontsize=12)
    plt.gca().set_xscale('log')
    filename = os.path.basename(image2)
    t = filename.split('.coadd')
    plt.title(t[0],fontsize=12)
    plt.show()
    plt.savefig(image_dir+t[0]+'-filter-ratio.png')
    if return_flag:
        return ave, std

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ="Run sextractor in two-image mode.  \n To run from within ipython:\n %run ~/github/HalphaImaging/uat_sextractor_2image.py --image1 pointing-1_R.coadd.fits --image2 pointing-1_ha4.coadd.fits --plot --imagedir './' ")
    #parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Run sextractor to create object catalogs')
    parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files')
    parser.add_argument('--image1',dest = 'image1', default = None,  help = 'image used to define apertures (R-band)')
    parser.add_argument('--image2',dest = 'image2', default = None,  help = 'image used to for measuring phot based on image1 (typically this is the Halpha image)')
    parser.add_argument('--plot',dest = 'plot', default = False, action = 'store_true', help = 'make diagnostic plots')
    parser.add_argument('--imagedir',dest = 'imagedir', default = '.', help = 'directory for saving plots')

    args = parser.parse_args()

    # get input files
    #print 'cp ' +args.d + '/default.* .'
    os.system('cp ' +args.d + '/default.* .')
    #files = sorted(glob.glob(args.filestring))

    #nfiles = len(files)
    i = 1
    run_sextractor(args.image1, args.image2, default_se_dir=args.d)
    if args.plot:
        make_plot(args.image1, args.image2, image_dir = args.imagedir)

    
