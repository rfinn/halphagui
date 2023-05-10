#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to run SExtractor in two image mode, to detect objects on the R-band image and
  apply the same apertures to the Halpha image.

PROCEDURE:
  - copy default setup files to current directory
  - Run sextractor on each image

EXAMPLE:

EXTRA NOTES:
- updated 8/14/19 to encase code in two functions.  this makes it easier to import into other programs, like the Halpha gui.

WRITTEN BY:
Rose Finn, 04 Jan 2017

'''
import os
from astropy.io import fits
import argparse
import numpy as np
from astropy.stats.sigma_clipping import sigma_clip


def run_sextractor(image1,image2, default_se_dir = '/Users/rfinn/github/halphagui/astromatic'):
    # get magnitude zeropoint for image 1 and image 2
    header1 = fits.getheader(image1)
    header2 = fits.getheader(image2)
    try:
        ZP1 = header1['PHOTZP']
        print('got ZP = ',ZP1)
        zp1flag = True
    except KeyError:
        print('no PHOTZP found in image 1 header.  Too bad :(')
        print('did you run getzp.py?')
        zp1flag = False
    try:
        ZP2 = header2['PHOTZP']
        print('got ZP = ',ZP2)
        zp2flag = True
    except KeyError:
        print('no PHOTZP found in image 2 header.  Too bad :(')
        print('did you run getzp.py?')
        zp2flag = False

    base = os.path.basename(image1)
    froot1 = os.path.splitext(base)[0]
    base = os.path.basename(image2)
    froot2 = os.path.splitext(base)[0]
    
    # check if output catalogs exist - if they do, don't rerun SE

    secatalog1 = froot1+'.cat'
    secatalog2 = froot2+'.cat'
    if os.path.exists(secatalog1) and os.path.exists(secatalog2):
        print("FOUND SE CATALOGS - NOT RERUNNING")
    else:
        print('RUNNING SOURCE EXTRACTOR')

        if zp1flag:
            os.system('sex ' + image1+','+image1 + ' -c default.sex.HDI -CATALOG_NAME ' + froot1 + '.cat -MAG_ZEROPOINT '+str(ZP1))
        else:
            os.system('sex ' + image1+','+image1 + ' -c default.sex.HDI -CATALOG_NAME ' + froot1 + '.cat')
        if zp2flag:
            os.system('sex ' + image1+','+image2 + ' -c default.sex.HDI -CATALOG_NAME ' + froot2 + '.cat -MAG_ZEROPOINT '+str(ZP2))
        else:
            os.system('sex ' + image1+','+image2 + ' -c default.sex.HDI -CATALOG_NAME ' + froot2 + '.cat')


def make_plot(image1, image2, return_flag = False, image_dir = './'):
    from matplotlib import pyplot as plt
    from scipy.stats import scoreatpercentile
    base = os.path.basename(image1)
    froot1 = os.path.splitext(base)[0]
    #cat1 = fits.getdata(froot1+'.cat')
    cat1 = fits.getdata(froot1+'.cat',2)
    print(froot1+'.cat')

    base = os.path.basename(image2)
    froot2 = os.path.splitext(base)[0]
    
    cat2 = fits.getdata(froot2+'.cat',2)
    plt.figure(figsize=(6,6))
    plt.subplots_adjust(bottom=.2,left=.15,right=.95,hspace=.5)
    plt.subplot(2,1,1)

    # cut out extreme outliers using scoreatpercentile
    xmax = scoreatpercentile(cat1.FLUX_AUTO,95)
    ymax = scoreatpercentile(cat2.FLUX_AUTO,95)

    keepflag = (cat1.FLUX_AUTO < xmax) & (cat1.FLUX_AUTO > 0) & (cat2.FLUX_AUTO < ymax) & (cat2.FLUX_AUTO > 0)  
    plt.plot(cat1.FLUX_AUTO[keepflag],cat2.FLUX_AUTO[keepflag],'k.',alpha=.4)
    c = np.polyfit(cat1.FLUX_AUTO[keepflag],cat2.FLUX_AUTO[keepflag],1,cov=True)
    print("results from polyfit = ",c)

    
    xline = np.linspace(0,xmax,100)
    plt.plot(xline,np.polyval(c[0],xline),ls='--')

    print()
    #plt.xlim(0,xmax)
    #plt.ylim(0,ymax)
    plt.xlabel("r-band FLUX_AUTO")
    plt.ylabel("NB FLUX_AUTO")

    filename = os.path.basename(image2)
    t = filename.replace('.fits','')    
    plt.title(t,fontsize=12)
    plt.text(0.05,.9,'$med ratio = %.4f (%.4f)$'%(c[0][0],np.sqrt(c[1][0][0])),transform=plt.gca().transAxes,fontsize=8)    
    
    plt.subplot(2,1,2)
    x = cat2.FLUX_AUTO
    y = cat2.FLUX_AUTO/cat1.FLUX_AUTO
    flag = keepflag & (cat1.FLAGS == 0)
    plt.plot(x[flag],y[flag],'k.',alpha=.4)
    #plt.axis([0,500,.02,.08])
    x1,x2 = plt.xlim()
    #plt.xlim(.2,x2)
    #  TODO this is not working well for BOK
    #  need to do sigma clipping first?
    data = y[flag]
    clipped_data = sigma_clip(data,cenfunc='median',stdfunc='mad_std',sigma=2)
    ave = np.ma.median(clipped_data)
    # use the MAD instead
    std = np.ma.std(clipped_data)

    
    plt.axhline(y=ave,ls='--')
    plt.ylim(-0.5*ave,3*ave)
    print('%.4f (%.4f)'%(ave,std))
    plt.ylabel('Flux(Halpha)/Flux(R)')
    plt.xlabel('Flux(Halpha) (ADU)')
    plt.text(0.05,.9,'$med ratio = %.4f (%.4f)$'%(ave,std),transform=plt.gca().transAxes,fontsize=8)
    #plt.gca().set_xscale('log')

    #plt.show()
    #plt.axis([0,5000,-.06,.06])
    plt.savefig(image_dir+'/'+t+'-filter-ratio.png')

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

    
