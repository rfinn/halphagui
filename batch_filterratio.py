#!/usr/bin/env python

'''
GOAL: 
- get filter ratio for all images in the current directory

- need to account for the different naming conventions of the different Halpha filters


USAGE:

Run from the directory with the coadds /media/rfinn/hdata/coadds/all-virgo-coadds

python ~/github/halphagui/batch_filterratio.py

'''

import os
import glob
from astropy.io import fits
import matplotlib
import time
import filterratio as runse
import multiprocessing as mp

matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")

# overwrite output files if they exist
overwrite = False

import argparse

image_results = []
def collect_results(result):

    global results
    image_results.append(result)


def subtract_images(rimage,himage,filter_ratio):
    r = fits.open(rimage)
    ha = fits.open(himage)
    halpha_cs = ha[0].data - filter_ratio*r[0].data
    hacoadd_cs_fname = himage.split('.fits')[0]+'-CS.fits'
    print("subtracting images: ",himage," -> ",hacoadd_cs_fname)
    fits.writeto(hacoadd_cs_fname,halpha_cs,header=ha[0].header,overwrite=True)
    r.close()
    ha.close()
    pass

parser = argparse.ArgumentParser(description ='Run filterratio.py for all images in the specified directory')
parser.add_argument('--plotdir',dest = 'plotdir', default ='plots-filterratio', help = 'directory for to store the plots in')
parser.add_argument('--prefix', dest = 'prefx', default = 'VF-',help = "prefix for coadds.  the default is VF-")

args = parser.parse_args()



plotdir = os.path.join(os.getcwd(),args.plotdir)
if not os.path.exists(plotdir):
    os.mkdir(plotdir)


def getoneratio(fname,instrument):
    # find the corresponding Halpha image
    if f == 'INT':
        # could end in Halpha or Ha6657
        t = rimage.split('-')
        #print(t)
        if rimage[10] == '+':
            dateobs = t[3]
            pointing = t[4]
        else:
            dateobs = t[4]
            pointing = t[5]
        if dateobs == '20190531': # ugh - new month too!
            newdate = '20190601'
        else:
            newdate = dateobs[:-2]
        hastring = 'VF-*INT-'+newdate+'*-'+pointing+'*-Halpha.fits'
        hfiles = glob.glob(hastring)
        if len(hfiles) < 1:
            hastring = 'VF-*INT-'+newdate+'*-'+pointing+'*-Ha6657.fits'
            hfiles = glob.glob(hastring)

        #print("found these ",hfiles)
        if len(hfiles) > 0:
            himage = hfiles[0]
        else:
            print()
            print('Warning: no halpha image found for ',rimage)
            print('moving to the next image...')
            return
    elif f == 'BOK':
        # should end in Ha4.fits
        testname = rimage.replace('-r.fits','-Ha4.fits')
        if os.path.exists(testname):
            himage = testname
        else:
            print()
            print('Warning: no halpha image found for ',rimage)
            print('moving to the next image...')
            return
    elif f == 'HDI':
        # should end in ha4.fits
        testname = rimage.replace('-r.fits','-ha4.fits').replace('-R.fits','-ha4.fits')

        print('testname = ',testname)
        if os.path.exists(testname):
            himage = testname

        else:
            # halpha image could have a different date
            t = rimage.split('-')

            if rimage[10] == '+':
                dateobs = t[3]
                pointing = t[4]
            else:
                dateobs = t[4]
                pointing = t[5]
            hastring = 'VF-*'+f+'-*'+dateobs[:-2]+'*-'+pointing+'*-ha4.fits'
            hfiles = glob.glob(hastring)

            #print("found these ",hfiles)
            if len(hfiles) > 0:
                himage = hfiles[0]
                print('halpha image = ',himage)                
            else:
                print()
                print('Warning: no halpha image found for ',rimage)
                print('moving to the next image...')
                return


    print()
    print('##########################################')        
    print('GETTING FILTER RATIO FOR: ',rimage,himage)
    print('##########################################')

    start_time = time.perf_counter()

    runse.run_sextractor(rimage, himage)
    ave, std = runse.make_plot(rimage, himage, return_flag = True, plotdir = plotdir)
    #print(ave,std)

    subtract_images(rimage,himage,ave)

    # add ratio to r-band image headers
    r,header = fits.getdata(rimage,header=True)
    header.set('FLTRATIO',ave)
    header.set('FLTR_ERR',std)
    header.set('HAIMAGE',os.path.basename(himage))
    fits.writeto(rimage,r,header=header,overwrite=True)        
    # clock time to get filter ratio
    end_time = time.perf_counter()
    print('\t total time = ',end_time - start_time)


    
########################################################
# loop through telescopes and get r-band images
########################################################

# INT
# r-shifted.fits
# Halpha or Ha6657

# BOK
# r.fits
# Ha4

# HDI
# -r.fits -R.fits
# ha4

# will need to add more for mosaic images

# telescope/instrument names
inames = ["INT","BOK","HDI"]
inames = ["BOK","HDI"]
inames = ["HDI"]
inames = ["BOK"]


for i,f in enumerate(inames):
    # get list of current directory
    # this will grab the coadds but not the weight images
    if f == 'INT':
        #print("match string = ",'VF-*'+f+'*r-shifted.fits')
        rimages = glob.glob('VF-*'+f+'*r-shifted.fits')
    elif f == 'BOK':
        rimages = glob.glob('VF-*'+f+'*r.fits')
    elif f == 'HDI':
        rimages1 = glob.glob('VF-*'+f+'*r.fits')
        rimages2 = glob.glob('VF-*'+f+'*R.fits')
        rimages = rimages1 + rimages2
    print(f"found {len(rimages)} rband images for {f}")

    # sort the file list
    rimages.sort()
    
    #for rimage in rimages: # loop through list
    image_pool = mp.Pool(mp.cpu_count())
    myresults = [image_pool.apply_async(getoneratio,args=(im,f),callback=collect_results) for im in rimages]
    
    image_pool.close()
    image_pool.join()
    image_results = [r.get() for r in myresults]




