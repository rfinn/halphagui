#!/usr/bin/env python

'''
GOAL: 
- go into each subdirectory and run getzp on coadded images
- delete any intermediate images to save space

Run this from, e.g. /home/rfinn/research/Virgo/gui-output-2019
- this directory has a subdirectory for each pointing


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib

import argparse

parser = argparse.ArgumentParser(description ='Run sextractor, scamp, and swarp to determine WCS solution and make mosaics')
parser.add_argument('--coaddir',dest = 'coaddir', default ='/home/rfinn/data/reduced/virgo-coadds-feb2019-int/', help = 'directory for coadds. Default is /home/rfinn/data/reduced/virgo-coadds/feb2019-int/')

args = parser.parse_args()



matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
telescope = 'INT'

# get list of current directory
imagedir = args.coaddir
flist1 = glob.glob(imagedir+'VF-*r-shifted.fits')
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True

# just doing this for rogue pointings 022 and 026
#flist1a = glob.glob('VF-*p022*r.fits')
#flist1b = glob.glob('VF-*p026*r.fits')
#flist1 = flist1a+flist1b
flist1.sort()
#print(flist1)


# shift r to match Halpha
i = 0
for rimage in flist1: # loop through list

    print()
    print('##########################################')
    print('##########################################')
    print('##########################################')
    print('##########################################')        
    print('WORKING ON IMAGE: ',rimage)
    print('##########################################')
    print('##########################################')        
    
    # read in r-band images
    # find matching Halpha image

    # grab other coadds

    rootname = rimage.split('-r')[0]
    rweightimage = rootname+'-r.weight.fits'
    rweightimage = rweightimage.replace('-shifted','')
    # last entry is the pointing name - match on this
    # because sometimes the Halpha coordinates are slightly different
    # or the UT date changed between Halpha and r images
    pointing = rootname.split('-')[-1]
    dirname = os.path.dirname(rimage)
    coadds = glob.glob(dirname+'/VF*'+pointing+'*.fits')
    #print('rootname = ',rootname)
    #print(coadds)
    #print(coadds)
    haimage = None
    for c in coadds:
        #print(c)
        if c.find('CS') > -1:
            print('skipping CS image')
            continue
        if (c.find('-Halpha.fits') > -1) & (c.find('weight') < 0):
            haimage = c
            hfilter = 'inthalpha'
            print('haimage = ',c)
        elif (c.find('-Ha6657') > -1) & (c.find('weight') < 0):
            haimage = c
            hfilter = 'intha6657'
            print('haimage = ',c)            
            #print('matching ha image: ',haimage)
    if haimage is not None:
        #print(rootname)
        pointing = os.path.basename(rootname).split('-')[4]
        prefix = 'v19'+pointing
        print('prefix = ',prefix)
        #command_string = 'python ~/github/HalphaImaging/python3/INT_align_images.py --image1 {} --image2 {} --weight2 {}'.format(haimage,rimage,rweightimage)
        command_string = 'python  ~/github/halphagui/halphamain.py --virgo --rimage {} --haimage {} --filter {} --psfdir /home/rfinn/data/reduced/psf-images/ --tabledir /home/rfinn/research/Virgo/tables-north/v1/ --prefix {} --auto'.format(rimage,haimage,hfilter,prefix)

        # check to see if shifted r-band image exists.  if 
        try:
            print('running : ',command_string)
            os.system(command_string)
        except:
            print('##########################################')
            print('WARNING: problem running align_images for ',rimage)
            print('##########################################')

    # just running on one directory for testing purposes
    #i += 1
    #if i > 3:
    #    break


