#!/usr/bin/env python

'''
GOAL: 

Run this from, e.g. /home/rfinn/research/Virgo/gui-output-2019
- the gui will create a cutout folder in this directory that has a subdirectory for each pointing


python ~/github/halphagui/batch_gui.py --coaddir /media/rfinn/hdata/coadds/BOK2021pipeline/ --bok --psfdir /media/rfinn/hdata/psf-images/


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib

import argparse

parser = argparse.ArgumentParser(description ='run the halpha gui in automatic mode.  Run this from the directory where you want the output data stored.  For example: /home/rfinn/research/Virgo/gui-output-NGC5846/')
parser.add_argument('--coaddir',dest = 'coaddir', default ='/home/rfinn/data/reduced/virgo-coadds-feb2019-int/', help = 'directory for coadds. Default is /home/rfinn/data/reduced/virgo-coadds/feb2019-int/')
parser.add_argument('--hdi',dest = 'hdi', default =False, action='store_true', help = 'set this for HDI data.  it will grab the filenames according to the HDI naming convention for the coadded images.')
parser.add_argument('--mosaic',dest = 'mosaic', default =False, action='store_true', help = "set this for Mosaic data.  it will grab the filenames according to Becky's naming convention for the coadded images.")
parser.add_argument('--bok',dest = 'bok', default =False, action='store_true', help = "set this for Bok 90prime data.  it will grab the filenames according to naming convention for the 90Prime images.")
parser.add_argument('--psfdir',dest = 'psfdir', default='/home/rfinn/data/reduced/psf-images/', help = "directory containing PSF images.  Default is /home/rfinn/data/reduced/psf-images/, which is for laptop.  When running on virgo vms, set to /mnt/qnap_home/rfinn/Halpha/reduced/psf-images/")
#parser.add_argument('--getgalsonly',dest = 'psfdir', default='/home/rfinn/data/reduced/psf-images/', help = "directory containing PSF images.  Default is /home/rfinn/data/reduced/psf-images/, which is for laptop.  When running on virgo vms, set to /mnt/qnap_home/rfinn/Halpha/reduced/psf-images/")

args = parser.parse_args()



matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
telescope = 'INT'

# get list of current directory
imagedir = args.coaddir
if args.hdi:
    # HDI has rband images taken in r and R, so grab both
    flista = glob.glob(imagedir+'VF-*r-noback-coadd.fits')
    flistb = glob.glob(imagedir+'VF-*R-noback-coadd.fits')
    flist1 = flista + flistb
elif args.mosaic:
    flist1 = glob.glob(imagedir+'n*R.fits')
elif args.bok:
    flist1 = glob.glob(imagedir+'VF-*r.fits')
else: # assume INT
    flist1 = glob.glob(imagedir+'VF-*r-shifted.fits')
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True

flist1.sort()


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

    if args.hdi:
        print()
        print('#######  HDI DATA #########')
        print()        
        if rimage.find('-R-') > -1:
            rfilter = 'R'
        elif rimage.find('-r-') > -1:
            rfilter = 'r'
        print(rimage)
        rootname = rimage.split('-'+rfilter+'-')[0] # should be VF-2018-03-16-HDI-p054
        print(rootname)        
        rweightimage = rimage.split('.fits')[0]+'.weight.fits'
        pointing = rootname.split('-')[-1]
        dirname = os.path.dirname(rimage)
        coadds = glob.glob(rootname+'*.fits')

    elif args.mosaic:
        print()
        print('#######  MOSAIC DATA #########')
        print()        
        if rimage.find('R.fits') > -1:
            rfilter = 'R'
        elif rimage.find('r.fits') > -1:
            rfilter = 'r'
        print(rimage)
        rootname = rimage.split(rfilter+'.')[0] # should be VF-2018-03-16-HDI-p054
        print(rootname)        
        rweightimage = None
        pointing = rootname.split('_')[-1].replace('R','')
        dirname = os.path.dirname(rimage)
        coadds = glob.glob(rootname+'*.fits')

    elif args.bok:
        rootname = rimage.split('-r')[0]
        rweightimage = rootname+'-r.weight.fits'

        # last entry is the pointing name - VFIDXXXX for the case of the Bok data
        pointing = rootname.split('-')[-1]
        dirname = os.path.dirname(rimage)
        coadds = glob.glob(dirname+'/VF*'+pointing+'*.fits')
    else:
        rootname = rimage.split('-r')[0]
        rweightimage = rootname+'-r.weight.fits'
        rweightimage = rweightimage.replace('-shifted','')

        # last entry is the pointing name - match on this
        # because sometimes the Halpha coordinates are slightly different
        # or the UT date changed between Halpha and r images
        pointing = rootname.split('-')[-1]
        print(f"\nPointing name = {pointing}\n")
        dirname = os.path.dirname(rimage)
        coadds = glob.glob(dirname+'/VF*'+pointing+'*.fits')
    #print('rootname = ',rootname)
    #print(coadds)
    #print(coadds)
    haimage = None
    print(coadds)
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
        elif (c.find('-ha4') > -1) & (c.find('weight') < 0):
            haimage = c
            hfilter = '4'
            print('haimage = ',c)            
            #print('matching ha image: ',haimage)
        
        elif (c.find('-Ha4') > -1) & (c.find('weight') < 0):
            haimage = c
            hfilter = '4'
            print('haimage = ',c)            
            
        elif (c.find('Ha.fits') > -1) & (c.find('weight') < 0):
            haimage = c
            hfilter = '4'
            print('haimage = ',c)            
            #print('matching ha image: ',haimage)
    if haimage is not None:
        #print(rootname)
        if args.hdi:
            prefix = os.path.basename(rootname)
            #prefix = None
        elif args.mosaic:
            prefix = os.path.basename(rootname)
        elif args.bok:
            prefix = os.path.basename(rootname)
            #prefix = None
        else:
            prefix = 'v19'+pointing
        print('prefix = ',prefix)
        #command_string = 'python ~/github/HalphaImaging/python3/INT_align_images.py --image1 {} --image2 {} --weight2 {}'.format(haimage,rimage,rweightimage)
        #command_string = 'python  ~/github/halphagui/halphamain.py --virgo --rimage {} --haimage {} --filter {} --psfdir {} --tabledir /home/rfinn/research/Virgo/tables-north/v1/ --prefix {} --auto'.format(rimage,haimage,hfilter,args.psfdir,prefix)

        command_string = 'python  ~/github/halphagui/halphamain.py --virgo --rimage {} --haimage {} --filter {} --psfdir {} --tabledir /home/rfinn/research/Virgo/tables-north/v2/ --prefix {} --auto'.format(rimage,haimage,hfilter,args.psfdir,prefix)

        # check to see if shifted r-band image exists.  if 
        try:
            print('running : ',command_string)
            os.system(command_string)
        except:
            print('##########################################')
            print('WARNING: problem running auto gui on ',rimage)
            print('##########################################')

    #just running on one directory for testing purposes
    #i += 1
    #if i > 0:
    #    break


