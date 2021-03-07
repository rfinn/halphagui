#!/usr/bin/env python

'''
GOAL: 
- create psf images for all the coadds

Run this from, e.g. /home/rfinn/data/reduced/psf-images/

Change the coadd_image_directory as needed before running


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib
import time
matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
telescope = 'INT'
working_dir = os.getcwd()

# overwrite output files if they exist
overwrite = True

import argparse

parser = argparse.ArgumentParser(description ='Run buildpsf.py for all images in the specified directory')
parser.add_argument('--coaddir',dest = 'coaddir', default ='/home/rfinn/data/reduced/virgo-coadds-feb2019-int/', help = 'directory for coadds. Default is /home/rfinn/data/reduced/virgo-coadds/feb2019-int/')
parser.add_argument('--int', dest = 'int', default = False,action='store_true', help = 'set this for INT data')

args = parser.parse_args()

coadd_image_directory = args.coaddir




filters = ['r','Halpha','Ha6657','ha4','R']
#saturate_level = [100,30,30]

for i,f in enumerate(filters):
    # get list of current directory
    # this will grab the coadds but not the weight images
    flist1 = glob.glob(coadd_image_directory+'VF-*-'+f+'.fits')
    print(i,f,len(flist1))
    # changing this to just do pointing 022 and 026
    #flistp20 = glob.glob(coadd_image_directory+'VF-*p017*-'+f+'.fits')
    #flistp26 = glob.glob(coadd_image_directory+'VF-*p072*-'+f+'.fits')
    #flist1 = flistp20+flistp26
    flist1.sort()
    for fimage in flist1: # loop through list

        print('##########################################')
        print('##########################################')        
        print('BUILDING PSF FOR: ',fimage)
        print('##########################################')
        print('##########################################')
        start_time = time.perf_counter()
        # adding saturation limit for normalized images
        # I estimated this from the r-band image for p001
        # this is in counts/sec
        if args.int:
            command_string = 'python ~/github/halphagui/buildpsf.py --image {} --int'.format(fimage)
        else:
            command_string = 'python ~/github/halphagui/buildpsf.py --image {} '.format(fimage)
        try:
            print('running : ',command_string)
            os.system(command_string)
        except:
            print('##########################################')
            print('WARNING: problem running buildpsf.py for ',fimage)
            print('##########################################')

        # clock time to run buildpsf
        end_time = time.perf_counter()
        print('\t total time = ',end_time - start_time)

        # just running on one directory for testing purposes
        #break
    #break



