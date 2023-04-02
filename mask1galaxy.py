#!/usr/bin/env python

"""
GOAL:
* to create a mask from the r-band image of a galaxy
* to reproject that mask onto the WISE pixel scale

USAGE:

~/github/halphagui/mask1galaxy.py VFID_dirname


NOTES:

when running on vms, the topdir is:
/mnt/astrophysics/muchogalfit-output

when running on grawp, the topdir is:
/mnt/astrophysics/rfinn/muchogalfit-output

this switch is handled automatically

"""
import os
import sys
#import glob
import numpy as np
import subprocess



def funpack_image(input,output,nhdu=1):
    command = 'funpack -O {} {}'.format(output,input)
    print(command)
    os.system(command)


if __name__ == '__main__':
    homedir = os.getenv("HOME")
    topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
    image_source_dir = '/mnt/astrophysics/virgofilaments-data/'
    try:
        os.chdir(topdir)
    except FileNotFoundError: # assuming that we are running on virgo vms
        topdir = '/mnt/astrophysics/muchogalfit-output/'
        image_source_dir = '/mnt/virgofilaments-data/'        
        os.chdir(topdir)
    # take as input the galaxy name
    galname = sys.argv[1]

    # move to muchogalfit-output directory
    output_dir = topdir+galname+'/'
    if not os.path.exists(output_dir):
        print('WARNING: {} does not exist\n Be sure to run setup_galfit.py first')
        os.chdir(topdir)
        sys.exit()
    
    os.chdir(output_dir)

    # read in input file to get galname, objname, ra, dec, and bandpass
    sourcefile = open(galname+'sourcelist','r')
    galaxies = sourcefile.readlines()
    if len(galaxies) > 1:
        # set the flag to have more than one galaxy in the galfit input file
        multiflag = True
    elif len(galaxies) == 1:
        multiflag = False
    else:
        print('Problem reading sourcelist for {}'.format(galname))
        print('Please check the setup directory')
        os.chdir(topdir)
        sys.exit()

    # parse information from file
    vfid, objname, ra, dec, bandpass = galaxies[0].rstrip().split()
    ra = float(ra)
    dec = float(dec)

    # set up path name for image directory
    # directory where galaxy images are
    data_dir = f"{image_source_dir}/{int(ra)}/{objname}/"





    # construct image name
    image = f"{objname}-custom-image-r.fits"
    # check if r-band image exists
    if not os.path.exists(image):
        # if not, copy from john's directories and funpack
        funpack_image(os.path.join(data_dir,image),os.path.join(output_dir,image.replace('.fz','')))
        # get information from the image, like the image size and position of gal
        image = image.replace('.fz','')
        



    # call maskwrapper.py
    cmd = f"{homedir}+'/github/halphagui/maskwrapper.py --image {image} --auto"
    os.system(cmd)

    # define the mask name created by maskwrapper.py

    mask = image.replace('.fits','-mask.fits')
    # reproject mask onto wise image
    
    # construct WISE W3 image name
    wimage = f"{objname}-custom-image-W3.fits"
    # check if r-band image exists
    if not os.path.exists(wimage):
        # if not, copy from john's directories and funpack
        funpack_image(os.path.join(data_dir,wimage),os.path.join(output_dir,wimage.replace('.fz','')))
        # get information from the image, like the image size and position of gal
        wimage = wimage.replace('.fz','')
        
    # call reproject_mask.py
    cmd = f"python {homedir}+'/github/halphagui/reproject_mask.py {mask} {wimage}"
    subprocess.call(cmd)
