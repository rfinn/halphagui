#!/usr/bin/env python

"""
GOAL:
* reproject the r-band mask to the wise pixel scale
* this is designed to work with the cutouts from John Moustakas

USEAGE:

reproject_mask.py maskfilename wiseimage


"""

import sys
from astropy.io import fits
from reproject import reproject_interp

maskfile = sys.argv[1]

reffile = sys.argv[2]


hmask = fits.open(maskfile)

href = fits.open(reffile)

# reproject r-band mask onto W3 header

wisemask,footprint = reproject_interp(hmask,href[0].header)

# all wise images have the same pixel scale, so we only need one wise mask
outname = maskfile.replace('r-mask','wise-mask')
fits.writeto(outname,wisemask,href[0].header,overwrite=True)




