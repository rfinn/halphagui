# halphagui
development of gui interface for H-alpha narrowband imaging

My goal is to create a gui interface that will handle the analysis process from making cutouts through surface photometry and radial profiles.  I haven't made a gui before, so this is going to be a learning experience.

My thought is to use qt designer to make the layout, and then to use ginga to handle the image display.  Sounds easy enough, right?  So far, not so much.  Here are some questions I have.

I'm guessing I set up a main window in designer, and then add the subframes that will contain:
* the R-band mosaic image, with the positions of NSA galaxies marked
 * a drop-down menu to select an NSA galaxy to work on (from the list of galaxies in the FOV)
* another image panel that shows the cutout of the galaxy in R and Halpha (tabs)
 * this one should have a slider to allow the user to adjust the size of the cutout image.
 * also shows ellipse used for photometry, with ability to hide
* same as above for continuum subtracted image
 * this one should have a slider that allows the user to adjust the filter ratio to use with the image subtaction.
 * also shows ellipse used for photometry, with ability to hide
* a subframe to show the mask and edit it if necessary
* a subframe of the R radial profiles
* a subframe of the CS radial profiles
* a subframe of the normalized R and CS profiles

So this is 6 subframes!

1 - coadded image, with dropdown menu to select galaxy
2 - Ha, R images with sliders to select size and filter ratio
3 - CS cutout image, with ellipse for measuring photometry (and ability to hide this)
4 - image frame showing mask, with ability to edit mask
5 - R radial profile
6 - CS radial profile
7 - normalized R and CS profile overplotted

Can make a grid with 2 rows and 5 columns, with the coadded image given a 2x2 slot.

so maybe a big frame on the left to hold coadd, then a 2x2 grid on the right

# Workflow #

```
pyuic5 -x halpha.v1.ui -o halphav1.py
python halphamain.py
```

# Requirements ##

## Python Modules (need to add versions) ##
* ginga
* astropy
* pyqt5
* photutils

## Other software ##
* Source Extractor

# Notes #

## some issues I ran into setting up ginga

when installing ginga using pip (anaconda 3.7), the util module is out of date.

* I downloaded the ginga source code from github, and then ran setup.py
* I then copied the util directory to my anaconda installation of ginga
  * cp -r ~/github/ginga/ginga/util ~/anaconda3/lib/python3.7/site-packages/ginga/.
* it's probably better to copy the entire package, but this is what I did for starters.

## some issue with getting qt to run

at first, I wasn't able to run any of the ginga examples, and this was because I hadn't defined an environment variable correctly.  Sorry - I can't remember the details now, but the instructions set to set the environment variable to qt, but it actually needed to be qt5.  

