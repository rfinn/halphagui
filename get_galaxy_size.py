
'''
INPUT
* image
* coordinates of galaxies (gximage,gyimage)

OUTPUT:
bbox_xmin,bbox_ymax,bbox_ymin,bbox_ymax

PROCEDURE:
* read in image
* run photutils to get segmentation map
* get value of segmentation image at the position of galaxies
* return size of box, scaled by a factor
# 

'''

from photutils import detect_threshold
from photutils import detect_sources
from photutils import source_properties
from astropy.visualization import simple_norm
from astropy.io import fits
import numpy as np
import os


def display_image(image,percent=99.5):
    norm = simple_norm(image,stretch='asinh',percent=percent)            
    plt.imshow(image,origin='upper',cmap='gray_r', norm=norm)
    pass


def getsegmentation(image):
    ''' 
    pass in image name, array of x coord, array of y coord.  
    Returns object size 
    '''
    # read in image
    data = fits.getdata(image)
    threshold = detect_threshold(data, nsigma=2.)
    # create segmentation map
    segm = detect_sources(data, threshold, npixels=25)

    # skipping deblending for now
    
    # create catalog
    cat = source_properties(data, segm)
    tbl = cat.to_table()    
    return segm.data, tbl

def getsegmentation_se(image):
    #print 'cp ' +args.d + '/default.* .'
    os.system('ln -s ~/github/halphagui/astromatic/default.param .')
    os.system('ln -s ~/github/halphagui/astromatic/default.se.objsize .')
    os.system('ln -s ~/github/halphagui/astromatic/default.conv .')
    os.system('ln -s ~/github/halphagui/astromatic/default.nnw .')            
    #files = sorted(glob.glob(args.filestring))
    # making output files have the 
    command = f'sex {image} -c default.se.objsize -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME {os.path.basename(image)}-segmentation.fits -CATALOG_NAME {os.path.basename(image)}_cat.fits -CATALOG_TYPE FITS_LDAC'
    print(command)
    os.system(command)
    segm = fits.getdata(f'{image}-segmentation.fits')
    tbl = fits.getdata(f'{image}_cat.fits',2)
    return segm,tbl

def getobjectsize(image,xobj,yobj,scale=1.75,plotflag=False,usese=True):
    if usese:
        segm, tbl = getsegmentation_se(image)
    else:
        segm, tbl = getsegmentation(image)
    # get catalog row for objects

    matchflag = np.ones(len(xobj),'bool')
    objids = segm[yobj,xobj]
    print(objids)
    if usese:
        xminkey= 'XMIN_IMAGE'
        xmaxkey= 'XMAX_IMAGE'
        yminkey= 'YMIN_IMAGE'
        ymaxkey= 'YMAX_IMAGE'
    else:
        xminkey= 'bbox_xmin'
        xmaxkey= 'bbox_xmax'
        yminkey= 'bbox_ymin'
        ymaxkey= 'bbox_ymax'
        
    # get max dimension of bbox_ymax - bbox_ymin or bbox_xmax-bbox_xmin
    ymin = tbl[yminkey][objids-1]
    ymax = tbl[ymaxkey][objids-1]
    dy = ymax-ymin
    
    xmin = tbl[xminkey][objids-1]
    xmax = tbl[xmaxkey][objids-1]
    dx = xmax-xmin
    print(xmin,xmax)
    print(dx)
    print(ymin,ymax)
    print(dy)
    objsize = np.zeros(len(dx),'i')
    for i in range(len(objsize)):
        if usese:
            objsize[i] = int(round(max(dx[i],dy[i])))
        else:
            objsize[i] = int(round(max(dx[i].value,dy[i].value)))
    # scale max dimension by a scalefactor
    objsize = objsize*scale

    # return scaled size
    return objsize


                         
