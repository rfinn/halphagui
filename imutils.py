#!/usr/bin/env python


def get_pixel_scale(imheader):
    ''' takes in image header and returns the pixel scale in arcsec  '''
    from astropy.wcs import WCS
    import astropy.units as u
    # get pixel scale from image header
    # convert from degrees/pix to arcsec/pix
    
    ## making more general - not all images have CD1_1 keyword
    #self.pscale = abs(float(self.image_header['CD1_1'])*3600)

    # found a better way to get the pixel scale
    
    image_wcs = WCS(imheader)        
    pscalex,pscaley = image_wcs.proj_plane_pixel_scales()
    
    pscale = pscalex.to(u.arcsec).value # convert deg/pix to arcsec/pix

    return pscale
