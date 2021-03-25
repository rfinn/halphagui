#!/usr/bin/env python

import os
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from astropy.io import ascii

wave_halpha = 6563. # angstrom

class filter_trace():
    def __init__(self,hafilter,filterpath=None,mintrans=10.):
        '''
        hafilter can be 4, 8, 12, 16

        filter path should point to filter_trace subdirectory of github/halphagui

        mintrans is the minimum transmission to use for selecting galaxies
        '''
        # read in filter trace
        # assume github directory is off main dir
        self.hafilter = hafilter
        if filterpath == None:
            self.filterpath = os.getenv('HOME')+'/github/halphagui/filter_traces/'
        else:
            self.filterpath = filterpath

        wavescale = 1
        if hafilter == 'inthalpha':
            filterfile = self.filterpath+'/wfc-int-197-halpha.txt'
            # wavelength is in nm so scale by 10
            wavescale=10
            pass
        elif hafilter == 'intha6657':
            filterfile = self.filterpath+'/wfc-int-227-ha6657.txt'
            wavescale=10
            pass
        else:
            filterfile = self.filterpath+'/ha'+str(hafilter)+'-sim04.txt'
        tab = ascii.read(filterfile)
        self.wave = tab['col1']*wavescale # wavelength in angstrom
        self.trans = tab['col2'] # transmission percent
        self.maxtrans = np.max(self.trans)
        # get wavelengths where transmission crosses 10 percent level
        ids = np.where(self.trans > 10.)
        # calculate min and max redshifts that correspond to transmission cut
        self.minz_trans10 = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans10 = (self.wave[ids[0][-1]]/wave_halpha -1.)
        ids = np.where(self.trans > 30.)
        self.minz_trans30 = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans30 = (self.wave[ids[0][-1]]/wave_halpha -1.)
        ids = np.where(self.trans > 50.)
        self.minz_trans50 = (self.wave[ids[0][0]]/wave_halpha -1.)
        self.maxz_trans50 = (self.wave[ids[0][-1]]/wave_halpha -1.)
        self.spline_fit()
    def spline_fit(self):
        
        # create spline fit
        self.spline_fit = interpolate.splrep(self.wave, self.trans)

    def get_transmission(self, wave):
        '''
        INTPUT:
        - wavelength, either individual value or array

        OUTPUT:
        - transmission (spline fit to transmission curve at each wavelength
        - fitflag - True if wavelength is within the filter trace, false otherwise
        '''
        # make sure that the wavelength is in the right range
        if len(wave) > 1: #array
            #print('input is an array')
            transmission = np.zeros(len(wave),'f')
            self.fitflag = (wave > np.min(self.wave)) & (wave < np.max(self.wave))
            #print(self.fitflag)
            if sum(self.fitflag) == 0:
                print('all wavelengths out of range')
                return None, False
            elif sum(self.fitflag) < len(wave):
                print('WARNING: some galaxies are outside the filter window')
                print('wavelength out of range.  needs be between %.1f and %.1f Angstrom'%( np.min(self.wave),np.max(self.wave)))
            #print(self.fitflag)
            #print(wave)
            transmission[self.fitflag] = interpolate.splev(wave[self.fitflag],self.spline_fit)
            return transmission, self.fitflag
        elif (wave < np.min(self.wave)) | (wave > np.max(self.wave)): # check wavelength for a is single value
            print('wavelength out of range.  needs be between %.1f and %.1f Angstrom'%( np.min(self.wave),np.max(self.wave)))
            return None, False
        else:
            # function = return transmission for a given wavelength
            return interpolate.splev(wave,self.spline_fit), True
    def get_trans_correction(self, redshift):
        # return ratio of max filter transmission
        # to transmission at that wavelength
        # wave can be a single value or an array
        wave = (redshift+1)*wave_halpha
        transmission, flag = self.get_transmission(wave)
        #self.test = transmission
        correction = np.zeros(len(transmission),'f')
        correction = self.maxtrans/transmission
        plt.figure()
        plt.plot(self.wave, self.trans/10,'k-')
        plt.hist(wave, bins=5)
        plt.xlim((self.minz_trans10+1)*wave_halpha-50,(self.maxz_trans10+1)*wave_halpha+50)
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Transmission %/10')
        titlestring = 'Halpha Filter = {}'.format(self.hafilter)
        plt.title(titlestring)
        #plt.show()
        plt.savefig('galaxies_in_filter.png')
        return correction



