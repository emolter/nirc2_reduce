#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import warnings

class Image:
    
    def __init__(self, fname):
        warnings.filterwarnings('ignore','The following header keyword is invalid or follows an unrecognized non-standard convention')
        warnings.filterwarnings('ignore','non-ASCII characters are present in the FITS file header and have been replaced')
        self.hdulist = fits.open(fname,ignore_missing_end=True)
        self.header = self.hdulist[0].header
        self.data = self.hdulist[0].data
        targ = self.header['OBJECT']
        try:
            self.target = targ.split()[0].strip(', \n')
        except:
            self.target = 'Unknown'
        try:
            self.filt = self.header['FWINAME'].strip(', \n')
        except:
            self.filt = 'Unknown'
        
    def plot(self):
        plt.imshow(self.data,origin='lower left')
        plt.show()