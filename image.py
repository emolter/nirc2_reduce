#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import warnings

class Image:
    
    def __init__(self, fname):
        warnings.filterwarnings('ignore','The following header keyword is invalid or follows an unrecognized non-standard convention')
        warnings.filterwarnings('ignore','non-ASCII characters are present in the FITS file header and have been replaced')
        warnings.filterwarnings('ignore','Header block contains null bytes instead of spaces for padding, and is not FITS-compliant')
        self.hdulist = fits.open(fname,ignore_missing_end=True)
        self.header = self.hdulist[0].header
        self.data = self.hdulist[0].data
        targ = self.header['OBJECT']
        try:
            self.target = targ.split()[0].strip(', \n').capitalize()
        except:
            self.target = 'Unknown'
        try:
            fwi = self.header['FWINAME'].strip(', \n')
            fwo = self.header['FWONAME'].strip(', \n')
            if fwo == 'clear' or fwo == 'PK50_1.5':
                self.filt = fwi
            elif fwi == 'clear' or fwi == 'PK50_1.5':
                self.filt = fwo
        except:
            self.filt = 'Unknown'
        
    def plot(self):
        plt.imshow(self.data,origin='lower left')
        plt.show()