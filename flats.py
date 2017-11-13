#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from .image import Image
#from image import Image

     
class Flats:
    
    def __init__(self,fnames_off,fnames_on):
        self.dummy_fits = Image(fnames_on[0]) #use this to hijack header info
        self.frames_off = np.asarray([Image(f).data for f in fnames_off])
        self.frames_on = np.asarray([Image(f).data for f in fnames_on])

        off = np.median(self.frames_off,axis=0)
        on = np.median(self.frames_on,axis=0)
        flat = on - off
        self.flat = flat/np.median(flat)
    
    def write(self,outfile):
        '''subtract lights off from lights on,
        divide by median value to center it around 1,
        write it out'''
        # change some header info and write to .fits
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header['OBJECT'] = 'DOME_FLAT_MASTER'
        hdulist_out[0].data = self.flat
        hdulist_out[0].writeto(outfile, overwrite=True)
        
    def plot(self):
        plt.imshow(self.flat, origin = 'lower left')
        plt.show()
    
    def make_badpx_map(self, outfile, tol, blocksize):
        '''Find pixels whose values are very far from the average of their neighbors
        tol is the fraction different they can be
        blocksize is the number of pixels in each direction over which to average'''
        badpx_map = np.ones(self.flat.shape)
        for i in range(0,self.flat.shape[0]+blocksize,blocksize):
            for j in range(0,self.flat.shape[1]+blocksize,blocksize):
                flatblock = self.flat[i:i+blocksize,j:j+blocksize]
                mapblock = badpx_map[i:i+blocksize,j:j+blocksize]
                med = np.median(flatblock)
                
                #if not within tolerance, set to NaN
                mapblock[np.where(flatblock/med > 1 + tol)] = 0
                mapblock[np.where(flatblock/med < 1 - tol)] = 0
                badpx_map[i:i+blocksize,j:j+blocksize] = mapblock
        self.badpx_map = badpx_map
        # change some header info and write to .fits
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header['OBJECT'] = 'BADPX_MAP'
        hdulist_out[0].data = self.badpx_map
        hdulist_out[0].writeto(outfile, overwrite=True)
        
