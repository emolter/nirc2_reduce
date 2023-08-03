#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from .image import Image
from scipy.interpolate import RectBivariateSpline
from scipy.signal import medfilt
import astroscrappy


class Nod:
    
    def __init__(self, skyf, imagef):
        '''Input two file names, a science image and a sky.'''
        if type(skyf) != str and type(skyf) != np.str_:
            if len(skyf) == 1:
                skyf = skyf[0]
            else:
                print('ERROR: passed 2 or more sky files into nod.py')
        if type(imagef) != str and type(imagef) != np.str_:
            if len(imagef) == 1:
                imagef = imagef[0]
            else:
                print('ERROR: passed 2 or more image files into nod.py')
            
        self.image = Image(imagef)
        self.data = self.image.data
        self.sky = Image(skyf).data
        self.subc = self.image.hdulist[0].header['NAXIS1']
        if self.subc == 256:
            self.data = self.data[4:-4,:]
            self.sky = self.sky[4:-4,:]
        self.target = self.image.hdulist[0].header['OBJECT']
        
        #handle case where sky is larger than image subc
        if self.sky.shape[0] != self.subc:
            ctr = int(self.sky.shape[0]/2)
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            self.sky = self.sky[ll:ul,ll:ul]
        
    def apply_sky(self):
        self.data = self.data - self.sky
        
    def apply_flat(self,fname):
        flat = Image(fname).data
        #handle subarrays
        if flat.shape[0] != self.subc:
            ctr = int(flat.shape[0]/2)
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            flat = flat[ll:ul,ll:ul]
        
        with np.errstate(divide='ignore',invalid='ignore'):
            self.data = self.data / flat
            
    def apply_badpx_map(self,fname):
        '''Replaces all bad pixels in input map with the median of the pixels
        around it.'''
        badpx_map = Image(fname).data
        #handle subarrays
        if badpx_map.shape[0] != self.subc:
            ctr = int(badpx_map.shape[0]/2)
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            badpx_map = badpx_map[ll:ul,ll:ul]
        bad_indices = np.where(badpx_map == 0)
        ## test case - kernel of 7 shown to remove all bad pixels including weird stripe in bottom left of detector
        #smoothed_badpx_map = medfilt(badpx_map,kernel_size = 7) 
        #fig, (ax0,ax1) = plt.subplots(1,2, figsize=(10,6))
        #ax0.imshow(badpx_map, origin = 'lower left')
        #ax1.imshow(smoothed_badpx_map, origin = 'lower left')
        #plt.show()
        
        smoothed = medfilt(self.data,kernel_size = 7)
        self.data[bad_indices] = smoothed[bad_indices]

    def dewarp(self):
        '''Using updated maps from Jessica Lu galactic center group (Service et al. 2016).
        Spline interpolation used to interpolate original image. The spline
        is then evaluated at new pixel locations computed from the map.
        '''
        warpx = fits.getdata('/Users/emolter/ned_laptop/Python/nirc2_reduce/nirc2_distort_X_post20150413_v1.fits')
        warpy = fits.getdata('/Users/emolter/ned_laptop/Python/nirc2_reduce/nirc2_distort_Y_post20150413_v1.fits')
        #handle subarrays
        if self.subc != 1024: #hardcode because dewarp arrays will always be full size of detector
            ctr = 512
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            warpx = warpx[ll:ul,ll:ul]
            warpy = warpy[ll:ul,ll:ul]
        #plt.imshow(warpx,origin='lower left')
        #plt.show()
        szx = self.data.shape[0]
        szy = self.data.shape[0]
        xx = np.linspace(0,szx-1,szx)
        yy = np.linspace(0,szy-1,szy)
        x,y = np.meshgrid(xx,yy)
        mapx = x - warpx #check if plus or minus!!
        mapy = y - warpy #check if plus or minus!!
        flatx = mapx.flatten()
        flaty = mapy.flatten()
        
        spline = RectBivariateSpline(xx,yy,self.data)
        self.data = spline.ev(flatx,flaty).reshape(self.data.shape).T

    def remove_cosmic_rays(self):
        '''Detects cosmic rays using the astroscrappy package.'''
        crmask, cleanarr = astroscrappy.detect_cosmics(self.data, cleantype='medmask')
        #fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize=(15,6))
        #ax0.imshow(frame,origin = 'lower left')
        #ax1.imshow(crmask,origin = 'lower left')
        #ax2.imshow(cleanarr,origin = 'lower left')
        #plt.show()
        self.data = cleanarr

    def per_second(self):
        '''Changes units to counts/second'''
        header = self.image.hdulist[0].header
        tint = header['ITIME']
        coadd = header['COADDS']
        self.data = self.data/(tint*coadd)
        
    def apply_photometry(self,flux_per):
        '''Simply multiplies each frame by flux density / (cts/s)'''
        self.data = self.data*flux_per
        
    def crop(self,bw):
        '''Applies crop to the final image. bw is border width'''
        szx,szy = self.data.shape[0],self.data.shape[1]
        self.data = self.data[bw:szx-bw,bw:szy-bw]
        
    def uranus_crop(self, bw):
        '''Custom crop for Uranus data. We use the right side of the
        NIRC2 detector to avoid the bad pixels in the lower left corner
        so just cut that out'''
        szx,szy = self.data.shape[0],self.data.shape[1]
        self.data = self.data[bw:szx-bw,2*bw:]        
        
    def write(self,outfile,png = False,png_file=''):
        hdulist_out = self.image.hdulist
        hdulist_out[0].header['OBJECT'] = self.target+'_REDUCED'
        hdulist_out[0].data = self.data
        hdulist_out[0].writeto(outfile, overwrite=True)
        if png:
            fig, ax = plt.subplots(1,1, figsize=(8,8))
            ax.imshow(self.data,origin='lower left')
            plt.tight_layout()
            fig.savefig(png_file,bbox='None')
            plt.close()
            
    def help(self):
        
        helpstr = '''
        Contains tasks for reducing data taken with a target/sky nod
        
        Functions (see DOCUMENTATION.py for use):
            apply_sky(self,fname)
            apply_flat(self,fname)
            apply_badpx_map(self,fname)
            dewarp(self)
            remove_cosmic_rays(self)
            per_second(self)
            apply_photometry(self,flux_per)
            crop(self,bw)
            uranus_crop(self, bw)
            write(self,outfile,png = False,png_file='')
             
        Attributes:
            image: Image object used to hijack header info
            data: the image data. this gets updated in every step
            sky: the sky data
            subc: size of image array
            target: name of planet
        '''
        
        print(helpstr)