#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from .image import Image
from .prettycolors import make_colormap, get_colormap
#from image import Image
from scipy.signal import medfilt, fftconvolve
from scipy.interpolate import RectBivariateSpline
import astroscrappy, sys
from image_registration.chi2_shifts import chi2_shift
from image_registration.fft_tools.shift import shiftnd, shift2d


class Observation:
    
    def __init__(self, fnames):
        '''
        Parameters
        ----------
        fnames : list or str, required. filenames of data frames. 
            if str, assume passing single data frame
        
        Attributes
        ----------
        dummy_fits : 
        frames : 
        
        '''
        if type(fnames) == str:
            fnames = [fnames,]
        self.dummy_fits = Image(fnames[0]) #used to hijack header info
        self.frames = np.asarray([Image(f).data for f in fnames])
        self.subc = self.dummy_fits.hdulist[0].header['NAXIS1']
        self.target = self.dummy_fits.hdulist[0].header['OBJECT']        
        
    
    def apply_flat(self,fname):
        flat = Image(fname).data
        #handle subarrays
        if flat.shape[0] != self.subc:
            ctr = int(flat.shape[0]/2)
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            flat = flat[ll:ul,ll:ul]
        
        for frame in self.frames:
            with np.errstate(divide='ignore',invalid='ignore'):
                frames_flat.append(frame/flat)
        self.frames = np.asarray(frames_flat)
        
        
    
    
    
class Nod(Observation):
    
    def __init__(self, data_fname, sky_fname):
        
        super().__init__(data_fname)
        
        
    
        
           
class Bxy3(Observation):
    
    def __init__(self,fnames):
        '''Input three file names as list. Must be in the conventional order
        where target is in position: top left | bottom right | top right'''
        self.dummy_fits = Image(fnames[0]) #used to hijack header info
        self.frames = np.asarray([Image(f).data for f in fnames])
        self.subc = self.dummy_fits.hdulist[0].header['NAXIS1']
        self.target = self.dummy_fits.hdulist[0].header['OBJECT']
        if self.frames[0].shape[0] != self.subc or self.frames[1].shape[1] != self.subc:
            #subarrays smaller than 512x512 on Keck are not square. chop out center
            frames_square = []
            for frame in self.frames:
                shp = np.min([frame.shape[0], frame.shape[1]])
                ctr = int(shp/2)
                ll = int(ctr - self.subc/2)
                ul = int(ctr + self.subc/2)
                frame = frame[ll:ul,ll:ul]
                frames_square.append(frame)
            self.frames = np.asarray(frames_square)
                              
    
    def make_sky(self, outfile):
        #make the master sky
        self.sky = np.median(self.frames, axis = 0)
        # change some header info and write to .fits
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header['OBJECT'] = 'SKY_FULL'
        hdulist_out[0].data = self.sky
        hdulist_out[0].writeto(outfile, overwrite=True)
        
    def apply_sky(self,fname):
        '''identify patches of "normal" sky in each of the three bxy3 frames. 
        comes up with a median average value of background for each frame. 
        Then you normalize the "full" sky frame we got from median of bxy3 
        to the value of the sky brightness in the single frame.
        Also normalizes everything to be around 1, so can directly multiply by data
        This ONLY works if bxy3 frames are input in the correct (usual) order,
        top left | bottom right | top right'''
        
        full_sky = Image(fname).data #could just as easily use self.sky, but this way could theoretically load a different sky background if desired
        if full_sky.shape[0] != self.subc:
            ctr = int(full_sky.shape[0]/2)
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            full_sky = full_sky[ll:ul,ll:ul]

        sz = int(self.subc/16) #just hardcoded this with a reasonable value
        c0 = (int(self.subc/2 + self.subc/4), int(self.subc/2 - self.subc/4))
        c1 = (int(self.subc/2 - self.subc/4), int(self.subc/2 + self.subc/4))
        c2 = (int(self.subc/2 + self.subc/4), int(self.subc/2 + self.subc/4))
        box0 = [c0[0]-sz,c0[0]+sz,c0[1]-sz,c0[1]+sz]
        box1 = [c1[0]-sz,c1[0]+sz,c1[1]-sz,c1[1]+sz]
        box2 = [c2[0]-sz,c2[0]+sz,c2[1]-sz,c2[1]+sz]
        #take bits of sky where planet is not
        skybits0 = [self.frames[0][box1[0]:box1[1],box1[2]:box1[3]],self.frames[0][box2[0]:box2[1],box2[2]:box2[3]]]
        skybits1 = [self.frames[1][box0[0]:box0[1],box0[2]:box0[3]],self.frames[1][box2[0]:box2[1],box2[2]:box2[3]]]
        skybits2 = [self.frames[2][box0[0]:box0[1],box0[2]:box0[3]],self.frames[2][box1[0]:box1[1],box1[2]:box1[3]]]
        
        frames_skysub = []
        fullmed = np.median(full_sky)
        if np.abs(fullmed) > 5.0: #minimum average counts hardcoded here
            for i in range(3):
                skybits = [skybits0,skybits1,skybits2][i]
                frame = self.frames[i]
                med0 = np.median(skybits[0])
                med1 = np.median(skybits[1])
                med = np.mean([med0, med1])
                norm = med/fullmed
                sky_norm = full_sky*norm
                frames_skysub.append(frame - sky_norm)
            self.frames = np.asarray(frames_skysub)
        else:
            print('Sky subtraction not applied (counts too low for good stats)')
                    
    def apply_flat(self,fname):
        flat = Image(fname).data
        #handle subarrays
        if flat.shape[0] != self.subc:
            ctr = int(flat.shape[0]/2)
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            flat = flat[ll:ul,ll:ul]
        frames_flat = []
        for frame in self.frames:
            with np.errstate(divide='ignore',invalid='ignore'):
                frames_flat.append(frame/flat)
        self.frames = np.asarray(frames_flat)

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
        #ax0.imshow(badpx_map, origin = 'lower')
        #ax1.imshow(smoothed_badpx_map, origin = 'lower')
        #plt.show()
        frames_badpx = []
        for frame in self.frames:
            smoothed = medfilt(frame,kernel_size = 7)
            frame[bad_indices] = smoothed[bad_indices]
            frames_badpx.append(frame)
        self.frames = np.asarray(frames_badpx)
            
    def dewarp(self):
        '''
        Using updated maps from Ghez, Lu galactic center group (Service et al. 2016)
        doi:10.1088/1538-3873/128/967/095004
        Spline interpolation used to interpolate original image. 
        The spline is then evaluated at new pixel locations computed from the map.
        
        Parameters
        ----------
        frame: 2-D np array containing the nirc2 data
        
        Returns
        -------
        dewarped frame
        
        To do
        -----
        relative import of fits files
        
        Notes
        -----
        From Service+16: 
            "These are lookup tables generated by evaluating the fits 
            from the previous section at the center of every pixel on the 
            NIRC2 detector and are the values that should be added to 
            raw NIRC2 positions to shift them to a distortion-free frame"
        So the position offsets in the fits files should be added. 
        That this works right can be checked by looking at difference maps 
            and comparing with the expected vectors in Service+16
        '''
        warpx = fits.getdata('/Users/emolter/Python/nirc2_reduce/nirc2_distort_X_post20150413_v1.fits')
        warpy = fits.getdata('/Users/emolter/Python/nirc2_reduce/nirc2_distort_Y_post20150413_v1.fits')
        #handle subarrays
        if self.subc != 1024: #hardcode because dewarp arrays will always be full size of detector
            ctr = 512
            ll = int(ctr - self.subc/2)
            ul = int(ctr + self.subc/2)
            warpx = warpx[ll:ul,ll:ul]
            warpy = warpy[ll:ul,ll:ul]
        #plt.imshow(warpx,origin='lower')
        #plt.show()
        szx = self.frames[0].shape[0]
        szy = self.frames[0].shape[0]
        xx = np.linspace(0,szx-1,szx)
        yy = np.linspace(0,szy-1,szy)
        x,y = np.meshgrid(xx,yy)
        mapx = x - warpy #check if plus or minus!!
        mapy = y + warpx #check if plus or minus!!
        flatx = mapx.flatten()
        flaty = mapy.flatten()
        frames_dewarp = []
        for frame in self.frames:
            frame_spline = RectBivariateSpline(xx,yy,frame)
            frame_dw = frame_spline.ev(flatx,flaty).reshape(frame.shape).T
            frames_dewarp.append(frame_dw)
            #plt.imshow(frame_dw, origin='lower')
            #plt.show()
        self.frames = np.asarray(frames_dewarp)
        
    def trim(self):
        '''Clips each image in bxy3 to its own quadrant.
        relies on frames input in correct order'''
        sz = int(self.subc/4)
        c0 = (int(self.subc/2 + self.subc/4), int(self.subc/2 - self.subc/4))
        c1 = (int(self.subc/2 - self.subc/4), int(self.subc/2 + self.subc/4))
        c2 = (int(self.subc/2 + self.subc/4), int(self.subc/2 + self.subc/4))
        box0 = [c0[0]-sz,c0[0]+sz,c0[1]-sz,c0[1]+sz]
        box1 = [c1[0]-sz,c1[0]+sz,c1[1]-sz,c1[1]+sz]
        box2 = [c2[0]-sz,c2[0]+sz,c2[1]-sz,c2[1]+sz]
        ftrim0 = self.frames[0][box0[0]:box0[1],box0[2]:box0[3]]
        ftrim1 = self.frames[1][box1[0]:box1[1],box1[2]:box1[3]]
        ftrim2 = self.frames[2][box2[0]:box2[1],box2[2]:box2[3]]
        self.frames = np.asarray([ftrim0,ftrim1,ftrim2])
        
    def remove_cosmic_rays(self):
        '''Detects cosmic rays using the astroscrappy package.'''
        frames_cosray = []
        for frame in self.frames:
            crmask, cleanarr = astroscrappy.detect_cosmics(frame, cleantype='medmask')
            #fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize=(15,6))
            #ax0.imshow(frame,origin = 'lower')
            #ax1.imshow(crmask,origin = 'lower')
            #ax2.imshow(cleanarr,origin = 'lower')
            #plt.show()
            frames_cosray.append(cleanarr)
        self.frames = np.asarray(frames_cosray)         

    def per_second(self):
        '''Changes units to counts/second'''
        header = self.dummy_fits.hdulist[0].header
        tint = header['ITIME']
        coadd = header['COADDS']
        persec_frames = []
        for frame in self.frames:
            persec_frames.append(frame/(tint*coadd))
        self.frames = np.asarray(persec_frames)
    
    def apply_photometry_frames(self,flux_per):
        '''Simply multiplies each frame by flux density / (cts/s)'''
        self.frames = np.asarray([frame*flux_per for frame in self.frames])
        
    def apply_photometry_final(self,flux_per):
        '''Simply multiplies final image by flux density / (cts/s)'''
        self.final = self.final*flux_per
            
    def stack(self):
        '''Cross-correlate the images applying sub-pixel shift.
        Shift found using DFT upsampling method as written by 
        Stack them on top of each other to increase SNR.'''
        shifted_data = [self.frames[0]]
        for frame in self.frames[1:]:
            [dx,dy,dxerr,dyerr] = chi2_shift(self.frames[0],frame)
            #error is nonzero only if you include per-pixel error of each image as an input. Should eventually do that, but no need for now.
            shifted = shift2d(frame,-1*dx,-1*dy)
            shifted_data.append(shifted)
            
        self.final = np.median(shifted_data,axis=0)
    
    def crop(self,bw):
        '''Applies crop to the final image. bw is border width'''
        szx,szy = self.final.shape[0],self.final.shape[1]
        self.final = self.final[bw:szx-bw,bw:szy-bw]
        
    def write_frames(self,outfiles):
        '''Write each individual frame to .fits at any step in the process.'''
        for i in range(len(self.frames)):
            frame = self.frames[i]
            outfile = outfiles[i]
            hdulist_out = self.dummy_fits.hdulist
            hdulist_out[0].data = frame
            hdulist_out[0].writeto(outfile, overwrite=True) 
            
    def plot_frames(self):
        fig, axes = plt.subplots(1,3, figsize = (9,5))
        for i in range(len(axes)):
            ax = axes[i]
            ax.imshow(self.frames[i], origin = 'lower')
            ax.set_title('Frame %d'%i)
            ax.set_xticks([])
            ax.set_yticks([])
        plt.show()
        
    def write_final(self,outfile,png = False,png_file=''):
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header['OBJECT'] = self.target+'_STACKED'
        hdulist_out[0].data = self.final
        hdulist_out[0].writeto(outfile, overwrite=True)
        if png:
            fig, ax = plt.subplots(1,1, figsize=(8,8))
            cmap = get_colormap(self.target)
            ax.imshow(self.final, cmap = cmap, origin='lower')
            plt.tight_layout()
            fig.savefig(png_file,bbox='None')
            plt.close()
            
    def help(self):
        
        helpstr = '''
        Contains tasks for reducing data taken with the Keck bxy3 dither
        
        Functions (see DOCUMENTATION.py for use):
            make_sky(self,outfile)
            apply_sky(self,fname)
            apply_flat(self,fname)
            apply_badpx_map(self,fname)
            dewarp(self)
            trim(self)
            remove_cosmic_rays(self)
            per_second(self)
            apply_photometry_frames(self,flux_per)
            apply_photometry_final(self,flux_per)
            stack(self)
            crop(self,bw)
            write_frames(self,outfiles)
            plot_frames(self)
            write_final(self,outfile,png = False,png_file='')
             
        Attributes:
            dummy_fits: Image object used to hijack header info
            frames: the three grids of image data. this gets updated in every step
            subc: size of image array
            target: name of planet
            sky: master sky map produced by make_sky
        '''
        print(helpstr)

