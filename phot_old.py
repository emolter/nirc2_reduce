#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from .bxy3 import Bxy3
from . import image
from .nod import Nod
#from bxy3 import Bxy3
#import image

'''Routines for handling photometry in bxy3 configuration'''

def get_filt_info(filt):
    '''Helper to get_standard_mag. Contains dictionary of filter info.
    key "vega" is value coverting Vega mags into erg s-1 cm-2 um-1
    key "lambda" is central wavelength in um
    key "dlambda" is bandpass width in um
    Still need to make these for the other filters'''
    filt_info = {'h':{'vega':1.15641e-06, 'lambda': 1.633,'dlambda':0.296},
                 'hcont':{'vega':1.29E-6, 'lambda': 1.5804,'dlambda':0.0232},
                 'feii':{'vega':1.070E-6, 'lambda': 1.6455,'dlambda':0.0256},
                 'j':{'vega':2.933E-6, 'lambda': 1.248,'dlambda':0.163},
                 'pabeta':{'vega':2.574E-6, 'lambda': 1.2903,'dlambda':0.0193},
                 'jcont':{'vega':3.252E-6, 'lambda': 1.2132,'dlambda':0.0198},
                 'h210':{'vega':4.37013e-07, 'lambda': 2.1281,'dlambda':0.0342},
                 'kcont':{'vega':3.43553e-07, 'lambda': 2.2706,'dlambda':0.0296},
                 'ch4_short':{'vega':1.216E-6, 'lambda': 1.5923,'dlambda':0.1257},
                 'ch4_s':{'vega':1.216E-6, 'lambda': 1.5923,'dlambda':0.1257},
                 'kp':{'vega':4.417e-07, 'lambda': 2.124,'dlambda':0.351},
                 'lp':{'vega':5.411e-08, 'lambda': 3.776,'dlambda':0.700},
                 'ms':{'vega':2.231e-08, 'lambda': 4.670,'dlambda':0.241},
                 }
    return filt_info[filt]['vega'],filt_info[filt]['lambda'],filt_info[filt]['dlambda']

def load_standard(infile):
    '''Helper to get_standard_mag.
    load standard star in Keck starlist format, set variables for mags'''
    with open(infile,'r') as f:
        line = f.readline()
        objname = line[0:15].strip(' ')
        ra = line[16:27]
        dec = line[28:39]
        epoch = line[40:46]
        l = line[46:].split()
        pairs = [val.split('=') for val in l]
        mags = {}
        for val in pairs:
            if val[0].endswith('mag'):
                mags[val[0][0]] = float(val[1])
        return mags

class bxy3Phot(Bxy3):
    
    def __init__(self,filt,fnames, load = False):
        #handles loading already-reduced images of photometric standard - just skip bxy3_reduce
        Bxy3.__init__(self,fnames)
        self.filt = filt
    
    def reduce(self,skyfname,flatfname,badpxfname,outfnames):
        
        self.apply_sky(skyfname)
        self.apply_flat(flatfname)
        self.apply_badpx_map(badpxfname)
        self.dewarp()
        self.trim()
        self.remove_cosmic_rays()
        # at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; need to check if those exist in flats or not
        self.per_second()
        self.write_frames(outfnames)
    
    def get_mag(self, infile):
        '''Takes standard star in Keck format, find its magnitude in the given
        filter taking into account its color.
        Does not work! Ignore for now!'''
        johnson = {'j':{'lambda':1.235},
                   'h':{'lambda':1.662},
                   'k':{'lambda':2.159},
                   'lp':{'lambda':3.776},
                   'ms':{'lambda':4.68}
                   } #lp and ms are not Johnson - they're 2MASS
        mags = load_standard(infile)
        try:
            #if observation is in regular j, h, or k filter, 
            johnson[self.filt.lower()]
            self.star_mag = mags[self.filt.lower()]
        except:
            vega, w, dw = get_filt_info(self.filt.lower())
            
            #find the two Johnson filters this one is between
            if w < johnson['j']['lambda']:
                sys.exit('ERROR: could not compute color because wavelength too short!')
            elif johnson['j']['lambda'] < w < johnson['h']['lambda']:
                filt_l = 'j'
                filt_u = 'h'
            elif johnson['h']['lambda'] < w < johnson['k']['lambda']:
                filt_l = 'h'
                filt_u = 'k'
            elif johnson['k']['lambda'] < w < johnson['lp']['lambda']:
                filt_l = 'k'
                filt_u = 'lp'
            elif johnson['lp']['lambda'] < w < johnson['ms']['lambda']:
                filt_l = 'lp'
                filt_u = 'ms'
            elif johnson['ms']['lambda'] < w:
                sys.exit('ERROR: could not compute color because wavelength too long!')
   
            #linearly interpolate between the filters
            color = mags[filt_l] - mags[filt_u]
            delta_w = w - johnson[filt_l]['lambda']
            delta_w_johnson = johnson[filt_u]['lambda'] - johnson[filt_l]['lambda']
            self.star_mag = mags[filt_l] + color*(delta_w/delta_w_johnson)
                
    def find_cts(self, dist = 50, plot = True):
        '''Get photons/s coming from the entire star
        dist specifies how far from star to integrate up the flux
        Plots will appear to ensure dist is the correct value'''
        star_cts = []
        for frame in self.frames:
            star_loc = np.where(frame == np.max(frame))
            
            #make circular aperture
            xx = np.arange(frame.shape[0]) - star_loc[0]
            yy = np.arange(frame.shape[1]) - star_loc[1]
            x, y = np.meshgrid(xx,yy)
            mask = np.zeros(frame.shape)
            mask[np.where(np.sqrt(x**2 + y**2) <= dist)] = 1.0
            
            #simple square aperture
            box = [star_loc[0] - dist, star_loc[0] + dist, star_loc[1] - dist, star_loc[1] + dist]
            #ct = np.sum(frame[box[0][0]:box[1][0],box[2][0]:box[3][0]])
            ct = np.sum(frame*mask)
            star_cts.append(ct)
            
            if plot:
                #plot a cut across the star to make sure box encapsulates all flux
                cut_x = frame[:,star_loc[1]]
                cut_y = frame[star_loc[0],:]
                fig, (ax0, ax1) = plt.subplots(2,1, figsize=(8,8))
                ax0.plot(cut_x.flatten())
                ax0.axvline(box[0], linestyle = ':', color = 'k')
                ax0.axvline(box[1], linestyle = ':', color = 'k')
                ax0.axhline(linestyle = ':', color = 'r')
                ax1.plot(cut_y.flatten())
                ax1.axvline(box[2], linestyle = ':', color = 'k', label = 'Bounding Box')
                ax1.axvline(box[3], linestyle = ':', color = 'k')
                ax1.axhline(linestyle = ':', color = 'r')
                ax0.set_xlabel('X position')
                ax0.set_ylabel('Counts')
                ax1.set_xlabel('Y position')
                ax1.set_ylabel('Counts')
                ax1.legend()
                plt.show()
        self.star_cts = star_cts
              
    def find_flux_conversion(self):
        '''Given star flux in cts/s and star mag (Vega) in a given filter, 
        calculate conversion between image counts and flux density units'''
        flux_vega, w, dw = get_filt_info(self.filt)
        flux_star = flux_vega * 10**(-(1/2.5)*self.star_mag)
        flux_per =  [flux_star/ct for ct in self.star_cts] #units erg s-1 cm-2 um-1 / cts s-1
        self.flux_per = np.asarray(flux_per)
        self.median_flux_per = np.median(self.flux_per)




class nodPhot(Nod):
    
    def __init__(self,filt,skyf,imagef, load = False):
        #load = True is for already-reduced images of photometric standard - just skip bxy3_reduce
        Nod.__init__(self,skyf,imagef)
        self.filt = filt
    
    def reduce(self, flatfname,badpxfname,outfname):
        self.apply_sky()
        self.apply_flat(flatfname)
        self.apply_badpx_map(badpxfname)
        self.dewarp()
        #self.crop(50)
        self.remove_cosmic_rays()
        # at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; need to check if those exist in flats or not
        self.per_second()
        self.write(outfname)      
                
    def find_cts(self, dist, plot = True):
        '''Get photons/s coming from the entire star
        dist specifies how far from star to integrate up the flux
        Plots will appear to ensure dist is the correct value'''
        
        star_loc = np.where(self.data == np.max(self.data))
        
        #make circular aperture
        xx = np.arange(self.data.shape[0]) - star_loc[0]
        yy = np.arange(self.data.shape[1]) - star_loc[1]
        x, y = np.meshgrid(xx,yy)
        mask = np.zeros(self.data.shape)
        mask[np.where(np.sqrt(x**2 + y**2) <= dist)] = 1.0
        
        box = [star_loc[0] - dist, star_loc[0] + dist, star_loc[1] - dist, star_loc[1] + dist]
        #simple square aperture
        #ct = np.sum(frame[box[0][0]:box[1][0],box[2][0]:box[3][0]])
        self.star_cts = np.sum(self.data*mask)
        
        if plot:
            #plot a cut across the star to make sure box encapsulates all flux
            cut_x = self.data[:,star_loc[1]]
            cut_y = self.data[star_loc[0],:]
            fig, (ax0, ax1) = plt.subplots(2,1, figsize=(8,8))
            ax0.plot(cut_x.flatten())
            ax0.axvline(box[0], linestyle = ':', color = 'k')
            ax0.axvline(box[1], linestyle = ':', color = 'k')
            ax0.axhline(linestyle = ':', color = 'r')
            ax1.plot(cut_y.flatten())
            ax1.axvline(box[2], linestyle = ':', color = 'k', label = 'Bounding Box')
            ax1.axvline(box[3], linestyle = ':', color = 'k')
            ax1.axhline(linestyle = ':', color = 'r')
            ax0.set_xlabel('X position')
            ax0.set_ylabel('Counts')
            ax1.set_xlabel('Y position')
            ax1.set_ylabel('Counts')
            ax1.legend()
            plt.show()
              
    def find_flux_conversion(self, mag):
        '''Given star flux in cts/s and star mag (Vega) in a given filter, 
        calculate conversion between image counts and flux density units'''
        self.star_mag = mag
        flux_vega, w, dw = get_filt_info(self.filt)
        flux_star = flux_vega * 10**(-(1/2.5)*self.star_mag)
        self.flux_per = flux_star/self.star_cts #units erg s-1 cm-2 um-1 / cts s-1


