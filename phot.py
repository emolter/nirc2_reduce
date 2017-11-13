#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from .bxy3 import Bxy3
from . import image
from .nod import Nod

'''This group of routines finds the counts from photometric standard star,
computes the conversion factor between counts and specific intensity
to be applied to final science images.'''

nearest_stand_filt = {'j':'j', 'he1a':'j', 'pagamma':'j', 'jcont':'j', 'pabeta':'j',
                      'h':'h', 'hcont':'h', 'ch4s':'h', 'feii':'h', 'ch4l':'h',
                      'k':'k', 'kp':'k', 'kcont':'k', 'ks':'k', 'he1b':'k', 'brgamma':'k', 'h210':'k', 'h221':'k', 'co':'k',
                      'lp':'l', 'lw':'l', 'bracont':'l', 'bra':'l',
                      'ms':'m'}

def get_filt_info(filt):
    with open('/Users/emolter/Python/nirc2_reduce/filter_passbands/vega_fluxes.txt','r') as f:
        f.readline() #header
        for line in f:
            l = line.split(',')
            fl = l[0].strip(', \n')
            if fl == filt:
                wl = float(l[1].strip(', \n'))
                vega_mag = float(l[2].strip(', \n'))
                return wl, vega_mag
                
def get_factor(keck_filt):
    '''Queries table of conversions between flux in catalog filters and 
    flux in Keck filters. table was produced to work for A star standards only,
    otherwise color correction may become important'''
    cat_filt = nearest_stand_filt[keck_filt]
    if cat_filt == 'h':
        catalog_filt = '2mass_h'
    elif cat_filt == 'j':
        catalog_filt = '2mass_j'
    elif cat_filt == 'k':
        catalog_filt = '2mass_ks'
    with open('/Users/emolter/Python/nirc2_reduce/filter_passbands/filter_conversions.txt','r') as f:
        f.readline() #header
        for line in f:
            l = line.split()
            fl1 = l[0].strip(', \n')
            fl2 = l[1].strip(', \n')
            factor = float(l[2].strip(', \n'))
            if fl1 == catalog_filt and fl2 == keck_filt:
                return factor   
        return 'Failed! Check filter_passbands/filter_conversions.txt to ensure proper filter conversion exists'     

class bxy3Phot(Bxy3):
    
    def __init__(self,fnames):
        #fnames can be EITHER non-reduced OR reduced images: if loading reduced images just skip bxy3_reduce
        Bxy3.__init__(self,fnames)
    
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
        
    def find_flux_conversion(self, star_mag, keck_filt):
        '''Given star flux in cts/s and star mag in a given filter, 
        calculate conversion between image counts and flux density units'''
        w, flux_vega = get_filt_info(keck_filt)
        factor = get_factor(keck_filt) #vega flux in keck filter / vega flux in catalog filter
        flux_star = flux_vega * (1/factor) * 10**(-(1/2.5)*star_mag)
        flux_per =  [flux_star/ct for ct in self.star_cts] #units erg s-1 cm-2 um-1 / cts s-1
        flux_per = np.asarray(flux_per)
        print('all values ', flux_per)
        self.flux_per = np.median(flux_per)
        print('Median flux per = ',self.flux_per)
        
        
class nodPhot(Nod):
    
    def __init__(self,skyf,imagef):
        
        Nod.__init__(self,skyf,imagef)
    
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
        
    def find_cts(self, dist = 50, plot = True):
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
            
    def find_flux_conversion(self, star_mag, keck_filt):
        '''Given star flux in cts/s and star mag in a given filter, 
        calculate conversion between image counts and flux density units'''
        w, flux_vega = get_filt_info(keck_filt)
        factor = get_factor(keck_filt) #vega flux in keck filter / vega flux in catalog filter
        flux_star = flux_vega * (1/factor) * 10**(-(1/2.5)*star_mag)
        self.flux_per =  flux_star/self.star_cts #units erg s-1 cm-2 um-1 / cts s-1
           