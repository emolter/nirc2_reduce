#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from .bxy3 import Bxy3
from . import image
from .nod import Nod

'''This group of routines finds the counts from photometric standard star,
computes the conversion factor between counts and specific intensity
to be applied to final science images.'''

nearest_stand_filt = {'j':'j', 'he1a':'j', 'he1_a':'j', 'pagamma':'j', 'jcont':'j', 'pabeta':'j',
                      'h':'h', 'hcont':'h', 'ch4s':'h', 'ch4_short':'h', 'feii':'h', 'ch4l':'h', 'ch4_long':'h',
                      'k':'k', 'kp':'k', 'kcont':'k', 'ks':'k', 'he1b':'k', 'brgamma':'k', 'br_gamma':'k', 'h210':'k', 'h221':'k', 'co':'k', 'h2o':'k', 'h2o_ice':'k',
                      'lp':'l', 'lw':'l', 'bracont':'l', 'bra':'l', 'br_alpha':'l', 'bra_cont':'l', 'pah':'l',
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


def airmass_correction(air_t, air_c, filt):
    '''Helper to I/F. Computes correction factor to photometry based on airmass.
       Multiply 
       air_t is airmass of science target.
       air_c is airmass of calibrator.
       filt options are j, h, k, l, m. Use nearest one for narrowband filters'''
    cdict = {'j': 0.102,
             'h': 0.059,
             'k': 0.088,
             'l': 0.093,
             'm': 0.220} #from https://www2.keck.hawaii.edu/realpublic/inst/nirc/exts.html
    if filt == None:
        return 1.0
    tau = cdict[filt]
    factor = np.exp(tau*air_t)/np.exp(tau*air_c)
    return factor
    
                
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
    elif cat_filt == 'l':
        catalog_filt = 'bessell_lp'
    elif cat_filt == 'm':
        catalog_filt = 'bessell_m'
    with open('/Users/emolter/Python/nirc2_reduce/filter_passbands/filter_conversions.txt','r') as f:
        f.readline() #header
        for line in f:
            l = line.split()
            fl1 = l[0].strip(', \n')
            fl2 = l[1].strip(', \n')
            factor = float(l[2].strip(', \n'))
            if fl1 == catalog_filt and fl2 == keck_filt:
                return factor   
        print('Failed! Check filter_passbands/filter_conversions.txt to ensure proper filter conversion exists')
        return


class bxy3Phot(Bxy3):
    
    def __init__(self,fnames):
        #fnames can be EITHER non-reduced OR reduced images: if loading reduced images just skip bxy3_reduce
        Bxy3.__init__(self,fnames)
    
    def reduce(self,skyfname,flatfname,badpxfname,outfnames):
        
        self.make_sky(skyfname)
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
        if plot:
            fig, axes = plt.subplots(3,2, figsize = (13,9))
        
        star_cts = []
        for i in range(len(self.frames)):
            frame = self.frames[i]
            star_loc = np.where(frame == np.max(frame))
            if len(star_loc[0]) > 1:
                print('WARN: bxy3Phot.find_cts: Found more than one maximum pixel in frame %d...be wary of unexpected results'%i)
                star_loc = (star_loc[0][0], star_loc[0][1])
            
            
            # make circular aperture
            xx = np.arange(frame.shape[0]) - star_loc[0]
            yy = np.arange(frame.shape[1]) - star_loc[1]
            x, y = np.meshgrid(xx,yy)
            mask = np.zeros(frame.shape)
            mask[np.where(np.sqrt(x**2 + y**2) <= dist)] = 1.0
            
            # simple square aperture
            box = [star_loc[0] - dist, star_loc[0] + dist, star_loc[1] - dist, star_loc[1] + dist]
            #ct = np.sum(frame[box[0][0]:box[1][0],box[2][0]:box[3][0]])
            
            #handle nonzero background, e.g. if sky brightens frame-to-frame
            bkgd = np.copy(frame)
            bkgd[mask != 0] = np.nan
            zeropt = np.nanmean(bkgd)
            corr_frame = frame-zeropt
            
            #put together count
            ct = np.sum(corr_frame*mask)
            star_cts.append(ct)
            
            if plot:
                #plot a cut across the star to make sure box encapsulates all flux
                cut_x = corr_frame[:,star_loc[1]]
                cut_y = corr_frame[star_loc[0],:]
                ax0 = axes[i][0]
                ax1 = axes[i][1]
                ax0.plot(cut_x.flatten())
                ax0.axvline(box[0], linestyle = ':', color = 'k')
                ax0.axvline(box[1], linestyle = ':', color = 'k')
                ax0.axhline(linestyle = ':', color = 'r')
                ax1.plot(cut_y.flatten())
                ax1.axvline(box[2], linestyle = ':', color = 'k')
                ax1.axvline(box[3], linestyle = ':', color = 'k')
                ax1.axhline(linestyle = ':', color = 'r')
                ax0.set_ylabel('Counts')
                if i == 2:
                    ax0.set_xlabel('X position')
                    ax1.set_xlabel('Y position')
                ax0.set_title('Frame %d'%i, loc = 'left', y = 0.83)
                
        self.star_cts = star_cts
        if plot:
            plt.show()
        
    def find_flux_conversion(self, star_mag, keck_filt, verbose = False):
        '''Given star flux in cts/s and star mag in a given filter, 
        calculate conversion between image counts and flux density units'''
        w, flux_vega = get_filt_info(keck_filt)
        #factor = get_factor(keck_filt) #vega flux in keck filter / vega flux in catalog filter #this is not needed
        flux_star = flux_vega * 10**(-(1/2.5)*star_mag)
        flux_per =  [flux_star/ct for ct in self.star_cts] #units erg s-1 cm-2 um-1 / cts s-1
        flux_per = np.asarray(flux_per)
        self.all_flux_per = flux_per
        self.flux_per = np.median(flux_per)
        if verbose:
            print('all values ', flux_per)
            print('Median flux per = ',self.flux_per)
            
    def help(self):
        
        helpstr = '''
        Contains tasks to reduce standard star observations taken with the
            keck bxy3 dither, and extract photometric information
        
        Functions (see DOCUMENTATION.py for use):
            reduce(self,skyfname,flatfname,badpxfname,outfnames)
            find_cts(self, dist = 50, plot = True)
            find_flux_conversion(self, star_mag, keck_filt, verbose = False)
             
        Attributes:
            *In addition to below list, contains all attributes of Bxy3 object
            star_cts: number of total cts s-1 from the star for each frame
            all_flux_per: flux per second for each frame, units erg s-1 cm-2 um-1 / cts s-1
            flux_per: flux per second averaged over the three frames, units erg s-1 cm-2 um-1 / cts s-1
        '''
        print(helpstr)
        
        
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
        
    def help(self):
        
        helpstr = '''
        Contains tasks to reduce standard star observations taken with the
            a target/sky nod, and extract photometric information
        
        Functions (see DOCUMENTATION.py for use):
            reduce(self, flatfname,badpxfname,outfname)
            find_cts(self, dist = 50, plot = True)
            find_flux_conversion(self, star_mag, keck_filt)
             
        Attributes:
            *In addition to below list, contains all attributes of Nod object
            star_cts: number of total cts s-1 from the star
            flux_per: flux per second, units erg s-1 cm-2 um-1 / cts s-1
        '''
        print(helpstr)
           