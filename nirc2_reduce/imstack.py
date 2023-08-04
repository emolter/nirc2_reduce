#!/usr/bin/env python
from nirc2_reduce import image, filt
import matplotlib.pyplot as plt
import numpy as np

'''Functions to compute filter ratios, extract features at given locations,
generate points to compare against RT codes.
Takes a list of images with I/F already calculated'''

#def get_lambda_eff(filt):
#    tablef = '/Users/emolter/Python/nirc2_reduce/filter_passbands/vega_fluxes.txt'
#    with open(tablef, 'r') as f:
#        f.readline() #header line
#        for line in f:
#            l = line.split(',')
#            ft = l[0].strip(', \n')
#            wl = float(l[1].strip(', \n'))
#            if ft == filt:
#                return wl

class Stack:
    
    def __init__(self, fnames, filts = None):
        '''fnames should be a list of filt_centered.fits files'''
        if not filts:
            filts = []
            for fn in fnames:
                key = '_'.join(fn.split('_')[:-1])
                filts.append(key)
            self.filts = filts
        else:
            self.filts = filts
        
        stack = {} #dict in format filter:CoordGrid object
        for i in range(len(fnames)):
            key = self.filts[i]
            stack[key] = image.Image(fnames[i])
        self.stack = stack
        
        '''Make images of filter ratios. For every possible combination of
        filters, compute filt1/filt2 over the image'''
        ratios = {}
        filter_objects = [filt.Filt(f) for f in self.filts]
        self.wls_eff = [f.wl_eff for f in filter_objects]
        for i in range(len(self.filts)):
            for j in range(len(self.filts)):
                if i != j:
                    im1 = self.stack[self.filts[i]].data
                    im2 = self.stack[self.filts[j]].data
                    ratio_key = self.filts[i]+'/'+self.filts[j]
                    ratioim = im1/im2
                    ratios[ratio_key] = ratioim         
        self.ratios = ratios
        print('Image stack object created with filters {}'.format(self.filts))
    
    def plot_one(self, f):
        print('Note if this is a stack of projected images, axes may be flipped!')
        im = self.stack[f]
        plt.imshow(im.data, origin = 'lower left')
        plt.show()
        
    def plot_ratio(self, filt1, filt2):
        print('Note if this is a stack of projected images, axes may be flipped!')
        ratioim = self.ratios[filt1+'/'+filt2]
        fig, ax = plt.subplots(1,1, figsize = (9,9))
        ax.imshow(ratioim, origin = 'lower left')
        plt.show()
        
    def write(self, fname, filt1, filt2):
        ratioim = self.ratios[filt1+'/'+filt2]
        hdulist_out = self.stack[filt1].hdulist
        hdulist_out[0].header['OBJECT'] = 'ratio_'+filt1+'/'+filt2
        hdulist_out[0].data = ratioim
        hdulist_out[0].writeto(fname, overwrite=True)
        
    def extract_point(self, x, y):
        '''Print the I/F value and central wavelength of each filter
        for a given x,y coordinate in the image'''
        fluxes = [self.stack[f].data[x,y] for f in self.filts]
        sort = np.argsort(self.wls_eff)
        print('filter    wl_eff (um)    I/F')
        for j in sort:
            print(self.filts[j] + '    ' + str(self.wls_eff[j]) + '    ' + str(fluxes[j]))
        
    def extract_feature(self, frac, which = 'all', showme = 0):
        '''Find mean flux of a feature by outlining that feature with a contour
        frac: fraction of max flux in initial box to contour as the region
        filt: only used because need to show an image for extraction region
        230, 20; 266, 75
        showme: which filter do you want displayed? int from 0 to len(self.filts)'''

        if which == 'all':
            which = list(self.filts)
            
        
        plt.imshow(self.stack[which[showme]].data, origin = 'lower left')
        plt.show()
        print('Define a box around the feature you want to track. Note x,y are reversed in image due to weird Python indexing!')
        pix_l = input('Enter lower left pixel x,y separated by a comma: ')
        pix_u = input('Enter upper right pixel x,y separated by a comma: ')
        
        p0x, p0y = int(pix_l.split(',')[0].strip(', \n')),int(pix_l.split(',')[1].strip(', \n'))
        p1x, p1y = int(pix_u.split(',')[0].strip(', \n')),int(pix_u.split(',')[1].strip(', \n'))      
        
        def make_contour(data, frac):
            rgn = np.copy(data)
            rgn[rgn < frac*np.max(rgn)] = 0.0
            rgn[rgn > 0.0] = 1.0
            rgn[rgn < 1.0] = np.nan
            return rgn
        
        print('Contour level: %f'%frac)
        print('filter    I/F')
        ifvals = []
        wlseff = []
        
        which.insert(0, which.pop(showme))
        for i in range(len(which)):
            #rstr = sorted(which(i))
            rstr = which[i]
            wleff = self.wls_eff[i]
            alldata = self.stack[rstr].data
            data = alldata[p0x:p1x,p0y:p1y]
            if i == 0:
                rgn = make_contour(data, frac)
                rgn_plot = np.copy(rgn)
                rgn_plot[np.isnan(rgn_plot)] = 0.0
                plt.imshow(data, origin = 'lower left')
                plt.contour(rgn_plot, [0.5], colors = 'white')
                plt.show()
            flux = np.nanmean(rgn * data)
            wlseff.append(wleff)
            ifvals.append(flux)
            print(rstr + '    ' + str(flux))
        return wlseff, ifvals
        
        
        
    def extract_ratios(self, x, y, which = 'all'):
        '''Print the filter ratios for a given x,y coordinate in the image
        optionally specify which ratios you want as list, 
        e.g. ['kp/h', 'kp/ch4s', 'kp/pabeta']'''
        if which == 'all':
            which = list(self.ratios.keys())
        print('ratio    I/F')
        for rstr in sorted(which):
            flux = self.ratios[rstr][x,y]
            print(rstr + '    ' + str(flux))
            
    def help(self):
        
        helpstr = '''
        Contains tasks for extracting information from a stack of images
             taken in different nirc2 filters. relies on naming conventions
             output by coordgrid.py to load things properly
        
        Functions (see DOCUMENTATION.py for use):
            plot_one(self, filt)
            plot_ratio(self, filt1, filt2)
            write(self, fname, filt1, filt2)
            extract_point(self, x, y)
            extract_feature(self, frac, which = 'all')
             
        Attributes:
            stack: list of images
            filts: list of filters for the images
            wls_eff: list of effective wavelengths for the filters
            ratios: all of the flux ratios as a dict, e.g. ratios['kp/h'] = float
        '''
        print(helpstr)
    
    
    
    
    