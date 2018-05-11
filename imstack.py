#!/usr/bin/env python
from nirc2_reduce import image
import matplotlib.pyplot as plt
import numpy as np

'''Functions to compute filter ratios, extract features at given locations,
generate points to compare against RT codes.
Takes a list of images with I/F already calculated'''

def get_lambda_eff(filt):
    tablef = '/Users/emolter/Python/nirc2_reduce/filter_passbands/vega_fluxes.txt'
    with open(tablef, 'r') as f:
        f.readline() #header line
        for line in f:
            l = line.split(',')
            ft = l[0].strip(', \n')
            wl = float(l[1].strip(', \n'))
            if ft == filt:
                return wl

class Stack:
    
    def __init__(self, fnames):
        '''fnames should be a list of filt_centered.fits files'''
        stack = {} #dict in format filter:CoordGrid object
        for fn in fnames:
            key = '_'.join(fn.split('_')[:-1])
            stack[key] = image.Image(fn)
        self.stack = stack
        
        '''Make images of filter ratios. For every possible combination of
        filters, compute filt1/filt2 over the image'''
        ratios = {}
        self.filts = list(self.stack.keys())
        self.wls_eff = [get_lambda_eff(filt) for filt in self.filts]
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
    
    def plot_one(self, filt):
        im = self.stack[filt]
        plt.imshow(im.data, origin = 'lower left')
        plt.show()
        
    def plot_ratio(self, filt1, filt2):
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
        fluxes = [self.stack[filt].data[x,y] for filt in self.filts]
        sort = np.argsort(self.wls_eff)
        print('filter    wl_eff (um)    I/F')
        for j in sort:
            print(self.filts[j] + '    ' + str(self.wls_eff[j]) + '    ' + str(fluxes[j]))
        
    def extract_feature(self, frac, which = 'all'):
        '''Find mean flux of a feature by outlining that feature with a contour
        frac: fraction of max flux in initial box to contour as the region
        filt: only used because need to show an image for extraction region
        230, 20; 266, 75'''

        if which == 'all':
            which = list(self.stack.keys())
        
        plt.imshow(self.stack[which[0]].data, origin = 'lower left')
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
        for i in range(len(sorted(which))):
            rstr = sorted(which)[i]
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
            print(rstr + '    ' + str(flux))
        
        
        
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
    
    
    
    
    