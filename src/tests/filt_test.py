#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from nirc2_reduce import filt


filter_curve_h = filt.Filt('H') # H band, taken from nirc2_reduce/filter_passbands
filter_curve_wfpc2 = filt.Filt('wfpc2_f555w', isfsps = True) # HST filter, taken from FSPS
curve_list = [filter_curve_h, filter_curve_wfpc2]
labels = ['H', 'F555W']
colors = ['red', 'blue']


spec_in = np.loadtxt('/Users/emolter/Python/nirc2_reduce/filter_passbands/sun_spectrum.txt', comments = '#') # wl in nm, irradiance in W m-2 nm-1
# Erandi - note that since convloving with a transmission curve is basically just a weighted average, the brightness units don't matter!  This works the same for I/F, or SI units, or cgs, or whatever
wl, flux = spec_in.T[0], spec_in.T[1]
# nevertheless for readability we will convert to um from nm
wl*=1e-3
flux*=1e3

fig, ax = plt.subplots(1,1, figsize = (7,5))
ax.semilogx(wl, flux, color = 'k', label = 'solar spectrum')
ax2 = ax.twinx()
for i in range(2):
    curve = curve_list[i]
    
    ### HERE IS THE EXCITING PART
    # evaluate the transmission curve at the spectrum's native wavelengths
    trans = curve.eval(wl)
    ax2.semilogx(wl, trans, label = labels[i], color = colors[i])
    
    # weighted average to get flux in filter
    flux_in_filter = np.sum(trans*flux)/(np.sum(trans))
    
    # the module also contains the effective wavelength of the filter, i.e., what you would get by convolving the transmission curve with a flat spectrum
    ax.semilogx(curve.wl_eff, flux_in_filter, label = 'Sun %s-band flux'%(labels[i]), color = colors[i], marker = 'o', markersize = 10, linestyle = '')
    
    
    

ax.set_xlim([0.3, 3])    
ax.set_xlabel(r'Wavelength ($\mu$m)') 
ax.set_ylabel(r'Flux (W m$^{-2}$ $\mu$m$^{-1}$)') 
ax2.set_ylabel('Transmission')
ax.legend()

plt.savefig('filt_test.png', bbox = None)
plt.show()