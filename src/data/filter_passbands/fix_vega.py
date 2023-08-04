#!/usr/bin/env python
import numpy as np
from nirc2_reduce.filt import Spectrum
import matplotlib.pyplot as plt

spec = Spectrum('vega_spectrum.txt', units = 'vega')
print(spec.flux[spec.wl > 2.5])

def jeans(w,T):
    '''wl in microns, T in kelvin, return erg s-1 cm-2 um-1'''
    c = 2.997e8 #m s-1
    kb = 1.38e-23 #J K-1
    factor = 1e5
    return factor * 2*c*kb*T / w**4
    
tvega = 11400
wls = np.arange(2.61,5.0,0.01)
fluxes = jeans(wls, tvega)

plt.plot(wls, fluxes)
plt.plot(spec.wl, spec.flux)
plt.xlim([2.0,3.5])
plt.ylim([0, 2e-6])
plt.show()

#so it works. let's add those values into the spectrum
wls *= 1e4 #back to angstroms
fluxes *= 1e-4 #back to angstroms

wl_out = np.concatenate((spec.wl*1e4, wls))
flux_out = np.concatenate((spec.flux*1e-4, fluxes))

newarr = np.asarray([wl_out,flux_out]).T
np.savetxt('vega_new.txt', newarr, delimiter = ' ')