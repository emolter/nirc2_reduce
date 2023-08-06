import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import fsps
import pysynphot


nearest_stand_filt = {'j':'j', 'he1a':'j', 'he1_a':'j', 'pagamma':'j', 'jcont':'j', 'pabeta':'j',
                      'h':'h', 'hcont':'h', 'ch4s':'h', 'ch4_short':'h', 'feii':'h', 'ch4l':'h', 'ch4_long':'h',
                      'k':'k', 'kp':'k', 'kcont':'k', 'ks':'k', 'he1b':'k', 'brgamma':'k', 'br_gamma':'k', 'h210':'k', 'h221':'k', 'co':'k', 'h2o':'k', 'h2o_ice':'k',
                      'lp':'l', 'lw':'l', 'bracont':'l', 'bra':'l', 'br_alpha':'l', 'bra_cont':'l', 'pah':'l',
                      'ms':'m'}


class Filt:
    
    def __init__(self, filt, isfsps = False):
        '''wl will be output in microns'''
        self.filtname = filt.lower()
        #try:
        if not isfsps:
            try:
                # Keck
                self.load_data()
            except:
                # pysynphot
                bp = pysynphot.ObsBandpass(self.filtname)
                self.wl = bp.binset * 1.0e-4 #convert to micron from angstrom
                trans = bp(bp.binset)
                self.trans = trans/np.max(trans)
                self.wl_eff = np.average(self.wl, weights = self.trans)
                self.interp_f = interp1d(self.wl,self.trans, bounds_error = False, fill_value = 0.0)
                
        else:
            # pysynphot replaces this - really only needed so legacy code keeps working
            ftr = fsps.get_filter(self.filtname)
            wl, trans = ftr.transmission
            self.wl = wl *1.0e-4 #convert to micron from angstrom
            self.trans = trans/np.max(trans)
            self.wl_eff = np.average(self.wl, weights = self.trans)
            self.interp_f = interp1d(self.wl,self.trans, bounds_error = False, fill_value = 0.0)
        #except:
        #    print('Could not find filter name in Keck NIRC2 filter table nor in pysynphot')
        #    print('Check spelling, else check that data for that filter exist')
    
    def load_data(self):
        wl = []
        trans = []
        infile = '/Users/emolter/Python/nirc2_reduce/filter_passbands/'+self.filtname+'.csv'
        with open(infile,'r') as f:
            f.readline()
            for line in f:
                l = line.split(',')
                try:
                    wl.append(float(l[0].strip(', \n')))
                    trans.append(float(l[1].strip(', \n')))
                except:
                    pass
        self.wl = np.asarray(wl)
        self.trans = np.asarray(trans)/100
        self.wl_eff = np.average(self.wl, weights = self.trans)
        self.interp_f = interp1d(self.wl,self.trans, bounds_error = False, fill_value = 0.0)
    
    def plot(self):
        plt.plot(self.wl,self.trans, color = 'r')
        plt.xlabel(r'Wavelength ($\mu$m)')
        plt.ylabel('Transmission')
        plt.show()
    
    def find_vega_flux(self):
        '''Returns the flux of Vega in that filter in units erg s-1 cm-2 um-1'''
        vega = Spectrum('filter_passbands/vega_spectrum.txt', 'vega')
        t = self.eval(vega.wl) #transmission at vega spectrum resolution
        f = np.sum(vega.flux * t)/np.sum(t)
        self.vega_flux = f
        return f
        
    def find_sun_flux(self):
        '''Returns the flux of the Sun in that filter in units erg s-1 cm-2 um-1'''
        sun = Spectrum('filter_passbands/sun_spectrum.txt', 'sun')
        t = self.eval(sun.wl) #transmission at vega spectrum resolution
        f = np.sum(sun.flux * t)/np.sum(t)
        self.sun_flux = f
        return f        
        
    def eval(self,w):
        out = []
        for val in w:
            if ( self.wl[0] < val ) & ( val < self.wl[-1] ):
                interpval = self.interp_f(val)
                if interpval > 0:
                    out.append(self.interp_f(val))
                else:
                    out.append(0)
            else:
                out.append(0)
        return np.asarray(out)
        
    def help(self):
        
        helpstr = '''
        Contains nirc2 filter transmission data and simple tasks for querying it
        
        Functions:
            load_data(self)
            plot(self)
            find_vega_flux(self)
            eval(self,w)
             
        Attributes:
            filtname: name of filter
            wl: list of wls in microns
            trans: normalized transmission in filter at each wl
            wl_eff: effective wl of filt
            interp_f: interpolation function takes in given wl and outputs trans at that wl
            vega_flux: flux of vega in the filter as calculated by find_vega_flux
            sun_flux: flux of sun in the filter as calculated by find_sun_flux
        '''
        print(helpstr)
            
class Spectrum:
    
    def __init__(self, infile, units):
        '''two column spectrum 
        units = vega: Angstrom, erg s-1 cm-2 A-1
        units = sun: nm, W m-2 nm-1'''
        spec = np.loadtxt(infile, comments = '#')
        if units == 'vega':
            self.wl = spec.T[0] * 1.0e-4 #from angstroms to microns
            self.flux = spec.T[1] * 10000.0 #from angstroms to microns
        elif units == 'sun':
            self.wl = spec.T[0] * 1.0e-3 #from nm to microns
            fx = spec.T[1] * 1000. #from nm to microns
            fx *= 1.0e7 #from W to erg s-1
            self.flux = fx*1.0e-4 #from m-2 to cm-2
        
    def plot(self):
        plt.plot(self.wl, self.flux)
        plt.xlim([0,5])
        plt.xlabel('Wavelength ($\mu$m)')
        plt.ylabel('Flux density (erg s$^{-1}$ cm$^{-2}$ $\mu$m$^{-1}$)')
        plt.show()
        
    def flux_in_filter(self,filt):
        filter_object = Filt(filt)
        trans = filter_object.eval(self.wl)
        flux = np.sum(trans*self.flux)/(np.sum(trans))
        return flux
        
    def help(self):
        
        helpstr = '''
        Simple wrapper for 2-column spectra, helps Filt object
        
        Functions:
            plot(self)
            flux_in_filter(self,filt)
             
        Attributes:
            wl: wavelengths
            flux: fluxes
        '''
        print(helpstr)


        
def make_vega_table():
    '''This function to make table of effective wavelengths and fluxes
    for each Keck filter'''
    filtlist = ['brgamma','feii','h','hcont','j','k','kcont','kp','ks','lp','ms','pah','br_alpha','bra_cont','h2o_ice','ch4_short','pabeta','he1_a']
    with open('filter_passbands/vega_fluxes.txt', 'w') as f:
        f.write('filter, lambda_eff (um), Vega flux density (erg s-1 cm-2 um-1)'+'\n')
        for filtname in filtlist:
            filt = Filt(filtname)
            vega_mag = str(filt.find_vega_flux())
            wl_eff = str(filt.wl_eff)
            fnstr = filtname + '          '
            f.write(fnstr[:10] + ', ' + wl_eff + ', ' + vega_mag + '\n')

def make_sun_table():
    filtlist = ['brgamma','feii','h','hcont','j','k','kcont','kp','ks','lp','ms','pah','br_alpha','bra_cont','h2o_ice','ch4_short','pabeta','he1_a']
    with open('filter_passbands/sun_fluxes.txt', 'w') as f:
        f.write('filter, lambda_eff (um), Solar flux density (erg s-1 cm-2 um-1)'+'\n')
        for filtname in filtlist:
            filt = Filt(filtname)
            sun_mag = str(filt.find_sun_flux())
            wl_eff = str(filt.wl_eff)
            fnstr = filtname + '          '
            f.write(fnstr[:10] + ', ' + wl_eff + ', ' + sun_mag + '\n')
            
            
def find_vega_narrow(wl_in):
    vega = Spectrum('filter_passbands/vega_spectrum.txt', 'vega')
    interp = interp1d(vega.wl, vega.flux)
    return interp(wl_in)
    
def find_sun_narrow(wl_in):
    sun = Spectrum('filter_passbands/sun_spectrum.txt', 'sun')
    interp = interp1d(sun.wl, sun.flux)
    return interp(wl_in)
        
def convert_filters(f1, f2):
    '''f1 is the input filter, usually 2mass or some large survey.
    will be read by FSPS. f2 is the output Keck filter.'''
    
    filt1 = Filt(f1, isfsps = True)
    v1 = filt1.find_vega_flux()
    print(v1)
    filt2 = Filt(f2)
    v2 = filt2.find_vega_flux()
    print(v2)
    return v2/v1
    
'''
We need a way to convert between filters.
Say we have a standard star flux in 2MASS k band, but we want its flux
in NIRC2 kp band. How do we go between those? Can we just assume our
standard star has the same relative flux in those filters as Vega?
if so, can use 2mass filter definitions in Python FSPS to get Vega mags
for those, then compare to the NIRC2 filter Vega mags I compute to get a
conversion factor.
'''

if __name__ == "__main__":
    sun = Spectrum('filter_passbands/sun_spectrum.txt', 'sun')
    sun.plot()
    make_sun_table()        
    make_vega_table()   
    print(find_sun_narrow(1.2903))
    print(convert_filters('2mass_ks', 'h2o_ice'))
    print(convert_filters('bessell_m', 'ms'))