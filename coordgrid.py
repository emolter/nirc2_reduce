#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from .get_ephem import get_ephemerides, naif_lookup
from nirc2_reduce.image import Image
from nirc2_reduce.phot import nearest_stand_filt
from datetime import datetime, timedelta
import warnings
from skimage import feature
from image_registration.chi2_shifts import chi2_shift
from image_registration.fft_tools.shift import shiftnd, shift2d
from scipy.interpolate import interp2d, RectBivariateSpline, NearestNDInterpolator, griddata
#from .fit_gaussian import fitgaussian
from astropy.modeling import models, fitting
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage.interpolation import zoom

#from mpl_toolkits import basemap
import pyproj


def lat_lon(x,y,ob_lon,ob_lat,pixscale_km,np_ang,req,rpol):
    '''Find latitude and longitude on planet given x,y pixel locations and
    planet equatorial and polar radius'''
    np_ang = -np_ang
    x1 = pixscale_km*(np.cos(np.radians(np_ang))*x - np.sin(np.radians(np_ang))*y)
    y1 = pixscale_km*(np.sin(np.radians(np_ang))*x + np.cos(np.radians(np_ang))*y)
    olrad = np.radians(ob_lat)
    
    #set up quadratic equation for ellipsoid
    r2 = (req/rpol)**2
    a = 1 + r2*(np.tan(olrad))**2 #second order
    b = 2*y1*r2*np.sin(olrad) / (np.cos(olrad)**2) #first order
    c = x1**2 + r2*y1**2 / (np.cos(olrad))**2 - req**2 #constant

    radical = b**2 - 4*a*c
    #will equal nan outside planet since radical < 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") #suppresses error for taking sqrt nan
        x3s1=(-b+np.sqrt(radical))/(2*a)
        x3s2=(-b-np.sqrt(radical))/(2*a)
    z3s1=(y1+x3s1*np.sin(olrad))/np.cos(olrad)
    z3s2=(y1+x3s2*np.sin(olrad))/np.cos(olrad)
    odotr1=x3s1*np.cos(olrad)+z3s1*np.sin(olrad)
    odotr2=x3s2*np.cos(olrad)+z3s2*np.sin(olrad)
    #the two solutions are front and rear intersections with planet
    #only want front intersection
    
    #tricky way of putting all the positive solutions into one array
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") #suppresses error for taking < nan
        odotr2[odotr2 < 0] = np.nan
        x3s2[odotr2 < 0] = np.nan
        z3s2[odotr2 < 0] = np.nan
        odotr1[odotr1 < 0] = odotr2[odotr1 < 0]
        x3s1[odotr1 < 0] = x3s2[odotr1 < 0]
        z3s1[odotr1 < 0] = z3s2[odotr1 < 0]
    
    odotr,x3,z3 = odotr1,x3s1,z3s1
    y3 = x1
    r = np.sqrt(x3**2 + y3**2 + z3**2)
    
    #lon_w = np.degrees(np.arctan(y3/x3)) + ob_lon
    lon_w = np.degrees(np.arctan2(x3,y3)-np.pi/2) + ob_lon 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") #suppresses error for taking < nan
        lon_w[lon_w < 0] += 360
        lon_w = lon_w%360
    lat_c = np.degrees(np.arcsin(z3/r))
    lat_g = np.degrees(np.arctan(r2*np.tan(np.radians(lat_c))))
    #plt.imshow(lon_w, origin = 'lower left')
    #plt.show()
    return lat_g, lat_c, lon_w

def surface_normal(lat_g, lon_w, ob_lon):
    '''Returns the normal vector to the surface of the planet.
    Take dot product with sub-obs or sub-sun vector to find cosine of emission angle'''
    nx = np.cos(np.radians(lat_g))*np.cos(np.radians(lon_w-ob_lon))
    ny = np.cos(np.radians(lat_g))*np.sin(np.radians(lon_w-ob_lon))
    nz = np.sin(np.radians(lat_g))
    return np.asarray([nx,ny,nz])

def emission_angle(ob_lat, surf_n):
    '''Return the cosine of the emission angle of surface wrt observer'''
    ob = np.asarray([np.cos(np.radians(ob_lat)),0,np.sin(np.radians(ob_lat))])
    return np.dot(surf_n.T, ob).T
    
#def sun_angle(ob_lon, ob_lat, sun_lon, sun_lat):
#    return
    
def get_filt_info(filt):
    '''Helper to I/F. Will find flux of sun in given filter'''
    with open('/Users/emolter/Python/nirc2_reduce/filter_passbands/sun_fluxes.txt','r') as f:
        f.readline() #header
        for line in f:
            l = line.split(',')
            fl = l[0].strip(', \n')
            if fl == filt:
                wl = float(l[1].strip(', \n'))
                sun_mag = float(l[2].strip(', \n'))
                return wl, sun_mag
                
def find_airmass(observatory, time, object):
    '''Use the power of astropy to retrieve airmass of standard
    star automatically... obvs not done yet'''
    return
                
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

class CoordGrid:
    
    def __init__(self, infile, lead_string = None, req = 24764, rpol = 24341, scope = 'keck'):
        '''Pull ephemeris data, calculate lat and lon. Pixscale in arcsec, req and rpol in km'''
        
        self.im = Image(infile)
        self.req = req
        self.rpol = rpol
        scope = str.lower(scope)

        if scope == 'keck':
            pixscale = 0.009942
        elif scope == 'lick':
            pixscale = 0.033
        elif scope == 'alma':
            pixscale = np.abs(self.im.header['CDELT1']) * 3600 #deg to arcsec
        elif scope == 'hst':
            pixscale = np.abs(self.im.header['PIXSCAL'])
        else:
            pixscale = 0.0
        self.pixscale_arcsec = pixscale
        if scope == 'vla' or scope == 'alma':
            self.data = self.im.data[0,0,:,:]
        else:
            self.data = self.im.data
        
        #pull and reformat header info
        if not scope == 'hst_wfc3':
            targ = self.im.header['OBJECT'].split('_')[0]
            targ = targ.split(' ')[0]
            self.target = targ
        else:
            targ = 'Neptune'
            self.target = 'Neptune'
        if scope == 'vla' or scope == 'alma':
            date = self.im.header['DATE-OBS']
            expstart = date.split('T')[1]
            date = date.split('T')[0]
        elif scope == 'hst':
            date = self.im.header['DATE-OBS']
            expstart = self.im.header['TIME-OBS']
        elif scope == 'keck':
            expstart = self.im.header['EXPSTART']
            date = self.im.header['DATE-OBS']
        imsize_x = self.data.shape[0]
        imsize_y = self.data.shape[1]
        if scope == 'lick':
            tstart = datetime.strptime(self.im.header['DATE-BEG'][:-7],'%Y-%m-%dT%H:%M')
        else:
            tstart = datetime.strptime(date+' '+expstart[:5],'%Y-%m-%d %H:%M')
        tend = tstart + timedelta(minutes=1)
        tstart = datetime.strftime(tstart, '%Y-%m-%d %H:%M')
        tend = datetime.strftime(tend, '%Y-%m-%d %H:%M')
        self.date_time = tstart
        
        #pull ephemeris data
        naif = naif_lookup(targ)
        if scope == 'keck':
            obscode = 568
        elif scope == 'vla':
            obscode = -5
        elif scope == 'lick':
            obscode = 662
        elif scope == 'alma':
            obscode = -7
        elif scope == 'hst_wfc3' or scope == 'hst_opal' or scope == 'hst_wfc2':
            obscode = '500@-48'
        else:
            obscode = input('Enter Horizons observatory code: ')
        
        ## check ephem inputs
        #print(naif)    
        #print(obscode)
        #print(tstart)
        #print(tend)
        
        ephem = get_ephemerides(naif, obscode, tstart, tend, '1 minutes')[0][0] #just the row for start time
        ephem = [val.strip(' ') for val in ephem]
        time = ephem[0]
        ra, dec = ephem[3], ephem[4]
        dra, ddec = float(ephem[5]), float(ephem[6])
        if not scope == 'hst_wfc3':
            az, el = float(ephem[7]), float(ephem[8])
            self.airmass, extinction = float(ephem[9]), float(ephem[10])
        apmag, sbrt = float(ephem[11]), float(ephem[12])
        self.ang_diam = float(ephem[15])
        self.ob_lon, self.ob_lat = float(ephem[16]), float(ephem[17])
        self.sun_lon, self.sun_lat = float(ephem[18]), float(ephem[19])
        #self.sun_ang = sun_angle(ob_lon, ob_lat, sun_lon, sun_lat)
        self.np_ang, self.np_dist = float(ephem[20]), float(ephem[21])
        self.sun_dist = float(ephem[22])*1.496e8 #from AU to km
        self.dist = float(ephem[24])*1.496e8 #from AU to km        
        
        self.pixscale_km = self.dist*np.radians(pixscale/3600)
        avg_circumference = 2*np.pi*((req + rpol)/2.0)
        self.deg_per_px = self.pixscale_km * (1/avg_circumference) * 360 #approximate conversion between degrees and pixels at sub-observer point
        
        if lead_string != None:
            #if you already did the edge detection and centering and are loading centered image
            self.centered = self.im.data
            self.lat_g = Image(lead_string+'_latg.fits').data
            self.lon_w = Image(lead_string+'_lone.fits').data
            self.model_planet = np.nan_to_num(self.lat_g * 0.0 + 1.0)
            try:
                self.err_x = Image(lead_string+'_errx.fits').data
                self.err_y = Image(lead_string+'_erry.fits').data
            except:
                pass
            try:
                self.projected = Image(lead_string+'_proj.fits').data
            except:
                pass
            try:
                self.mu = Image(lead_string+'_mu.fits').data
            except:
                self.surf_n = surface_normal(self.lat_g, self.lon_w, self.ob_lon)
                self.mu = emission_angle(self.ob_lat, self.surf_n)
            try:
                self.mu_projected = Image(lead_string+'_mu_proj.fits').data
            except:
                pass
        
        else:
            xcen, ycen = int(imsize_x/2), int(imsize_y/2) #pixels at center of planet
            xx = np.arange(imsize_x) - xcen
            yy = np.arange(imsize_y) - ycen
            x,y = np.meshgrid(xx,yy)
            self.lat_g, self.lat_c, self.lon_w = lat_lon(x,y,self.ob_lon,self.ob_lat,self.pixscale_km,self.np_ang,req,rpol)
            self.surf_n = surface_normal(self.lat_g, self.lon_w, self.ob_lon)
            self.mu = emission_angle(self.ob_lat, self.surf_n)

    def ioverf(self, filt, flux_per, stand_airmass):
        '''Compute I/F ratio given an image in cts s-1 and a conversion between
        cts s-1 and erg s-1 cm-2 um-1 sr-1'''
        wl, sun_flux_earth = get_filt_info(filt)
        sun_flux = sun_flux_earth * (1/np.pi)*(1.496e8/self.sun_dist)**2 #factor of pi because F = pi*B
        sr_per_px = np.radians(self.pixscale_arcsec/3600)**2
        sun_flux_density = sun_flux * sr_per_px # from erg s-1 cm-2 um-1 sr-1 to erg s-1 cm-2 um-1 px-1
        print('Sun flux density ', sun_flux_density)
        
        #photometry correction
        airmass_filt = nearest_stand_filt[filt]
        air_corr = airmass_correction(self.airmass, stand_airmass, airmass_filt)
        print('Airmass correction ',air_corr)
        self.data = self.data * flux_per * air_corr / sun_flux_density
        if hasattr(self, 'centered'):
            self.centered = self.centered * flux_per * air_corr / sun_flux_density
    
    #def minnaert(self, k):
    #    '''Applies Minnaert correction to remove effects of limb darkening
    #    k is an empirically determined constant between 0 and 1, 0.7 for Io'''
    #    self.data = self.data * self.sun_ang**k * self.mu**(k - 1)
    #    if hasattr(self, 'centered'):
    #        self.centered = self.centered * mu0**k * self.mu**(k - 1)
            
    def write_photonly(self, outstr):
        '''If you want to just run ioverf and then write'''
        hdulist_out = self.im.hdulist
        #centered data
        hdulist_out[0].header['OBJECT'] = self.target+'_calibrated'
        hdulist_out[0].data = self.data
        hdulist_out[0].writeto(outstr, overwrite=True)

    def edge_detect(self, low_thresh = 0.01, high_thresh = 0.05, sigma = 5, plot = True):
        '''Uses skimage canny algorithm to find edges of planet, correlates
        that with edges of model, '''
        
        self.model_planet = np.nan_to_num(self.lat_g * 0.0 + 1.0)
        edges = feature.canny(self.data/np.max(self.data), sigma=sigma, low_threshold = low_thresh, high_threshold = high_thresh)
        model_edges = feature.canny(self.model_planet, sigma=sigma, low_threshold = low_thresh, high_threshold = high_thresh)
    
        [dx,dy,dxerr,dyerr] = chi2_shift(model_edges,edges)
        self.x_shift = -dx #need if we want to shift another filter the same amount
        self.y_shift = -dy
        print('Pixel shift X, Y = ', self.x_shift, self.y_shift)
        
        #error in position on surface is approximately equal to projected error times emission angle - in reality there is some asymmetry that would be important if error bars are large
        #These are lat/lon error in the x-hat and y-hat directions; their magnitude is correct for lat-lon space but their direction is not
        self.err_x = dxerr/(self.deg_per_px*self.mu)
        self.err_y = dyerr/(self.deg_per_px*self.mu)
        print('    Lat/Lon error at sub-obs point in x-hat direction = '+str(dxerr/self.deg_per_px))
        print('    Lat/Lon error at sub-obs point in y-hat direction = '+str(dyerr/self.deg_per_px))   
        
        self.centered = shift2d(self.data,-1*dx,-1*dy)
        self.edges = shift2d(edges,-1*dx,-1*dy)

        if plot:
            fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(10, 5))
            
            ax0.imshow(self.data, origin = 'lower left')
            ax0.set_title('Image')
            
            ax1.imshow(edges, origin = 'lower left')
            ax1.set_title('Canny filter, $\sigma=$%d'%sigma)
            
            ax2.imshow(self.edges, origin = 'lower left', alpha = 0.5)
            ax2.imshow(model_edges, origin = 'lower left', alpha = 0.5)
            ax2.set_title('Overlay model and data')
            
            plt.show()
    
    
    def edge_detect_error(self, niter, perturb_l, perturb_dist, low_thresh = 0.002, dist = 0.02, sigma = 5, doplot = True):
        '''Perturbs parameters of Canny algorithm to produce a variety of 
        edge detection solutions. Finds most probable one, and takes 
        standard deviation of those to produce an error
        perturb_l, perturb_dist        Factor by which low_thresh and distance between low and high threshold are changed. must be >= 1
        sigmavals                      List of sigma values to use. Usually some subset of [3,5,7]
        niter                          Number of iterations for each low, high threshold value'''
        
        #set up the model planet
        self.model_planet = np.nan_to_num(self.lat_g * 0.0 + 1.0)
        model_edges = feature.canny(self.model_planet, sigma=sigma, low_threshold = low_thresh, high_threshold = low_thresh + dist)

        if doplot:
            inp = False
            while not inp:
            
                #set up arrays of values
                l_vals = np.arange(low_thresh/perturb_l, low_thresh*perturb_l, (low_thresh*perturb_l - low_thresh/perturb_l)/niter)
                dist_vals = np.arange(dist/perturb_dist, dist*perturb_dist, (dist*perturb_dist - dist/perturb_dist)/niter)
                
                #check that lowest, middle, and highest params give you what is expected
                l_bounds = [l_vals[0], low_thresh, l_vals[-1]]
                d_bounds = [dist_vals[0], dist, dist_vals[-1]]
                fig, axes = plt.subplots(3, 3, figsize=(8, 12), sharex = True, sharey = True)
                
                for i in range(3):
                    #do the edges 
                    l = l_bounds[i]
                    d = d_bounds[i]
                    edges = feature.canny(self.data/np.max(self.data), sigma=sigma, low_threshold = l, high_threshold = l + d)
                    [dx,dy,dxerr,dyerr] = chi2_shift(model_edges,edges)
                    shift_edges = shift2d(edges,-1*dx,-1*dy)
                    
                    #plot things
                    axarr = axes[i]
                    axarr[0].imshow(self.data, origin = 'lower left')
                    axarr[1].imshow(edges, origin = 'lower left')
                    axarr[2].imshow(shift_edges, origin = 'lower left', alpha = 0.5)
                    axarr[2].imshow(model_edges, origin = 'lower left', alpha = 0.5) 
                    
                    #cosmetics
                    if i == 0:
                        axarr[0].set_title('Image')
                        axarr[1].set_title('Edges, $\sigma=$%d'%sigma)
                        axarr[2].set_title('Overlay model and data')
                    
                    for ax in axarr:
                        ax.set_xticklabels([])
                        ax.set_yticklabels([])
                        ax.set_ylim([0, edges.shape[1]])
                        ax.set_xlim([0, edges.shape[0]])
                        
                    axarr[0].set_ylabel('l = '+str(l)[:6]+', d = '+str(d)[:6])
                
                plt.subplots_adjust(wspace = 0, hspace = 0)
                plt.show()          
            
                yn = input('Are you okay with these? (y/n): ')
                if yn.lower().strip() == 'y' or yn.lower().strip() == 'yes':
                    inp = True
                else:
                    low_thresh = float(input('New value of low_thresh: '))
                    dist = float(input('New value of dist: '))
                    perturb_l = float(input('New value of perturb_l: '))
                    perturb_dist = float(input('New value of perturb_dist: '))
        else:
            l_vals = np.arange(low_thresh/perturb_l, low_thresh*perturb_l, (low_thresh*perturb_l - low_thresh/perturb_l)/niter)
            dist_vals = np.arange(dist/perturb_dist, dist*perturb_dist, (dist*perturb_dist - dist/perturb_dist)/niter)
        
        #now iterate over the parameter space and determine x,y for each
        dx_vals, dy_vals = [], []
        dxerr_vals , dyerr_vals = [], []
        for lt in l_vals:
            for d in dist_vals:
                ht = lt + d
                edges = feature.canny(self.data/np.max(self.data), sigma=sigma, low_threshold = lt, high_threshold = ht)
                [dx,dy,dxerr,dyerr] = chi2_shift(model_edges,edges)
                dx_vals.append(dx)
                dy_vals.append(dy)
                dxerr_vals.append(dxerr)
                dyerr_vals.append(dyerr)
                #print(d, dx, dy)
                plotevery = False
                if plotevery:
                    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(10, 5))
                    
                    ax0.imshow(self.data, origin = 'lower left')
                    ax0.set_title('Image')
                    
                    ax1.imshow(edges, origin = 'lower left')
                    ax1.set_title('Canny filter, $\sigma=$%d'%sigma)
                    
                    ax2.imshow(shift2d(edges,-1*dx,-1*dy), origin = 'lower left', alpha = 0.5)
                    ax2.imshow(model_edges, origin = 'lower left', alpha = 0.5)
                    ax2.set_title('Overlay model and data') 
                    
                    plt.show()                   
        dx_vals = np.asarray(dx_vals)
        dy_vals = np.asarray(dy_vals)
        dxerr_vals = np.asarray(dxerr_vals)
        dyerr_vals = np.asarray(dyerr_vals)
        print('Best X, Y = ', np.mean(dx_vals), np.mean(dy_vals))
        print('Sigma X, Y = ', np.std(dx_vals), np.std(dy_vals))
        print('Typical shift error X, Y = ', np.mean(dxerr_vals), np.mean(dyerr_vals))
        print('Spread in shift error X, Y = ', np.std(dxerr_vals), np.std(dyerr_vals))
        print('Total error X, Y = ', np.sqrt(np.std(dx_vals)**2 + np.mean(dxerr_vals)**2), np.sqrt(np.std(dy_vals)**2 + np.mean(dyerr_vals)**2))
        
        if doplot:
            #histogram of the samples
            fig, (ax0, ax1) = plt.subplots(2,1, figsize = (6,8))
            
            ax0.hist(dx_vals, bins = niter)
            ax1.hist(dy_vals, bins = niter)
            
            ax0.set_xlabel('Shift in X')
            ax0.set_ylabel('Number')
            ax1.set_xlabel('Shift in Y')
            ax1.set_ylabel('Number')
            plt.show()
        
    
    def manual_shift(self,dx,dy):
        self.centered = shift2d(self.data,dx,dy)
        
    def plot_latlon(self):
        '''Make pretty plot of lat_g and lon_w overlaid on planet'''
        fig, (ax0, ax1) = plt.subplots(1,2, figsize = (12,6))
        
        #little circle around planet - now does not depend on self.edges existing
        planetedge = np.copy(self.lat_g)
        nans = np.isnan(planetedge)
        planetedge[np.invert(nans)] = 100
        planetedge[nans] = 0
        
        #latitudes
        ax0.imshow(self.centered, origin = 'lower left')
        levels_lat = np.arange(-90,105,15)
        label_levels_lat = np.arange(-90,60,30)
        ctr_lat = ax0.contour(self.lat_g, levels_lat, colors='white', linewidths=2)
        ax0.clabel(ctr_lat, label_levels_lat, inline=1, inline_spacing = 2, fontsize=16, fmt='%d')
        ax0.contour(planetedge, colors = 'white', linewidths = 1)
        #ax0.set_title('Latitudes', fontsize = 18)
        ax0.get_xaxis().set_ticks([])
        ax0.axes.get_yaxis().set_ticks([])
        
        #longitudes
        ax1.imshow(self.centered, origin = 'lower left')
        #hack here to avoid discontinuity in contours - split longs in half
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") #suppresses error for taking < nan
            lon_w1 = np.copy(self.lon_w)
            lon_w1[lon_w1 >= 180] = np.nan
            lon_w2 = np.copy(self.lon_w)
            lon_w2[lon_w2 < 180] = np.nan
        
        levels_lon = range(0,360,30)
        levels_lon_hack = [1] + list(levels_lon[1:]) #make contour at zero actually 1 - otherwise won't plot it since it's at the edge
        ctr_lon1 = ax1.contour(lon_w1, levels_lon_hack, colors='white', linewidths=2)
        ctr_lon2 = ax1.contour(lon_w2, levels_lon_hack, colors='white', linewidths=2)
                
        fmt = {}
        vals = np.arange(0,360,30)
        for l, v in zip(levels_lon_hack, vals):
            fmt[l] = str(int(v)) #make it so the labels say the right things despite hack
        ax1.clabel(ctr_lon1, levels_lon_hack, fmt = fmt, inline=1, inline_spacing = 2, fontsize=16)
        ax1.clabel(ctr_lon2, levels_lon_hack, fmt = fmt, inline=1, inline_spacing = 2, fontsize=16)
        ax1.contour(planetedge, colors = 'white', linewidths = 1)
        #ax1.set_title('Longitudes', fontsize = 18)
        ax1.get_xaxis().set_ticks([])
        ax1.axes.get_yaxis().set_ticks([])        
                
        plt.tight_layout()
        plt.savefig('lat_lon_overlay.png')
        plt.show()
        
    def write(self, lead_string):
        '''Tertiary data products'''
        hdulist_out = self.im.hdulist
        #centered data
        hdulist_out[0].header['OBJECT'] = self.target+'_CENTERED'
        hdulist_out[0].data = self.centered
        hdulist_out[0].writeto(lead_string + '_centered.fits', overwrite=True)
        #latitudes
        hdulist_out[0].header['OBJECT'] = self.target+'_LATITUDES'
        hdulist_out[0].data = self.lat_g
        hdulist_out[0].writeto(lead_string + '_latg.fits', overwrite=True)
        #longitudes
        hdulist_out[0].header['OBJECT'] = self.target+'_LONGITUDES'
        hdulist_out[0].data = self.lon_w
        hdulist_out[0].writeto(lead_string + '_lone.fits', overwrite=True)
        #errors only exist if edge_detect was run. if manual shift, just ignore
        try:
            #error in x*mu
            hdulist_out[0].header['OBJECT'] = self.target+'_XERR'
            hdulist_out[0].data = self.err_x
            hdulist_out[0].writeto(lead_string + '_errx.fits', overwrite=True)
            #error in y*mu
            hdulist_out[0].header['OBJECT'] = self.target+'_YERR'
            hdulist_out[0].data = self.err_y
            hdulist_out[0].writeto(lead_string + '_erry.fits', overwrite=True)
        except:
            pass

    def bootstrap_func(self, order = 2):
        '''Takes a navigated image, plots flux as function of emission angle,
        fits (to nth order) to the minimum flux vs emission angle curve.
        returns the fit coefficients to IoverF(mu)'''
    
        onplanet = np.copy(self.centered)
        onplanet *= self.model_planet
    
        vals = onplanet.flatten()
        mus = self.mu.flatten()
    
        vals_2 = vals[vals > 0] #remove zeros on outside of image
        mus_2 = mus[vals > 0]
        
        np.savetxt('flux_vs_mu.txt', np.asarray([mus_2, vals_2]))
    
        '''It looks like the minimum value at each emission angle follows a fairly
        regular distribution. Try to make a function fit the bottom of it.'''
        bins = np.arange(0,1.01,0.01)
        x = np.digitize(mus_2, bins)
        mins = [np.min(vals_2[np.where(x == i)]) if len(vals_2[np.where(x == i)]) > 0 else 0.0 for i in range(bins.shape[0])]
        mins = np.asarray(mins)
        bins = bins[np.where(mins > 0)]
        mins = mins[np.where(mins > 0)]
        
        z = np.polyfit(bins,mins, order)
        func = np.poly1d(z)
    
        plt.semilogy(mus_2, vals_2, linestyle = '', marker = '.', color = 'k', markersize = 1)
        plt.semilogy(bins, func(bins), color = 'r')
        #plt.semilogy(bins, z[0]*bins**2 + z[1]*bins*1.1 + z[2])
        plt.xlabel(r'Emission angle $\mu$')
        plt.ylabel('I/F')
        plt.show()
    
        print('Polynomial parameters ... + A2*x**2 + A1*x + A2 = ',z)
        return z
        
    def locate_feature(self, outfile=None):
        plt.imshow(self.centered, origin = 'lower left')
        plt.show()
        print('Define a box around the feature you want to track. Note x,y are reversed in image due to weird Python indexing!')
        pix_l = input('Enter lower left pixel x,y separated by a comma: ')
        pix_u = input('Enter upper right pixel x,y separated by a comma: ')
        
        p0x, p0y = int(pix_l.split(',')[0].strip(', \n')),int(pix_l.split(',')[1].strip(', \n'))
        p1x, p1y = int(pix_u.split(',')[0].strip(', \n')),int(pix_u.split(',')[1].strip(', \n'))
        region = self.centered[p0x:p1x,p0y:p1y]        

        #Brightest spot in feature
        maxloc = np.where(self.centered == np.max(region))
        maxlat, maxlon = self.lat_g[maxloc], self.lon_w[maxloc]
        
        #Gaussian fit
        A0 = np.sum(region)
        x0, y0 = (p1x - p0x)/2, (p1y - p0y)/2
        x_std0, y_std0 = region.shape[0]/2, region.shape[1]/2
        theta0 = 0.0
        g_init = models.Gaussian2D(A0, x0, y0, x_std0, y_std0, theta0)
        xx, yy = np.mgrid[:region.shape[0], :region.shape[1]]
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, xx, yy, region)
        
        xfg, yfg = g.x_mean + p0x, g.y_mean + p0y
        latfg, lonfg = self.lat_g[int(round(xfg)),int(round(yfg))], self.lon_w[int(round(xfg)), int(round(yfg))]
        
        '''   #estimate lat and lon errors
        frac = 0.5 #can probably constrain the center much more, but depends on morphology of storm over time
        fracmax = np.where(g_overlay/np.max(g_overlay) > frac)
        inlat, inlon = self.lat_g[fracmax], self.lon_w[fracmax]
        minlat, maxlat = np.min(inlat), np.max(inlat)
        minlon, maxlon = np.min(inlon), np.max(inlon) #wrapping issues still present!
        
        lat_errl, lat_erru = np.abs(latf - minlat), np.abs(maxlat - latf) 
        lon_wrrl, lon_wrru = np.abs(lonf - minlon), np.abs(maxlon - lonf)'''

        #Contour method after Martin, de Pater, Marcus 2012
        def ctr_region(rgn, frac):
            rgn = np.copy(rgn)
            rgn[rgn < frac*np.max(rgn)] = 0.0
            rgn[rgn > 0.0] = 1.0
            xf, yf = center_of_mass(rgn)
            return xf, yf
            
        def interp_latlon_atpt(xd,yd):
            x, y = int(round(xd)), int(round(yd))
            xcoords, ycoords = np.arange(x-3, x+4), np.arange(y-3, y+4)
            latgrid = self.lat_g[x-3:x+4,y-3:y+4]
            longrid = self.lon_w[x-3:x+4,y-3:y+4]
            #RectBivariateSpline and interp2d both fail if lat/lon grid has NaNs, i.e. if near edge of planet
            if np.any(np.isnan(latgrid[2:5,2:5])):
                #Check to see if you are way too close
                print('    WARNING! DANGER! ERROR! HELP! OMG!')
                print('            Trying to find pixel locations very near edge of planet')
                print('            Lat-lon errors are probably wrong!!')
            latgrid[np.isnan(latgrid)] = 0.0
            longrid[np.isnan(longrid)] = 0.0
            interpf_lat = RectBivariateSpline(xcoords, ycoords, latgrid)
            interpf_lon = RectBivariateSpline(xcoords, ycoords, longrid)
            return interpf_lat.ev(xd, yd), interpf_lon.ev(xd, yd)
            #interpf_lat = interp2d(xcoords, ycoords, latgrid, kind = 'cubic')
            #interpf_lon = interp2d(xcoords, ycoords, longrid, kind = 'cubic')
            #return interpf_lat(xd, yd), interpf_lon(xd, yd)
            
        def ctr_fit(rgn, frac):
            (px, py) = ctr_region(rgn, frac)
            xfc, yfc = px + p0x, py + p0y
            #latfc, lonfc = self.lat_g[int(round(xfc)),int(round(yfc))], self.lon_w[int(round(xfc)), int(round(yfc))] #old way -  fails if error is <~ 1 pixel
            latfc, lonfc = interp_latlon_atpt(xfc, yfc)
            return xfc, yfc, latfc, lonfc
        
        poslist = []    
        for val in np.arange(0.68, 0.96, 0.01):
            xfc, yfc, latfc, lonfc = ctr_fit(region, val)
            poslist.append([xfc, yfc, latfc, lonfc])
        poslist = np.asarray(poslist)
        medposlist = np.median(poslist, axis = 0)
        meanposlist = np.mean(poslist, axis = 0)
        minposlist = np.min(poslist, axis = 0)
        maxposlist = np.max(poslist, axis = 0)
        sigmaposlist = np.std(poslist, axis = 0)            
        ## there will be longitude wrapping issues here later - will need to fix at some point

        print('Error in X,Y,lat,lon', sigmaposlist)
        print('Best X,Y,lat,lon (contour method)', medposlist)
        
        #write to file
        if outfile != None:
            print('Writing detailed outputs to file.')
            rows = ['#      ','X_pixel', 'Y_pixel', 'Latitude', 'Longitude']
            columns = ['Brightest', 'Contour_Median','Contour_Mean','Contour_Sigma','Contour_Min','Contour_Max','Gaussian_Center','Gaussian_Sigma','Edge_Detect_Err']
            with open(outfile, 'w') as f:
                f.write('#Finding center of extended feature with a few techniques.\n')
                f.write('#Note that edge detection error in lat, lon is actually edge detection lat/lon err in x-hat, y-hat\n')
                f.write('   '.join(rows)+'\n')
                f.write(columns[0]+'   '+str(maxloc[0])+'   '+str(maxloc[1])+'   '+str(maxlat[0])+'   '+str(maxlon[0])+'\n')
                f.write(columns[1]+'   '+str(medposlist[0])+'   '+str(medposlist[1])+'   '+str(medposlist[2])+'   '+str(medposlist[3])+'\n')
                f.write(columns[2]+'   '+str(meanposlist[0])+'   '+str(meanposlist[1])+'   '+str(meanposlist[2])+'   '+str(meanposlist[3])+'\n')
                f.write(columns[3]+'   '+str(sigmaposlist[0])+'   '+str(sigmaposlist[1])+'   '+str(sigmaposlist[2])+'   '+str(sigmaposlist[3])+'\n')
                f.write(columns[4]+'   '+str(minposlist[0])+'   '+str(minposlist[1])+'   '+str(minposlist[2])+'   '+str(minposlist[3])+'\n')
                f.write(columns[5]+'   '+str(maxposlist[0])+'   '+str(maxposlist[1])+'   '+str(maxposlist[2])+'   '+str(maxposlist[3])+'\n')
                f.write(columns[6]+'   '+str(xfg)+'   '+str(yfg)+'   '+str(latfg)+'   '+str(lonfg)+'\n')
                f.write(columns[7]+'   '+str(g.x_stddev + 0)+'   '+str(g.y_stddev + 0)+'   '+str(-999)+'   '+str(-999)+'\n')
                f.write(columns[8]+'   '+str(-999)+'   '+str(-999)+'   '+str(self.err_x[int(round(medposlist[0])),int(round(medposlist[1]))])+'   '+str(self.err_y[int(round(medposlist[0])),int(round(medposlist[1]))])+'\n')
        
        #replace Gaussian onto grid we had before for plotting
        eval_g = g(xx, yy)
        g_overlay = np.zeros(self.centered.shape)
        g_overlay[p0x:p1x,p0y:p1y] = eval_g
        #plot things up                                     
        fig, (ax0) = plt.subplots(1,1, figsize = (8,5))
        ax0.imshow(region, origin = 'lower left')
        levels = [0.68, 0.815, 0.95]
        cs = ax0.contour(region/np.max(region), levels, colors = 'white')
        cs_g = ax0.contour(eval_g/np.max(eval_g), [0.68], colors = 'red')
        plt.show()                                 
        
           
    def project(self, outstem = 'h', pixsz = None, interp = 'cubic'):
        '''Project the data onto a flat x-y grid.
        pixsz is in arcsec. if pixsz = None, translates the pixel scale
        of the image to a distance at the sub-observer point.
        interp asks whether to regrid using a nearest neighbor, linear, or cubic'''
        
        #determine the number of pixels in resampled image
        if pixsz == None:
            pixsz = self.pixscale_arcsec
        npix_per_degree = (1/self.deg_per_px) * (self.pixscale_arcsec / pixsz) # (old pixel / degree lat) * (arcsec / old pixel) / (arcsec / new pixel) = new pixel / degree lat
        npix = int(npix_per_degree * 180) + 1 #(new pixel / degree lat) * (degree lat / planet) = new pixel / planet
        print('New image will be %d by %d pixels'%(2*npix + 1, npix))
        print('Pixel scale %f km = %f pixels per degree'%(self.pixscale_km, npix_per_degree))
        
        #create new lon-lat grid
        extra_wrap_dist = 180
        newlon, newlat = np.arange(-extra_wrap_dist,360 + extra_wrap_dist, 1/npix_per_degree), np.arange(-90,90, 1/npix_per_degree)
        gridlon, gridlat = np.meshgrid(newlon, newlat)
        nans = np.isnan(self.lon_w.flatten())
        def input_helper(arr, nans):
            '''removing large region of NaNs speeds things up significantly'''
            return arr.flatten()[np.logical_not(nans)]
        inlon, inlat, indat = input_helper(self.lon_w, nans), input_helper(self.lat_g, nans), input_helper(self.centered, nans)

        #fix wrapping by adding dummy copies of small lons at > 360 lon
        inlon_near0 = inlon[inlon < extra_wrap_dist]
        inlon_near0 += 360
        inlon_near360 = inlon[inlon > 360 - extra_wrap_dist]
        inlon_near360 -= 360
        inlon_n = np.concatenate((inlon_near360, inlon, inlon_near0))
        inlat_n = np.concatenate((inlat[inlon > 360 - extra_wrap_dist], inlat, inlat[inlon < extra_wrap_dist]))
        indat_n = np.concatenate((indat[inlon > 360 - extra_wrap_dist], indat, indat[inlon < extra_wrap_dist]))

        #do the regridding
        datsort = griddata((inlon_n, inlat_n), indat_n, (gridlon, gridlat), method = interp)
        
        #trim extra data we got from wrapping
        wrap_i_l = len(gridlon[0][gridlon[0] < 0]) - 1
        wrap_i_u = len(gridlon[0][gridlon[0] >= 360])
        datsort = datsort[:,wrap_i_l:-wrap_i_u]
        gridlon = gridlon[:,wrap_i_l:-wrap_i_u]
        gridlat = gridlat[:,wrap_i_l:-wrap_i_u]
        
        # make far side of planet into NaNs
        snorm = surface_normal(gridlat, gridlon, self.ob_lon)
        #emang = emission_angle(self.ob_lat, snorm).T
        emang = emission_angle(self.ob_lat, snorm)
        farside = np.where(emang < 0.0)
        datsort[farside] = np.nan
        self.projected = datsort
        self.mu_projected = emang
        
        #write data to fits file    
        hdulist_out = self.im.hdulist
        ## projected data
        hdulist_out[0].header['OBJECT'] = self.target+'_projected'
        hdulist_out[0].data = datsort
        hdulist_out[0].writeto(outstem + '_proj.fits', overwrite=True)
        ## emission angles
        hdulist_out[0].header['OBJECT'] = self.target+'_mu_proj'
        hdulist_out[0].data = emang
        hdulist_out[0].writeto(outstem + '_mu_proj.fits', overwrite=True)
        print('Writing files %s'%outstem + '_proj.fits and %s'%outstem + '_mu_proj.fits')
        
    def plot_projected(self, outfname, ctrlon = 180, lat_limits = [-90, 90], lon_limits = [0, 360], cbarlabel = 'I/F'):
        '''Once projection has been run, plot it using this function'''  
        
        #apply center longitude to everything
        npix = self.projected.shape[1]
        npix_per_degree = 1.0 / self.deg_per_px
        print(npix_per_degree)
        offset = (ctrlon + 180)%360
        offsetpix = np.round(offset*npix_per_degree)
        uoffsetpix = npix - offsetpix
        newim = np.copy(self.projected)
        lefthalf = self.projected[:,:offsetpix]
        righthalf = self.projected[:,offsetpix:]
        newim[:,uoffsetpix:] = lefthalf #switch left and right halves
        newim[:,:uoffsetpix] = righthalf
        
        #extent = [ctrlon - 180, ctrlon + 180, -90, 90]
        extent = [ctrlon + 180, ctrlon - 180, -90, 90]
        parallels = np.arange(lat_limits[0],lat_limits[1] + 30, 30.)
        #meridians = np.arange(lon_limits[0],lon_limits[1] + 60, 60.)
        meridians = np.arange(lon_limits[1], lon_limits[0] - 60, -60.)
          
        #plot it
        fs = 14 #fontsize for plots
        fig, ax0 = plt.subplots(1,1, figsize = (10,7))
        
        cim = ax0.imshow(np.fliplr(newim), origin = 'lower left', cmap = 'gray', extent = extent)
        for loc in parallels:
            ax0.axhline(loc, color = 'cyan', linestyle = ':')
        for loc in meridians:
            ax0.axvline(loc, color = 'cyan', linestyle = ':')

        ax0.set_xlabel('Longitude (W)', fontsize = fs)
        ax0.set_ylabel('Latitude', fontsize = fs)
        ax0.set_ylim(lat_limits)
        ax0.set_xlim(lon_limits[::-1])
        ax0.set_title(self.date_time, fontsize = fs + 2)
        ax0.tick_params(which = 'both', labelsize = fs - 2)
        
        #plot the colorbar
        divider = make_axes_locatable(ax0)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(cim, cax = cax, orientation = 'vertical')
        cbar.set_label(cbarlabel, fontsize = fs)
        cax.tick_params(which = 'both', labelsize = fs - 2)
        
        plt.savefig(outfname, bbox = None)
        plt.show()
        
    def feature_size_projected(self):
        '''Find a feature, tell how big it is in lat-lon space'''
        
        if not hasattr(self, 'projected'):
            print('Must run on projected data. Run coords.project() first. Returning')
            return
        
        plt.imshow(self.projected, origin = 'lower left')
        plt.show()
        print('Define a box around the feature you want to track. Note x,y are reversed in image due to weird Python indexing!')
        pix_l = input('Enter lower left pixel x,y separated by a comma: ')
        pix_u = input('Enter upper right pixel x,y separated by a comma: ')
        
        p0x, p0y = int(pix_l.split(',')[0].strip(', \n')),int(pix_l.split(',')[1].strip(', \n'))
        p1x, p1y = int(pix_u.split(',')[0].strip(', \n')),int(pix_u.split(',')[1].strip(', \n'))
        region = self.projected[p0x:p1x,p0y:p1y]
        
        def build_contour(rgn, frac):
            rgn = np.copy(rgn)
            rgn[rgn < frac*np.max(rgn)] = 0.0
            rgn[rgn > 0.0] = 1.0
            return rgn
            
        level = 0.5
        fwhm_2d = build_contour(region, level)
        
        # find the longest rays in the x and y direction across the amorphous fwhm region
        collapse_x = np.sum(fwhm_2d, axis = 1)
        collapse_y = np.sum(fwhm_2d, axis = 0)
        fwhm_x, wherefwhm_x = np.max(collapse_x), np.argmax(collapse_x)
        fwhm_y, wherefwhm_y = np.max(collapse_y), np.argmax(collapse_y)
        
        # convert to lat-lon
        deg_lat, deg_lon = fwhm_x * self.deg_per_px, fwhm_y * self.deg_per_px
        km_lat, km_lon = fwhm_x * self.pixscale_km, fwhm_y * self.pixscale_km
        print('%f degrees lat, %f degrees lon'%(deg_lat, deg_lon))
        print('%f km in zonal direction, %f km in meridional direction'%(km_lat, km_lon))
        
        # for plotting, find min point and max point for each of these rays
        ray_x = np.where(fwhm_2d[wherefwhm_x, :] == 1)[0]
        ray_y = np.where(fwhm_2d[:, wherefwhm_y] == 1)[0]
        
        plt.imshow(region, origin = 'lower left', cmap = 'gray')
        plt.contour(fwhm_2d, levels = [level], colors = ['red'])
        plt.plot(ray_x, np.full(ray_x.shape, wherefwhm_x), color = 'r', lw = 2)
        plt.plot(np.full(ray_y.shape, wherefwhm_y), ray_y, color = 'r', lw = 2)
        plt.title('%s'%self.date_time[:-6])
        plt.savefig('storm_size.png')
        plt.show()
        
        
        
    def help(self):
        
        helpstr = '''
        Contains tasks for image navigation, image projection, calculating I/F, and
        other things that require Horizons ephemeris data to compute
        
        Functions (see DOCUMENTATION.py for use):
            ioverf(self, filt, flux_per, stand_airmass)
            write_photonly(self, outstr)
            edge_detect(self, low_thresh = 0.01, high_thresh = 0.05, sigma = 5, plot = True)
            edge_detect_error(self, niter, perturb_l, perturb_dist, low_thresh = 0.002, dist = 0.02, sigma = 5, doplot = True)
            manual_shift(self,dx,dy)
            plot_latlon(self)
            write(self, lead_string)
            bootstrap_func(self, order = 2)
            locate_feature(self, outfile=None)
            project(self, outstem = 'h', pixsz = None, interp = 'cubic')
            plot_projected(self, outfname, ctrlon = 180, lat_limits = [-90, 90], lon_limits = [0, 360], cbarlabel = 'I/F')
             
        Attributes:
            im: image object for infile
            req: equatorial radius
            rpol: polar radius
            data: the image data
            pixscale_arcsec: pixel scale of image in arcsec
            target: name of planet
            date_time: datetime object for start time
            airmass:
            ang_diam:
            ob_lon: sub observer longitude
            ob_lat: sub observer latitude
            sun_lon: sub solar longitude
            sun_lat: sub solar latitude
            np_ang: north pole position angle (CCW, or east, wrt celestial north)
            np_dist: angular dist of np from sub observer point
            sun_dist: distance from target to sun
            dist: distance from target to observer       
            pixscale_km: pixel scale of image in km
            deg_per_px: lat/lon degrees on planet per pixel at sub obs point
            lat_g: grid of planetographic latitudes on planet
            lat_c: grid of planetocentric latitudes on planet
            lon_w: grid of west longitudes on planet (west for IAU standard, as output by Horizons)
            err_x: lat/lon error from navigation in x direction, in px
            err_y: lat/lon error from navigation in y direction, in px
            centered: centered data after edge detection
            model_planet: ellipse with planet's projected shape based on ephemeris
            projected: image on regular lat/lon grid after projection
            mu: grid of emission angles on planet
            surf_n: grid of normal vectors to the surface of planet
            x_shift: pixel shift in x direction from edge detect
            y_shift: pixel shift in x direction from edge detect
            edges: map of edges found by edge detect
        '''
        print(helpstr)
        
        
        
        
        
        ''' 
    def change_projection(self):
        not written yet but would use basemap, or Charles Goullaud replacement, to reproject arbitrarily
        
        '''        
        