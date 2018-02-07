Sample workflow:

## start Python 3 in top-level directory, where raw and reduced are subfolders 

## setup ##
from nirc2_reduce import sort_rawfiles, bxy3, phot, flats
import matplotlib.pyplot as plt
date = '2017sep27'
outdir = 'reduced/'+date+'/'
filt_name = 'h'
target_name = 'Neptune'

## flats ##
domeflatoff, domeflaton = sort_rawfiles.get_flats(filt_name, date)
flat = flats.Flats(domeflatoff, domeflaton)
flat.write(outdir+'flat_master_'+filt_name+'.fits')
tol = 0.07 #how far from median value can a pixel be before it's bad?
blocksize = 4 #must be divisible into 1024. how big of a box to take the median?
flat.make_badpx_map(outdir+'badpx_map_'+filt_name+'.fits',tol,blocksize)

## check that these things look nice ##
plt.imshow(flat.badpx_map, origin = 'lower left')
plt.show()
flat.plot()

## if no flats were taken, get from a nearby date and just modify apply_badpx_map and apply_flat

## bxy3 ##
flat_filt = 'kp'
filenames_all = sort_rawfiles.find_object(target_name,date)
fnames = sort_rawfiles.find_filter(filenames_all, filt_name)
print(fnames)
#ensure these files are a bxy3 before continuing
obs = bxy3.Bxy3(fnames)
obs.make_sky(outdir+'sky_'+filt_name+'.fits')
obs.apply_sky(outdir+'sky_'+filt_name+'.fits')
obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
obs.dewarp()
obs.trim()
obs.remove_cosmic_rays() # at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; possibly failing to find them in flats with bad pixel search. Perhaps multi-layer search, e.g. large blocksize first, small blocksize second
obs.per_second()
#to check any step use:
obs.plot_frames()


## multiple bxy3 at once ##
from nirc2_reduce import multi_reduce
filt_list = ['lp','ms','pah','kcont','h2o_ice','br_alpha','bra_cont','hcont']
flatfilt_list = ['kp' for val in range(len(filt_list))]
date = '2017may27'
target_name = 'Io'
multi_reduce.multiBxy3(date, target_name, filt_list, flatfilt_list)
multi_reduce.multiBxy3Plot(date, filt_list)

## nod ##
from nirc2_reduce import nod
obs = nod.Nod(skyf, imagef)
obs.apply_sky()
obs.apply_flat(outdir+'flat_master_'+filt_name+'.fits')
obs.apply_badpx_map(outdir+'badpx_map_'+filt_name+'.fits')
obs.dewarp()
obs.remove_cosmic_rays()
obs.per_second()
obs.crop(200)
obs.write(outdir+'stacked_nophot_'+filt_name+'.fits', png = True, png_file = outdir+target_name+'_'+filt_name+'.png')

## now if there is no photometry, do the following ##
obs.write_frames([outdir+'frame0_nophot_'+filt_name+'.fits',outdir+'frame1_nophot_'+filt_name+'.fits',outdir+'frame2_nophot_'+filt_name+'.fits'])
obs.stack()
obs.crop(50)
#check it
plt.imshow(obs.final, origin = 'lower left')
plt.show()
obs.write_final(outdir+'stacked_nophot_'+filt_name+'.fits')#,png=True, png_file = outdir+target_name+'_'+filt_name+'.png')




## if there is photometry ## 
from nirc2_reduce import sort_rawfiles, bxy3, phot, flats
import matplotlib.pyplot as plt
standard_name = 'HD1160' 
star_mag = 7.040
keck_filt = 'kp'
date = '2015dec25'
outdir = 'reduced/'+date+'/'
sky_name = 'raw/'+date+'/n0220.fits'
image_name = 'raw/'+date+'/n0217.fits'

standard_files = sort_rawfiles.find_object(standard_name,date)
stand_fnames = sort_rawfiles.find_filter(standard_files, keck_filt)
#check that these are the ones you want
stand = phot.bxy3Phot(stand_fnames)  OR stand = phot.nodPhot(sky_name, image_name)
stand.reduce(outdir+'photsky_'+keck_filt+'.fits', outdir+'flat_master_'+keck_filt+'.fits', outdir+'badpx_map_'+keck_filt+'.fits', [outdir+'phot0_'+keck_filt+'.fits', outdir+'phot1_'+keck_filt+'.fits', outdir+'phot2_'+keck_filt+'.fits'])
OR 
stand.reduce(outdir+'flat_master_'+keck_filt+'.fits', outdir+'badpx_map_'+keck_filt+'.fits', outdir+'phot_'+keck_filt+'.fits')

stand.find_cts(50, plot = True) #first number is dist, number of pixel radius

stand.find_flux_conversion(star_mag, keck_filt)
print(stand.flux_per)


## applying photometry ##
# set up obs as reduced frames
obs.apply_photometry_frames(stand.flux_per)
obs.write_frames([outdir+'frame0_'+filt_name+'.fits',outdir+'frame1_'+filt_name+'.fits',outdir+'frame2_'+filt_name+'.fits'])
obs.stack()
obs.crop(50)
obs.apply_photometry_final(stand.flux_per)
obs.write_final(outdir+'stacked_'+filt_name+'.fits',png=False, png_file = '')




## multiple photometry at once ##
from nirc2_reduce import multi_reduce
filt_list = ['lp','ms','pah','kcont','h2o_ice','br_alpha','bra_cont']
mag_list = [7.31, 7.32, 7.31, 7.294, 7.294, 7.31, 7.31] #remember, these are mags in the standard filter whose wavelength matches most closely with the Keck filter
#in this case [bessell_lp, bessell_m, bessell_lp, 2mass_ks, 2mass_ks, bessell_lp, bessell_lp]
flatfilt_list = ['kp' for val in range(len(filt_list))]
date = '2017may27'
target_name = 'HD106965'
multi_reduce.multiPhot(date, target_name, filt_list, mag_list, flatfilt_list, doplots = True)

#apply multiple photometry to data at once
from nirc2_reduce import multi_reduce
airmass = 1.029
date = '2017jul31'
multi_reduce.multiApplyPhot(date, airmass)




## tertiary data products ##
# start Python 3 in directory containing stacked image
from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_calibrated.fits')
## apply final phot at this step if not yet applied
#coords.ioverf(filt, flux_per, stand_airmass)
coords.edge_detect(low_thresh = 0.01, high_thresh = 0.05, sigma = 5)
#for kp, coords.edge_detect(low_thresh = 0.0001, high_thresh = 0.001, sigma = 5)
#or
coords.manual_shift(xshift,yshift)
coords.plot_latlon()
coords.write('h')

# to manipulate latitudes and longitudes, use coords.lat_g, coords.lon_e
plt.imshow(coords.lat_g, origin = 'lower left')
plt.show()

# to determine locations of bright features at first glance (no deprojection yet)
coords.locate_feature(outfile = 'locate.txt')

# to get the coefficients of a photometry bootstrapping fit (IoverF_min vs mu)
coords.bootstrap_func(order = 2)

# to project onto flat
coords.project(outstem = 'h') #, pixsz = 0.009942)
coords.plot_projected(outstem = 'h')

# if you already did the edge detection and ran write_latlon, then want to reopen:
# loads image into self.centered
# this will also look for a projected image and load that into self.projected
from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.plot_latlon()





# extracting I/F values from a stack of images
from nirc2_reduce import imstack
import matplotlib.pyplot as plt
stack = imstack.Stack(['h_centered.fits', 'kp_centered.fits', 'ch4_short_centered.fits', 'pabeta_centered.fits'])
stack.extract_point(200,200) #look in locate.txt to find center of storm
stack.plot_one('kp')
stack.plot_ratio('kp','h')
stack.write('kp_over_h.fits', 'kp', 'h')




test
from nirc2_reduce import coordgrid, imstack
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.project(outstem = 'h')
coords = coordgrid.CoordGrid('kp_centered.fits', lead_string = 'kp')
coords.project(outstem = 'kp')
coords = coordgrid.CoordGrid('ch4_short_centered.fits', lead_string = 'ch4_short')
coords.project(outstem = 'ch4_short')
coords = coordgrid.CoordGrid('pabeta_centered.fits', lead_string = 'pabeta')
coords.project(outstem = 'pabeta')
stack = imstack.Stack(['h_proj.fits', 'kp_proj.fits', 'ch4_short_proj.fits', 'pabeta_proj.fits'])
stack.write('kp_over_h_proj.fits', 'kp', 'h')
stack.write('h_over_ch4_short_proj.fits', 'h', 'ch4_short')
stack.write('h_over_pabeta_proj.fits', 'h', 'pabeta')

test
from nirc2_reduce import multi_reduce
filt_list = ['h', 'kp', 'ch4_short', 'pabeta', 'he1_a']
flatfilt_list = ['h', 'kp', 'h', 'j', 'j']
mag_list = [7.013, 7.04, 7.013, 6.983, 6.983]
target_name = 'HD1160'
airmass = 1.044
date = '2017jul25'
dist_list = [100, 100, 100, 100, 100]
multi_reduce.multiPhot(date, target_name, filt_list, mag_list, flatfilt_list, doplots = True, dist_list = dist_list)
#multi_reduce.multiApplyPhot(date, airmass)

from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('stacked_nophot_h.fits')
coords.edge_detect(low_thresh = 0.01, high_thresh = 0.05, sigma = 5)
coords.write('h')

from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.project(outstem = 'h')
coords.plot_projected('h_proj.png', ctrlon = 180, lat_limits = [-60, 60], lon_limits = [0, 360])
