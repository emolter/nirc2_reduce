Sample workflow:

## start Python 3 in top-level directory, where raw and reduced are subfolders 

## setup ##
from nirc2_reduce import sort_rawfiles, bxy3, phot, flats
import matplotlib.pyplot as plt
date = '2017sep27'
outdir = 'reduced/'+date+'/'
filt_name = 'h'
target_name = 'Neptune'
filenames_all = sort_rawfiles.find_object('Neptune',date)

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
fnames = sort_rawfiles.find_filter(filenames_all, filt_name)
print(fnames)
#ensure these files are a bxy3 before continuing
obs = bxy3.Bxy3(fnames)
obs.make_sky(outdir+'sky_'+filt_name+'.fits')
obs.apply_sky(outdir+'sky_'+filt_name+'.fits')
obs.apply_flat(outdir+'flat_master_'+filt_name+'.fits')
obs.apply_badpx_map(outdir+'badpx_map_'+filt_name+'.fits')
obs.dewarp()
obs.trim()
obs.remove_cosmic_rays() # at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; possibly failing to find them in flats with bad pixel search. Perhaps multi-layer search, e.g. large blocksize first, small blocksize second
obs.per_second()
#to check anything, plot obs.frames
plt.imshow(obs.frames[0], origin = 'lower left')
plt.show()

## now if there is no photometry, do the following ##
obs.write_frames([outdir+'frame0_nophot_'+filt_name+'.fits',outdir+'frame1_nophot_'+filt_name+'.fits',outdir+'frame2_nophot_'+filt_name+'.fits'])
obs.stack()
obs.crop(50)
#check it
plt.imshow(obs.final, origin = 'lower left')
plt.show()
obs.write_final(outdir+'stacked_nophot_'+filt_name+'.fits',png=True, png_file = outdir+target_name+'_'+filt_name+'.png')


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
stand.reduce(outdir+'sky_'+keck_filt+'.fits', outdir+'flat_master_'+keck_filt+'.fits', outdir+'badpx_map_'+keck_filt+'.fits', [outdir+'phot0_'+keck_filt+'.fits', outdir+'phot1_'+keck_filt+'.fits', outdir+'phot2_'+keck_filt+'.fits'])
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







## tertiary data products ##
# start Python 3 in directory containing stacked image
from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('stacked_nophot_h.fits')
#apply final phot at this step
coords.ioverf(filt, flux_per, stand_airmass)
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

# if you already did the edge detection and ran write_latlon, then want to reopen:
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.plot_latlon()





# extracting I/F values from a stack of images
from nirc2_reduce import imstack
import matplotlib.pyplot as plt
stack = imstack.Stack(['h_centered.fits', 'kp_centered.fits', 'ch4s_centered.fits', 'pabeta_centered.fits'])
stack.extract_point(200,200)
plt.imshow(stack.ratios['kp/h'], origin = 'lower left')
plt.show()







test
from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.locate_feature(outfile = 'locate.txt')