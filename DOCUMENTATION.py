'''This contains sample workflows for reducing nirc2 planetary data.

File structure setup:
    -For data reduction
        -Create a "raw" directory with subdirectories of date, e.g.
        raw/2017sep03/n0001.fits
        -Create a directory "reduced" with subdirectories of date, e.g.
        reduced/2017sep03/
        -Start Python in directory that contains "raw" and "reduced" dirs
    -For navigation and analysis navigate to reduced/DATE/ and then start Python

For list of attributes in a class, call class.help()
    Works for Bxy3, bxy3Phot, CoordGrid, Filt, Flat, Image, Stack, Nod, nodPhot

Workflows in this document:
1. Reducing data that were taken with the Keck bxy3 dither pattern
2. Multiple standard bxy3 reductions at once, e.g for Twilight Zone data
3. Data reduction for target/sky nod, e.g. Uranus data
4. Multiple standard nod reductions at once
5. Photometry, bxy3 or nod
6. Multiple photometry at once
7. Navigating an image
8. Converting navigated image to I/F
9. Image projection onto a regular lat/lon grid
10. Extracting information from a stack of images taken in different filters
11. Extracting full-sized individual reduced frames with photometry, e.g. for moon photometry
999. Miscellaneous coordGrid tasks


NOTES FROM ERIN:
dependencies: anaconda, image_registration, astroscrappy (need pip install), pyproj (REMOVED)

tutorial on setup is needed! -- include that if multiple planets, need to re-make reduced/date directory

hardcoded path to service et al. warp solutions - make relative path - problem in bxy3 and nod
hardcoded path to naif IDs in get_ephem
'''

##################################################################################
### 1: Reducing data that were taken with the Keck bxy3 dither pattern
##################################################################################

# setup
# remember, bxy3 assumes top left, bottom right, top right in that order
from nirc2_reduce import sort_rawfiles, bxy3, flats
import matplotlib.pyplot as plt
date = '2017sep27'
input_dir = 'raw/'+date+'/'
outdir = 'reduced/'+date+'/'
filt_name = 'h'
flat_filt = 'h'
target_name = 'Neptune'

# make flats and badpx map
domeflatoff, domeflaton = sort_rawfiles.get_flats(input_dir, filt_name)
# domeflatoff, domeflaton are lists of filenames
flat = flats.Flats(domeflatoff, domeflaton)
flat.write(outdir+'flat_master_'+filt_name+'.fits')
tol = 0.07 #how far from median value can a pixel be before it's bad?
blocksize = 4 #must be divisible into 1024. how big of a box to take the median?
flat.make_badpx_map(outdir+'badpx_map_'+filt_name+'.fits',tol,blocksize)

# check that these things look nice
plt.imshow(flat.badpx_map, origin = 'lower left')
plt.show()
flat.plot()

# if no flats were taken, get from a nearby date and just modify apply_badpx_map and apply_flat

# bxy3
filenames_all = sort_rawfiles.find_object(input_dir, target_name)
fnames = sort_rawfiles.find_filter(filenames_all, filt_name)
fnames = np.sort(fnames)
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

# writing and stacking frames
obs.write_frames([outdir+'frame0_nophot_'+filt_name+'.fits',outdir+'frame1_nophot_'+filt_name+'.fits',outdir+'frame2_nophot_'+filt_name+'.fits'])
obs.stack()
obs.crop(50)
#check it
plt.imshow(obs.final, origin = 'lower left')
plt.show()
obs.write_final(outdir+'stacked_nophot_'+filt_name+'.fits')#,png=True, png_file = outdir+target_name+'_'+filt_name+'.png')


##################################################################################
### 2: Multiple standard bxy3 reductions at once, e.g for Twilight Zone data
##################################################################################

#set up 'reduced/date' directory, put flats in it
from nirc2_reduce import multi_reduce
filt_list = ['lp','ms','pah','kcont','h2o_ice','br_alpha','bra_cont','hcont','jcont']
flatfilt_list = ['kp' for val in range(len(filt_list))]
date = '2017may27'
target_name = 'Io'
multi_reduce.multiBxy3(date, target_name, filt_list, flatfilt_list)
multi_reduce.multiBxy3Plot(date, filt_list)


##################################################################################
### 3: Data reduction for target/sky nod, e.g. Uranus data
##################################################################################

# setup
from nirc2_reduce import nod, sort_rawfiles
import matplotlib.pyplot as plt
date = '2018jan20'
outdir = 'reduced/'+ date + '/'
filt_name = 'ch4_short'
flat_filt = 'h'
target_name = 'Uranus'

# get relevant files from raw
imagefiles_all = sort_rawfiles.find_object(target_name,date)
skyfiles_all = sort_rawfiles.find_object('sky',date)
imagef = sort_rawfiles.find_filter(imagefiles_all, filt_name)
skyf = sort_rawfiles.find_filter(skyfiles_all, filt_name)
print(imagef, skyf) #ensure this is an on-off pair

# nod
obs = nod.Nod(skyf, imagef)
obs.apply_sky()
obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
obs.dewarp()
obs.remove_cosmic_rays()
obs.per_second()
plt.imshow(obs.data, origin = 'lower left')
plt.show()
obs.write(outdir+'stacked_nophot_'+filt_name+'.fits', png = False) #, png = True, png_file = outdir+target_name+'_'+filt_name+'.png')


##################################################################################
### 4: Multiple standard nod reductions at once
##################################################################################

from nirc2_reduce import multi_reduce
filt_list = ['h', 'kp', 'ch4_short', 'pabeta', 'he1_a', 'ch4_long']
flatfilt_list = ['h', 'kp', 'h', 'j', 'j', 'h']
date = '2018jan20'
target_name = 'Uranus'
multi_reduce.multiNod(date, target_name, filt_list, flatfilt_list)
multi_reduce.multiNodPlot(date, filt_list)


##################################################################################
### 5: Photometry, bxy3 or nod
##################################################################################

# setup
from nirc2_reduce import sort_rawfiles, bxy3, phot, flats
import matplotlib.pyplot as plt
standard_name = 'HD1160' 
star_mag = 7.040
keck_filt = 'kp'
date = '2015dec25'
outdir = 'reduced/'+date+'/'
sky_name = 'raw/'+date+'/n0220.fits'
image_name = 'raw/'+date+'/n0217.fits'

# get relevant files from raw
standard_files = sort_rawfiles.find_object(standard_name,date)
stand_fnames = sort_rawfiles.find_filter(standard_files, keck_filt)
print(stand_fnames) #check that these are the ones you want

# bxy3
stand = phot.bxy3Phot(stand_fnames)
stand.reduce(outdir+'photsky_'+keck_filt+'.fits', outdir+'flat_master_'+keck_filt+'.fits', outdir+'badpx_map_'+keck_filt+'.fits', [outdir+'phot0_'+keck_filt+'.fits', outdir+'phot1_'+keck_filt+'.fits', outdir+'phot2_'+keck_filt+'.fits'])
# OR nod 
stand = phot.nodPhot(sky_name, image_name)
stand.reduce(outdir+'flat_master_'+keck_filt+'.fits', outdir+'badpx_map_'+keck_filt+'.fits', outdir+'phot_'+keck_filt+'.fits')

# get photometry
stand.find_cts(50, plot = True) #first number is radius in pixels within which to add up flux
stand.find_flux_conversion(star_mag, keck_filt)
print(stand.flux_per) #flux per second in CGS units, i.e. erg s-1 cm-2 um-1 / cts s-1
# remember to put this number into the photometry_history spreadsheet on Ned's Google Drive

# apply photometry onto the already reduced obs object
# note this is NOT necessary - will calculate I/F when doing image navigation in coordGrid starting with flux_per values - see Workflow 8
obs.apply_photometry_frames(stand.flux_per)
obs.write_frames([outdir+'frame0_'+filt_name+'.fits',outdir+'frame1_'+filt_name+'.fits',outdir+'frame2_'+filt_name+'.fits'])
obs.stack()
obs.crop(50)
obs.apply_photometry_final(stand.flux_per)
obs.write_final(outdir+'stacked_'+filt_name+'.fits',png=False, png_file = '')


##################################################################################
### 6: Multiple photometry at once
##################################################################################

from nirc2_reduce import multi_reduce
filt_list = ['lp','ms','pah','kcont','h2o_ice','br_alpha','bra_cont']
mag_list = [7.31, 7.32, 7.31, 7.294, 7.294, 7.31, 7.31] #remember, these are mags in the standard filter whose wavelength matches most closely with the Keck filter
#in this case [bessell_lp, bessell_m, bessell_lp, 2mass_ks, 2mass_ks, bessell_lp, bessell_lp]
dist_list = [100, 100, 100, 100, 100] # radius in px around which to add up the flux, uses 50 for all if not specified.
flatfilt_list = ['kp' for val in range(len(filt_list))]
date = '2017may27'
target_name = 'HD106965'
multi_reduce.multiPhot(date, target_name, filt_list, mag_list, flatfilt_list, doplots = True, dist_list = dist_list)
# remember to put numbers into the photometry_history spreadsheet on Ned's Google Drive

# apply multiple photometry to data at once
# note this is NOT necessary - will calculate I/F when doing image navigation in coordGrid starting with flux_per values - see Workflow 8
airmass = 1.029
multi_reduce.multiApplyPhot(date, airmass)


##################################################################################
### 7: Navigating an image
##################################################################################

# start in directory containing calibrated images
from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_calibrated.fits', req = 10000, rpol = 9000, scope = 'keck')
coords.edge_detect(low_thresh = 0.01, high_thresh = 0.05, sigma = 5)
#for kp, coords.edge_detect(low_thresh = 0.0001, high_thresh = 0.001, sigma = 5) seemed to work

# if want to determine error on that edge detection:
perturb_l = 3 #Factor by which low_thresh is changed (both + and - direction). must be >= 1
perturb_dist = 3 #Factor by which distance between low and high threshold is changed (both + and - direction). must be >= 1
niter = 10 #number of steps between lower and upper factors for each, i.e. if niter = 10 will run 100 iterations
coords.edge_detect_error(niter, perturb_l, perturb_dist, low_thresh = 0.002, dist = 0.02, sigma = 5, doplot = True)

# plot and write
coords.plot_latlon()
coords.write('h')

# to manipulate latitudes and longitudes, use coords.lat_g, coords.lon_e
plt.imshow(coords.lat_g, origin = 'lower left')
plt.show()

# to manipulate emission angles, use coords.mu
plt.imshow(coords.mu, origin = 'lower left')
plt.show()

# if you already did the edge detection and ran write(), then want to reopen:
# loads image into self.centered
# this will also look for a projected image and load that into self.projected
from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.plot_latlon()


##################################################################################
### 8: Converting navigated image to I/F
##################################################################################

from nirc2_reduce import coordgrid
filt = 'h'
flux_per = 9.3e-17 #erg s-1 cm-2 um-1 / cts s-1
stand_airmass = 1.10
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
coords.ioverf(filt, flux_per, stand_airmass)


##################################################################################
### 9: Image projection onto a regular lat/lon grid
##################################################################################

from nirc2_reduce import coordgrid
coords = coordgrid.CoordGrid('h_centered.fits', lead_string = 'h')
# to project onto flat
coords.project(outstem = 'h') #, pixsz = 0.009942)
coords.plot_projected('h_proj.png', ctrlon = 180, lat_limits = [-60, 60], lon_limits = [0, 360])

# running project() saves it to file automatically
# use coords.projected to access it within the class


##################################################################################
### 10: Extracting information from a stack of images taken in different filters
##################################################################################

from nirc2_reduce import imstack
import matplotlib.pyplot as plt
stack = imstack.Stack(['h_centered.fits', 'kp_centered.fits', 'ch4_short_centered.fits', 'pabeta_centered.fits'])
stack.extract_point(200,200) #look in locate.txt, created by coords.locate_feature (see workflow 999) to find center of a feature
stack.extract_feature(0.7) #make region around feature at some fractional contour level, print the mean in that contour
stack.plot_one('kp')
stack.plot_ratio('kp','h')
stack.write('kp_over_h.fits', 'kp', 'h')


##################################################################################
### 11: Extracting full-sized reduced frames, with photometry
###     Useful for, e.g., validating photometry on moons that are outside the FOV
###     of the bxy3 center but in one of the frames (tested on Proteus)
##################################################################################

from nirc2_reduce import sort_rawfiles, bxy3
date = '2017jul25'
outdir = 'reduced/'+date+'/'
filt_name = 'ch4_short'
target_name = 'Neptune'
flat_filt = 'h'
flux_per = 2.81E-16

filenames_all = sort_rawfiles.find_object(target_name,date)
fnames = sort_rawfiles.find_filter(filenames_all, filt_name)
print(fnames)
obs = bxy3.Bxy3(fnames)
obs.make_sky(outdir+'sky_'+filt_name+'.fits')
obs.apply_sky(outdir+'sky_'+filt_name+'.fits')
obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
obs.dewarp()
obs.remove_cosmic_rays()
obs.per_second()
obs.apply_photometry_frames(flux_per)
obs.plot_frames()
obs.write_frames([outdir+'frame0_full_phot_'+filt_name+'.fits',outdir+'frame1_full_phot_'+filt_name+'.fits',outdir+'frame2_full_phot_'+filt_name+'.fits'])


##################################################################################
### 999: Miscellaneous coordGrid tasks
##################################################################################

# to fftshift by a custom amout
coords.manual_shift(xshift,yshift)

# to determine locations of bright features at first glance (no projection yet)
# should really do this on projected images
coords.locate_feature(outfile = 'locate.txt')

# to get the coefficients of a photometry bootstrapping fit (IoverF_min vs mu)
# note never actually used for any papers - may not be the best way to do this
coords.bootstrap_func(order = 2)

# to write photometry without doing navigation and also writing lat/lon/mu
coords.write_photonly('h_phot.fits')



