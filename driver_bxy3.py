#!/usr/bin/env python
#import bxy3, flats, phot
from . import bxy3, flats, phot

'''Lists methods to call for reduction of nirc2 data'''

# first handle flats       
flat = flats.Flats(['test_flats/off0.fits','test_flats/off1.fits','test_flats/off2.fits','test_flats/off3.fits','test_flats/off4.fits'],
            ['test_flats/on0.fits','test_flats/on1.fits','test_flats/on2.fits','test_flats/on3.fits','test_flats/on4.fits'])

flat.make('out/flat_master.fits')
#flat.plot()
tol = 0.07 #how far from median value can a pixel be before it's bad?
blocksize = 4 #must be divisible into 1024. how big of a box to take the median?
flat.make_badpx_map('out/badpx_map.fits',tol,blocksize)

#handle bxy3 observation
obs = bxy3.Bxy3(['test_data/bxy3_1.fits','test_data/bxy3_2.fits','test_data/bxy3_3.fits'])
obs.make_sky('out/sky.fits')
obs.apply_sky('out/sky.fits')
obs.apply_flat('out/flat_master.fits')
obs.apply_badpx_map('out/badpx_map.fits')
obs.dewarp()
obs.trim()
obs.remove_cosmic_rays()
# at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; need to check if those exist in flats or not
obs.per_second()

#do photometry
filt = 'h'        
stand = phot.Phot(filt,['test_data/standard_1.fits','test_data/standard_2.fits','test_data/standard_3.fits'])
stand.bxy3_reduce('out/sky.fits','out/flat_master.fits','out/badpx_map.fits',['out/phot1.fits','out/phot2.fits','out/phot3.fits'])
stand.get_mag('test_data/starlist_standard.txt')
#print(stand.star_mag)
stand.find_cts(50, plot = False)
#print(obs.star_cts)
stand.find_flux_conversion()
#print(obs.flux_per) #Imke consistently gets between 7e-17 and 9e-17 in H

obs.apply_photometry(stand.median_flux_per)

#write out final images
obs.write_frames(['out/frame1.fits','out/frame2.fits','out/frame3.fits'])
obs.stack()
obs.crop(50)
obs.write_final('out/final.fits',png=True, png_file = 'out/final.png')