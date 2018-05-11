#!/usr/bin/env python

from . import sort_rawfiles, bxy3, phot, flats, image, coordgrid, nod
import matplotlib.pyplot as plt
import numpy as np

def multiBxy3(date, target_name, filt_list, flatfilt_list):
    '''Do bxy3 for many filters without human intervention'''
    outdir = 'reduced/'+date+'/'
    filenames_all = sort_rawfiles.find_object(target_name,date)
    
    nproblems = 0
    for i in range(len(filt_list)):
        filt_name = filt_list[i]
        flat_filt = flatfilt_list[i]
        fnames = sort_rawfiles.find_filter(filenames_all, filt_name)
        print('Starting filter %s (%d/%d)'%(filt_name, i+1, len(filt_list)))
        
        #check that it is a bxy3
        if len(fnames) != 3:
            nproblems += 1
            print('ERROR: found %d files for filt %s (3 required)'%(len(fnames),filt_name))
            print('Skipping this filter')
        else:
            try:
                #do full data reduction for that filter
                obs = bxy3.Bxy3(fnames)
                obs.make_sky(outdir+'sky_'+filt_name+'.fits')
                obs.apply_sky(outdir+'sky_'+filt_name+'.fits')
                if flat_filt != None:
                    obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
                    obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
                obs.dewarp()
                obs.trim()
                obs.remove_cosmic_rays() # at this step there remain a few spots where pixel value is much lower than Neptune pixels. Doubt this is the fault of cosmic ray program; possibly failing to find them in flats with bad pixel search. Perhaps multi-layer search, e.g. large blocksize first, small blocksize second
                obs.per_second()
                obs.write_frames([outdir+'frame0_nophot_'+filt_name+'.fits',outdir+'frame1_nophot_'+filt_name+'.fits',outdir+'frame2_nophot_'+filt_name+'.fits'])
                obs.stack()
                obs.crop(50)
                obs.write_final(outdir+'stacked_nophot_'+filt_name+'.fits')
                print('Finished filter %s'%filt_name)
            except:
                nproblems += 1
                print('ERROR: something went wrong internally with filter %s'%filt_name)
                print('Ensure flat and bad pixel map files exist in outdir')
                print('Else recommend trying manual reduction to find bug')
            
    print('Finished! Errors detected in %d filters... skipped these'%nproblems)

    
def multiBxy3Plot(date, filt_list):
    '''Plot results of multiBxy3'''
    outdir = 'reduced/'+date+'/'
    gridsz = int(np.sqrt(len(filt_list))) + 1
    
    fig, axes = plt.subplots(gridsz, gridsz, figsize = (9,9))
    fig.subplots_adjust(hspace = 0.2, wspace = 0.01)
    axes = axes.flatten()
    for j in range(len(filt_list)):
        ax = axes[j]
        filt_name = filt_list[j]
        fname = outdir+'stacked_nophot_'+filt_name+'.fits'
        img = image.Image(fname).data
        ax.imshow(img, origin = 'lower left')
        ax.set_title(filt_name)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        
    plt.show()


def multiNod(date, target_name, filt_list, flatfilt_list):
    '''Do nod for many filters without human intervention'''
    outdir = 'reduced/'+date+'/'
    imagefiles_all = sort_rawfiles.find_object(target_name,date)
    skyfiles_all = sort_rawfiles.find_object('sky',date)
    
    nproblems = 0
    for i in range(len(filt_list)):
        filt_name = filt_list[i]
        flat_filt = flatfilt_list[i]
        print('Starting filter %s (%d/%d)'%(filt_name, i+1, len(filt_list)))
        try:
            imagef = sort_rawfiles.find_filter(imagefiles_all, filt_name)
        except:
            nproblems += 1
            print('ERROR: Did not find any images in filter %s'%filt_name)
            print('Skipping this filter')
            continue
        try:
            skyf = sort_rawfiles.find_filter(skyfiles_all, filt_name)
        except:
            nproblems += 1
            print('ERROR: Did not find any sky frames in filter %s'%filt_name)
            print('Skipping this filter')
            continue
        #check that it is a bxy3
        if len(imagef) != 1:
            nproblems += 1
            print('ERROR: found %d image files for filt %s (1 required)'%(len(imagef),filt_name))
            print('Skipping this filter')
        elif len(skyf) != 1:
            nproblems += 1
            print('ERROR: found %d sky files for filt %s (1 required)'%(len(skyf),filt_name))
            print('Skipping this filter')            
        else:
            try:
                #do full data reduction for that filter
                obs = nod.Nod(skyf[0], imagef[0])
                obs.apply_sky()
                obs.apply_flat(outdir+'flat_master_'+flat_filt+'.fits')
                obs.apply_badpx_map(outdir+'badpx_map_'+flat_filt+'.fits')
                obs.dewarp()
                obs.remove_cosmic_rays()
                obs.per_second()
                obs.write(outdir+'stacked_nophot_'+filt_name+'.fits', png = False)
                print('Finished filter %s'%filt_name)
            except:
                nproblems += 1
                print('ERROR: something went wrong internally with filter %s'%filt_name)
                print('Ensure flat and bad pixel map files exist in outdir')
                print('Else recommend trying manual reduction to find bug')
            
    print('Finished! Errors detected in %d filters... skipped these'%nproblems)    
    
def multiNodPlot(date, filt_list):
    '''Plot results of multiNod'''
    outdir = 'reduced/'+date+'/'
    gridsz = int(np.sqrt(len(filt_list))) + 1
    
    fig, axes = plt.subplots(gridsz, gridsz, figsize = (9,9))
    fig.subplots_adjust(hspace = 0.2, wspace = 0.01)
    axes = axes.flatten()
    for j in range(len(filt_list)):
        ax = axes[j]
        filt_name = filt_list[j]
        fname = outdir+'stacked_nophot_'+filt_name+'.fits'
        img = image.Image(fname).data
        ax.imshow(img, origin = 'lower left')
        ax.set_title(filt_name)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        
    plt.show()    
    
    
    
    
        
        
#####################
# handle photometry #
#####################
        
def make_table(outfile, cts_dict):
    with open(outfile,'w') as f:
        f.write('### Conversion between counts per second and flux density: units erg s-1 cm-2 um-1 / cts s-1\n')
        f.write('###FILTER    MEDFLUX    FLUX0    FLUX1    FLUX2\n')
        for key, val in cts_dict.items():
            flux0 = str(cts_dict[key]['flux0'])
            flux1 = str(cts_dict[key]['flux1'])
            flux2 = str(cts_dict[key]['flux2'])
            meanflux = str(cts_dict[key]['meanflux'])
            outstr = '    '.join([key, meanflux, flux0, flux1, flux2])
            f.write(outstr+'\n')
    
def multiPhot(date, target_name, filt_list, mag_list, flatfilt_list, doplots = True, dist_list = None):
    '''Do bxy3 for many filters without human intervention'''
    outdir = 'reduced/'+date+'/'
    filenames_all = sort_rawfiles.find_object(target_name,date)
    
    cts_dict = {}
    nproblems = 0
    for i in range(len(filt_list)):
        star_mag = mag_list[i]
        filt_name = filt_list[i]
        cts_dict[filt_name] = {}
        flat_filt = flatfilt_list[i]
        fnames = sort_rawfiles.find_filter(filenames_all, filt_name)
        print('Starting filter %s (%d/%d)'%(filt_name, i+1, len(filt_list)))
        
        #check that it is a bxy3
        if len(fnames) != 3:
            nproblems += 1
            print('ERROR: found %d files for filt %s (3 required)'%(len(fnames),filt_name))
            print('Skipping this filter')
        else:
            try:
                #do full data reduction for that filter
                obs = phot.bxy3Phot(fnames)
                obs.reduce(outdir+'photsky_'+filt_name+'.fits', outdir+'flat_master_'+flat_filt+'.fits', outdir+'badpx_map_'+flat_filt+'.fits', [outdir+'phot0_'+filt_name+'.fits', outdir+'phot1_'+filt_name+'.fits', outdir+'phot2_'+filt_name+'.fits'])
                
                #do photometry
                dist = 50
                if dist_list != None:
                    dist = dist_list[i]
                obs.find_cts(dist, plot = doplots)
                obs.find_flux_conversion(star_mag, filt_name)
                cts_dict[filt_name]['flux0'] = obs.all_flux_per[0]
                cts_dict[filt_name]['flux1'] = obs.all_flux_per[1]
                cts_dict[filt_name]['flux2'] = obs.all_flux_per[2]
                cts_dict[filt_name]['meanflux'] = obs.flux_per
                print('Finished filter %s'%filt_name)
            except:
                nproblems += 1
                print('ERROR: something went wrong internally with filter %s'%filt_name)
                print('Ensure flat and bad pixel map files exist in outdir')
                print('Else recommend trying manual reduction to find bug')
            
    make_table(outdir+'phot_table.txt', cts_dict)
    print('Made file reduced/%s/phot_table.txt'%date)
    print('Finished! Errors detected in %d filters... skipped these'%nproblems)

    
def multiApplyPhot(date, stand_airmass):
    '''Read phot_table and apply to the correct images of Io.
    outputs a file called FILT_calibrated.fits''' 
    outdir = 'reduced/'+date+'/'
    phot_data = np.genfromtxt(outdir+'phot_table.txt', comments = '#', dtype = [('mystring', '|S10'), ('myfloat', '<f8'), ('myfloat2', '<f8'), ('myfloat3', '<f8'), ('myfloat4', '<f8')])
    filt_list = [val[0].decode("utf-8") for val in phot_data]
    flux_per_list = [val[1] for val in phot_data]
    
    for i in range(len(filt_list)):
        filt = filt_list[i]
        print('Filter %s'%filt)
        flux_per = flux_per_list[i]
        coords = coordgrid.CoordGrid(outdir + 'stacked_nophot_%s.fits'%filt)
        coords.ioverf(filt, flux_per, stand_airmass)
        coords.write_photonly(outdir + filt + '_calibrated.fits')
        
     
     
     
     
     
     
        