from . import sort_rawfiles, observation, flats, image
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import warnings
from astropy.time import Time


# lookup table to see which broadband flat filter
# corresponds best to the narrow band filter used for
# a given observation
nearest_stand_filt = {'j':'j', 'he1a':'j', 'he1_a':'j', 'pagamma':'j', 'jcont':'j', 'pabeta':'j',
                      'h':'h', 'hcont':'h', 'ch4s':'h', 'ch4_short':'h', 'feii':'h', 'ch4l':'h', 'ch4_long':'h', 'h + clear':'h',
                      'k':'k', 'kp':'k', 'kcont':'k', 'ks':'k', 'he1b':'k', 'brgamma':'k', 'br_gamma':'k', 'h210':'k', 'h221':'k', 'co':'k', 'h2o':'k', 'h2o_ice':'k',
                      'lp':'l', 'lw':'l', 'bracont':'l', 'bra':'l', 'br_alpha':'l', 'bra_cont':'l', 'pah':'l',
                      'ms':'m'}


class MultiReduce:
    '''
    Description
    -----------
    superclass to multiBxy3 and multiNod
    '''
    
    def __init__(self, rawdir, dfits_kws = ['OBJECT', 'DATE-OBS', 'FILTER', 'FLSPECTR', 'TARGNAME', 'AXESTAT']):
        '''
        Parameters
        ----------
        rawdir : str, required.
            directory containing the raw fits files
        '''
        self.rawdir = rawdir
        self.tab = sort_rawfiles.dfits_fitsort(os.path.join(self.rawdir, '*.fits'), fits_kws = dfits_kws)
        
        
    def process_flats(self, flatdir, date_kw = 'DATE-OBS', filter_kw = 'FILTER', tol=0.07, blocksize=6, **kwargs):
        '''
        Description
        -----------
        check for new dome flats in rawdir, 
        make new flats for each filter present and put into flatdir
        
        Parameters
        ----------
        flatdir : str, required.
            directory to put into and/or retrieve flats
            if retrieving flats without making them, ensure
            filenames match defaults; 
            {date}_badpx_map_{filt}.fits for bad pixel maps in given filter
            {date}_flat_master_{filt}.fits for flats in given filter
        **kwargs : dict, optional.
            keyword arguments to sort_rawfiles.get_flats()
        '''
        # sort by filter
        filts, tabs_by_filt = sort_rawfiles.split_by_kw(self.tab, filter_kw)
        
        counter = 0
        for i, tab in enumerate(tabs_by_filt):
            flatoff, flaton = sort_rawfiles.get_flats(self.tab, **kwargs)
            date = self.tab[date_kw].data[0]
            
            # if there are new flats in that filter, make master flat and badpx map
            # and put them into flatdir
            if (len(flatoff) > 0) and (len(flaton) > 0):
                counter += 1
                short_name = nearest_stand_filt[filts[i].lower()]
                badpx_outpath = os.path.join(flatdir, f'{date}_badpx_map_{short_name}.fits')
                flat_outpath = os.path.join(flatdir, f'{date}_flat_master_{short_name}.fits')
                flatobj = flats.Flats(flatoff, flaton)
                flatobj.write(flat_outpath)
                flatobj.make_badpx_map(badpx_outpath, tol, blocksize)
                print(f'wrote files {flat_outpath} and {badpx_outpath}')
        if counter == 0:
            warnings.warn('No flats processed! (none found by sort_rawfiles.get_flats)')
            
        return
        
        
    def _find_flats(self,flatdir, date, filt):
        '''
        Description
        -----------
        find nearest-in-time flat and badpx map
        
        Parameters
        ----------
        flatdir : str, required.
            directory to search for flats
        date : str, required.
        filt : str, required.
        
        Returns
        -------
        str
            filename of nearest-in-time flat
        str
            filename of nearest-in-time badpx map
        '''
        
        # find all the flats and bad pixel maps
        all_flats = glob.glob(f'{flatdir}/*flat*{filt}*.fits')
        all_badpx = glob.glob(f'{flatdir}/*badpx*{filt}*.fits')
        if (len(all_flats) < 1):
            raise FileNotFoundError(f'No flats in directory {flatdir} match wildcard *flat*{filt}*.fits')
        if (len(all_badpx) < 1):
            raise FileNotFoundError(f'No badpx maps in directory {flatdir} match wildcard *badpx*{filt}*.fits')
        flat_dates = [f.split('/')[-1].split('_')[0] for f in all_flats]
        badpx_dates = [f.split('/')[-1].split('_')[0] for f in all_badpx]
        
        # figure out which flat is nearest in time
        time_want = Time(date, format='isot').jd
        deltas = Time(np.array(flat_dates), format='isot').jd - time_want
        nearest_idx = np.argmin(np.abs(deltas))
        nearest = flat_dates[nearest_idx]

        # ensure that flat has an associated bad pixel map
        if nearest not in badpx_dates:
            raise ValueError(f'Nearest-in-time flat ({nearest}) has no associated bad pixel map (unmatched wildcard {nearest}*badpx*{filt}.fits)')
        
        return all_flats[nearest_idx], all_badpx[nearest_idx]
        
        
class MultiBxy3(MultiReduce):
    '''
    Description
    -----------
    run multiple filters and objects in a single date of observing
    
    Example
    -------
    obs = MultiBxy3(rawdir)
    obs.process_flats(flatdir) # if new flats were taken
    obs.run(outdir, flatdir)
    '''
    
    def __init__(self,rawdir):
        
        super().__init__(rawdir)
        
        
    def run(self, outdir, flatdir, filts_want = None, object_kw='OBJECT', filter_kw='FILTER', date_kw = 'DATE-OBS', flatposkw='AXESTAT', flatposarg='not'):
        '''
        Parameters
        ----------
        outdir : str, required.
            directory to save outputs to
        flatdir : str, required.
            directory in which to look for flats
        filts_want : list-like or None, optional. default None.
            list of filter names to run
            if None, runs all filters found in rawdir
        
        Writes
        ------
        
        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        # ignore all files with dome in flat position
        isinflatpos, flatpostabs = sort_rawfiles.split_by_kw(self.tab, flatposkw)
        tab = flatpostabs[np.argwhere([not(s.startswith(flatposarg)) for s in isinflatpos])[0,0]]
        date = tab[date_kw].data[0]
        
        # split by object
        targets, target_tabs = sort_rawfiles.split_by_kw(tab, object_kw)
        
        # loop over all targets
        for i, targ in enumerate(targets):
            print(f'Starting object {targ} ({i+1} of {len(targets)})')
            targ_tab = target_tabs[i]
            filts, filt_tabs = sort_rawfiles.split_by_kw(targ_tab, filter_kw)
            
            # loop over all filters
            for i, filt_name in enumerate(filts):
                if (filts_want is not None) and (filt_name not in filts_want):
                    continue
                    
                print(f'Starting filter {filt_name} ({i+1} of {len(filts)})')
                failures = 0
                filt_tab = filt_tabs[i]
                fnames = np.sort(filt_tab['FILENAME'].data)

                # check right number of frames for bxy3
                if len(fnames) < 3:
                    warnings.warn(f'Fewer than 3 frames found for object {targ}, filter {filt_name}; skipping filter {filt_name}!')
                    failures += 1
                    continue
                if len(fnames) > 3:
                    warnings.warn(f'More than 3 frames found for object {targ}, filter {filt_name}; using last three!')
                    fnames = fnames[-3:]
                
                # find corresponding wideband flatfield and badpx map
                try:
                    flat_filt = nearest_stand_filt[filt_name.lower()]
                    flat_fname, badpx_fname = self._find_flats(flatdir, date, flat_filt)
                except (ValueError, FileNotFoundError):
                    warnings.warn(f'Could not find flat for filter {flat_filt} as requested by {date} {targ} {filt_name}; skipping filter {filt_name}!')
                    failures += 1
                    continue
                
                # run all steps of bxy3
                filt_str = filt_name.replace(' ', '')
                filt_str = filt_str.replace('+', '')
                targ_str = targ.split(' ')[0]
                obs = observation.Bxy3(fnames)
                obs.make_sky(os.path.join(outdir, f'{date}_{targ_str}_sky_{filt_str}.fits'))
                obs.apply_sky(os.path.join(outdir, f'{date}_{targ_str}_sky_{filt_str}.fits'))
                obs.apply_flat(flat_fname)
                obs.apply_badpx_map(badpx_fname)
                obs.remove_cosmic_rays()
                obs.dewarp()
                obs.rotate()
                obs.per_second()
                obs.write_frames(
                            [os.path.join(outdir, f'{date}_{targ_str}_0_nophot_{filt_str}.fits'),
                            os.path.join(outdir, f'{date}_{targ_str}_1_nophot_{filt_str}.fits'),
                            os.path.join(outdir, f'{date}_{targ_str}_2_nophot_{filt_str}.fits')],
                            )
                obs.trim()
                obs.stack()
                obs.crop_final(50)
                obs.plot_final(show=False, png_file=os.path.join(outdir, f'{date}_{targ_str}_{filt_str}.png'))
                obs.write_final(os.path.join(outdir, f'{date}_{targ_str}_stacked_nophot_{filt_str}.fits'))
            print(f'Object {targ_str} finished with {failures} failed filters')
                
        return
        
    
class MultiNod(MultiReduce):
    
    def __init__(self, rawdir):
        
        raise NotImplementedError()
        
        super().__init__(rawdir)
        
        
    def run(self):
        
        return
