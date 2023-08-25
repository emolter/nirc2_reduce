from . import sort_rawfiles, observation, flats, image
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import warnings
from astropy.time import Time
import importlib
import yaml
import nirc2_reduce.data.header_kw_dicts as inst_info

'''
To do: instead of nearest_stand_filt, use header's effective wavelength!
'''
standard_wleff = {'j': 1.248, 'h': 1.633, 'kp': 2.124, 'lp':3.776, 'ms':4.674,}

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

class MultiReduce:
    """
    
    
    superclass to multiBxy3 and multiNod
    """

    def __init__(self,rawdir,instrument):
        """
        Parameters
        ----------
        rawdir : str, required.
            directory containing the raw fits files
        instrument : str, required.
            name of instrument used, e.g. nirc2. 
            will look for file data/{instrument}.yaml
            in order to scrub header keywords
        
        Attributes
        ----------
        rawdir : str
        instrument : str
        header_kw_dict : dict
        tab : astropy.table.Table
        """
        self.rawdir = rawdir
        self.instrument = instrument.lower()
        with importlib.resources.open_binary(inst_info, f"{self.instrument}.yaml") as file:
            yaml_bytes = file.read()
            self.header_kw_dict = yaml.safe_load(yaml_bytes)
        self.tab = sort_rawfiles.dfits_fitsort(
            os.path.join(self.rawdir, "*.fits"), fits_kws=self.header_kw_dict['dfits_kws'])
        

    def process_flats(
        self,
        flatdir,
        **badpx_kwargs):
        """
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
        **badpx_kwargs : dict, optional.
            keyword arguments to flats.Flats.make_badpx_map()
        """
        date_kw = self.header_kw_dict['date']
        filter_kw = self.header_kw_dict['filter']
        subc_kw = self.header_kw_dict['subc']
        wl_kw = self.header_kw_dict['wl']
        standard_filts = np.array(list(standard_wleff.keys()))
        standard_wls = np.array([float(standard_wleff[key]) for key in standard_filts])
        
    
        # sort by filter
        filts, tabs_by_filt = sort_rawfiles.split_by_kw(self.tab, filter_kw)

        counter = 0
        for i, tab in enumerate(tabs_by_filt):
            
            # sort by subarray
            subcs, tabs_by_subc = sort_rawfiles.split_by_kw(tab, subc_kw)
            for j, tab in enumerate(tabs_by_subc):
            
                flatoff, flaton = sort_rawfiles.get_flats(tab, self.instrument)
                # if there are new flats in that filter, make master flat and badpx map
                # and put them into flatdir
                if (len(flatoff) > 0) and (len(flaton) > 0):
                    counter += 1
                    date = tab[date_kw].data[0]
                    wl_eff = float(tab[wl_kw].data[0])
                    standard_idx, _ = find_nearest(standard_wls, wl_eff)
                    nearest_stand_filt = standard_filts[standard_idx]
                    short_name = nearest_stand_filt
                    badpx_outpath = os.path.join(
                        flatdir, f"{date}_badpx_map_{subcs[j]}_{short_name}.fits"
                    )
                    flat_outpath = os.path.join(
                        flatdir, f"{date}_flat_master_{subcs[j]}_{short_name}.fits"
                    )
                    flatobj = flats.Flats(flatoff, flaton)
                    flatobj.write(flat_outpath)
                    flatobj.make_badpx_map(badpx_outpath, **badpx_kwargs)
                    print(f"wrote files {flat_outpath} and {badpx_outpath}")
        if counter == 0:
            warnings.warn("No flats processed! (none found by sort_rawfiles.get_flats)")

        return

    def _find_flats(self, flatdir, date, filt, subc):
        """
        find nearest-in-time flat and badpx map
        searches flatdir for nearest-in-time flats
        and badpx maps. assumes filenames in flats/ match the wildcards
        DATE-OBS*flat*{filt}*.fits and 
        DATE-OBS*badpx*{filt}*.fits
        
        Parameters
        ----------
        flatdir : str, required.
            directory to search for flats
        date : str, required.
        filt : str, required.
        subc : str, required.
        
        Returns
        -------
        str
            filename of nearest-in-time flat
        str
            filename of nearest-in-time badpx map
        """

        # find all the flats and bad pixel maps
        all_flats = glob.glob(f"{flatdir}/*flat*{subc}*{filt}*.fits")
        all_badpx = glob.glob(f"{flatdir}/*badpx*{subc}*{filt}*.fits")
        if len(all_flats) < 1:
            raise FileNotFoundError(
                f"No flats in directory {flatdir} match wildcard *flat*{subc}*{filt}*.fits"
            )
        if len(all_badpx) < 1:
            raise FileNotFoundError(
                f"No badpx maps in directory {flatdir} match wildcard *badpx*{subc}*{filt}*.fits"
            )
        flat_dates = [f.split("/")[-1].split("_")[0] for f in all_flats]
        badpx_dates = [f.split("/")[-1].split("_")[0] for f in all_badpx]

        # figure out which flat is nearest in time
        time_want = Time(date, format="isot").jd
        deltas = Time(np.array(flat_dates), format="isot").jd - time_want
        nearest_idx = np.argmin(np.abs(deltas))
        nearest = flat_dates[nearest_idx]

        # ensure that flat has an associated bad pixel map
        if nearest not in badpx_dates:
            raise ValueError(
                f"Nearest-in-time flat ({nearest}) has no associated bad pixel map (unmatched wildcard {nearest}*badpx*{filt}.fits)"
            )

        return all_flats[nearest_idx], all_badpx[nearest_idx]


class MultiBxy3(MultiReduce):
    """
    run multiple filters and objects in a single date of observing
    
    Example
    -------
    obs = MultiBxy3(rawdir)
    obs.process_flats(flatdir) # if new flats were taken
    obs.run(outdir, flatdir)
    """

    def __init__(self, rawdir, instrument):

        super().__init__(rawdir, instrument)

    def run(
        self,
        outdir,
        flatdir,
        filts_want=None,
        show=False):
        """
        Parameters
        ----------
        outdir : str, required.
            directory to save outputs to
        flatdir : str, required.
            directory in which to look for flats
        filts_want : list-like or None, optional. default None.
            list of filter names to run
            if None, runs all filters found in rawdir
        show : bool, optional. default False
            show final quicklook images while running?
        
        Writes
        ------
        
        """
        object_kw = self.header_kw_dict['object']
        date_kw = self.header_kw_dict['date']
        time_kw = self.header_kw_dict['time']
        wl_kw = self.header_kw_dict['wl']
        subc_kw = self.header_kw_dict['subc']
        filter_kw = self.header_kw_dict['filter']
        flatposkw = self.header_kw_dict['isdome']['kw']
        flatposarg = self.header_kw_dict['isdome']['yesdome']
        trackingarg = self.header_kw_dict['isdome']['nodome']
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        # read in standard filters    
        standard_filts = np.array(list(standard_wleff.keys()))
        standard_wls = np.array([float(standard_wleff[key]) for key in standard_filts])

        # select only files where scope is tracking, i.e., the science frames
        isinflatpos, flatpostabs = sort_rawfiles.split_by_kw(self.tab, flatposkw)
        tab = flatpostabs[np.argwhere(isinflatpos == trackingarg)[0,0]]
        date = tab[date_kw].data[0]

        # split by object
        targets, target_tabs = sort_rawfiles.split_by_kw(tab, object_kw)

        # loop over all targets
        for i, targ in enumerate(targets):
            print(f"Starting object {targ} ({i+1} of {len(targets)})")
            targ_tab = target_tabs[i]
            filts, filt_tabs = sort_rawfiles.split_by_kw(targ_tab, filter_kw)

            # loop over all filters
            failures = 0
            for i, filt_name in enumerate(filts):
                filt_str = filt_name.replace(" ", "")
                filt_str = filt_str.replace("+", "")
                filt_str = filt_str.replace("clear", "")
                targ_str = targ.split(" ")[0]
                if (filts_want is not None) and (filt_name not in filts_want):
                    continue

                print(f"Starting filter {filt_name} ({i+1} of {len(filts)})")
                filt_tab = filt_tabs[i]
                fnames_i = np.argsort(filt_tab["FILENAME"].data)
                fnames = filt_tab["FILENAME"].data[fnames_i]
                time_strs = filt_tab[time_kw].data[fnames_i]
                time_strs = [s[:5].replace(":", "")+'UT' for s in time_strs]

                # check right number of frames for bxy3
                if len(fnames) < 3:
                    warnings.warn(
                        f"Fewer than 3 frames found for object {targ}, filter {filt_name}; skipping filter {filt_name}!", stacklevel=2
                    )
                    failures += 1
                    continue
                if len(fnames) > 3:
                    warnings.warn(
                        f"More than 3 frames found for object {targ}, filter {filt_name}; using last three!", stacklevel=2
                    )
                    fnames = fnames[-3:]

                # find corresponding wideband flatfield and badpx map
                wl_eff = float(tab[wl_kw].data[0])
                standard_idx, _ = find_nearest(standard_wls, wl_eff)
                flat_filt = standard_filts[standard_idx]
                subc = filt_tab[subc_kw].data[0]
                try:
                    flat_fname, badpx_fname = self._find_flats(flatdir, date, flat_filt, subc)
                    print(f'Applying flat {flat_fname}, badpx {badpx_fname}')
                except (ValueError, FileNotFoundError):
                    warnings.warn(
                        f"Could not find flat for filter {flat_filt} as requested by {date} {targ} {filt_name} {subc}; skipping filter {filt_name}!", stacklevel=2
                    )
                    failures += 1
                    continue

                # run all steps of bxy3
                obs = observation.Bxy3(fnames, self.instrument)
                obs.make_sky(
                    os.path.join(outdir, f"{date}_{targ_str}_sky_{filt_str}.fits")
                )
                obs.apply_sky(
                    os.path.join(outdir, f"{date}_{targ_str}_sky_{filt_str}.fits")
                )
                obs.apply_flat(flat_fname)
                obs.apply_badpx_map(badpx_fname)
                obs.remove_cosmic_rays()
                obs.dewarp()
                obs.rotate()
                obs.per_second()
                obs.write_frames(
                    [
                        os.path.join(
                            outdir, f"{date}_{targ_str}_0_nophot_{filt_str}.fits"
                        ),
                        os.path.join(
                            outdir, f"{date}_{targ_str}_1_nophot_{filt_str}.fits"
                        ),
                        os.path.join(
                            outdir, f"{date}_{targ_str}_2_nophot_{filt_str}.fits"
                        ),
                    ],
                )
                obs.trim()
                obs.stack()
                obs.write_final(
                    os.path.join(
                        outdir, f"{date}_{targ_str}_stacked_nophot_{filt_str}.fits"
                    )
                )
                obs.com_crop_final(100)
                obs.plot_final(
                    show=show,
                    png_file=os.path.join(outdir, f"{targ_str}_{filt_str}_{time_strs[0]}.png"),
                )
            print(f"Object {targ_str} finished with {failures} failed filters")

        return


class MultiNod(MultiReduce):
    def __init__(self, rawdir):

        raise NotImplementedError()

        super().__init__(rawdir)

    def run(self):

        return
