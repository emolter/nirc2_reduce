import warnings
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table
import nirc2_reduce.data.header_kw_dicts as inst_info
import importlib
import yaml

warnings.filterwarnings(
    "ignore",
    "The following header keyword is invalid or follows an unrecognized non-standard convention",
)
warnings.filterwarnings(
    "ignore",
    "non-ASCII characters are present in the FITS file header and have been replaced",
)


def dfits_fitsort(
    input_wildcard,
    fits_kws=[ "OBJECT", "DATE-OBS", "FILTER"]
    ):
    """
    Python implementation of the dfits | fitsort bash script workflow
    Searches headers of all fits in input_dir
    to find the requested keywords
    makes a table of filename, kw0, kw1, ...
    
    Parameters
    ----------
    input_wildcard : str, required. path to filenames
        wildcards are allowed, as long as glob.glob()
        can read it.
        example 'raw/2023jul14/\*.fits'
    fits_kws : list, optional. default ['OBJECT', 'DATE-OBS', 'FILTER']
        which header keywords to include in the output table
    
    Returns
    -------
    astropy.table.Table
        filenames and values corresponding to fits_kws
    """
    fnames = glob.glob(input_wildcard)
    fnames = np.sort(fnames)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdrs = [fits.getheader(f, 0, ignore_missing_end=True) for f in fnames]
        vals = [
            [fnames[i],] + [hdr[key] for key in fits_kws] for i, hdr in enumerate(hdrs)
        ]

    names = ["FILENAME",] + fits_kws
    dtype = [str,] * len(names)
    tab = Table(np.array(vals), names=names, dtype=dtype)

    return tab


def split_by_kw(tab, kw):
    """
    Parameters
    ----------
    tab : astropy Table, required.
        table of fits header params output by dfits_fitsort()
    kw : str, required.
        table header column to split fnames by
    
    Returns
    -------
    list
        unique values of tab[kw]
    list
        Astropy tables corresponding to those values
    """
    tab.add_index(kw)
    unq = np.unique(tab[kw])
    split_tables = [tab.loc[kw, val] for val in unq]
    return unq.data, split_tables


def get_flats(
    tab,
    instrument,
    ignore_objects=["FLAT", "BADPX", "FLAT_MASTER", "DOME_FLAT_MASTER", "BADPX_MAP"],
):
    """
    scrub table from dfits_fitsort to find domeflaton, domeflatoff filenames
    default params are for NIRC2
    
    Parameters
    ----------
    tab : astropy Table, required.
        table of fits header params output by dfits_fitsort()
        must have FILENAME kw
    instrument : str, required.
        name of instrument used, e.g. nirc2. 
        will look for file data/{instrument}.yaml
        in order to scrub header keywords
    ignore_objects : list, optional.
        special values in header[object] to ignore
    
    Returns
    -------
    list 
        filenames for dome flat off
    list
        filenames for dome flat on
    """
    instrument = instrument.lower()
    with importlib.resources.open_binary(inst_info, f"{instrument}.yaml") as file:
        yaml_bytes = file.read()
        header_kw_dict = yaml.safe_load(yaml_bytes)

    # ignore pre-existing bad pixel maps and master flats
    rm_bool = np.sum(np.array([tab[header_kw_dict['object']] == arg for arg in ignore_objects]), axis=0)
    rm_i = list(np.argwhere(rm_bool).flatten())
    tab.remove_rows(rm_i)

    # find dome position
    targnames, tabs = split_by_kw(tab, header_kw_dict['isdome']['kw'])
    dometab = tabs[np.argwhere([s.startswith(header_kw_dict['isdome']['yesdome'].strip()) for s in targnames])[0, 0]]

    # find ons and offs
    onoff, dometabs = split_by_kw(dometab, header_kw_dict['lamps']['kw'])
    ons = dometabs[np.argwhere(onoff == header_kw_dict['lamps']['lampon'])[0, 0]]
    offs = dometabs[np.argwhere(onoff == header_kw_dict['lamps']['lampoff'])[0, 0]]
    flatoff = offs["FILENAME"].data
    flaton = ons["FILENAME"].data

    return flatoff, flaton
