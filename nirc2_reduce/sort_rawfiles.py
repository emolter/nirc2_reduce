import warnings
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table

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
    fits_kws=[
        "OBJECT",
        "DATE-OBS",
        "FILTER",
        "FLIMAGIN",
        "FLSPECTR",
        "CURRINST",
        "TARGNAME",
        "AXESTAT",
    ],
):
    """
    Description
    -----------
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
    fits_kws : list, optional. default ['OBJECT', 'FILTER']
        which keywords to include in the output table
    
    Returns
    -------
    astropy.table.Table of filenames and values corresponding to fits_kws
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
    list of kw vals, list of Astropy tables
    """
    tab.add_index(kw)
    unq = np.unique(tab[kw])
    split_tables = [tab.loc[kw, val] for val in unq]
    return unq.data, split_tables


def get_flats(
    tab,
    isdome_kw="TARGNAME",
    isdome_arg="DOME",
    lamps_kw="FLSPECTR",
    lampson_arg="on",
    lampsoff_arg="off",
    remove_kw="OBJECT",
    remove_args=["FLAT_MASTER", "DOME_FLAT_MASTER", "BADPX_MAP"],
):
    """
    Description
    -----------
    scrub table from dfits_fitsort to find domeflaton, domeflatoff filenames
    default params are for NIRC2
    
    Parameters
    ----------
    tab : astropy Table, required.
        table of fits header params output by dfits_fitsort()
        must have FILENAME kw
    isdome_kw : str, optional.
        the values in column isdome_kw that start with isdome_arg
        should select images in dome flat position
    isdome_arg : str, optional.
    lamps_kw : str, optional.
        the values in column lamps_kw that equal lampson_arg
        should select images with domeflat lights on, and
        the values in column lamps_kw that equal lampsoff_arg
        should select images with domeflat lights off
    lampson_arg : str, optional.
    lampsoff_arg : str, optional.
    remove_kw : str_optional.
        the values in column remove_kw that equal any of remove_args
        will be removed from both returned lists
    remove_args : str, optional.
    
    Returns
    -------
    list 
        filenames for dome flat off
    list
        filenames for dome flat on
    """
    # ignore pre-existing bad pixel maps and master flats
    rm_bool = np.sum(np.array([tab[remove_kw] == arg for arg in remove_args]), axis=0)
    rm_i = list(np.argwhere(rm_bool).flatten())
    tab.remove_rows(rm_i)

    # find dome position
    targnames, tabs = split_by_kw(tab, isdome_kw)
    dometab = tabs[np.argwhere([s.startswith(isdome_arg) for s in targnames])[0, 0]]

    # find ons and offs
    onoff, dometabs = split_by_kw(dometab, lamps_kw)
    ons = dometabs[np.argwhere(onoff == lampson_arg)[0, 0]]
    offs = dometabs[np.argwhere(onoff == lampsoff_arg)[0, 0]]
    flatoff = offs["FILENAME"].data
    flaton = ons["FILENAME"].data

    return flatoff, flaton
