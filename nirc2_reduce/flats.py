import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from .image import Image


class Flats:
    """
    Description
    -----------
    class for dome flats
    """

    def __init__(self, fnames_off, fnames_on):
        """
        Description
        -----------
        subtract lights off from lights on,
        divide by median value to center it around 1
        
        Parameters
        ----------
        fnames_off : list, required. fits files for dome flat OFF
        fnames_on : list, required. fits files for dome flat ON
        
        Attributes
        ----------
        dummy_fits : Image object used to hijack header info
        frames_off : list of Image objects corresponding to fnames_off
        frames_on : list of Image objects corresponding to fnames_on
        flat : the flat-field
        """
        self.dummy_fits = Image(fnames_on[0])  # use this to hijack header info
        self.frames_off = np.asarray([Image(f).data for f in fnames_off])
        self.frames_on = np.asarray([Image(f).data for f in fnames_on])

        off = np.median(self.frames_off, axis=0)
        on = np.median(self.frames_on, axis=0)
        flat = on - off
        self.flat = flat / np.median(flat)

    def write(self, outfile):
        """
        Parameters
        ----------
        outfile : str, required. output fits filename
        """
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header["OBJECT"] = "FLAT_MASTER"
        hdulist_out[0].header["TARGNAME"] = "FLAT_MASTER"
        hdulist_out[0].data = self.flat
        hdulist_out[0].writeto(outfile, overwrite=True)

    def plot(self):
        plt.imshow(self.flat, origin="lower")
        plt.show()

    def make_badpx_map(self, outfile, tol, blocksize):
        """
        Description
        -----------
        Find pixels whose values are very far from the average of their neighbors
        Bad pixel is defined as 
        abs(pixel value / median of nearby pixels - 1) > tol
        
        Parameters
        ----------
        outfile : str, required. fits filename to write to
        tol : float, required. fractional tolerance. 
        blocksize : int, required. number of pixels 
            over which to average in each direction
        """
        badpx_map = np.ones(self.flat.shape)
        for i in range(0, self.flat.shape[0] + blocksize, blocksize):
            for j in range(0, self.flat.shape[1] + blocksize, blocksize):
                flatblock = self.flat[i : i + blocksize, j : j + blocksize]
                mapblock = badpx_map[i : i + blocksize, j : j + blocksize]
                med = np.median(flatblock)

                # if not within tolerance, set to NaN
                mapblock[np.where(flatblock / med > 1 + tol)] = 0
                mapblock[np.where(flatblock / med < 1 - tol)] = 0
                badpx_map[i : i + blocksize, j : j + blocksize] = mapblock
        self.badpx_map = badpx_map
        # change some header info and write to .fits
        hdulist_out = self.dummy_fits.hdulist
        hdulist_out[0].header["OBJECT"] = "BADPX_MAP"
        hdulist_out[0].header["TARGNAME"] = "BADPX_MAP"
        hdulist_out[0].data = self.badpx_map
        hdulist_out[0].writeto(outfile, overwrite=True)
