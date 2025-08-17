import warnings

import matplotlib.pyplot as plt
from astropy.io import fits


class Image:
    """
    Thin wrapper to astropy.io.fits.

    To do
    -----
    make agnostic to specific .fits header keywords
    expand plot, write functionality
    """

    def __init__(self, fname):
        """
        Provide convenience function for opening fits files.

        Parameters
        ----------
        fname : str, required. input .fits filename

        Attributes
        ----------
        hdulist: .fits format hdulist structure
        header: image header, dictionary-like
        data: the image data
        target: name of planet
        filt: filter in which image was taken
        """
        warnings.filterwarnings(
            "ignore",
            "The following header keyword is invalid or follows an unrecognized "
            "non-standard convention",
        )
        warnings.filterwarnings(
            "ignore",
            "non-ASCII characters are present in the FITS file header and have been replaced",
        )
        warnings.filterwarnings(
            "ignore",
            "Header block contains null bytes instead of spaces for padding, "
            "and is not FITS-compliant",
        )
        self.hdulist = fits.open(fname, ignore_missing_end=True)
        self.header = self.hdulist[0].header
        self.data = self.hdulist[0].data
        try:
            targ = self.header["OBJECT"]
            self.target = targ.split()[0].strip(", \n").capitalize()
        except KeyError:
            self.target = "Unknown"
        try:
            fwi = self.header["FWINAME"].strip(", \n")
            fwo = self.header["FWONAME"].strip(", \n")
            if fwo == "clear" or fwo == "PK50_1.5":
                self.filt = fwi
            elif fwi == "clear" or fwi == "PK50_1.5":
                self.filt = fwo
            else:
                self.filt = "Unknown"
        except KeyError:
            self.filt = "Unknown"

    def plot(self):
        """Make a simple plot of the data."""
        plt.imshow(self.data, origin="lower left")
        plt.show()

    def write(self, fname):
        """
        Write the image to fits.

        Parameters
        ----------
        fname : str
            Filename to write to.
        """
        self.hdulist[0].header = self.header
        self.hdulist[0].data = self.data
        self.hdulist.writeto(fname, overwrite=True)
