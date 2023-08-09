Description
***********

nirc2_reduce is a pure Python package for infrared imaging data reduction. Originally designed for the `Keck Twilight Zone <https://www2.keck.hawaii.edu/inst/tda/TwilightZone.html#>`_ program on Keck's NIRC2 instrument, it has been packaged and made generic enough to process any type of imaging data.

Features include

* A pure Python implementation of the dfits | fitsort shell scripts for file management
* Standard flatfielding and bad pixel removal
* Cosmic ray removal using `astroscrappy <https://astroscrappy.readthedocs.io/en/latest/api/astroscrappy.detect_cosmics.html#astroscrappy.detect_cosmics>`_
* Distortion and astrometry corrections, e.g., the ones for NIRC2 `found here <https://github.com/jluastro/nirc2_distortion/wiki>`_
* Sub-pixel cross-correlation of images for stacking using DFT upsampling with `image_registration <https://image-registration.readthedocs.io/en/latest/>`_
* (limited) support for automatic processing of multiple targets and filters at once

At the time of writing, the package has only been tested on NIRC2 and the OSIRIS imager on Keck 1; see :ref:`Add a New Instrument` for instructions on adding a new instrument.