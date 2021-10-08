# nirc2_reduce
Software for calibration, imaging, and analysis of infrared images of solar system objects.  Features include:
- Dark subtraction, flat-fielding, cosmic ray removal, dewarping, stacking
- Flux calibration and conversion to I/F
- Retrieval of ephemerides from JPL Horizons to produce model planet ellipsoids at arbitrary viewing geometries
- Automatic planet detection for centering model on image
- Projection onto longitude-latitude grid
- Flexible infrastructure for applying new filter passbands, dewarp solutions, etc. for use cases on new telescopes

Has been applied successfully to NIRC2, OSIRIS, and Lick ShARCS data. 

# usage
See DOCUMENTATION.py for sample workflows.

# caveats
This is not an official release version and some customization will be required. This software was built with Python 3.5; there may be compatibility issues with other versions of Python and packages may be deprecated. The many package dependencies are as-yet undocumented. Your development and pull requests are welcomed!

# cite
If you find this software useful for your research, please cite this page as a software citation - most journals will allow that nowadays.
