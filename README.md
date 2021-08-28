# nirc2_reduce
Software for calibration, imaging, and analysis of infrared images of solar system objects.  Features include:
- Dark subtraction, flat-fielding, cosmic ray removal, dewarping, stacking
- Flux calibration and conversion to I/F
- Retrieval of ephemerides from JPL Horizons to produce model planet ellipsoids at arbitrary viewing geometries
- Automatic planet detection for centering model on image
- Projection onto longitude-latitude grid
- Flexible infrastructure for applying new filter passbands, dewarp solutions, etc. for use cases on new telescopes

Has been applied successfully to NIRC2, OSIRIS, and Lick ShARCS data. 

See DOCUMENTATION.py for sample workflows.

