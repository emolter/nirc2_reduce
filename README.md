[![Documentation Status](https://readthedocs.org/projects/nirc2-reduce/badge/?version=latest)](https://nirc2-reduce.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6584662.svg)](https://doi.org/10.5281/zenodo.6584662)


# Infrared Imaging Data Pipeline
The nirc2\_reduce contains software for calibration, imaging, and analysis of infrared images of solar system objects.  Features include:
- Dark subtraction, flat-fielding, cosmic ray removal
- Distortion correction and application of astrometry solutions
- Easy reduction of the famous NIRC2 bxy3 dither pattern (among other supported dithers)
- Flux calibration and photometry
- Flexible infrastructure for applying new dither patterns, filter passbands, distortion and  solutions, etc. for use cases on new telescopes

This package has been applied successfully to NIRC2, OSIRIS, and Lick ShARCS data.

Example images of Solar System planets made using this pipeline are available at the [Keck Twilight Zone website](https://www2.keck.hawaii.edu/inst/tda/TwilightZone.html#). 

# Usage
The readthedocs page is in progress. For now, see DOCUMENTATION.py for sample workflows.

# Caveats
This is not yet an official release version, and some customization will be required. This software was originally written in 2016 when I was a first-year grad student and had no idea how to package code. At the time of writing, only the flats and observation modules have been rewritten and nicely packaged; I'm working on the rest. Your development and pull requests are welcomed!

Please note that some functionality that used to be in nirc2\_reduce is now located in the [pylanetary](https://github.com/emolter/pylanetary) package. If you got here following the Zenodo link from one of my publications prior to 2023, please look through the GitHub history for an older version.

# Cite
If you find this software useful for your research, please cite it using the Zenodo DOI above.

# Publications
This pipeline has been used in the following refereed publications:
- [Chavez et al. 2023, Evolution of Neptune at near-infrared wavelengths from 1994 through 2022, Icarus](https://doi.org/10.1016/j.icarus.2023.115667)
- [Chavez et al. 2023, Drift rates of major Neptunian features between 2018 and 2021, Icarus](https://doi.org/10.1016/j.icarus.2023.115604)
- [de Kleer et al. 2019, Io's Volcanic Activity from Time Domain Adaptive Optics Observations: 2013–2018, Astronomical Journal](https://doi.org/10.3847/1538-3881/ab2380)
- [Molter et al. 2019, Analysis of Neptune's 2017 bright equatorial storm, Icarus](https://doi.org/10.1016/j.icarus.2018.11.018)
