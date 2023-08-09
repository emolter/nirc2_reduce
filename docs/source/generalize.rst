Add a New Instrument
********************

Take the following steps to add support for a new observatory and/or instrument:

* Make a .yaml file for the instrument. In the ``nirc2_reduce/data/header_kw_dicts/`` directory, copy the nirc2.yaml file and give it a short, descriptive, all-lowercase name corresponding to the new instrument.
* Update the new .yaml file with the header keywords the code will scrub. Most of the keywords should be self-explanatory.
* If you'd like to apply distortion corrections, put new x, y distortion maps into the ``nirc2_reduce/data/`` directory; these can be called with optional keywords to ``observation.Observation.dewarp(warpx='', warpy='')``
* If needed, define new dither patterns in ``nirc2_reduce.observation``. They can most likely at least partially inherit from the base ``Observation`` class.
* If desired, define new custom colormaps for your targets of interest in ``nirc2_reduce.prettycolors``
* Test the package out on some of your data, and submit any issues you find with the codebase in the process.
* Make a pull request for your additions so that the community can benefit! See :ref:`Contributing`.
