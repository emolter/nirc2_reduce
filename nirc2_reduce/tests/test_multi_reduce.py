import glob
import os

import numpy as np
from astropy.io import fits
from pytest import fixture

from nirc2_reduce import multi_reduce


@fixture
def datadir(request):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data")
    return path


@fixture
def rawdir(request):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "raw")
    return path


@fixture
def reddir(request):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "reduced")
    return path


@fixture
def flatdir(request):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "flats")
    return path


def test_multi_bxy3(datadir, rawdir, reddir, flatdir):
    obs = multi_reduce.MultiBxy3(rawdir, "nirc2_pre_oct23")
    obs.process_flats(flatdir)  # looks for raw flat frames in rawdir, puts into flatdir
    obs.run(reddir, flatdir)

    # test major output files were written correctly
    stack_test = fits.open(os.path.join(reddir, "2017-07-25_Neptune_stacked_nophot_H.fits"))[0].data
    stack_expected = fits.open(os.path.join(datadir, "bxy3_stack_expected.fits"))[0].data
    assert np.allclose(stack_test, stack_expected, rtol=1e-3)

    flat_test = fits.open(os.path.join(flatdir, "2017-07-25_flat_master_1024_h.fits"))[0].data
    flat_expected = fits.open(os.path.join(datadir, "flat_expected.fits"))[0].data
    assert np.allclose(flat_test, flat_expected, rtol=1e-3)

    badpx_test = fits.open(os.path.join(flatdir, "2017-07-25_badpx_map_1024_h.fits"))[0].data
    badpx_expected = fits.open(os.path.join(datadir, "badpx_map_expected.fits"))[0].data
    assert np.allclose(badpx_test, badpx_expected, rtol=1e-3)

    ## cleanup: remove test files that were made
    fnames = (
        glob.glob(flatdir + "/*.fits")
        + glob.glob(reddir + "/*.fits")
        + glob.glob(reddir + "/*.png")
    )
    for f in fnames:
        os.remove(f)
