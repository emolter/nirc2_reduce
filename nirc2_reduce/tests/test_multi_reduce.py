import pytest
from pytest import fixture
import importlib
from .. import multi_reduce
import numpy as np
import os
import glob
from astropy.io import fits
import warnings


@fixture
def datadir(request, tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data")
    return path


@fixture
def rawdir(request, tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "raw")
    return path


@fixture
def reddir(request, tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "reduced")
    return path


@fixture
def flatdir(request, tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "flats")
    return path


def test_MultiBxy3(datadir, rawdir, reddir, flatdir):

    obs = multi_reduce.MultiBxy3(rawdir, "nirc2")
    obs.process_flats(flatdir)  # looks for raw flat frames in rawdir, puts into flatdir
    obs.run(reddir, flatdir)

    # test major output files were written correctly
    stack_test = fits.open(
        os.path.join(reddir, "2017-07-25_Neptune_stacked_nophot_H.fits")
    )[0].data
    stack_expected = fits.open(os.path.join(datadir, "bxy3_stack_expected.fits"))[
        0
    ].data
    assert np.allclose(stack_test, stack_expected, rtol=1e-3)

    flat_test = fits.open(os.path.join(flatdir, "2017-07-25_flat_master_1024_h.fits"))[
        0
    ].data
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
