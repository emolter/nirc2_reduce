import pytest
from pytest import fixture
import importlib
from .. import flats
import numpy as np
from distutils import dir_util
import os, shutil
from astropy.io import fits

# until packaging is better and rootdir is found, can use
# pytest --rootdir=/Users/emolter/Python/nirc2_reduce test_flats.py


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


def test_datafiles_working(rawdir, datadir):

    assert os.path.isfile(os.path.join(rawdir, "off0.fits"))
    assert os.path.isfile(os.path.join(datadir, "flat_expected.fits"))


def test_flats(datadir, rawdir):
    """
    fnames required in data/raw:
        f'off{i}.fits') for i in range(5)
        f'on{i}.fits') for i in range(5)
    fnames required in data/
        flat_expected.fits
        badpx_map_expected.fits
    """

    # make the flat
    domeflatoff = [os.path.join(rawdir, f"off{i}.fits") for i in range(5)]
    domeflaton = [os.path.join(rawdir, f"on{i}.fits") for i in range(5)]
    flat = flats.Flats(domeflatoff, domeflaton)
    assert (not np.any(np.isnan(flat.flat)))
    flat.write(os.path.join(datadir, "flat_test.fits"))

    # make the badpx map
    tol = 0.07  # how far from median value can a pixel be before it's bad?
    blocksize = 6  # must be divisible into 1024. how big of a box to take the median?
    flat.make_badpx_map(os.path.join(datadir, "badpx_map_test.fits"), tol, blocksize)

    # load expected flat and test
    flat_expected = fits.open(os.path.join(datadir, "flat_expected.fits"))[0].data
    flat_test = fits.open(os.path.join(datadir, "flat_test.fits"))[0].data
    assert np.allclose(flat_expected, flat_test, rtol=1e-3)

    # load expected badpx map and test
    flat_expected = fits.open(os.path.join(datadir, "badpx_map_expected.fits"))[0].data
    flat_test = fits.open(os.path.join(datadir, "badpx_map_test.fits"))[0].data
    assert np.allclose(flat_expected, flat_test, rtol=1e-3)

    # cleanup: remove test fits
    os.remove(os.path.join(datadir, "badpx_map_test.fits"))
    os.remove(os.path.join(datadir, "flat_test.fits"))
