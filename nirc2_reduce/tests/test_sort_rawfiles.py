import pytest
from pytest import fixture
import importlib
from .. import sort_rawfiles
import numpy as np

# from distutils import dir_util
# import shutil
import os
from astropy.io import fits
from astropy import table
import warnings


@fixture
def rawdir(request, tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "raw")
    return path


def test_dfits_fitsort(rawdir):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        input_wildcard = os.path.join(rawdir, "o*.fits")
        tab = sort_rawfiles.dfits_fitsort(input_wildcard,
                ["OBJECT", "DATE-OBS", "FILTER", "AXESTAT", "FLSPECTR", "SHRNAME"])
        flatoff, flaton = sort_rawfiles.get_flats(tab, "nirc2_pre_oct23")

    offs_expected = np.array([os.path.join(rawdir, f"off{i}.fits") for i in range(5)])
    ons_expected = np.array([os.path.join(rawdir, f"on{i}.fits") for i in range(5)])

    assert type(tab) == table.Table
    assert np.all(np.sort(flatoff) == offs_expected)
    assert np.all(np.sort(flaton) == ons_expected)
