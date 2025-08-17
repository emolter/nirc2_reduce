import os
import warnings

import numpy as np
from astropy import table
from pytest import fixture

from nirc2_reduce import sort_rawfiles


@fixture
def rawdir(request):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, "nirc2_reduce", "tests", "data", "raw")
    return path


def test_dfits_fitsort(rawdir):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        input_wildcard = os.path.join(rawdir, "o*.fits")
        tab = sort_rawfiles.dfits_fitsort(
            input_wildcard, ["OBJECT", "DATE-OBS", "FILTER", "AXESTAT", "FLSPECTR", "SHRNAME"]
        )
        flatoff, flaton = sort_rawfiles.get_flats(tab, "nirc2_pre_oct23")

    offs_expected = np.array([os.path.join(rawdir, f"off{i}.fits") for i in range(5)])
    ons_expected = np.array([os.path.join(rawdir, f"on{i}.fits") for i in range(5)])

    assert isinstance(tab, table.Table)
    assert np.all(np.sort(flatoff) == offs_expected)
    assert np.all(np.sort(flaton) == ons_expected)
