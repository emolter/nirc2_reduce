import pytest
from pytest import fixture
import importlib
from .. import sort_rawfiles
import numpy as np
#from distutils import dir_util
#import shutil
import os
from astropy.io import fits
from astropy import table
import warnings


@fixture
def datadir(request,tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, 'nirc2_reduce', 'tests', 'data')
    return path


def test_dfits_fitsort(datadir):
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        input_wildcard = os.path.join(datadir, 'o*.fits')
        tab = sort_rawfiles.dfits_fitsort(
                input_wildcard, 
                )
        flatoff, flaton = sort_rawfiles.get_flats(tab)

    offs_expected = np.array([os.path.join(datadir, f'off{i}.fits') for i in range(5)])
    ons_expected = np.array([os.path.join(datadir, f'on{i}.fits') for i in range(5)])
    
    assert type(tab) == table.Table
    assert np.all(np.sort(flatoff) == offs_expected)
    assert np.all(np.sort(flaton) == ons_expected)