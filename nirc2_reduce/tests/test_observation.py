import pytest
from pytest import fixture
import importlib
from .. import observation
import numpy as np
from distutils import dir_util
import os, shutil

# until packaging is better and rootdir is found, can use
# pytest --rootdir=/Users/emolter/Python/nirc2_reduce test_observation.py

@fixture
def datadir(request,tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, 'nirc2_reduce', 'tests', 'data')
    return path
    
    
def test_datafiles_working(datadir):
    
    assert os.path.isfile(os.path.join(datadir, 'bxy3_1.fits'))


def test_bxy3():
    
    assert True
    
    
def test_nod():
    
    assert True