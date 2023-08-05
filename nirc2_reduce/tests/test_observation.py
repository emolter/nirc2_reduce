import pytest
from pytest import fixture
import importlib
from .. import observation
import numpy as np
#from distutils import dir_util
#import shutil
import os
from astropy.io import fits


@fixture
def datadir(request,tmpdir):
    rootdir = request.config.rootdir
    path = os.path.join(rootdir, 'nirc2_reduce', 'tests', 'data')
    return path
    
    
def test_datafiles_working(datadir):
    
    assert os.path.isfile(os.path.join(datadir, 'bxy3_1.fits'))


def test_bxy3(datadir):
    '''
    fnames required in data/: 
        sky_expected.fits
        flat_expected.fits
        badpx_map_expected.fits
        bxy3_stack_expected.fits
    '''
    
    # run the bxy3 workflow
    fnames = [os.path.join(datadir, 'bxy3_1.fits'), 
                os.path.join(datadir, 'bxy3_2.fits'), 
                os.path.join(datadir, 'bxy3_3.fits')]
    obs = observation.Bxy3(fnames)
    obs.make_sky(os.path.join(datadir, 'sky_test.fits'))
    obs.apply_sky(os.path.join(datadir, 'sky_expected.fits'))
    obs.apply_flat(os.path.join(datadir, 'flat_expected.fits'))
    obs.apply_badpx_map(os.path.join(datadir, 'badpx_map_expected.fits'))
    obs.dewarp()
    obs.remove_cosmic_rays()
    obs.per_second()
    obs.trim()
    obs.stack()
    obs.crop_final(50)
    obs.write_final(os.path.join(datadir, 'bxy3_stack_test.fits'))
    
    # test the sky frame
    sky_test = fits.open(os.path.join(datadir, 'sky_test.fits'))[0].data
    sky_expected = fits.open(os.path.join(datadir, 'sky_expected.fits'))[0].data
    assert np.allclose(sky_test, sky_expected, rtol=1e-3)
    
    # test the final frame
    stack_test = fits.open(os.path.join(datadir, 'bxy3_stack_test.fits'))[0].data
    stack_expected = fits.open(os.path.join(datadir, 'bxy3_stack_expected.fits'))[0].data
    assert np.allclose(stack_test, stack_expected, rtol=1e-3)
    
    ## cleanup: remove test fits
    os.remove(os.path.join(datadir, 'sky_test.fits'))
    os.remove(os.path.join(datadir, 'bxy3_stack_test.fits'))
    
    
def test_nod(datadir):
    '''
    fnames required in data/: 
        sky_expected.fits
        flat_expected.fits
        badpx_map_expected.fits
        nod_expected.fits
    '''
    obs = observation.Nod(os.path.join(datadir, 'bxy3_2.fits'), 
        os.path.join(datadir, 'sky_expected.fits'))
    obs.apply_sky()
    obs.apply_flat(os.path.join(datadir, 'flat_expected.fits'))
    obs.apply_badpx_map(os.path.join(datadir, 'badpx_map_expected.fits'))
    obs.dewarp()
    obs.remove_cosmic_rays()
    obs.per_second()
    obs.uranus_crop(100)
    obs.write_frames([os.path.join(datadir, 'nod_test.fits')])
    
    # test the final frame
    nod_test = fits.open(os.path.join(datadir, 'nod_test.fits'))[0].data
    nod_expected = fits.open(os.path.join(datadir, 'nod_expected.fits'))[0].data
    assert np.allclose(nod_test, nod_expected, rtol=1e-3)
    
    ## cleanup: remove test fits
    os.remove(os.path.join(datadir, 'nod_test.fits'))   
    