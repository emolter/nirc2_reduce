#!/usr/bin/env python

'''Routines for sorting raw files for nirc2 data. Designed to feed into
bxy3 and driver_bxy3'''

import os, warnings
import numpy as np
from astropy.io import fits

warnings.filterwarnings('ignore','The following header keyword is invalid or follows an unrecognized non-standard convention')
warnings.filterwarnings('ignore','non-ASCII characters are present in the FITS file header and have been replaced')

def find_object(obj_name,date):
    '''Run this one first - it searches directory for given date to find
    images of a given target, e.g. domeflatoff or neptune'''
    cwd = os.getcwd()
    fnames = os.listdir(cwd+'/raw/'+date)
    fnames = np.asarray([cwd+'/raw/'+date+'/'+name for name in fnames])
    good = np.asarray([s.endswith('.fits') for s in fnames])
    fnames = fnames[good]
    obj = np.asarray([fits.getheader(f, 0, ignore_missing_end=True)['OBJECT'].split(' ')[0].lower() for f in fnames])
    return fnames[obj == obj_name.split(' ')[0].lower()]

def find_filter(fnames, filt_name):
    '''Take in filenames, check filter of all files, return file names that are for
    that filter'''
    filt = np.asarray([fits.getheader(f, 0, ignore_missing_end=True)['FWINAME'].split(' ')[0].lower() for f in fnames])
    return fnames[filt == filt_name.split(' ')[0].lower()]

def get_flats(filt_name, date):
    '''return domeflatoff, domeflaton filenames for filter'''
    flatoff = find_filter(find_object('domeflatoff',date),filt_name)
    flaton = find_filter(find_object('domeflaton',date),filt_name)
    return flatoff, flaton
    
#print(fits.getheader('raw/2017jul25/n0001.fits', 0, ignore_missing_end=True)['OBJECT'])

## example workflow
#date = '2017jul25' #must be this format
#domeflatoff, domeflaton = get_flats('h', date)
#neptune_files = find_object('neptune',date)
#neptune_files_h = find_filter(neptune_files, 'h')
#standard_files = find_object('HD1160',date)
#standard_files_h = find_filter(standard_files, 'h')
#then run steps in driver_bxy3 for these



