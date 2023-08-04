#!/usr/bin/env python

from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

VERSION = '0.1.0'
DESCRIPTION = "data reduction pipelines for infrared images from research-grade telescopes"

CLASSIFIERS = list(filter(None, map(str.strip,
                                    """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Programming Language :: Python
Programming Language :: Python :: 3.6
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3.11
Topic :: Software Development :: Libraries :: Python Modules
""".splitlines())))

install_requires = [
    'numpy',
    'matplotlib',
    'scipy',
    'astropy',
    'astroquery',
    'astroscrappy',
    'pysynphot',
    'image_registration @ git+https://github.com/keflavich/image_registration.git',
    'sphinx',
    'nbsphinx',
    'pytest',
    ]

setup(
    name="nirc2_reduce",
    version=VERSION,
    description=DESCRIPTION,
    long_description=readme(),
    long_description_content_type="text/markdown",
    classifiers=CLASSIFIERS,
    author="Ned Molter",
    author_email="emolter@berkeley.edu",
    url="https://github.com/emolter/nirc2_reduce",
    python_requires='>=3.6',
    license="GPL3",
    keywords='planetary astronomy keck infrared detector',
    packages=find_packages(),
    py_modules=['nirc2_reduce'],
    platforms=['any'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=install_requires,
)
