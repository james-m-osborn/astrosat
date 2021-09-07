#!/usr/bin/env python

from distutils.core import setup


setup(
    name='astrosat',
    author='James Osborn',
    author_email='james.osborn@durham.ac.uk, matthew.townson@durham.ac.uk',
    url='https://github.com/james-m-osborn/astrosat',
    packages=['astrosat',
              'astrosat.data',
              ],
    package_dir={'astrosat': 'astrosat'},
    package_data={'astrosat': ['data/bsc.dat']},
    description='A tool for forecasting satellite transits in astronomical observations',
    version='0.1',
    install_requires=[
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
)
