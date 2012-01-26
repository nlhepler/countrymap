#!/usr/bin/env python

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from countryplot import __version__ as _cplot_version

setup(name='countryplot',
      version=_cplot_version,
      description='Facilities to assist in plotting data per-country',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/veg/countryplot',
      license='GNU GPL version 3',
      packages=['countryplot'],
      package_dir={'countryplot': 'lib/countryplot'},
      package_data={'countryplot': ['data/tmwb-0.3/TM_WORLD_BORDERS-0.3.*']},
      data_files=[('/usr/local/bin', [
            'bin/plotgcnet'
      ])],
      requires=['matplotlib', 'mpl_toolkits.basemap', 'numpy']
     )
