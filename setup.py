#!/usr/bin/env python

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from countrymap import __version__ as _cmap_version

setup(name='countrymap',
      version=_cmap_version,
      description='Facilities to assist in plotting data per-country',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/veg/countrymap',
      license='GNU GPL version 3',
      packages=['countrymap'],
      package_dir={'countrymap': 'lib/countrymap'},
      package_data={'countrymap': ['data/tmwb-0.3/TM_WORLD_BORDERS-0.3.*']},
      data_files=[('/usr/local/bin', [
            'bin/plotgcnet'
      ])],
      requires=['matplotlib', 'mpl_toolkits.basemap', 'numpy', 'pydot']
     )
