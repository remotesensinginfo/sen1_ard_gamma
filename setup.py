#!/usr/bin/env python
"""
Setup script for sen1_ard_gamma. Use like this for Unix:

$ python setup.py install

"""
# Copyright 2019 Pete Bunting
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Purpose:  Install the sen1_ard_gamma software
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.


from distutils.core import setup
import os

setup(name='sen1_ard_gamma',
    version='0.2.0',
    description='Software for processing ESA Sentinel-1 data to generate ARD products',
    author='Pete Bunting',
    author_email='pfb@aber.ac.uk',
    scripts=['bin/sen1_grd_ard.py'],
    packages=['sen1_ard_gamma'],
    package_dir={'sen1_ard_gamma': 'sen1_ard_gamma'},
    data_files=[(os.path.join('share','sen1_ard_gamma'),
                [os.path.join('share','sen1_ard_gamma', 'loggingconfig.json')])],
    license='LICENSE.txt',
    url='https://bitbucket.org/petebunting/sen1_ard_gamma',
    classifiers=['Intended Audience :: Developers',
    	  'Intended Audience :: Remote Sensing Scientists',
    	  'Intended Audience :: Atmospheric Scientists',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.2',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6'])
