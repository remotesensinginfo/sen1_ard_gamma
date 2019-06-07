#!/usr/bin/env python
"""
sen1_ard_gamma - this file is needed to ensure it can be imported

See other source files for details
"""
# This file is part of 'sen1_ard_gamma'
# A set of tools to produce Sentinel-1 ARD using Gamma.
#
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
#
# Purpose:  Setup variables and imports across the whole module
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

from distutils.version import LooseVersion
import os
import logging
import logging.config
import json

SEN1_ARD_GAMMA_VERSION_MAJOR = 0
SEN1_ARD_GAMMA_VERSION_MINOR = 1
SEN1_ARD_GAMMA_VERSION_PATCH = 1

SEN1_ARD_GAMMA_VERSION = str(SEN1_ARD_GAMMA_VERSION_MAJOR) + "."  + str(SEN1_ARD_GAMMA_VERSION_MINOR) + "." + str(SEN1_ARD_GAMMA_VERSION_PATCH)
SEN1_ARD_GAMMA_VERSION_OBJ = LooseVersion(SEN1_ARD_GAMMA_VERSION)

SEN1_ARD_GAMMA_COPYRIGHT_YEAR = "2019"
SEN1_ARD_GAMMA_COPYRIGHT_NAMES = "Pete Bunting"

SEN1_ARD_GAMMA_SUPPORT_EMAIL = "pfb@aber.ac.uk"

SEN1_POLS = ['VV', 'VH']

log_default_level=logging.INFO

# Get install prefix
install_prefix = __file__[:__file__.find('lib')]

log_config_path = os.path.join(install_prefix, "share","sen1_ard_gamma", "loggingconfig.json")

log_config_value = os.getenv('S1ARD_LOG_CFG', None)
if log_config_value:
    log_config_path = log_config_value
if os.path.exists(log_config_path):
    with open(log_config_path, 'rt') as f:
        config = json.load(f)
    logging.config.dictConfig(config)
else:
    print('Warning: did not find a logging configuration file.')
    logging.basicConfig(level=log_default_level)
