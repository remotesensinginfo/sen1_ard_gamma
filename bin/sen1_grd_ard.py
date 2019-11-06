#!/usr/bin/env python
"""
sen1_ard_gamma - tools for Sentinel-1 GRD processing using Gamma

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
# Purpose:  Command line user interface processing Sentinel-1 GRD data products
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import sys
sys.path.append("..")

import argparse
import os
import os.path
import subprocess

import sen1_ard_gamma
import sen1_ard_gamma.sen1_grd_ard_tools
import sen1_ard_gamma.sen1_ard_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Provide the input SAFE file.")
    parser.add_argument("-o", "--output", type=str, required=True, help='''Provide the output directory for 
                                                                           final products.''')
    parser.add_argument("-t", "--tmpdir", type=str, required=True, help='''Provide the tmp directory for 
                                                                           intermediate products.''')
    parser.add_argument("-d", "--dem", type=str, required=True, help="Provide GDAL readable DEM file.")
    parser.add_argument("-r", "--resolution", type=float, required=True, help="Provide output image resolution.")
    parser.add_argument("-p", "--projepsg", type=int, default=None, help='''Provide EPSG code for output projection 
                                                                            (Default=None).''')
    parser.add_argument("--pol", type=str, nargs='+', default=None, choices=sen1_ard_gamma.SEN1_POLS,
                        help='''Specify the polarisations to be processed. If None then all processed.
                                Default is None.''')
    parser.add_argument("-f", "--format", type=str, default="GTIFF", help="Provide GDAL format for final output files.")
    parser.add_argument("--nostats", action='store_true', default=False,
                        help="Specifies that the image statistics and pyramids should be build for all output images.")
    parser.add_argument("--keepfiles", action='store_true', default=False,
                        help="Specifies intermediate files within the --tmpdir should be kept (i.e., for debugging).")
    parser.add_argument("--nodemcheck", action='store_true', default=False,
                        help="Specifies that the DEM should not be checked to ensure the minimum value is 1 with a "
                             " no data value of 0 (zero) which is used by Gamma.")
    parser.add_argument("--zip", action='store_true', default=False, help="Specifies that the input SAFE file is a "
                                                                          "zip file and needs extracting.")
    parser.add_argument("--intimgs", action='store_true', default=False,
                        help="Specifies that the output images should have a gain (x1000) applied and integerised.")

    args = parser.parse_args()

    unzip_tmp_dir_created = False
    unzip_dir = ""
    if args.zip:
        uid_val = sen1_ard_gamma.sen1_ard_utils.uidGenerator()
        base_file_name = os.path.splitext(os.path.basename(args.input))[0]
        unzip_dir = os.path.join(args.tmpdir, "{}_{}".format(base_file_name, uid_val))
        if not os.path.exists(unzip_dir):
            os.makedirs(unzip_dir)
            unzip_tmp_dir_created = True
        cwd = os.getcwd()
        os.chdir(unzip_dir)
        print(args.input)
        cmd = "unzip {}".format(args.input)
        subprocess.call(cmd, shell=True)
        input_safe = os.path.join(unzip_dir, "{}.SAFE".format(base_file_name))
        os.chdir(cwd)
    else:
        input_safe = args.input

    # Run analysis
    sen1_ard_gamma.sen1_grd_ard_tools.run_sen1_grd_ard_analysis(input_safe, args.output, args.tmpdir, args.dem,
                                                                args.resolution, args.projepsg, args.pol, args.format,
                                                                args.nostats, args.keepfiles, args.nodemcheck,
                                                                args.intimgs)
    if args.zip and unzip_tmp_dir_created:
        import shutil
        shutil.rmtree(unzip_dir)
