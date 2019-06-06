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
import os.path
import glob

import osgeo.gdal as gdal

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

    args = parser.parse_args()

    scn_metadata_info = sen1_ard_gamma.sen1_ard_utils.retrieve_sentinel1_metadata(args.input)

    scn_safe_files = sen1_ard_gamma.sen1_ard_utils.find_sen1_ard_files(args.input)

    scn_basename = sen1_ard_gamma.sen1_ard_utils.create_sentinel1_basename(scn_metadata_info)
    print("Basename for scene: {}".format(scn_basename))

    polarisations = args.pol
    if polarisations is None:
        polarisations = scn_metadata_info['product_polarisations']
    else:
        for pol in polarisations:
            if pol not in scn_metadata_info['product_polarisations']:
                raise Exception("Polarisation {} is not within the scene provided.")

    c_uid = 'cde44f'#sen1_ard_gamma.sen1_ard_utils.uidGenerator()
    c_tmp_dir = os.path.join(args.tmpdir, '{}_tmp_{}'.format(scn_basename, c_uid))
    c_tmp_dir_created = False
    if not os.path.exists(c_tmp_dir):
        os.mkdir(c_tmp_dir)
        c_tmp_dir_created = True

    c_out_dir = os.path.join(args.tmpdir, '{}_{}'.format(scn_basename, c_uid))
    c_out_dir_created = False
    if not os.path.exists(c_out_dir):
        os.mkdir(c_out_dir)
        c_out_dir_created = True

    demparfile = os.path.join(c_tmp_dir, scn_basename+'_gamma_scn_dem.dem_par')

    out_scns = dict()
    first = True
    for pol in polarisations:
        pol_lower = pol.lower()
        c_scn_basename = scn_basename + '_' + pol_lower
        """
        sen1_ard_gamma.sen1_grd_ard_tools.run_gamma_grd_ard_processing(scn_safe_files['measure_'+pol_lower],
                                                                       scn_safe_files['annotation_'+pol_lower],
                                                                       scn_safe_files['calibration_'+pol_lower],
                                                                       scn_safe_files['noise_'+pol_lower],
                                                                       args.dem, demparfile, c_scn_basename, c_out_dir,
                                                                       c_tmp_dir, args.resolution, -args.resolution,
                                                                       args.projepsg, use_dem_file=(not first),
                                                                       check_in_dem_filename=True,
                                                                       dem_resample_method=gdal.GRA_CubicSpline)
        """
        out_scns[pol_lower] = dict()
        out_files = glob.glob(os.path.join(c_out_dir, "{}*.tif".format(c_scn_basename)))
        for img_file in out_files:
            img_file_name = os.path.basename(img_file)
            if 'inc' in img_file_name:
                out_scns[pol_lower]['inc'] = img_file
            elif 'lsmap' in img_file_name:
                out_scns[pol_lower]['lsmap'] = img_file
            elif 'pix' in img_file_name:
                out_scns[pol_lower]['pix'] = img_file
            elif 'pwr' in img_file_name:
                out_scns[pol_lower]['pwr'] = img_file
        first = False

    print(out_scns)
    fnl_out_imgs = sen1_ard_gamma.sen1_grd_ard_tools.create_pol_stacked_products(out_scns, scn_basename, c_out_dir,
                                                                                 args.output, args.format)

    if not args.nostats:
        print("Calc image stats.")
