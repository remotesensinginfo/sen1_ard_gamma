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
# Purpose:  Command line tool to edit DEM to be used in Gamma for geocoding.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 08/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import rsgislib
import rsgislib.imagecalc
import rsgislib.imageutils

import argparse
import glob
import os.path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Provide the input DEM file. (or string "
                                                                        "for glob to find input files, will need"
                                                                       "quoting to include required '*')")
    parser.add_argument("-o", "--output", type=str, required=True, help="Provide the output DEM file (or directory "
                                                                        "for output files)")
    parser.add_argument("-f", "--format", type=str, default="GTIFF", help="Provide GDAL format for final output files.")
    parser.add_argument("--genbatchcmds", action='store_true', default=False, help="Create list of commands for batch"
                                                                                   "processing of DEM tiles.")
    parser.add_argument("--outfile", type=str, default="cmds_dem_lst.sh", help="Output file which contains the list of "
                                                                               "commands to be processed")
    args = parser.parse_args()

    rsgis_utils = rsgislib.RSGISPyUtils()

    if args.genbatchcmds:
        input_img_files = glob.glob(args.input)
        out_file_ext = rsgis_utils.getFileExtension(args.format)
        cmds = []
        for img in input_img_files:
            print(img)
            img_basename = os.path.splitext(os.path.basename(img))[0]
            out_img = os.path.join(args.output, "{}{}".format(img_basename, out_file_ext))
            cmd = "prep_dem_gamma.py -i {0} -o {1} -f {2}".format(img, out_img, args.format)
            cmds.append(cmd)
        rsgis_utils.writeList2File(cmds, args.outfile)
    else:
        no_data_val = rsgis_utils.getImageNoDataValue(args.input, 1)

        exp = "(b1==0)||(b1=={})?1:b1".format(no_data_val)
        rsgislib.imagecalc.imageMath(args.input, args.output, exp, args.format, rsgislib.TYPE_16INT)
        rsgislib.imageutils.popImageStats(args.output, usenodataval=True, nodataval=0, calcpyramids=True)
