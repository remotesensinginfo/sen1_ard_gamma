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
# Purpose:  Command line tool to produce the commands for sentinel-1 processing
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 09/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import argparse
import glob
import os.path

import sen1_ard_gamma.sen1_ard_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Provide string for glob to find input SAFE "
                                                                       "files, will need quoting to include at least "
                                                                       "one required '*'")
    parser.add_argument("-o", "--output", type=str, required=True, help="Provide an output directory for output files "
                                                                        "of the processing.")
    parser.add_argument("-t", "--tmpdir", type=str, required=True, help="Provide an tmp directory for output files "
                                                                        "of the processing.")
    parser.add_argument("-d", "--dem", type=str, required=True, help="Provide GDAL readable DEM file.")
    parser.add_argument("-f", "--format", type=str, default="GTIFF", help="Provide GDAL format for final output files.")
    parser.add_argument("-r", "--resolution", type=float, required=True, help="Provide output image resolution.")
    parser.add_argument("--outfile", type=str, default="cmds_dem_lst.sh", help="Output file which contains the list of "
                                                                               "commands.")
    parser.add_argument("--cmd", type=str, default="SEN1_GRD_ARD", choices=["SEN1_GRD_ARD"],
                        help="Specify the tool for which commands are to be generated.")
    parser.add_argument("--zip", action='store_true', default=False, help="Specifies that the input SAFE file is a "
                                                                          "zip file and needs extracting.")
    args = parser.parse_args()

    if args.cmd == "SEN1_GRD_ARD":
        input_img_files = glob.glob(args.input)
        print("{} input files.".format(len(input_img_files)))
        out_file_ext = sen1_ard_gamma.sen1_ard_utils.getFileExtension(args.format)
        cmds = []
        out_cmd = True
        for img in input_img_files:
            out_cmd = True
            if args.zip:
                if os.path.isfile(img):
                    out_cmd = True
                else:
                    out_cmd = False

            cmd = "python sen1_grd_ard.py -i {0} -o {1} -t {2} -d {3} -f {4} -r {5} ".format(img, args.output,
                                                                                             args.tmpdir, args.dem,
                                                                                             args.format,
                                                                                             args.resolution)
            if args.zip:
                cmd = cmd + " --zip "
            if out_cmd:
                cmds.append(cmd)

        sen1_ard_gamma.sen1_ard_utils.write_list_to_file(cmds, args.outfile)
        print("Complete.")
    else:
        raise Exception("Command tool specified is not recognised.")
