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
# Purpose:  Pre-processing DEMs for use in Gamma.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import os.path
import subprocess

import osgeo.osr as osr
import osgeo.gdal as gdal

import math
import numpy

import sen1_ard_gamma.sen1_ard_utils

import logging

logger = logging.getLogger(__name__)

gdal.UseExceptions()

def subset_reproj_utm_dem_file(in_dem_file, bbox_wgs84, out_res_x, out_res_y, out_proj_epsg, out_dem_file,
                               out_dem_par_file, tmp_dir, no_dem_check=False, eResampleAlg=gdal.GRA_CubicSpline):
    """
    Function to subset and reproject a DEM provided in a GDAL compatiable format (e.g., GTIFF) into a format
    compatiable with Gamma.

    :param in_dem_file: file path to input image file.
    :param bbox_wgs84: the bbox (MinX, MaxX, MinY, MaxY) in WGS84 (lat/long)
    :param out_res_x: the output resolution (in metres) in the x axis.
    :param out_res_y: the output resolution (in metres) in the y axis.
    :param out_proj_epsg: EPSG code for the output UTM zone
    :param out_dem_file: output DEM binary file
    :param out_dem_par_file: output DEM parameter file
    :param tmp_dir: directory for temporary outputs.
    :param no_dem_check: If True the DEM is not checked for no data value and minimum value for Gamma compatibility
    :param eResampleAlg: specify the gdal resampling method for the DEM processing (if DEM being subsetted and
                         reprojected). Default: gdal.GRA_CubicSpline
    """
    logger.debug("Starting to run subset_reproj_utm_dem_file")
    in_dem_file_ds = gdal.Open(in_dem_file, gdal.GA_ReadOnly)
    if in_dem_file_ds == None:
        raise Exception('Could not open raster image: \'' + in_dem_file + '\'')
    in_dem_band = in_dem_file_ds.GetRasterBand(1)
    if in_dem_band == None:
        raise Exception('Could not open DEM raster band: \'' + in_dem_file + '\'')

    no_data_val = in_dem_band.GetNoDataValue()

    logger.debug("Calculating UTM BBox.")
    in_spat_ref = osr.SpatialReference()
    in_spat_ref.ImportFromEPSG(4326)

    out_spat_ref = osr.SpatialReference()
    out_spat_ref.ImportFromEPSG(out_proj_epsg)

    utm_zone = out_spat_ref.GetUTMZone()

    tlX = bbox_wgs84[0]
    tlY = bbox_wgs84[3]
    brX = bbox_wgs84[1]
    brY = bbox_wgs84[2]

    out_tlX = 0.0
    out_tlY = 0.0
    out_brX = 0.0
    out_brY = 0.0

    out_tlX, out_tlY = sen1_ard_gamma.sen1_ard_utils.reproj_point(in_spat_ref, out_spat_ref, tlX, tlY)
    out_brX, out_brY = sen1_ard_gamma.sen1_ard_utils.reproj_point(in_spat_ref, out_spat_ref, brX, brY)
    logger.debug("Calculated UTM BBox: [{}, {}, {}, {}]".format(out_tlX, out_brX, out_brY, out_tlY))

    # Snap TL to whole pixel so outputs on whole pixel grid.
    n_pxl_tl_x = math.floor(out_tlX / out_res_x)
    n_pxl_tl_y = math.ceil(out_tlY / ((-1) * out_res_y))
    logger.debug("N Pixels to TL: [{}, {}]".format(n_pxl_tl_x, n_pxl_tl_y))

    out_tlX = (n_pxl_tl_x * out_res_x) - (out_res_x / 2)
    out_tlY = (n_pxl_tl_y * ((-1) * out_res_y)) + (((-1) * out_res_y) / 2)
    logger.debug("Calculated UTM BBox (whole num grid): [{}, {}, {}, {}]".format(out_tlX, out_brX, out_brY, out_tlY))

    x_pxls = math.ceil((out_brX - out_tlX) / out_res_x)
    y_pxls = math.ceil((out_tlY - out_brY) / ((-1) * out_res_y))

    logger.debug("Output DEM image size (pixels): [{} x {}]".format(x_pxls, y_pxls))

    out_dem_tmp_file = os.path.join(tmp_dir ,os.path.splitext(os.path.basename(out_dem_file))[0] + '_tmpdem.dem')
    logger.debug("Output DEM tmp file: {}".format(out_dem_tmp_file))
    out_dem_file_ds = gdal.GetDriverByName('ENVI').Create(out_dem_tmp_file, x_pxls, y_pxls, 1, gdal.GDT_Int16)
    if out_dem_file_ds == None:
        raise Exception('Could not create output raster image: \'' + out_dem_tmp_file + '\'')
    out_dem_file_ds.SetGeoTransform((out_tlX, out_res_x, 0, out_tlY, 0, out_res_y))
    out_dem_file_ds.SetProjection(out_spat_ref.ExportToWkt())
    out_dem_band = out_dem_file_ds.GetRasterBand(1)
    if out_dem_band == None:
        raise Exception('Could not open DEM raster band: \'' + out_dem_tmp_file + '\'')
    out_dem_band.SetNoDataValue(no_data_val)
    logger.debug("Created output DEM tmp file: {}".format(out_dem_tmp_file))

    logger.debug("Start gdal warp for: {}".format(out_dem_tmp_file))
    wrpOpts = gdal.WarpOptions(resampleAlg=eResampleAlg, srcNodata=no_data_val, dstNodata=no_data_val,
                               multithread=False, callback=gdal.TermProgress)

    gdal.Warp(out_dem_file_ds, in_dem_file_ds, options=wrpOpts)
    logger.debug("Finished gdal warp for: {}".format(out_dem_tmp_file))
    in_dem_file_ds = None
    out_dem_file_ds = None

    if not no_dem_check:
        logger.debug("Going to check DEM if no data value is not zero.")
        if no_data_val != 0:
            out_dem_tmp_chk_file = os.path.join(tmp_dir,
                                                os.path.splitext(os.path.basename(out_dem_file))[0] + '_tmpdem_chk.dem')
            logger.debug("Checking DEM as no data value is not zero; output image: {}".format(out_dem_tmp_chk_file))
            check_dem_for_gamma_compat(out_dem_tmp_file, out_dem_tmp_chk_file, no_data_val)
            out_dem_tmp_file = out_dem_tmp_chk_file

    logger.debug("Start byte swap to create: {}".format(out_dem_file))
    cmd = 'swap_bytes {0} {1} 2'.format(out_dem_tmp_file, out_dem_file)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error('Cmd Failed: {}'.format(cmd))
        raise e
    logger.debug("Finished byte swap creating: {}".format(out_dem_file))

    logger.debug("Creating dem par file.")
    out_par_file = open(out_dem_par_file, 'w')
    imageparameters = dict()
    imageparameters['title'] = os.path.basename(in_dem_file)
    imageparameters['datatype'] = 'INTEGER*2'
    imageparameters['width'] = x_pxls
    imageparameters['nlines'] = y_pxls
    imageparameters['corner_north'] = out_tlY
    imageparameters['corner_east'] = out_tlX
    imageparameters['post_north'] = out_res_y
    imageparameters['post_east'] = out_res_x
    imageparameters['zone'] = utm_zone
    imageparameters['center_longitude'] = (utm_zone * 6) - 183.0

    out = '''Gamma DIFF&GEO DEM/MAP parameter file
title: {title}
DEM_projection:     UTM
data_format:        {datatype}
DEM_hgt_offset:          0.00000
DEM_scale:               1.00000
width:               {width}
nlines:              {nlines}
corner_north:  {corner_north}   m
corner_east:   {corner_east}   m
post_north:    {post_north}   m
post_east:     {post_east}   m

ellipsoid_name: WGS 84
ellipsoid_ra:        6378137.000   m
ellipsoid_reciprocal_flattening:  298.2572236

datum_name: WGS 1984
datum_shift_dx:              0.000   m
datum_shift_dy:              0.000   m
datum_shift_dz:              0.000   m
datum_scale_m:         0.00000e+00
datum_rotation_alpha:  0.00000e+00   arc-sec
datum_rotation_beta:   0.00000e+00   arc-sec
datum_rotation_gamma:  0.00000e+00   arc-sec
datum_country_list Global Definition, WGS84, World

projection_name: UTM
projection_zone:                 {zone}
false_easting:           500000.000   m
false_northing:          0000000.000   m
projection_k0:            0.9996000
center_longitude:         {center_longitude}   decimal degrees
center_latitude:          0.0000000   decimal degrees\n'''.format(**imageparameters)

    # false_northing:          10000000.000   m  # For Southern Hemisphere
    out_par_file.write(out)
    out_par_file.close()
    logger.debug("Created dem par file.")
    logger.debug("Finished running subset_reproj_utm_dem_file")


def check_dem_for_gamma_compat(input_dem_img, output_dem_img, no_data_val):
    """
    DEM is checked for no data value (will be 0) and minimum value of 1 for Gamma compatibility.
    The DEM will be outputted in ENVI format.

    :param input_dem_img: Input DEM image
    :param output_img: Output DEM image
    :param no_data_val: The current no data value for the DEM.

    """
    logger.debug("Starting to check DEM for Gamma compatibility.")
    img_ds = gdal.Open(input_dem_img)
    if img_ds is None:
        raise Exception("Could not open image: {}".format(input_dem_img))

    # Get Header information.
    geotransform = img_ds.GetGeoTransform()
    x_pxls = img_ds.RasterXSize
    y_pxls = img_ds.RasterYSize
    proj_str = img_ds.GetProjection()

    out_img_ds = gdal.GetDriverByName('ENVI').Create(output_dem_img, x_pxls, y_pxls, 1, gdal.GDT_Int16)
    if out_img_ds == None:
        raise Exception("Could not create DEM image output raster: '{}'.".format(output_dem_img))
    out_img_ds.SetGeoTransform(geotransform)
    out_img_ds.SetProjection(proj_str)

    img_band = img_ds.GetRasterBand(1)
    if img_band is None:
        raise Exception("Could not open image band 1 from {}".format(input_dem_img))
    val_arr = img_band.ReadAsArray()

    out_img_band = out_img_ds.GetRasterBand(1)
    if out_img_band == None:
        raise Exception("Could not open image band 1 from {}".format(output_dem_img))
    out_img_band.SetNoDataValue(0)
    logger.debug("Created output file and read input data.")

    val_arr[val_arr == 0] = 1
    val_arr[val_arr == no_data_val] = 0
    val_arr[val_arr < 0] = 0

    out_img_band.WriteArray(val_arr)

    img_ds = None
    out_img_ds = None
    logger.debug("Written revised DEM to output image file - finished Gamma DEM compatibility function .")
