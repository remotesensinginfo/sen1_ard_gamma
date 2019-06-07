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
# Purpose:  Setup variables and imports across the whole module
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import os
import os.path
import glob
import shutil

import math

import subprocess

import osgeo.gdal as gdal

import sen1_ard_gamma.sen1_ard_utils
import sen1_ard_gamma.gamma_dem_prep
import sen1_ard_gamma.calc_img_stats

import logging
logger = logging.getLogger(__name__)

gdal.UseExceptions()

def get_sen1_latlong_bbox(ann_sen1_xml):
    """

    :param ann_sen1_xml:
    :return:

    """
    import xml.etree.ElementTree as ET
    tree = ET.parse(ann_sen1_xml)
    root = tree.getroot()

    geo_loc_grid_tag = root.find('geolocationGrid')
    if geo_loc_grid_tag is None:
        raise Exception("Could not find the geolocationGrid tag within product tag in Annotations XML file.")

    geo_loc_grid_pt_lst_tag = geo_loc_grid_tag.find('geolocationGridPointList')
    if geo_loc_grid_pt_lst_tag is None:
        raise Exception(
            "Could not find the geolocationGridPointList tag within geolocationGrid tag in Annotations XML file.")

    min_lat = 0.0
    max_lat = 0.0
    min_lon = 0.0
    max_lon = 0.0
    first = True
    for geo_loc_grid_pt_tag in geo_loc_grid_pt_lst_tag:
        if geo_loc_grid_pt_tag.tag == 'geolocationGridPoint':
            lat_val_tag = geo_loc_grid_pt_tag.find('latitude')
            if lat_val_tag is None:
                raise Exception(
                    "Could not find the latitude tag within geolocationGridPoint tag in Annotations XML file.")
            lon_val_tag = geo_loc_grid_pt_tag.find('longitude')
            if lon_val_tag is None:
                raise Exception(
                    "Could not find the longitude tag within geolocationGridPoint tag in Annotations XML file.")
            lat_str_val = lat_val_tag.text.strip()
            lon_str_val = lon_val_tag.text.strip()

            lat_val = sen1_ard_gamma.sen1_ard_utils.str2float(lat_str_val)
            lon_val = sen1_ard_gamma.sen1_ard_utils.str2float(lon_str_val)

            if first:
                min_lat = lat_val
                max_lat = lat_val
                min_lon = lon_val
                max_lon = lon_val
                first = False
            else:
                if lat_val < min_lat:
                    min_lat = lat_val
                if lat_val > max_lat:
                    max_lat = lat_val

                if lon_val < min_lon:
                    min_lon = lon_val
                if lon_val > max_lon:
                    max_lon = lon_val

    return [min_lon, max_lon, min_lat, max_lat]


def exe_gamma_grd_ard_processing(sen1img, ann_sen1_xml, cal_sen1_xml, nos_sen1_xml, in_dem_file, demparfile,
                                 outbasename, out_dir, tmp_dir, out_res_x, out_res_y, out_proj_epsg=None,
                                 use_dem_file=True, check_in_dem_filename=False, no_dem_check=False,
                                 dem_resample_method=gdal.GRA_CubicSpline):
    """
    Function to run through a Gamma processing chain to geocode and topographically correct Sentinel-1 GRD data.

    To control the number of threads used you can use the OMP_NUM_THREADS environmental variable.

    If out_proj_epsg is None then the UTM zone for the scene centre point is identified used as the output projection.

    :param sen1img: input tiff file from SAFE file for the polarisation to be processed.
    :param ann_sen1_xml: the annotation XML file for the polarisation being processed.
    :param cal_sen1_xml: the calibration XML file for the polarisation being processed.
    :param nos_sen1_xml: the noise XML file for the polarisation being prcessed.
    :param in_dem_file: the DEM file which is to be used for processing. Either:
                        * GDAL readable file which will be subset/reprojected and converted to gamma binary
                        * Gamma binary DEM file which corresponds to the gamma dem par file.
    :param demparfile: The Gamma dem par file. If the file exists and use_dem_file is True then the in_dem_file
                       will be used. If check_in_dem_filename is True then the file name for the in_dem_file will be
                       derived from the demparfile name by removing extension and replacing with extension .dem.
    :param outbasename: The base name for the output file. All output files will be outbasename + process_stage.ext
    :param out_dir: A directory where final geotiff files generated by the processing chain will be stored.
    :param tmp_dir: A directory where temporary output files generated by the processing chain will be stored.
    :param out_res_x: A float with the output ARD resolution in the x axis.
    :param out_res_y: A float within the output ARD resolution in the y axis.
    :param out_proj_epsg: An int for the EPSG code for the output projection. If None then UTM will be used with the
                          zone automatically derived from the centre point of the scene.
    :param use_dem_file: Boolean. If True the DEM file is expected to be available within Gamma format.
    :param check_in_dem_filename: If True then the DEM file name will be derived from the demparfile (by removing
                                  extension and adding .dem).
    :param no_dem_check: If True the DEM is not checked for no data value and minimum value for Gamma compatibility
    :param dem_resample_method: specify the gdal resampling method for the DEM processing (if DEM being subsetted
                                and reprojected). Default: gdal.GRA_CubicSpline

    """
    logger.info("Starting processing of {}.".format(sen1img))

    # -------------------------- STEP 0: Debug record function inputs -------------------------- #
    logger.debug("Annotation XML file: {}.".format(ann_sen1_xml))
    logger.debug("Calibration XML file: {}.".format(cal_sen1_xml))
    logger.debug("Noise XML file: {}.".format(nos_sen1_xml))
    logger.debug("Input DEM file: {}.".format(in_dem_file))
    logger.debug("DEM Par file: {}.".format(demparfile))
    logger.debug("Output file basename: {}.".format(outbasename))
    logger.debug("Output directory: {}.".format(out_dir))
    logger.debug("Temporary directory: {}.".format(tmp_dir))
    logger.debug("Output resolution X axis: {}.".format(out_res_x))
    logger.debug("Output resolution Y axis: {}.".format(out_res_y))
    logger.debug("Output EPSG projection: {}.".format(out_proj_epsg))
    logger.debug("use_dem_file: {}.".format(use_dem_file))
    logger.debug("check_in_dem_filename: {}.".format(check_in_dem_filename))
    logger.debug("dem_resample_method: {}.".format(dem_resample_method))
    # -------------------------- STEP 0: End -------------------------- #

    # -------------------------- STEP 1: Import DEM -------------------------- #
    logger.info("Import DEM")
    if os.path.exists(demparfile) and use_dem_file:
        logger.debug("Expecting DEM to already be in Gamma format...")
        if check_in_dem_filename:
            dem_test_file = os.path.splitext(demparfile)[0] + '.dem'
            logger.debug("Expecting DEM to have file name: '{}'".format(dem_test_file))
            if os.path.exists(dem_test_file):
                demfile = dem_test_file
        else:
            logger.debug("Expecting DEM to have given file name: '{}'".format(in_dem_file))
            demfile = in_dem_file
    else:
        logger.debug("Expecting DEM to be in GDAL readable format: '{}'".format(in_dem_file))
        dem_reproj_file = os.path.splitext(demparfile)[0] + '.dem'
        logger.debug("Output DEM will have file name: '{}'".format(dem_reproj_file))

        sen1_bbox_wgs84 = get_sen1_latlong_bbox(ann_sen1_xml)
        logger.debug("WGS84 BBOX of input Sentinel 1 scene: [{}, {}, {}, {}]".format(sen1_bbox_wgs84[0],
                                                                                     sen1_bbox_wgs84[1],
                                                                                     sen1_bbox_wgs84[2],
                                                                                     sen1_bbox_wgs84[3]))

        if out_proj_epsg is None:
            logger.debug("Finding UTM zone for scene from centre point")
            centre_lat = sen1_bbox_wgs84[2] + ((sen1_bbox_wgs84[3] - sen1_bbox_wgs84[2]) / 2)
            centre_lon = sen1_bbox_wgs84[0] + ((sen1_bbox_wgs84[1] - sen1_bbox_wgs84[0]) / 2)
            logger.debug("WGS84 Sentinel 1 scene centre point: [{}, {}]".format(centre_lat, centre_lon))
            utm_zone = sen1_ard_gamma.sen1_ard_utils.latlon_to_zone_number(centre_lat, centre_lon)
            logger.debug("UTM Zone for Scene: {}".format(utm_zone))
            out_proj_epsg = sen1_ard_gamma.sen1_ard_utils.epsg_for_UTM(utm_zone, 'N')

        logger.debug("EPSG Code for output: {}".format(out_proj_epsg))

        sen1_ard_gamma.gamma_dem_prep.subset_reproj_utm_dem_file(in_dem_file, sen1_bbox_wgs84, out_res_x, out_res_y,
                                                                 out_proj_epsg, dem_reproj_file, demparfile, tmp_dir,
                                                                 no_dem_check, eResampleAlg=dem_resample_method)

        demfile = dem_reproj_file
    logger.debug("Finished Importing DEM.")
    # -------------------------- STEP 1: End -------------------------- #

    # -------------------------- STEP 2: Import to Gamma -------------------------- #
    logger.info("Import Sentinel-1 Scene to Gamma.")
    sen1_img_import = os.path.join(tmp_dir, outbasename + '_import')
    sen1_par_import = os.path.join(tmp_dir, outbasename + '_import.par')
    sen1_nesz_import = os.path.join(tmp_dir, outbasename + '_import_nesz')

    cmd = 'par_S1_GRD {0} {1} {2} {3} {4} {5} - - 1 - {6}'.format(sen1img, ann_sen1_xml, cal_sen1_xml, nos_sen1_xml,
                                                                  sen1_par_import, sen1_img_import, sen1_nesz_import)

    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.debug("Finished importing Sentinel-1 Scene to Gamma.")
    # -------------------------- STEP 2: End -------------------------- #

    # ------ STEP 3: Read resolution and inc. angle info from image and dem parameter file ------------ #
    pardict = sen1_ard_gamma.sen1_ard_utils.gamma_par_dict(sen1_par_import)
    mli_rgres = float(pardict['range_pixel_spacing'][0])
    mli_azres = float(pardict['azimuth_pixel_spacing'][0])
    center_lat = float(pardict['center_latitude'][0])
    inc = float(pardict['incidence_angle'][0])

    demdict = sen1_ard_gamma.sen1_ard_utils.gamma_par_dict(demparfile)
    proj = demdict['DEM_projection'][0]
    if proj == 'EQA':
        dem_xres = float(demdict['post_lon'][0])
        dem_yres = abs(float(demdict['post_lat'][0]))
    else:
        dem_xres = float(demdict['post_east'][0])
        dem_yres = abs(float(demdict['post_north'][0]))
    # -------------------------- STEP 3: End -------------------------- #

    # -------------------------- STEP 4: Calculate DEM oversampling factors -------------------------- #
    ovr_lon = abs(dem_xres / out_res_y)
    ovr_lat = abs(dem_yres / out_res_x)
    # -------------------------- STEP 4: End -------------------------- #

    # -------------------------- STEP 5: Calculate multilooking factors -------------------------- #
    if proj == 'EQA':
        deg_in_m = 111319 * math.cos(center_lat * 3.141 / 180)
        geo_xres_m = deg_in_m * out_res_y
        geo_yres_m = 110946 * abs(out_res_x)
        ml_rg = int(math.ceil(geo_xres_m / mli_rgres))
        ml_az = int(math.ceil(geo_yres_m / mli_azres))
        if ((ml_rg * mli_rgres) / geo_xres_m) > 1.25:
            ml_rg = int(math.floor(geo_xres_m / mli_rgres))
        if ((ml_az * mli_azres) / geo_yres_m) > 1.25:
            ml_az = int(math.floor(geo_yres_m / mli_azres))
    else:
        ml_rg = int(math.ceil(out_res_y / mli_rgres))
        ml_az = int(math.ceil(out_res_x / mli_azres))
        if ((ml_rg * mli_rgres) / out_res_y) > 1.25:
            ml_rg = int(math.floor(out_res_y / mli_rgres))
        if ((ml_az * mli_azres) / out_res_x) > 1.25:
            ml_az = int(math.floor(out_res_x / mli_azres))

            # In case calculation results in ml-factors of 0
    if (ml_rg < 1):
        ml_rg = 1
    if (ml_az < 1):
        ml_az = 1
    # -------------------------- STEP 5: End -------------------------- #

    # -------------------------- STEP 6: Multilooking -------------------------- #
    logger.info("Multilooking: {0} x {1}".format(ml_rg, ml_az))

    sen1_img_mli = os.path.join(tmp_dir, outbasename + '_import_mli')
    sen1_par_mli = os.path.join(tmp_dir, outbasename + '_import_mli.par')

    cmd = 'multi_look_MLI {0} {1} {2} {3} {4} {5}'.format(sen1_img_import, sen1_par_import, sen1_img_mli, sen1_par_mli,
                                                          ml_rg, ml_az)

    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    logger.debug("Finished Multilooking: {0} x {1}".format(ml_rg, ml_az))
    # -------------------------- STEP 6: End -------------------------- #

    # -------------------------- STEP 7: Read width and lines from mlipar -------------------------- #
    mlidict = sen1_ard_gamma.sen1_ard_utils.gamma_par_dict(sen1_par_mli)
    mli_width = int(mlidict['range_samples'][0])
    mli_lines = int(mlidict['azimuth_lines'][0])
    # -------------------------- STEP 7: End -------------------------- #

    # -------------------------- STEP 8: Subtract noise equivalent sigma nought from MLI -------------------------- #
    logger.info("Multilook S1 GRD noise image...")
    sen1_nesz_mli = os.path.join(tmp_dir, outbasename + '_import_nesz_mli')
    sen1_nesz_par_mli = os.path.join(tmp_dir, outbasename + '_import_nesz_mli.par')
    cmd = 'multi_look_MLI {0} {1} {2} {3} {4} {5}'.format(sen1_nesz_import, sen1_par_import, sen1_nesz_mli,
                                                          sen1_nesz_par_mli, ml_rg, ml_az)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.debug("Finished multilook S1 GRD noise image...")

    logger.info('Subtract NESZ from MLI')
    sen1_img_mli_cor = os.path.join(tmp_dir, outbasename + '_import_mli_cor')

    cmd = 'float_math {0} {1} {2} {3} 1 - - - - - 1'.format(sen1_img_mli, sen1_nesz_mli, sen1_img_mli_cor, mli_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.info('Finished Subtract NESZ from MLI')

    logger.info("Remove pwr values below -40 dB")
    sen1_img_mli_cor2 = os.path.join(tmp_dir, outbasename + '_import_mli_cor2')
    cmd = 'replace_values {0} 0.0001 0.0001 {1} {2} 2 2'.format(sen1_img_mli_cor, sen1_img_mli_cor2, mli_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_mli = sen1_img_mli_cor2
    logger.debug("Finished remove pwr values below -40 dB")
    # -------------------------- STEP 8: End -------------------------- #

    # -------------------------- STEP 9: Calibration -------------------------- #
    logger.info("Running calibration")
    sen1_img_cmli = os.path.join(tmp_dir, outbasename + '_import_cmli')
    sen1_img_pix_ell = os.path.join(tmp_dir, outbasename + '_import_pix_ell')
    cmd = 'radcal_MLI {0} {1} - {2} - 0 0 1 0 - {3}'.format(sen1_img_mli, sen1_par_mli, sen1_img_cmli, sen1_img_pix_ell)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_cmli = sen1_img_mli
    logger.debug("Finished calibration")
    # -------------------------- STEP 9: End -------------------------- #

    # -------------------------- STEP 10: Compute Geocoding Lookup Table -------------------------- #
    logger.info('Compute Geocoding Lookup Table')

    sen1_img_inc = os.path.join(tmp_dir, outbasename + '_import_inc')
    sen1_img_pix = os.path.join(tmp_dir, outbasename + '_import_pix')
    sen1_img_lsmap = os.path.join(tmp_dir, outbasename + '_import_lsmap')
    sen1_img_geo2rdc = os.path.join(tmp_dir, outbasename + '_import_geo2rdc')
    sen1_img_simsar = os.path.join(tmp_dir, outbasename + '_import_simsar')
    sen1_img_gcdem = os.path.join(tmp_dir, outbasename + '_import_gcdem')
    sen1_par_gcdem = os.path.join(tmp_dir, outbasename + '_import_gcdem.par')

    cmd = 'gc_map {0} - {1} {2} {3} {4} {5} {6} {7} {8} - - {9} - {10} {11} - 2'.format(sen1_par_mli, demparfile,
                                                                                        demfile, sen1_par_gcdem,
                                                                                        sen1_img_gcdem,
                                                                                        sen1_img_geo2rdc, ovr_lon,
                                                                                        ovr_lat, sen1_img_simsar,
                                                                                        sen1_img_inc, sen1_img_pix,
                                                                                        sen1_img_lsmap)

    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    logger.debug('Finished Computing Geocoding Lookup Table')
    # -------------------------- STEP 10: End -------------------------- #

    # -------------------------- STEP 11: Read width and lines from sen1_par_gcdem -------------------------- #
    gcdemdict = sen1_ard_gamma.sen1_ard_utils.gamma_par_dict(sen1_par_gcdem)
    dem_width = int(gcdemdict['width'][0])
    # -------------------------- STEP 11: End -------------------------- #

    # -------------------------- STEP 12: Estimate pixel scattering area based on DEM -------------------------- #
    logger.info('Estimate pixel scattering area based on DEM')

    sen1_img_pixdem = os.path.join(tmp_dir, outbasename + '_import_pix_dem')

    cmd = 'pixel_area {0} {1} {2} {3} {4} {5} {6}'.format(sen1_par_mli, sen1_par_gcdem, sen1_img_gcdem,
                                                          sen1_img_geo2rdc, sen1_img_lsmap, sen1_img_inc,
                                                          sen1_img_pixdem)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    logger.debug('Finished estimating pixel scattering area based on DEM')
    # -------------------------- STEP 12: End -------------------------- #

    # -------------------------- STEP 13: Improve Geocoding Lookup Table -------------------------- #
    logger.info('Check Accuracy of Geocoding Lookup Table')

    sen1_par_diff = os.path.join(tmp_dir, outbasename + '_import_diff.par')
    cmd = 'create_diff_par {0} - {1} 1 0'.format(sen1_par_mli, sen1_par_diff)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_offs = os.path.join(tmp_dir, outbasename + '_import_offs')
    sen1_img_ccp = os.path.join(tmp_dir, outbasename + '_import_ccp')
    sen1_img_offsets = os.path.join(tmp_dir, outbasename + '_import_offsets')
    cmd = 'offset_pwrm {0} {1} {2} {3} {4} 128 128 {5} - 64 64 0.2'.format(sen1_img_cmli, sen1_img_pixdem,
                                                                           sen1_par_diff, sen1_img_offs, sen1_img_ccp,
                                                                           sen1_img_offsets)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_coffs = os.path.join(tmp_dir, outbasename + '_import_coffs')
    sen1_img_coffsets = os.path.join(tmp_dir, outbasename + '_import_coffsets')
    cmd = 'offset_fitm {0} {1} {2} {3} {4} 0.2 - 0'.format(sen1_img_offs, sen1_img_ccp, sen1_par_diff, sen1_img_coffs,
                                                           sen1_img_coffsets)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    logger.debug('Finished Checking Accuracy of Geocoding Lookup Table')
    # -------------------------- STEP 13: End -------------------------- #

    # -------------------------- STEP 14: Perform Fine Registration -------------------------- #
    logger.info('Improve Geocoding Lookup Table')

    sen1_img_geo2rdc_fine = os.path.join(tmp_dir, outbasename + '_import_geo2rdc_fine')
    cmd = 'gc_map_fine {0} {1} {2} {3} 0'.format(sen1_img_geo2rdc, dem_width, sen1_par_diff, sen1_img_geo2rdc_fine)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.info('Finished improving Geocoding Lookup Table')
    sen1_img_geo2rdc = sen1_img_geo2rdc_fine

    logger.info('Update Pixel Area map')
    cmd = 'pixel_area {0} {1} {2} {3} {4} {5} {6}'.format(sen1_par_mli, sen1_par_gcdem, sen1_img_gcdem,
                                                          sen1_img_geo2rdc, sen1_img_lsmap, sen1_img_inc,
                                                          sen1_img_pixdem)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.debug('Finished updating Pixel Area map')
    # -------------------------- STEP 14: End -------------------------- #

    # -------------------------- STEP 15: Pixel Area Correction -------------------------- #
    logger.info('Pixel Area Correction')

    cmd = 'float_math {0} {1} {2} {3} 3 - - 1 1 - '.format(sen1_img_pix_ell, sen1_img_pixdem, sen1_img_pix, mli_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_s0 = os.path.join(tmp_dir, outbasename + '_import_s0')
    cmd = 'float_math {0} {1} {2} {3}  2 - - 1 1 - '.format(sen1_img_cmli, sen1_img_pix, sen1_img_s0, mli_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.debug('Finished Pixel Area Correction')
    # -------------------------- STEP 15: End -------------------------- #

    # -------------------------- STEP 16: Geocode and pixel area correction for GRDs -------------------------- #
    logger.info('Geocode to projection: {}'.format(proj))
    sen1_img_s0geo = os.path.join(tmp_dir, outbasename + '_import_s0geo')
    cmd = 'geocode_back {0} {1} {2} {3} {4} - 7 0'.format(sen1_img_s0, mli_width, sen1_img_geo2rdc, sen1_img_s0geo,
                                                          dem_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_pixgeo = os.path.join(tmp_dir, outbasename + '_import_pixgeo')
    cmd = 'geocode_back {0} {1} {2} {3} {4} - 1 0 '.format(sen1_img_pix, mli_width, sen1_img_geo2rdc, sen1_img_pixgeo,
                                                           dem_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    sen1_img_pix = sen1_img_pixgeo
    logger.debug('Finished geocode to projection: {}'.format(proj))
    # -------------------------- STEP 16: End -------------------------- #

    # -------------------------- STEP 17: Sigma to Gamma conversion -------------------------- #
    logger.info('Convert Sigma to gamma.')
    sen1_img_g0geo = os.path.join(tmp_dir, outbasename + '_import_g0geo')
    cmd = 'sigma2gamma {0} {1} {2} {3} '.format(sen1_img_s0geo, sen1_img_inc, sen1_img_g0geo, dem_width)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.debug('Finished converting Sigma to gamma.')

    # -------------------------- STEP 17: End -------------------------- #

    # -------------------------- STEP 18: Export to geotiff -------------------------- #
    logger.info('Convert file outputs to GeoTiff')

    sen1_img_pwrout = os.path.join(out_dir, outbasename + '_final_pwr.tif')
    sen1_img_incout = os.path.join(out_dir, outbasename + '_final_inc.tif')
    sen1_img_lsout = os.path.join(out_dir, outbasename + '_final_lsmap.tif')
    sen1_img_pixout = os.path.join(out_dir, outbasename + '_final_pix.tif')

    cmd = 'data2geotiff {0} {1} 2 {2}'.format(sen1_par_gcdem, sen1_img_g0geo, sen1_img_pwrout)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    cmd = 'data2geotiff {0} {1} 2 {2}'.format(sen1_par_gcdem, sen1_img_inc, sen1_img_incout)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    cmd = 'data2geotiff {0} {1} 2 {2}'.format(sen1_par_gcdem, sen1_img_pix, sen1_img_pixout)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e

    cmd = 'data2geotiff {0} {1} 2 {2}'.format(sen1_par_gcdem, sen1_img_lsmap, sen1_img_lsout)
    try:
        logger.debug("Running following command using subprocess '{}'".format(cmd))
        cmd = sen1_ard_gamma.sen1_ard_utils.preappend_cmd(cmd)
        subprocess.call(cmd, shell=True)
    except OSError as e:
        logger.error("Command failed: '{}'".format(cmd))
        raise e
    logger.debug('Finished converting file outputs to GeoTiff')
    # -------------------------- STEP 18: End -------------------------- #

    logger.info("Completed processing of {}.".format(sen1img))


def create_pol_stacked_products(out_scns_dict, out_base_name, tmp_dir, out_dir, gdal_format='GTIFF'):
    """
    A function to create the final stacked ARD image products.

    :param out_scns_dict: dictionary with polarisation at the first level and second level of (inc, lsmap, pix, pwr)
    :param out_base_name: output image file base name
    :param tmp_dir: tmp directory to sort tempory images to
    :param out_dir: output directory to export final images to
    :param gdal_format: the output GDAL image format (Default: GTIFF)

    """
    logger.info("Starting to stack images to create final products")
    out_img_list = dict()
    img_ext = sen1_ard_gamma.sen1_ard_utils.get_file_extension(gdal_format)
    scn_keys = list(out_scns_dict.keys())
    gdal_opts = ""
    if gdal_format == 'GTIFF':
        gdal_opts = "-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_NEEDED "
    if len(scn_keys) == 2:
        if ('vv' in scn_keys) and ('vh' in scn_keys):
            # ------------ Process inc Product ------------- #
            out_img = os.path.join(out_dir, "{}_inc{}".format(out_base_name, img_ext))
            cmd = "gdal_translate -of {0} {1} {2} {3}".format(gdal_format, gdal_opts, out_scns_dict['vv']['inc'],
                                                              out_img)
            try:
                logger.debug("Running following command using subprocess '{}'".format(cmd))
                subprocess.call(cmd, shell=True)
            except OSError as e:
                logger.error('Cmd Failed: {}'.format(cmd))
                raise e
            sen1_ard_gamma.sen1_ard_utils.set_band_names(out_img, ['inc'])
            out_img_list['inc'] = out_img
            # ------------ Processed inc Product ------------- #

            # ------------ Process pix Product ------------- #
            out_vrt_tmp = os.path.join(tmp_dir, "{}_pix_tmp.vrt".format(out_base_name))
            cmd = "gdalbuildvrt -separate {} {} {}".format(out_vrt_tmp, out_scns_dict['vv']['pix'],
                                                           out_scns_dict['vh']['pix'])
            try:
                logger.debug("Running following command using subprocess '{}'".format(cmd))
                subprocess.call(cmd, shell=True)
            except OSError as e:
                logger.error('Cmd Failed: {}'.format(cmd))
                raise e
            out_img = os.path.join(out_dir, "{}_pix{}".format(out_base_name, img_ext))
            cmd = "gdal_translate -of {0} {1} {2} {3}".format(gdal_format, gdal_opts, out_vrt_tmp, out_img)
            try:
                logger.debug("Running following command using subprocess '{}'".format(cmd))
                subprocess.call(cmd, shell=True)
            except OSError as e:
                logger.error('Cmd Failed: {}'.format(cmd))
                raise e
            sen1_ard_gamma.sen1_ard_utils.set_band_names(out_img, ['vv', 'vh'])
            out_img_list['pix'] = out_img
            # ------------ Processed pix Product ------------- #

            # ------------ Process pwr Product ------------- #
            out_ratio_tmp_file = os.path.join(tmp_dir, "{}_pwr_ratio_tmp.tif".format(out_base_name))
            sen1_ard_gamma.sen1_ard_utils.calc_ratio_img(out_scns_dict['vv']['pwr'], out_scns_dict['vh']['pwr'],
                                                         out_ratio_tmp_file, gdal_format)

            out_vrt_tmp = os.path.join(tmp_dir, "{}_pwr_tmp.vrt".format(out_base_name))
            cmd = "gdalbuildvrt -separate {} {} {} {}".format(out_vrt_tmp, out_scns_dict['vv']['pwr'],
                                                              out_scns_dict['vh']['pwr'], out_ratio_tmp_file)
            try:
                logger.debug("Running following command using subprocess '{}'".format(cmd))
                subprocess.call(cmd, shell=True)
            except OSError as e:
                logger.error('Cmd Failed: {}'.format(cmd))
                raise e
            out_img = os.path.join(out_dir, "{}_pwr{}".format(out_base_name, img_ext))
            cmd = "gdal_translate -of {0} {1} {2} {3}".format(gdal_format, gdal_opts, out_vrt_tmp, out_img)
            try:
                logger.debug("Running following command using subprocess '{}'".format(cmd))
                subprocess.call(cmd, shell=True)
            except OSError as e:
                logger.error('Cmd Failed: {}'.format(cmd))
                raise e
            sen1_ard_gamma.sen1_ard_utils.set_band_names(out_img, ['vv', 'vh', 'vv/vh'])
            out_img_list['pwr'] = out_img
            # ------------ Processed pwr Product ------------- #
        else:
            raise Exception("Do not recognise the polarisations provided")
    elif len(scn_keys) == 1:
        pol = scn_keys[0]
        for prod in ['inc', 'pix', 'pwr']:
            out_img = os.path.join(out_dir, "{}_{}{}".format(out_base_name, prod, img_ext))
            cmd = "gdal_translate -of {0} {1} {2} {3}".format(gdal_format, gdal_opts, out_scns_dict[pol][prod], out_img)
            try:
                logger.debug("Running following command using subprocess '{}'".format(cmd))
                subprocess.call(cmd, shell=True)
            except OSError as e:
                logger.error('Cmd Failed: {}'.format(cmd))
                raise e
            if prod in ['inc', 'lsmap']:
                sen1_ard_gamma.sen1_ard_utils.set_band_names(out_img, [prod])
            else:
                sen1_ard_gamma.sen1_ard_utils.set_band_names(out_img, [pol])
            out_img_list[prod] = out_img
    else:
        raise Exception("Only know how to stack images with 2 bands or copy single band images.")
    if 'pwr' in out_img_list:
        logger.debug("Calculate dB image from power.")
        out_dB_img = os.path.join(out_dir, "{}_dB{}".format(out_base_name, img_ext))
        sen1_ard_gamma.sen1_ard_utils.convert_to_dB(out_img_list['pwr'], out_dB_img, gdal_format)
        out_img_list['dB'] = out_dB_img
    return out_img_list


def run_sen1_grd_ard_analysis(input_safe_file, output_dir, tmp_dir, dem_img_file, out_img_res, out_proj_epsg,
                              polarisations, gdal_format, calc_no_stats, keep_files, no_dem_check):
    """
    High level function which runs the analysis of a Sentinel-1 GRD scene to an ARD product.

    :param input_safe_file:
    :param output_dir:
    :param tmp_dir:
    :param dem_img_file:
    :param out_img_res:
    :param out_proj_epsg:
    :param polarisations:
    :param gdal_format:
    :param calc_no_stats:
    :param keep_files:
    :param no_dem_check: If True the DEM is not checked for no data value and minimum value for Gamma compatibility
    """
    scn_metadata_info = sen1_ard_gamma.sen1_ard_utils.retrieve_sentinel1_metadata(input_safe_file)

    scn_safe_files = sen1_ard_gamma.sen1_ard_utils.find_sen1_ard_files(input_safe_file)

    scn_basename = sen1_ard_gamma.sen1_ard_utils.create_sentinel1_basename(scn_metadata_info)
    print("Basename for scene: {}".format(scn_basename))

    if polarisations is None:
        polarisations = scn_metadata_info['product_polarisations']
    else:
        for pol in polarisations:
            if pol not in scn_metadata_info['product_polarisations']:
                raise Exception("Polarisation {} is not within the scene provided.")

    c_uid = sen1_ard_gamma.sen1_ard_utils.uidGenerator()
    c_tmp_dir = os.path.join(tmp_dir, '{}_tmp_{}'.format(scn_basename, c_uid))
    c_tmp_dir_created = False
    if not os.path.exists(c_tmp_dir):
        os.mkdir(c_tmp_dir)
        c_tmp_dir_created = True

    c_out_dir = os.path.join(tmp_dir, '{}_{}'.format(scn_basename, c_uid))
    c_out_dir_created = False
    if not os.path.exists(c_out_dir):
        os.mkdir(c_out_dir)
        c_out_dir_created = True

    demparfile = os.path.join(c_tmp_dir, scn_basename + '_gamma_scn_dem.dem_par')

    out_scns = dict()
    first = True
    for pol in polarisations:
        pol_lower = pol.lower()
        c_scn_basename = scn_basename + '_' + pol_lower
        logger.info("Processing {} Polarisation.".format(pol))

        sen1_ard_gamma.sen1_grd_ard_tools.exe_gamma_grd_ard_processing(scn_safe_files['measure_' + pol_lower],
                                                                       scn_safe_files['annotation_' + pol_lower],
                                                                       scn_safe_files['calibration_' + pol_lower],
                                                                       scn_safe_files['noise_' + pol_lower],
                                                                       dem_img_file, demparfile, c_scn_basename, c_out_dir,
                                                                       c_tmp_dir, out_img_res, -out_img_res,
                                                                       out_proj_epsg, use_dem_file=(not first),
                                                                       check_in_dem_filename=True, no_dem_check=False,
                                                                       dem_resample_method=gdal.GRA_CubicSpline)

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

    logger.info("Creating final image outputs")
    # Create final output images with bands stacked where appropriate and ratio image (vv/vh) calculated.
    fnl_out_imgs = sen1_ard_gamma.sen1_grd_ard_tools.create_pol_stacked_products(out_scns, scn_basename, c_out_dir,
                                                                                 output_dir, gdal_format)

    if not calc_no_stats:
        logger.info("Calculating image statistics and pyramids for final outputs")
        if 'inc' in fnl_out_imgs:
            sen1_ard_gamma.calc_img_stats.run_calc_img_stats_pyramids(fnl_out_imgs['inc'], no_data_val=0.0)
        if 'pix' in fnl_out_imgs:
            sen1_ard_gamma.calc_img_stats.run_calc_img_stats_pyramids(fnl_out_imgs['pix'], no_data_val=0.0)
        if 'pwr' in fnl_out_imgs:
            sen1_ard_gamma.calc_img_stats.run_calc_img_stats_pyramids(fnl_out_imgs['pwr'], no_data_val=0.0)
        if 'dB' in fnl_out_imgs:
            sen1_ard_gamma.calc_img_stats.run_calc_img_stats_pyramids(fnl_out_imgs['dB'], no_data_val=999)

    if not keep_files:
        if c_tmp_dir_created:
            shutil.rmtree(c_tmp_dir)
        if c_out_dir_created:
            shutil.rmtree(c_out_dir)

    logger.info("Completed all data processing stages; outputs in {}".format(output_dir))
