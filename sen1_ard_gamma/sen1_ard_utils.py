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
# Purpose:  Provide utility functions for the module
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import logging
import os
import os.path
import glob
import datetime

import numpy

import osgeo.ogr as ogr
import osgeo.gdal as gdal

logger = logging.getLogger(__name__)

gdal.UseExceptions()


def preappend_cmd(cmd):
    """
    For using docker or singularity the command needs to be pre-appended and potentially the
    file paths need to be updated to use a local mount within a docker image.

    For the command to be pre-appended the S1ARD_PAP_CMD environmental variable needs defined.

    For the file path to be updated the S1ARD_PAP_PATH environmental variable needs to be defined
    with local_path:mount_path

    :param cmd: string with command to be run.
    :return: string with pre-appended command.
    """
    pap_cmd = os.getenv('S1ARD_PAP_CMD', None)
    if pap_cmd is None:
        logger.debug("S1ARD_PAP_CMD is not defined.")
        return cmd
    logger.debug("S1ARD_PAP_CMD is defined: '{}'.".format(pap_cmd))

    out_cmd = cmd
    pap_path = os.getenv('S1ARD_PAP_PATH', None)
    if pap_cmd is not None:
        logger.debug("S1ARD_PAP_PATH is defined: '{}'.".format(pap_path))
        pap_paths = pap_path.split(':')
        lcl_path = pap_paths[0]
        rmt_path = pap_paths[1]
        logger.debug("Local Path to be replaced: '{}'.".format(lcl_path))
        logger.debug("Remote Path to be populated: '{}'.".format(rmt_path))
        out_cmd = cmd.replace(lcl_path, rmt_path)
    logger.debug("Command file paths updated: '{}'.".format(out_cmd))

    fnl_out_cmd = "{} {}".format(pap_cmd, out_cmd)
    logger.debug("Pre-appended command outputted: '{}'.".format(fnl_out_cmd))
    return fnl_out_cmd


def uidGenerator(size=6):
    """
    A function which will generate a 'random' string of the specified length based on the UUID
    """
    import uuid
    randomStr = str(uuid.uuid4())
    randomStr = randomStr.replace("-","")
    return randomStr[0:size]


def gamma_par_dict(par_file):
    """
    Parse a Gamma parameter (par) file to create a dictionary

    :param parfile: input par file
    :return:
    """
    gamma_dict = dict()
    lines = open(par_file,"r").readlines()
    for l in lines:
        l = l.replace('\n','')
        if len(l.split(':')) > 1 and len(l.split(':')) < 3:
            key,values = l.split(':')
            values = values.split()
            gamma_dict[key] = values
    return gamma_dict


def str2float(str_val, err_val=None):
    """
    Convert a string to a float value, returning an error
    value if an error occurs. If no error value is provided
    then an exception is thrown.

    :param str_val: string variable containing float value.
    :param err_val: value to be returned if error occurs. If None then exception returned. Default None.
    :return: float
    """
    str_val = str(str_val).strip()
    out_flt = 0.0
    try:
        out_flt = float(str_val)
    except ValueError:
        if not err_val is None:
            out_flt = float(err_val)
        else:
            raise Exception("Could not convert string to float: \'" + str_val + '\'.')
    return out_flt


def metres_to_degrees(latitude, xsize, ysize):
    """

    :param latitude: latitude is degrees
    :param xsize: (numpy array) value of x pixel sizes (m)
    :param ysize: (numpy array) value of y pixel sizes (m)
    :return: (lonsize, latsize) (numpy array(s)) value(s) of x and y pixel sizes (degrees)

    Example::

    lonsize, latsize = metres_to_degrees(52,1.0,1.0)

    """
    # Set up parameters for ellipse; Semi-major and semi-minor for WGS-84 ellipse
    ellipse = [6378137.0, 6356752.314245]

    radlat = numpy.deg2rad(latitude)

    Rsq = (ellipse[0] * numpy.cos(radlat)) ** 2 + (ellipse[1] * numpy.sin(radlat)) ** 2
    Mlat = (ellipse[0] * ellipse[1]) ** 2 / (Rsq ** 1.5)
    Nlon = ellipse[0] ** 2 / numpy.sqrt(Rsq)
    lonsize = xsize / (numpy.pi / 180 * numpy.cos(radlat) * Nlon)
    latsize = ysize / (numpy.pi / 180 * Mlat)

    return lonsize, latsize


def reproj_point(in_proj_osr_obj, out_proj_osr_obj, x, y):
    """
    Reproject a point from 'in_proj_osr_obj' to 'out_proj_osr_obj' where they are gdal
    osgeo.osr.SpatialReference objects.

    :param in_proj_osr_obj: osgeo.osr.SpatialReference representation of projection of input point
    :param out_proj_osr_obj: osgeo.osr.SpatialReference representation of projection of output point
    :param x: x coordinate value (float)
    :param y: y coordinate value (float)
    :return: x, y
    """
    if in_proj_osr_obj.EPSGTreatsAsLatLong():
        wktPt = 'POINT(%s %s)' % (y, x)
    else:
        wktPt = 'POINT(%s %s)' % (x, y)
    point = ogr.CreateGeometryFromWkt(wktPt)
    point.AssignSpatialReference(in_proj_osr_obj)
    point.TransformTo(out_proj_osr_obj)
    if out_proj_osr_obj.EPSGTreatsAsLatLong():
        outX = point.GetY()
        outY = point.GetX()
    else:
        outX = point.GetX()
        outY = point.GetY()
    return outX, outY


def epsg_for_UTM(zone, hemisphere):
    """
Return EPSG code for given UTM zone and hemisphere using WGS84 datum.
:param zone: UTM zone
:param hemisphere: hemisphere either 'N' or 'S'

:return: corresponding EPSG code

"""
    if hemisphere not in ['N', 'S']:
        raise Exception('Invalid hemisphere ("N" or "S").')

    if zone < 0 or zone > 60:
        raise Exception('UTM zone outside valid range.')

    if hemisphere == 'N':
        ns = 600
    else:
        ns = 700

    if zone == 0:
        zone = 61

    return int(32000 + ns + zone)


def latlon_to_zone_number(latitude, longitude):
    """
Find the UTM zone number for a give latitude and longitude. If the input is a numpy array, just use the
first element user responsibility to make sure that all points are in one zone

:param latitude: float
:param longitude: float

:return: int

"""
    if 56 <= latitude < 64 and 3 <= longitude < 12:
        return 32

    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude < 9:
            return 31
        elif longitude < 21:
            return 33
        elif longitude < 33:
            return 35
        elif longitude < 42:
            return 37

    return int((longitude + 180) / 6) + 1


def findFile(dirPath, fileSearch, raise_exp=True):
    """
    Search for a single file with a path using glob. Therefore, the file
    path returned is a true path. Within the fileSearch provide the file
    name with '*' as wildcard(s).

    :return: string

    """
    files = glob.glob(os.path.join(dirPath, fileSearch))
    rtn_file = None
    if len(files) != 1:
        if raise_exp:
            raise Exception('Could not find a single file ('+fileSearch+'); found ' + str(len(files)) + ' files.')
    else:
        rtn_file = files[0]
    return rtn_file


def retrieve_sentinel1_metadata(input_safe_file):
    """
    The function parses the input manifest.safe to create a base name for the image.

    :param input_safe_file: the SAFE file directory (*.SAFE)
    :return: a dict with sentinel-1 metadata

    """
    import xml.etree.ElementTree as ET

    manifest_file = findFile(input_safe_file, 'manifest*.safe')

    tree = ET.parse(manifest_file)
    root = tree.getroot()

    metadata_section_tag = root.find('metadataSection')
    if metadata_section_tag is None:
        raise Exception("Could not find the metadataSection tag within product tag in {}.".format(manifest_file))

    scn_metadata_info = dict()

    for metadata_obj_tag in metadata_section_tag:
        if metadata_obj_tag.tag == 'metadataObject':
            if metadata_obj_tag.attrib['ID'] == 'platform':
                platform_tag = metadata_obj_tag.find("metadataWrap").find("xmlData").find(
                    "{http://www.esa.int/safe/sentinel-1.0}platform")
                scn_metadata_info['nssdc_identifier'] = platform_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0}nssdcIdentifier").text.strip()
                scn_metadata_info['family_name'] = platform_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0}familyName").text.strip()
                scn_metadata_info['sensor_number'] = platform_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0}number").text.strip()
                scn_metadata_info['mode'] = platform_tag.find("{http://www.esa.int/safe/sentinel-1.0}instrument").find(
                    "{http://www.esa.int/safe/sentinel-1.0}extension").find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}instrumentMode").find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}mode").text.strip()
            elif metadata_obj_tag.attrib['ID'] == 'measurementOrbitReference':
                orbit_ref_tag = metadata_obj_tag.find("metadataWrap").find("xmlData").find(
                    "{http://www.esa.int/safe/sentinel-1.0}orbitReference")
                scn_metadata_info['cycle_number'] = orbit_ref_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0}cycleNumber").text.strip()
                scn_metadata_info['pass'] = orbit_ref_tag.find("{http://www.esa.int/safe/sentinel-1.0}extension").find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1}orbitProperties").find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1}pass").text.strip()
                for orbit_tag in orbit_ref_tag:
                    if (orbit_tag.tag == '{http://www.esa.int/safe/sentinel-1.0}orbitNumber') and orbit_tag.attrib['type'] == 'start':
                        scn_metadata_info['orbit_number'] = orbit_tag.text.strip()
                    elif (orbit_tag.tag == '{http://www.esa.int/safe/sentinel-1.0}relativeOrbitNumber') and orbit_tag.attrib['type'] == 'start':
                        scn_metadata_info['rel_orbit_number'] = orbit_tag.text.strip()
            elif metadata_obj_tag.attrib['ID'] == 'generalProductInformation':
                product_info_tag = metadata_obj_tag.find("metadataWrap").find("xmlData").find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}standAloneProductInformation")

                scn_metadata_info['instrument_config_id'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}instrumentConfigurationID").text.strip()
                scn_metadata_info['mission_data_take_id'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}missionDataTakeID").text.strip()
                scn_metadata_info['product_class'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}productClass").text.strip()
                scn_metadata_info['product_class_dscp'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}productClassDescription").text.strip()
                scn_metadata_info['product_comp'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}productComposition").text.strip()
                scn_metadata_info['product_type'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}productType").text.strip()
                scn_metadata_info['product_timelines_cat'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}productTimelinessCategory").text.strip()
                scn_metadata_info['slice_prod_flag'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}sliceProductFlag").text.strip()
                scn_metadata_info['seg_start_time'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}segmentStartTime").text.strip()
                scn_metadata_info['slice_number'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}sliceNumber").text.strip()
                scn_metadata_info['total_slices'] = product_info_tag.find(
                    "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}totalSlices").text.strip()
                scn_metadata_info['product_polarisations'] = []
                for prod_tag in product_info_tag:
                    if prod_tag.tag == "{http://www.esa.int/safe/sentinel-1.0/sentinel-1/sar/level-1}transmitterReceiverPolarisation":
                        scn_metadata_info['product_polarisations'].append(prod_tag.text.strip())
            elif metadata_obj_tag.attrib['ID'] == 'acquisitionPeriod':
                acq_period_tag = metadata_obj_tag.find("metadataWrap").find("xmlData").find(
                    "{http://www.esa.int/safe/sentinel-1.0}acquisitionPeriod")
                scn_metadata_info['start_time'] = acq_period_tag.find("{http://www.esa.int/safe/sentinel-1.0}startTime").text.strip()
                scn_metadata_info['stop_time'] = acq_period_tag.find("{http://www.esa.int/safe/sentinel-1.0}stopTime").text.strip()
            elif metadata_obj_tag.attrib['ID'] == 'measurementFrameSet':
                foot_print_tag = metadata_obj_tag.find("metadataWrap").find("xmlData").find(
                    "{http://www.esa.int/safe/sentinel-1.0}frameSet").find(
                    "{http://www.esa.int/safe/sentinel-1.0}frame").find(
                    "{http://www.esa.int/safe/sentinel-1.0}footPrint")
                coords_tag = foot_print_tag.find("{http://www.opengis.net/gml}coordinates")
                coords_str = coords_tag.text.strip()
                coords_str_arr = coords_str.split(" ")

                scn_metadata_info['footprint'] = []
                lat_pt_sum = 0.0
                lon_pt_sum = 0.0
                n = 0.0
                for coord in coords_str_arr:
                    pt_coords = coord.split(",")
                    lat_pt = float(pt_coords[0])
                    lon_pt = float(pt_coords[1])
                    scn_metadata_info['footprint'].append([lat_pt, lon_pt])
                    lat_pt_sum += lat_pt
                    lon_pt_sum += lon_pt
                    n += 1.0
                scn_metadata_info['centre_pt'] = [(lat_pt_sum/n), (lon_pt_sum/n)]
    return scn_metadata_info


def create_sentinel1_basename(scn_metadata_info):
    """
    The function uses the metadata from manifest.safe to create a base name for the image.

    :param scn_metadata_info: metadata dict from manifest.safe; use function retrieve_sentinel1_metadata.
    :return: a string base name for outputs related to the input SAFE file.

    """
    base_name = 'sen1' + scn_metadata_info['sensor_number'].lower()

    start_time_obj = datetime.datetime.strptime(scn_metadata_info['start_time'], "%Y-%m-%dT%H:%M:%S.%f")
    date_str = start_time_obj.strftime("%Y%m%d")
    base_name = base_name + '_' + date_str

    lat_pt = scn_metadata_info['centre_pt'][0]
    lon_pt = scn_metadata_info['centre_pt'][1]
    east_west = 'e'
    if lon_pt < 0:
        east_west = 'w'
    north_south = 'n'
    if lat_pt < 0:
        north_south = 's'
    pos = "lat" + north_south + str(round(lat_pt, 1)).replace('.', '').replace('-', '') + "lon" + east_west + str(
        round(lon_pt, 1)).replace('.', '').replace('-', '')

    base_name = base_name + '_' + pos

    base_name = base_name + '_' + scn_metadata_info['mode'].lower()
    if scn_metadata_info['pass'].upper() == 'ASCENDING':
        base_name = base_name + '_asc'
    else:
        base_name = base_name + '_dsc'

    base_name = base_name + '_orb' + scn_metadata_info['orbit_number'] + '_cyc' + scn_metadata_info[
        'cycle_number'] + '_slc' + scn_metadata_info['slice_number']

    return base_name


def find_sen1_ard_files(input_safe_file):
    """
    Find the files which are needed for ARD processing of data from within the SAFE file.

    Looking for:
    1. measurement TIFF files
    2. annotation XML files
    3. calibration XML files
    4. noise XML files

    :param input_safe_file:
    :return: dict of files

    """
    safe_files = dict()

    safe_files['measure_vv'] = findFile(os.path.join(input_safe_file, 'measurement'), '*vv*.tiff', raise_exp=False)
    safe_files['measure_vh'] = findFile(os.path.join(input_safe_file, 'measurement'), '*vh*.tiff', raise_exp=False)

    safe_files['annotation_vv'] = findFile(os.path.join(input_safe_file, 'annotation'), '*vv*.xml', raise_exp=False)
    safe_files['annotation_vh'] = findFile(os.path.join(input_safe_file, 'annotation'), '*vh*.xml', raise_exp=False)

    safe_files['calibration_vv'] = findFile(os.path.join(input_safe_file, 'annotation', 'calibration'), 'calibration*vv*.xml', raise_exp=False)
    safe_files['calibration_vh'] = findFile(os.path.join(input_safe_file, 'annotation', 'calibration'), 'calibration*vh*.xml', raise_exp=False)

    safe_files['noise_vv'] = findFile(os.path.join(input_safe_file, 'annotation', 'calibration'), 'noise*vv*.xml', raise_exp=False)
    safe_files['noise_vh'] = findFile(os.path.join(input_safe_file, 'annotation', 'calibration'), 'noise*vh*.xml', raise_exp=False)

    return safe_files


def get_file_extension(gdalformat):
    """
    A function to get the extension for a given file format
    (NOTE, currently only KEA, GTIFF, HFA, PCI and ENVI are supported).

    :return: string

    """
    ext = ".NA"
    if gdalformat.lower() == "kea":
        ext = ".kea"
    elif gdalformat.lower() == "gtiff":
        ext = ".tif"
    elif gdalformat.lower() == "hfa":
        ext = ".img"
    elif gdalformat.lower() == "envi":
        ext = ".env"
    elif gdalformat.lower() == "pcidsk":
        ext = ".pix"
    else:
        raise Exception("The extension for the gdalformat specified is unknown.")
    return ext


def set_band_names(input_image, band_names, feedback=False):
    """
    A utility function to set band names.

    :param input_image: is the input image
    :param band_names: is a list of band names
    :param feedback: is a boolean specifying whether feedback will be printed to the console
                     (True: Printed / False (default) Not Printed)

    """
    dataset = gdal.Open(input_image, gdal.GA_Update)
    if dataset is None:
        raise Exception("Could not open input image: {}".format(input_image))
    for i in range(len(band_names)):
        band = i + 1
        band_name = band_names[i]
        img_band = dataset.GetRasterBand(band)
        # Check the image band is available
        if not img_band is None:
            if feedback:
                print('Setting Band {0} to "{1}"'.format(band, band_name))
            img_band.SetDescription(band_name)
        else:
            raise Exception("Could not open the image band: ", band)


def calc_ratio_img(vv_img, vh_img, out_img, gdal_format):
    """
    A function which calculates the ratio image of the VV/HV polarisations.

    :param vv_img: GDAL image with VV polarisation
    :param vh_img: GDAL image with VH polarisation
    :param out_img: Output image file.
    :param gdal_format: Output image file format.

    """
    vv_img_ds = gdal.Open(vv_img)
    if vv_img_ds is None:
        raise Exception("Could not open image: {}".format(vv_img))
    vv_img_band = vv_img_ds.GetRasterBand(1)
    if vv_img_band is None:
        raise Exception("Could not open image band {}".format(vv_img))
    vv_val_arr = vv_img_band.ReadAsArray()
    # Get Header information.
    vv_geotransform = vv_img_ds.GetGeoTransform()
    vv_x_pxls = vv_img_ds.RasterXSize
    vv_y_pxls = vv_img_ds.RasterYSize
    vv_proj_str = vv_img_ds.GetProjection()
    vv_img_ds = None

    vh_img_ds = gdal.Open(vh_img)
    if vh_img_ds is None:
        raise Exception("Could not open image: {}".format(vh_img))
    vh_img_band = vh_img_ds.GetRasterBand(1)
    if vh_img_band is None:
        raise Exception("Could not open image band {}".format(vh_img))
    vh_val_arr = vh_img_band.ReadAsArray()
    # Get Header information.
    vh_geotransform = vh_img_ds.GetGeoTransform()
    vh_x_pxls = vh_img_ds.RasterXSize
    vh_y_pxls = vh_img_ds.RasterYSize
    vh_img_ds = None

    if (vv_x_pxls != vh_x_pxls) and (vv_y_pxls != vh_y_pxls):
        raise Exception("The VV and HV images have different numbers of pixels.")

    if (vv_geotransform[0] != vh_geotransform[0]) and (vv_geotransform[3] != vh_geotransform[3]):
        raise Exception("The Top-Left corner coordinate is not the same for the VV and HV images.")

    if (vv_geotransform[1] != vh_geotransform[1]) and (vv_geotransform[5] != vh_geotransform[5]):
        raise Exception("The image pixel resolution is not the same for the VV and HV images.")

    if (vv_geotransform[2] != vh_geotransform[2]) and (vv_geotransform[4] != vh_geotransform[4]):
        raise Exception("The image rotation is not the same for the VV and HV images.")

    ratio_img_arr = numpy.where(
        (numpy.isfinite(vv_val_arr) & (vv_val_arr > 0.0) & numpy.isfinite(vh_val_arr) & (vh_val_arr > 0.0)),
        vv_val_arr / vh_val_arr, 0.0)

    ratio_img_arr[numpy.isnan(ratio_img_arr)] = 0.0
    ratio_img_arr[numpy.isinf(ratio_img_arr)] = 0.0

    out_ratios_file_ds = gdal.GetDriverByName(gdal_format).Create(out_img, vv_x_pxls, vv_y_pxls, 1, gdal.GDT_Float32)
    if out_ratios_file_ds == None:
        raise Exception('Could not create ratio image output raster: \'' + out_img + '\'')
    out_ratios_file_ds.SetGeoTransform(vv_geotransform)
    out_ratios_file_ds.SetProjection(vv_proj_str)
    out_ratio_band = out_ratios_file_ds.GetRasterBand(1)
    if out_ratio_band == None:
        raise Exception('Could not open ratio image band: \'' + out_img + '\'')
    out_ratio_band.SetNoDataValue(0.0)
    out_ratio_band.WriteArray(ratio_img_arr)
    logger.debug("Created ratio output image file: {}".format(out_img))


def convert_to_dB(input_img, output_img, gdal_format):
    """
    Convert power image to decibels (dB) by applying 10 x log10(pwr)

    :param input_img: Input power image
    :param output_img: Output dB image
    :param gdal_format: GDAL image format for output image

    """
    img_ds = gdal.Open(input_img)
    if img_ds is None:
        raise Exception("Could not open image: {}".format(input_img))

    # Get Header information.
    geotransform = img_ds.GetGeoTransform()
    x_pxls = img_ds.RasterXSize
    y_pxls = img_ds.RasterYSize
    proj_str = img_ds.GetProjection()
    n_bands = img_ds.RasterCount

    co = []
    if gdal_format == 'GTIFF':
        co = ["TILED=YES", "COMPRESS=LZW", "BIGTIFF=IF_NEEDED"]

    out_img_ds = gdal.GetDriverByName(gdal_format).Create(output_img, x_pxls, y_pxls, n_bands, gdal.GDT_Float32,
                                                          options=co)
    if out_img_ds == None:
        raise Exception("Could not create dB image output raster: '{}'.".format(output_img))
    out_img_ds.SetGeoTransform(geotransform)
    out_img_ds.SetProjection(proj_str)

    for band in range(n_bands):
        band = band + 1
        img_band = img_ds.GetRasterBand(band)
        if img_band is None:
            raise Exception("Could not open image band {} from {}".format(band, input_img))
        val_arr = img_band.ReadAsArray()

        if band == 1:
            msk_arr = numpy.where(((val_arr > 0.0) & numpy.isfinite(val_arr)), 1, 0)
        else:
            msk_arr = numpy.where(((val_arr > 0.0) & numpy.isfinite(val_arr) & (msk_arr == 1)), 1, 0)

    for band in range(n_bands):
        band = band + 1
        img_band = img_ds.GetRasterBand(band)
        if img_band is None:
            raise Exception("Could not open image band {} from {}".format(band, input_img))
        val_arr = img_band.ReadAsArray()

        dB_img_arr = numpy.where(((val_arr > 0.0) & numpy.isfinite(val_arr) & (msk_arr == 1)),
                                 10 * numpy.log10(val_arr), 999)
        dB_img_arr[numpy.isnan(dB_img_arr)] = 999
        dB_img_arr[numpy.isinf(dB_img_arr)] = 999

        out_img_band = out_img_ds.GetRasterBand(band)
        if out_img_band == None:
            raise Exception("Could not open image band {} from {}".format(band, output_img))
        out_img_band.SetNoDataValue(999)
        out_img_band.WriteArray(dB_img_arr)

    img_ds = None
    out_img_ds = None


def write_list_to_file(data_lst, out_file):
    """
    Write a list a text file, one line per item.

    :param dataList: List of items to be written
    :param outFile:

    """
    logger.debug("Creating output file: {}".format(out_file))
    f = open(out_file, 'w')
    logger.debug("Created output file: {}".format(out_file))
    for item in data_lst:
        f.write(str(item) + '\n')
    f.flush()
    f.close()
    logger.debug("Finished writing and close output file: {}".format(out_file))
