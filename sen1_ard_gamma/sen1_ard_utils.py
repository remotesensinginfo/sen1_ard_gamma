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
import os.path
import glob
import datetime

import numpy

import osgeo.ogr as ogr

logger = logging.getLogger(__name__)


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
    :return: x, y. (note if returning long, lat you might need to invert)
    """
    wktPt = 'POINT(%s %s)' % (x, y)
    point = ogr.CreateGeometryFromWkt(wktPt)
    point.AssignSpatialReference(in_proj_osr_obj)
    point.TransformTo(out_proj_osr_obj)
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

    base_name = base_name + '_' + date_str + '_' + scn_metadata_info['mode'].lower()
    if scn_metadata_info['pass'].upper() == 'ASCENDING':
        base_name = base_name + '_asc'
    else:
        base_name = base_name + '_dsc'

    base_name = base_name + '_orb' + scn_metadata_info['orbit_number'] + '_cyc' + scn_metadata_info[
        'cycle_number'] + '_mdt' + scn_metadata_info['mission_data_take_id']

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

    safe_files['calibration_vv'] = findFile(os.path.join(input_safe_file, 'calibration'), 'calibration*vv*.xml', raise_exp=False)
    safe_files['calibration_vh'] = findFile(os.path.join(input_safe_file, 'calibration'), 'calibration*vh*.xml', raise_exp=False)

    safe_files['noise_vv'] = findFile(os.path.join(input_safe_file, 'calibration'), 'noise*vv*.xml', raise_exp=False)
    safe_files['noise_vh'] = findFile(os.path.join(input_safe_file, 'calibration'), 'noise*vh*.xml', raise_exp=False)

    return safe_files
