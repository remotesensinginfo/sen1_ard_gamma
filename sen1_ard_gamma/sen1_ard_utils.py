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
