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
# Purpose:  This file allows image stats and pyramids to be calculated and
#           has been written using the implementation in the RIOS library as
#           a guide (https://bitbucket.org/chchrsc/rios/src/default/rios/calcstats.py)
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 06/06/2019
# Version: 1.0
#
# History:
# Version 1.0 - Created.

import numpy
from osgeo import gdal
gdal.UseExceptions()

import logging
logger = logging.getLogger(__name__)

def add_image_pyramids(img_ds):
    """
    A function which uses gdal.Dataset.BuildOverviews() to calculate image pyramids.

    :param ds: gdal dataset object for the image to be processed.

    """
    minoverviewdim = 33
    levels = [4, 8, 16, 32, 64, 128, 256, 512]

    if img_ds.RasterXSize < img_ds.RasterYSize:
        mindim = img_ds.RasterXSize
    else:
        mindim = img_ds.RasterYSize

    nOverviews = 0
    for i in levels:
        if (mindim // i) > minoverviewdim:
            nOverviews = nOverviews + 1

    img_ds.BuildOverviews("AVERAGE", levels[:nOverviews])


def find_or_create_col(rat_obj, usage, name, dtype):
    """
    A function which returns index of an existing column with the given usage.
    If a column does not exist then one is created using the supplied name and dtype.
    :param rat_obj: An object for the GDAL band RAT.
    :param usage: Specifies the purpose of the column (e.g., gdal.GFU_PixelCount)
    :param name: Specifies the name of a new column if created (e.g., 'Histogram')
    :param dtype: Specifies the datatype of a new column if created (e.g., gdal.GFT_Real)
    :return: tuple with index and a boolean specifying if it is a new column or not

    """
    ncols = rat_obj.GetColumnCount()
    for col in range(ncols):
        if rat_obj.GetUsageOfCol(col) == usage:
            return col, False

    # got here so can't exist
    rat_obj.CreateColumn(name, dtype, usage)
    # new one will be last col
    return ncols, True


def calc_img_band_stats(img_ds, no_data_val=None):
    """
    Calculates the image band stats and sorts them in image band metadata.
    The GDAL function gdal.Band.ComputeStatistics() is used to calculate the
    mean, stddev, min and max and gdal.Band.GetHistogram() the histogram from
    which the median and mode are derived.

    :param img_ds:
    :param no_data_val:
    :return:
    """

    # flush the cache. The ensures that any unwritten data is
    # written to file so we get the right stats. It also
    # makes sure any metdata is written on HFA. This means
    # the LAYER_TYPE setting will be picked up by rat.SetLinearBinning()
    img_ds.FlushCache()

    gdalLargeIntTypes = set([gdal.GDT_Int16, gdal.GDT_UInt16, gdal.GDT_Int32, gdal.GDT_UInt32])
    gdalFloatTypes = set([gdal.GDT_Float32, gdal.GDT_Float64])

    for n_band in range(img_ds.RasterCount):
        band = img_ds.GetRasterBand(n_band + 1)

        # fill in the metadata
        tmp_meta = band.GetMetadata()

        if no_data_val is not None:
            band.SetNoDataValue(no_data_val)
            tmp_meta["STATISTICS_EXCLUDEDVALUES"] = repr(no_data_val)  # doesn't seem to do anything

        try:
            (min_val, max_val, mean_val, stddev_val) = band.ComputeStatistics(False)
        except RuntimeError as e:
            if str(e).endswith('Failed to compute statistics, no valid pixels found in sampling.'):
                min_val = no_data_val
                max_val = no_data_val
                mean_val = no_data_val
                stddev_val = 0
            else:
                raise e

        tmp_meta["STATISTICS_MINIMUM"] = repr(min_val)
        tmp_meta["STATISTICS_MAXIMUM"] = repr(max_val)
        tmp_meta["STATISTICS_MEAN"] = repr(mean_val)
        tmp_meta["STATISTICS_STDDEV"] = repr(stddev_val)
        # because we did at full res - these are the default anyway
        tmp_meta["STATISTICS_SKIPFACTORX"] = "1"
        tmp_meta["STATISTICS_SKIPFACTORY"] = "1"

        # create a histogram so we can do the mode and median
        if band.DataType == gdal.GDT_Byte:
            # if byte data use 256 bins and the whole range
            hist_min = 0
            hist_max = 255
            hist_step = 1.0
            hist_calc_min = -0.5
            hist_calc_max = 255.5
            hist_n_bins = 256
            tmp_meta["STATISTICS_HISTOBINFUNCTION"] = 'direct'
        elif band.DataType in gdalLargeIntTypes:
            hist_range = int(numpy.ceil(max_val) - numpy.floor(min_val)) + 1
            hist_min = min_val
            hist_max = max_val
            if hist_range <= 256:
                hist_n_bins = hist_range
                hist_step = 1.0
                tmp_meta["STATISTICS_HISTOBINFUNCTION"] = 'direct'
                hist_calc_min = hist_min - 0.5
                hist_calc_max = hist_max + 0.5
            else:
                hist_n_bins = 256
                tmp_meta["STATISTICS_HISTOBINFUNCTION"] = 'linear'
                hist_calc_min = hist_min
                hist_calc_max = hist_max
                hist_step = float(hist_calc_max - hist_calc_min) / hist_n_bins
        elif band.DataType in gdalFloatTypes:
            hist_n_bins = 256
            hist_min = min_val
            hist_max = max_val
            tmp_meta["STATISTICS_HISTOBINFUNCTION"] = 'linear'
            hist_calc_min = hist_min
            hist_calc_max = hist_max
            hist_step = float(hist_calc_max - hist_calc_min) / hist_n_bins

        # get histogram and force GDAL to recalculate it
        hist = band.GetHistogram(hist_calc_min, hist_calc_max, hist_n_bins, False, False)

        # Check if GDAL's histogram code overflowed. This is not a fool-proof test,
        # as some overflows will not result in negative counts.
        histogram_overflow = (min(hist) < 0)

        # we may use this ratObj reference for the colours below also
        # may be None if format does not support RATs
        rat_obj = band.GetDefaultRAT()

        if not histogram_overflow:
            # comes back as a list for some reason
            hist = numpy.array(hist)

            # Note that we have explicitly set histstep in each datatype case
            # above. In principle, this can be calculated, as it is done in the
            # float case, but for some of the others we need it to be exactly
            # equal to 1, so we set it explicitly there, to avoid rounding
            # error problems.

            # do the mode - bin with the highest count
            mode_bin = numpy.argmax(hist)
            mode_val = mode_bin * hist_step + hist_min
            if band.DataType == gdal.GDT_Float32 or band.DataType == gdal.GDT_Float64:
                tmp_meta["STATISTICS_MODE"] = repr(mode_val)
            else:
                tmp_meta["STATISTICS_MODE"] = repr(int(round(mode_val)))

            if rat_obj is not None:
                hist_indx, hist_new = find_or_create_col(rat_obj, gdal.GFU_PixelCount, "Histogram", gdal.GFT_Real)
                # write the hist in a single go
                rat_obj.SetRowCount(hist_n_bins)
                rat_obj.WriteArray(hist, hist_indx)

                rat_obj.SetLinearBinning(hist_min, (hist_calc_max - hist_calc_min) / hist_n_bins)

                # The HFA driver still honours the STATISTICS_HISTOBINVALUES
                # metadata item. If we are recalculating the histogram the old
                # values will be copied across with the metadata so clobber it
                if "STATISTICS_HISTOBINVALUES" in tmp_meta:
                    del tmp_meta["STATISTICS_HISTOBINVALUES"]
            else:
                # old method
                tmp_meta["STATISTICS_HISTOBINVALUES"] = '|'.join(map(repr, hist)) + '|'

                tmp_meta["STATISTICS_HISTOMIN"] = repr(hist_min)
                tmp_meta["STATISTICS_HISTOMAX"] = repr(hist_max)
                tmp_meta["STATISTICS_HISTONUMBINS"] = repr(hist_n_bins)

            # estimate the median - bin with the middle number
            middle_num = hist.sum() / 2
            gt_middle = hist.cumsum() >= middle_num
            median_bin = gt_middle.nonzero()[0][0]
            median_val = median_bin * hist_step + hist_min
            if band.DataType == gdal.GDT_Float32 or band.DataType == gdal.GDT_Float64:
                tmp_meta["STATISTICS_MEDIAN"] = repr(median_val)
            else:
                tmp_meta["STATISTICS_MEDIAN"] = repr(int(round(median_val)))

        # set the data
        band.SetMetadata(tmp_meta)

        if (rat_obj is not None) and (not rat_obj.ChangesAreWrittenToFile()):
            # For drivers that require the in memory thing
            band.SetDefaultRAT(rat_obj)


def run_calc_img_stats_pyramids(input_img, no_data_val=None):
    """
    A function which runs functions to calculate the image statistics and pyramids using GDAL.

    Note. these functions are only for continuous data. If categorical data then see either rsgislib or rios.

    :param input_img: file path to GDAL compatible image file.
    :param no_data_val: the no data value which should be ignored. Can be None; Default.

    """
    logger.info("Computing image statistics and pyramids for {}".format(input_img))
    img_ds = gdal.Open(input_img, gdal.GA_Update)
    if img_ds is None:
        raise Exception("Could not open image: {}".format(input_img))

    calc_img_band_stats(img_ds, no_data_val)

    add_image_pyramids(img_ds)

    img_ds.FlushCache()
    img_ds = None
    logger.debug("Finished computing image statistics and pyramids for {}".format(input_img))
