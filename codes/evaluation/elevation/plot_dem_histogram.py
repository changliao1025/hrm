import os, sys
from osgeo import ogr, osr, gdal
import numpy as np
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.visual.map.vector.map_vector_polygon_file import map_vector_polygon_file
from pyearth.visual.histogram.histogram_plot import histogram_plot

sWorkspace = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation'

sFilename_mpas  = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation/dem_subset_mpas.geojson'
sFilename_drt_8th = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation/dem_subset_drt_8th.geojson'
sFilename_drt = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation/dem_subset_drt.geojson'

#read all the mean elevation data

pDriver = ogr.GetDriverByName('GeoJSON')

pDataset_mpas = pDriver.Open(sFilename_mpas, 0)
pLayer_mpas = pDataset_mpas.GetLayer()
pLayer_mpas.ResetReading()
aFeature_mpas = pLayer_mpas.GetNextFeature()
aElevation_mpas = list()
while aFeature_mpas:
    aElevation_mpas.append(aFeature_mpas.GetField('mean'))
    aFeature_mpas = pLayer_mpas.GetNextFeature()
pLayer_mpas.ResetReading()
pDataset_mpas = None

#read drt data
pDataset_drt = pDriver.Open(sFilename_drt_8th, 0)
pLayer_drt = pDataset_drt.GetLayer()
pLayer_drt.ResetReading()
aFeature_drt = pLayer_drt.GetNextFeature()
aElevation_drt_8th = list()
while aFeature_drt:
    aElevation_drt_8th.append(aFeature_drt.GetField('mean'))
    aFeature_drt = pLayer_drt.GetNextFeature()
pLayer_drt.ResetReading()
pDataset_drt = None

#read drt data
pDataset_drt = pDriver.Open(sFilename_drt, 0)
pLayer_drt = pDataset_drt.GetLayer()
pLayer_drt.ResetReading()
aFeature_drt = pLayer_drt.GetNextFeature()
aElevation_drt = list()
while aFeature_drt:
    aElevation_drt.append(aFeature_drt.GetField('mean'))
    aFeature_drt = pLayer_drt.GetNextFeature()
pLayer_drt.ResetReading()
pDataset_drt = None

#convert to numpy array
aElevation_mpas = np.array(aElevation_mpas)
aElevation_drt_8th = np.array(aElevation_drt_8th)
aElevation_drt = np.array(aElevation_drt)

#call pyearth for histogram plot
sFilename_output_in = os.path.join(sWorkspace, 'elevation_histogram.png')
aData_all_in = [aElevation_drt_8th, aElevation_drt, aElevation_mpas]

aLabel_legend = list()
aLabel_legend.append('Structured Latlon mesh-based (1/8th)')
aLabel_legend.append('Structured Latlon mesh-based (1/16th)')
aLabel_legend.append('Unstructured MPAS mesh-based')
aOrder_in = [ 2, 0, 1 ]
histogram_plot(aData_all_in,
                   sFilename_output_in = sFilename_output_in,
                   iSize_x_in=12,
                   iSize_y_in=4,
                   ncolumn_in=None,
                   iFlag_scientific_notation_in=None,
                   iFlag_normalize_in=  1,
                   iFlag_log_in=None,
                   sHisttype_in='barstacked',
                   aColor_in=['red', 'orange', 'blue'],
                   aPlot_order_in=aOrder_in,
                   iDPI_in=None,
                   dMin_x_in=0,
                   dMax_x_in=3000,
                   dSpace_x_in=50,
                   sLabel_x_in='Surface elevation (m)',
                   sLabel_y_in=None,
                   sFont_in=None,
                   aLocation_legend_in=None,
                   sLocation_legend_in=None,
                   sTitle_in='Histogram of elevation (mean)',
                   aLabel_legend_in=aLabel_legend)