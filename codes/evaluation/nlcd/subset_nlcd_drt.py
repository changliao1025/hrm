import os, sys
import time
from pathlib import Path
from os.path import realpath
import numpy as np
#this is the pre-processing of the data
#basically the input needs at least (1) boundary, (2) dem
from pyearth.toolbox.conversion.convert_vector_to_geojson import convert_vector_to_geojson


sPath = '/qfs/people/liao313/workspace/python/subset/'
sys.path.append(sPath)
from subset.subset_raster import subset_raster
sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'
sWorkspace_output_latlon = os.path.join(sWorkspace_out, 'latlon')
sFilename_raster_in = '/compyfs/liao313/00raw/nlcd/Annual_NLCD_LndCov_2023_CU_C1V0.tif'


sFilename_geojson_in = os.path.join(sWorkspace_output_latlon, 'latlon_conus_8th.geojson')
#example 1, call without saving the output
#time the process


tStart = time.time()
sFilename_vector_out = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/nlcd/nlcd_subset_drt_8th.shp'
vData = subset_raster(sFilename_raster_in,
                      sFilename_geojson_in,
                      iFlag_save_clipped_raster_in = None,
                  iFlag_save_in_vector_in = 1,
                       iFlag_stat_in = None,
                  sFilename_vector_out= sFilename_vector_out)
tEnd = time.time()
print('Elapsed time is: ', tEnd - tStart)

#convert shp to geojson
sFilename_geojson_out = sFilename_vector_out.replace('.shp', '.geojson')
convert_vector_to_geojson(sFilename_vector_out, sFilename_geojson_out)



tEnd = time.time()
print('Elapsed time is: ', tEnd - tStart)