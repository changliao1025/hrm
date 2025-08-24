import sys
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

sFilename_raster_in = '/compyfs/liao313/00raw/hydrology/hydroshed/hydroshed/dem/na_dem_3s.tif'
#sFilename_geojson_in = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501004/pyflowline/mpas.geojson'
sFilename_geojson_in = '/compyfs/liao313/04model/pyflowline/conus/pyflowline20241201019/mpas.geojson'
#example 1, call without saving the output
#time the process


#example 2, call and save the output in folder
sFolder_raster_out = '//compyfs/liao313/04model/subset/yukon'
tStart = time.time()
sFilename_vector_out = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation/dem_subset_mpas2.shp'
vData = subset_raster(sFilename_raster_in, sFilename_geojson_in, iFlag_save_clipped_raster_in = None,
                  iFlag_save_in_vector_in = 1,
                  iFlag_stat_in = 1,
                  sFilename_vector_out= sFilename_vector_out)
tEnd = time.time()
#print('Elapsed time is: ', tEnd - tStart)

#convert shp to geojson
sFilename_geojson_out = sFilename_vector_out.replace('.shp', '.geojson')
convert_vector_to_geojson(sFilename_vector_out, sFilename_geojson_out)


#calculate the elevation profile using result from example 1
nElevation_profile = 11
tStart = time.time()
print('Number of cell:', len(vData))
#for pData in vData:
#    #create a 11 element numpy array to store the result
#    aElevation_profile = np.zeros(nElevation_profile)
#    #remove the missing value
#    aData = pData[np.where(pData != -9999)]
#    #call the numpy percentile function
#    aElevation_profile[0] = np.min(aData)
#    aElevation_profile[nElevation_profile-1] = np.max(aData)
#    for i in range(1, nElevation_profile-1):
#        aElevation_profile[i] = np.percentile(aData, i*10)

    #print(aElevation_profile)

tEnd = time.time()
print('Elapsed time is: ', tEnd - tStart)