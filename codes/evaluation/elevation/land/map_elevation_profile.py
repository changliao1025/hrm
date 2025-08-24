import os, sys, stat
import glob
import numpy as np
from pathlib import Path
from os.path import realpath
from osgeo import osr, gdal, ogr
import elevation
from tinyr import RTree
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.visual.plot_xy_data import plot_xy_data
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.visual.map.raster.map_raster_file import map_raster_file
sPath = '/qfs/people/liao313/workspace/python/subset/'
sys.path.append(sPath)

sWorkspace_output = '/compyfs/liao313/04model/subset/mississippi/elevation/land'
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)
    os.chmod(sWorkspace_output, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

nElevation_profile = 11
sFilename_mesh_geojson = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250201003/pyflowline/mpas.geojson'

sFilename = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501002/pyflowline/00000001/flowline_edge.geojson'

aFlowline, pProjection_geojson = read_flowline_geojson(sFilename)

nEdge = len(aFlowline)

#find all the subnfolder under the workspace using glob function

#Find all subfolders (directories) under sWorkspace_output
aFolders = [f for f in glob.glob(os.path.join(sWorkspace_output, '**/'), recursive=True) if os.path.abspath(f) != os.path.abspath(sWorkspace_output)]

for sFolder in aFolders:
    #extract the folder name
    sCellID = os.path.basename(os.path.normpath(sFolder))
    sFilename_dem_box = os.path.join(sFolder, 'dem_box.tif')
    sFilename_dem_clip = os.path.join(sFolder, 'dem_clip.tif')
    aExtent0 = gdal_get_raster_extent(sFilename_dem_box)
    aExtent1 = gdal_get_raster_extent(sFilename_dem_clip)
    #get the largest extent encoompass both extent
    #to have an apple to apple comparison, we need to use the same extent for both rasters(mpas and latlon)
    aExtent = [min(aExtent0[0], aExtent1[0]), max(aExtent0[1], aExtent1[1]),
               min(aExtent0[2], aExtent1[2]), max(aExtent0[3], aExtent1[3])]



    #read the dem to get max and min elevation
    dum = gdal_read_geotiff_file(sFilename_dem_box)
    data_box_dummy = dum['dataOut']
    nan_value = dum['missingValue']
    aData_box_index = np.where(data_box_dummy != nan_value)
    aData_box = data_box_dummy[aData_box_index]
    dMin0 = np.min(aData_box)
    dMax0 = np.max(aData_box)

    dum = gdal_read_geotiff_file(sFilename_dem_clip)
    data_clip_dummy = dum['dataOut']
    nan_value = dum['missingValue']
    aData_clip_index = np.where(data_clip_dummy != nan_value)
    aData_clip = data_clip_dummy[aData_clip_index]
    dMin1 = np.min(aData_clip)
    dMax1 = np.max(aData_clip)
    dMin = min(dMin0, dMin1)
    dMax = max(dMax0, dMax1)

    sColormap_in = 'terrain'

    sFilename_output_in = os.path.join(sFolder, 'dem_clip.png')
    map_raster_file(sFilename_dem_clip, sFilename_output_in=sFilename_output_in,
                    iFlag_colorbar_in= 1, sColormap_in = sColormap_in, iFlag_zebra_in= 1,
                    dData_max_in= dMax, dData_min_in= dMin, sUnit_in= 'Unit: m', aExtent_in=aExtent)



    sFilename_output_in = os.path.join(sFolder, 'dem_box.png')
    map_raster_file(sFilename_dem_box, sFilename_output_in=sFilename_output_in,
                    iFlag_colorbar_in= 1, sColormap_in = sColormap_in, iFlag_zebra_in= 1,
                    dData_max_in= dMax, dData_min_in= dMin, sUnit_in= 'Unit: m', aExtent_in=aExtent)







