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

# Find all subfolders (directories) under sWorkspace_output
aFolders = [f for f in glob.glob(os.path.join(sWorkspace_output, '**/'), recursive=True) if os.path.abspath(f) != os.path.abspath(sWorkspace_output)]
aElevation_profile_clip = np.zeros(nElevation_profile)
aElevation_profile_box = np.zeros(nElevation_profile)
for sFolder in aFolders:
    #extract the folder name
    sCellID = os.path.basename(os.path.normpath(sFolder))
    sFilename_dem_clip = os.path.join(sFolder, 'dem_clip.tif')
    #read the dem to get max and min elevation
    dum = gdal_read_geotiff_file(sFilename_dem_clip)
    data = dum['dataOut']
    nan_value = dum['missingValue']
    dMin0 = np.min(data[data != nan_value])
    dMax0 = np.max(data[data != nan_value])
    aData_clip = data[data != nan_value]
    aElevation_profile_clip[0] = np.min(aData_clip)
    aElevation_profile_clip[nElevation_profile-1] = np.max(aData_clip)
    for i in range(1, nElevation_profile-1):
        aElevation_profile_clip[i] = np.percentile(aData_clip, i*10)

    sFilename_dem_box = os.path.join(sFolder, 'dem_box.tif')
    dum = gdal_read_geotiff_file(sFilename_dem_box)
    data = dum['dataOut']
    nan_value = dum['missingValue']
    dMin1 = np.min(data[data != nan_value])
    dMax1 = np.max(data[data != nan_value])
    #update min and max using dMin0 and dMin1
    dMin = np.min([dMin0, dMin1])
    dMax = np.max([dMax0, dMax1])

    aData_box = data[data != nan_value]
    aElevation_profile_box[0] = np.min(aData_box)
    aElevation_profile_box[nElevation_profile-1] = np.max(aData_box)
    for i in range(1, nElevation_profile-1):
        aElevation_profile_box[i] = np.percentile(aData_box, i*10)

    x = np.linspace(0, 100, nElevation_profile)
    #create a list of x values
    aX_all = [x, x]
    #create a list of y values
    aY_all = [aElevation_profile_box, aElevation_profile_clip]
    sFilename_out = os.path.join(sFolder, 'elevation_profile.png')
    plot_xy_data(aX_all,
         aY_all,
         sFilename_out, dMin_x_in=0, dMax_x_in=100,
         dMin_y_in=dMin, dMax_y_in=dMax, dSpace_x_in=10,
         aColor_in=['red', 'blue'], aLinestyle_in=['-', '-'],
         aLabel_legend_in=['Structured mesh-based', 'Unstructured mesh-based'],
         sLabel_x_in='Area fraction (%)', sLabel_y_in='Elevation (m)')







