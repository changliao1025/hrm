import os, sys, stat
import random
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

sWorkspace_output = '/compyfs/liao313/04model/subset/mississippi/elevation/river'
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)
    os.chmod(sWorkspace_output, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

nElevation_profile = 11
sFilename_mesh_geojson = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250201003/pyflowline/mpas.geojson'

sFilename = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501002/pyflowline/00000001/flowline_edge.geojson'

aFlowline, pProjection_geojson = read_flowline_geojson(sFilename)

nEdge = len(aFlowline)




#now we need to find which mesh cell include the start point
index_cell = RTree(max_cap=5, min_cap=2)
pDriver = ogr.GetDriverByName('GeoJSON')
pDataset_mesh = pDriver.Open(sFilename_mesh_geojson, 0)
pSpatialRef = osr.SpatialReference()
pSpatialRef.ImportFromEPSG(4326)

if pDataset_mesh is None:
    print('Failed to open the mesh file')
    sys.exit(1)
pLayer_mesh = pDataset_mesh.GetLayer(0)
nCell = pLayer_mesh.GetFeatureCount()
for i in range(nCell):
    lID = i
    pFeature = pLayer_mesh.GetFeature(lID)
    pGeometry = pFeature.GetGeometryRef()
    #get bounding box of the mesh cell
    pEnvelope = pGeometry.GetEnvelope()
    left = pEnvelope[0]
    bottom = pEnvelope[2]
    right = pEnvelope[1]
    top = pEnvelope[3]
    pBound= (left, bottom, right, top)
    index_cell.insert(lID, pBound)  #

#pick a random edge using the random number generator
print('The number of edges is: ', nEdge)
aElevation_profile_clip = np.zeros(nElevation_profile)
aElevation_profile_box = np.zeros(nElevation_profile)
#for iEdge in range(0, 400):
for iEdge in range(530,  nEdge, 10):
    print('The edge number is: ', iEdge)
    oFlowline = aFlowline[iEdge]
    #get the start and end point of the edge
    pVertex_start = oFlowline.pVertex_start
    pVertex_end = oFlowline.pVertex_end
    #create a point from the start point
    pPoint_start = ogr.Geometry(ogr.wkbPoint)
    pPoint_start.AddPoint(pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree)
    #create a point from the end point
    pPoint_end = ogr.Geometry(ogr.wkbPoint)
    pPoint_end.AddPoint(pVertex_end.dLongitude_degree, pVertex_end.dLatitude_degree)
    pBound_edge = oFlowline.pBound
    #use the bounding box of the start point to find the mesh cell
    aIntersect = list(index_cell.search(pBound_edge))
    nIntersect = len(aIntersect)
    for i in range(nIntersect):
        lID = aIntersect[i]
        pFeature = pLayer_mesh.GetFeature(lID)
        lCellID = pFeature.GetField('cellid')
        sCellID = "{:09d}".format(lCellID)
        if sCellID != '000234146':
            continue
        #check whether the start point is in the mesh cell
        pGeometry = pFeature.GetGeometryRef()
        if pGeometry.Contains(pPoint_start):
            sFolder = os.path.join(sWorkspace_output, sCellID)
            sFilename_dem_box = os.path.join(sFolder, 'dem_box.tif')
            sFilename_dem_clip = os.path.join(sFolder, 'dem_clip.tif')
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
            dMin = 100
            dMax = 200

            aElevation_profile_clip[0] = np.min(aData_clip)
            aElevation_profile_clip[nElevation_profile-1] = np.max(aData_clip)
            for i in range(1, nElevation_profile-1):
                aElevation_profile_clip[i] = np.percentile(aData_clip, i*10)


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
                 aLabel_legend_in=[ 'Structured Latlon mesh-based', 'Unstructured MPAS mesh-based'],
                 aLabel_tag_in = ['(e)'],
                 sLabel_x_in='Area fraction (%)', sLabel_y_in='Elevation (m)')
