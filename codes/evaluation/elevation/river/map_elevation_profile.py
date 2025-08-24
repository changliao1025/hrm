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
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
sPath_hrm = '/qfs/people/liao313/workspace/python/hrm/'
sys.path.append(sPath_hrm)
from codes.shared.map_raster_file import map_raster_file


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
    #flush the output
    sys.stdout.flush()
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

            aExtent0 = gdal_get_raster_extent(sFilename_dem_box)
            aExtent1 = gdal_get_raster_extent(sFilename_dem_clip)

            sFilename_box = os.path.join(sFolder, 'box.geojson')
            sFilename_mesh = os.path.join(sFolder, 'mesh_cell.geojson')
            sFilename_box_cross_section = os.path.join(sFolder, 'dem_box_cross_section.geojson')
            sFilename_clip_cross_section = os.path.join(sFolder, 'dem_clip_cross_section.geojson')

            if os.path.exists(sFilename_box_cross_section) and os.path.exists(sFilename_clip_cross_section):
                aFilename_vector_box = [sFilename_box, sFilename_mesh, sFilename_box_cross_section]
                aFilename_vector_clip = [sFilename_box, sFilename_mesh, sFilename_clip_cross_section]
                aColor_vector = ['red', 'blue', 'black']
            else:
                aFilename_vector_box = [sFilename_box, sFilename_mesh]
                aFilename_vector_clip = [sFilename_box, sFilename_mesh]
                aColor_vector = ['red', 'blue']

            aExtent2 = gdal_get_vector_extent(sFilename_box)
            aExtent3 = gdal_get_vector_extent(sFilename_mesh)
            #get the largest extent encoompass all 4 extents
            aExtent = [ min(aExtent0[0], aExtent1[0], aExtent2[0], aExtent3[0]),
                        max(aExtent0[1], aExtent1[1], aExtent2[1], aExtent3[1]),
                        min(aExtent0[2], aExtent1[2], aExtent2[2], aExtent3[2]),
                        max(aExtent0[3], aExtent1[3], aExtent2[3], aExtent3[3])]

            #aExtent = [min(aExtent0[0], aExtent1[0]), max(aExtent0[1], aExtent1[1]),
            #           min(aExtent0[2], aExtent1[2]), max(aExtent0[3], aExtent1[3])]

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
            aLabel_legend = ['(b)', 'Unstructured MPAS mesh-based']

            sFilename_output_in = os.path.join(sFolder, 'dem_clip.png')
            map_raster_file(sFilename_dem_clip, sFilename_output_in=sFilename_output_in,
                            iFlag_colorbar_in= 1, sColormap_in = sColormap_in, iFlag_zebra_in= 0,
                            dData_max_in= dMax, dData_min_in= dMin, sUnit_in= 'Unit: m', aExtent_in=aExtent,
                            aFilename_vector_in=aFilename_vector_clip,
                            aLabel_legend_in = aLabel_legend,
                            aColor_vector_in=aColor_vector)

            aLabel_legend = ['(a)', 'Structured Latlon mesh-based']
            sFilename_output_in = os.path.join(sFolder, 'dem_box.png')
            map_raster_file(sFilename_dem_box, sFilename_output_in=sFilename_output_in,
                            iFlag_colorbar_in= 1, sColormap_in = sColormap_in, iFlag_zebra_in= 0,
                            dData_max_in= dMax, dData_min_in= dMin, sUnit_in= 'Unit: m', aExtent_in=aExtent,
                            aFilename_vector_in=aFilename_vector_box,
                            aLabel_legend_in = aLabel_legend,
                            aColor_vector_in=aColor_vector)







