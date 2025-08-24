import os, sys, stat
import shutil
import numpy as np
from contextlib import redirect_stdout, redirect_stderr
from osgeo import osr, gdal, ogr
import elevation
from tinyr import RTree
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.gis.geometry.create_box_from_longitude_latitude import create_box_from_longitude_latitude, find_nearest_resolution
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import clip_raster_by_polygon_file
from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
import matplotlib.pyplot as plt
from pyearth.visual.plot_xy_data import plot_xy_data
sPath = '/qfs/people/liao313/workspace/python/subset/'
sys.path.append(sPath)

sPath = '/qfs/people/liao313/workspace/python/hrm/'
sys.path.append(sPath)
from codes.shared.get_cross_section_profile import get_raster_cross_section_profile

sWorkspace_output = '/compyfs/liao313/04model/subset/mississippi/elevation/river'
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)
    os.chmod(sWorkspace_output, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

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

index_edge = RTree(max_cap=5, min_cap=2)
for i in range(nEdge):
    oFlowline = aFlowline[i]
    pVertex_start = oFlowline.pVertex_start
    pVertex_end = oFlowline.pVertex_end
    #get the bounding box of the edge
    pBound_edge = oFlowline.pBound
    #expand the bounding box by 0.001 degree
    pBound_edge = (pBound_edge[0] - 0.001, pBound_edge[1] - 0.001,
                   pBound_edge[2] + 0.001, pBound_edge[3] + 0.001)
    #insert the edge into the index
    index_edge.insert(i, pBound_edge)

#pick a random edge using the random number generator
print('The number of edges is: ', nEdge)


iFlag_download_data = 0
iFlag_cross_section = 1
if iFlag_download_data == 1:
    sFolder_cache = elevation.CACHE_DIR
    if os.path.exists(sFolder_cache):
        shutil.rmtree(sFolder_cache)


for iEdge in range(5450,  nEdge, 10):
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
    dBearing0 = oFlowline.calculate_bearing_angle()
    aIntersect = list(index_cell.search(pBound_edge))
    nIntersect = len(aIntersect)
    for i in range(nIntersect):
        lID = aIntersect[i]
        pFeature = pLayer_mesh.GetFeature(lID)
        lCellID = pFeature.GetField('cellid')
        dArea = pFeature.GetField('area')
        #check whether the start point is in the mesh cell
        sCellID = "{:09d}".format(lCellID)
        if sCellID != '000234146':
            continue
        pGeometry = pFeature.GetGeometryRef()
        if pGeometry.Contains(pPoint_start):
            sFolder = os.path.join(sWorkspace_output, sCellID)
            if not os.path.exists(sFolder):
                os.makedirs(sFolder)
                os.chmod(sFolder, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)

            sFilename_dem_full = os.path.join(sFolder, 'dem_full.tif')
            sFilename_dem_clip = os.path.join(sFolder, 'dem_clip.tif')
            sFilename_dem_box = os.path.join(sFolder, 'dem_box.tif')
            sFilename_mesh_cell = os.path.join(sFolder, 'mesh_cell.geojson')
            sFilename_output_in = os.path.join(sFolder, 'box.geojson')

            if iFlag_download_data == 1:
                if os.path.exists(sFilename_dem_full):
                    os.remove(sFilename_dem_full)
                    pass
                else:
                    pass
                #now use the mesh cell to download the elevation data
                bounds = pGeometry.GetEnvelope()
                minx, maxx, miny, maxy = bounds
                bounds = (minx, miny, maxx, maxy)
                # Suppress output for elevation.clip
                elevation.clip(bounds=bounds, output=sFilename_dem_full)
                if os.path.exists(sFilename_dem_clip):
                    os.remove(sFilename_dem_clip)
                    pass
                else:
                    pass
                #call pyearth function to extract the elevation data using the geometry
                if os.path.exists(sFilename_mesh_cell):
                    pass
                else:
                    pDataset_cell = pDriver.CreateDataSource(sFilename_mesh_cell)
                    pLayer_cell = pDataset_cell.CreateLayer('mesh_cell', geom_type=ogr.wkbPolygon, srs=pSpatialRef)
                    #create a field for the mesh cell ID
                    sField = ogr.FieldDefn('cellid', ogr.OFTInteger)
                    pLayer_cell.CreateField(sField)
                    pFeature_cell = ogr.Feature(pLayer_cell.GetLayerDefn())
                    pFeature_cell.SetGeometry(pGeometry)
                    pFeature_cell.SetField('cellid', lCellID)
                    pLayer_cell.CreateFeature(pFeature_cell)
                    pFeature_cell.Destroy()
                    pDataset_cell = None
                #now we can use the pyearth function to extract the elevation data
                clip_raster_by_polygon_file(sFilename_dem_full, sFilename_mesh_cell, sFilename_dem_clip)
                dResolution_appro = np.sqrt(dArea)
                #find the nearest resolution
                dResolution_degree = find_nearest_resolution(dResolution_appro, pVertex_start.dLatitude_degree)
                #check file exist
                if os.path.isfile(sFilename_output_in):
                    os.remove(sFilename_output_in)
                    pass
                pBox_start, aCoordinates_out = create_box_from_longitude_latitude(pVertex_start.dLongitude_degree,
                                                                                      pVertex_start.dLatitude_degree,
                                                                 dResolution_degree,
                                                                   dResolution_degree,
                                                                   sFilename_output_in=sFilename_output_in)
                minx, maxx, miny, maxy = pBox_start
                bounds = (minx, miny, maxx, maxy)
                if os.path.exists(sFilename_dem_box):
                    os.remove(sFilename_dem_box)
                # Suppress output for elevation.clip
                elevation.clip(bounds=bounds, output=sFilename_dem_box)

            if iFlag_cross_section == 1:
                aIntersect_point = list(index_edge.search_surrounding([pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree]))
                nIntersect_point = len(aIntersect_point)
                oFlowline_start = None
                if nIntersect_point > 0:
                    for j in range(nIntersect_point):
                        lID = aIntersect_point[j]
                        oFlowline = aFlowline[lID]
                        #check whether it is the start or the end
                        if oFlowline.pVertex_start == pVertex_start:
                            print("Found matching start vertex.")
                        elif oFlowline.pVertex_end == pVertex_start:
                            #this is the one we need
                            oFlowline_start = oFlowline
                            break
                if oFlowline_start is None:
                    print('No matching flowline found for the start point.')
                    continue
                #calcualte bearing angle using the start and end point
                dBearing1 = oFlowline_start.calculate_bearing_angle()
                dBearing = (dBearing0 + dBearing1) / 2.0
                #now take the perpendicular bearing angle
                dBearing = dBearing + 90.0
                if dBearing >= 360.0:
                    dBearing = dBearing - 360.0
                sFilename_raster1 = os.path.join(sFolder, 'dem_clip.tif')
                aExtent = gdal_get_raster_extent(sFilename_raster1)
                #get center
                center_lon = (aExtent[0] + aExtent[1]) / 2
                center_lat = (aExtent[2] + aExtent[3]) / 2
                center_point_wgs84 = (center_lon, center_lat)
                sFilename_raster0 = os.path.join(sFolder, 'dem_box.tif')
                aExtent = gdal_get_raster_extent(sFilename_raster0)
                #get center
                center_lon = (aExtent[0] + aExtent[1]) / 2
                center_lat = (aExtent[2] + aExtent[3]) / 2
                center_point_wgs84_0 = (center_lon, center_lat)
                pixel_coords, geographic_coords, wgs84_coords, values, distances = get_raster_cross_section_profile(sFilename_raster0, center_point_wgs84_0, dBearing)
                #save geographic_coords as a geojson file
                sFilename_geojson = os.path.join(sFolder, 'dem_box_cross_section.geojson'    )
                if os.path.exists(sFilename_geojson):
                    os.remove(sFilename_geojson)
                pDataset = pDriver.CreateDataSource(sFilename_geojson)
                pLayer = pDataset.CreateLayer('cross_section', geom_type=ogr.wkbLineString, srs=pSpatialRef)
                pLayer.CreateField(ogr.FieldDefn('distance', ogr.OFTReal))
                pLayer.CreateField(ogr.FieldDefn('elev', ogr.OFTReal))
                pFeature = ogr.Feature(pLayer.GetLayerDefn())
                pLine = ogr.Geometry(ogr.wkbLineString)
                for k in range(len(geographic_coords)):
                    lon, lat = geographic_coords[k]
                    pLine.AddPoint(lon, lat)
                    pFeature.SetGeometry(pLine)
                    pFeature.SetField('distance', distances[k])
                    pFeature.SetField('elev', float(values[k]))
                    pLayer.CreateFeature(pFeature)
                pFeature.Destroy()
                pDataset = None

                ymin=100
                ymax=200
                sFilename_out = os.path.join(sFolder, 'dem_box_cross_section.png')
                plot_xy_data([distances],
                 [values],
                 sFilename_out, dMin_x_in=0, dMax_x_in=None,
                 dMin_y_in=ymin, dMax_y_in=ymax,
                 aColor_in=['red'], aLinestyle_in=['-'],
                    aLabel_tag_in = ['(c)', 'Structured Latlon mesh-based'],
                    sTitle_in='Cross-Section Profile',
                 sLabel_x_in='Distance (unit: m)', sLabel_y_in='Surface elevation (unit: m)')


                pixel_coords, geographic_coords, wgs84_coords, values, distances = get_raster_cross_section_profile(sFilename_raster1, center_point_wgs84, dBearing)
                sFilename_geojson = os.path.join(sFolder, 'dem_clip_cross_section.geojson'    )
                if os.path.exists(sFilename_geojson):
                    os.remove(sFilename_geojson)
                #save geographic_coords as a geojson file
                pDataset = pDriver.CreateDataSource(sFilename_geojson)
                pLayer = pDataset.CreateLayer('cross_section', geom_type=ogr.wkbLineString, srs=pSpatialRef)
                pLayer.CreateField(ogr.FieldDefn('distance', ogr.OFTReal))
                pLayer.CreateField(ogr.FieldDefn('elev', ogr.OFTReal))
                pFeature = ogr.Feature(pLayer.GetLayerDefn())
                pLine = ogr.Geometry(ogr.wkbLineString)
                for k in range(len(geographic_coords)):
                    lon, lat = geographic_coords[k]
                    pLine.AddPoint(lon, lat)
                    pFeature.SetGeometry(pLine)
                    pFeature.SetField('distance', distances[k])
                    pFeature.SetField('elev', float(values[k]))
                    pLayer.CreateFeature(pFeature)
                pFeature.Destroy()
                pDataset = None

                sFilename_out = os.path.join(sFolder, 'dem_clip_cross_section.png')
                plot_xy_data([distances],
                 [values],
                 sFilename_out, dMin_x_in=0, dMax_x_in=None,
                 dMin_y_in=ymin, dMax_y_in=ymax, 
                 aColor_in=['blue'], aLinestyle_in=['-'],
                    aLabel_tag_in = ['(d)','UNstructured MPAS mesh-based'],
                    sTitle_in='Cross-Section Profile',
                 sLabel_x_in='Distance (unit: m)', sLabel_y_in='Surface elevation (unit: m)')

elevation.clean()

