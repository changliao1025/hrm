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
sPath = '/qfs/people/liao313/workspace/python/subset/'
sys.path.append(sPath)


sWorkspace_output = '/compyfs/liao313/04model/subset/mississippi/elevation/land'
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

sFolder_cache = elevation.CACHE_DIR
if os.path.exists(sFolder_cache):
    shutil.rmtree(sFolder_cache)

iCount = 0
#for i in range(nCell):
for i, pFeature in enumerate(pLayer_mesh):

    fid = pFeature.GetFID()
    #pFeature = pLayer_mesh.GetFeature(fid)
    pGeometry = pFeature.GetGeometryRef()
    #get bounding box of the mesh cell

    lCellID = pFeature.GetField('cellid')
    dArea = pFeature.GetField('area')
    sCellID = "{:09d}".format(lCellID)
    #check whether the start point is in the mesh cell
    dArea_km = dArea / 1.E6
    if dArea_km > 1.E3:
        iCount += 1
        sFolder = os.path.join(sWorkspace_output, sCellID)
        if not os.path.exists(sFolder):
            os.makedirs(sFolder)
            os.chmod(sFolder, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        sFilename_dem_full = os.path.join(sFolder, 'dem_full.tif')
        if os.path.exists(sFilename_dem_full):
            os.remove(sFilename_dem_full)
            pass
        else:
            pass
        #now use the mesh cell to download the elevation data
        bounds = pGeometry.GetEnvelope()
        minx, maxx, miny, maxy = bounds
        bounds = (minx, miny, maxx, maxy)
        dLongitude_degree = (maxx + minx) / 2.0
        dLatitude_degree = (maxy + miny) / 2.0

        # Suppress output for elevation.clip
        elevation.clip(bounds=bounds, output=sFilename_dem_full)
        sFilename_dem_clip = os.path.join(sFolder, 'dem_clip.tif')
        if os.path.exists(sFilename_dem_clip):
            os.remove(sFilename_dem_clip)
            pass
        else:
            pass
        #call pyearth function to extract the elevation data using the geometry
        sFilename_mesh_cell = os.path.join(sFolder, 'mesh_cell.geojson')
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
        dResolution_degree = find_nearest_resolution(dResolution_appro, dLatitude_degree)
        sFilename_output_in = os.path.join(sFolder, 'box.geojson')
        #check file exist
        if os.path.isfile(sFilename_output_in):
            os.remove(sFilename_output_in)

        pBox_start, aCoordinates_out = create_box_from_longitude_latitude(dLongitude_degree,
                                                                              dLatitude_degree,
                                                         dResolution_degree,
                                                           dResolution_degree,
                                                           sFilename_output_in=sFilename_output_in)

        minx, maxx, miny, maxy = pBox_start
        bounds = (minx, miny, maxx, maxy)
        sFilename_dem_box = os.path.join(sFolder, 'dem_box.tif')
        if os.path.exists(sFilename_dem_box):
            os.remove(sFilename_dem_box)
        # Suppress output for elevation.clip
        elevation.clip(bounds=bounds, output=sFilename_dem_box)

elevation.clean()

