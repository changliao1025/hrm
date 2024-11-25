import sys, os, stat
import shutil
from shutil import copy2
import json
from osgeo import ogr, osr, gdal
from pyearth.system.define_global_variables import *
from pyearth.gis.spatialref.conversion_between_degree_and_meter import degree_to_meter
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyflowline.mesh.tin.create_tin_mesh import create_tin_mesh

dLongitude_center = -77
dLatitude_center = 41

iNumber_cell_x = 9
iNumber_cell_y = 9
dResolution_x = 0.5
dResolution_y = 0.5
#convert the resolution from degree to km
dResolution = degree_to_meter(dLatitude_center, dResolution_x)

dLongitude_left_in = dLongitude_center - 0.5*iNumber_cell_x*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*iNumber_cell_y*dResolution_y
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = dResolution_x

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'

iFlag_save_mesh = 1
sWorkspace_output_tin = os.path.join(sWorkspace_out, 'tin')
if not os.path.exists(sWorkspace_output_tin):
    os.makedirs(sWorkspace_output_tin)

sFilename_mesh = os.path.join(sWorkspace_output_tin, 'tin.geojson')

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_output_tin, 'boundary.geojson')

pDriver = ogr.GetDriverByName('GeoJSON')
if os.path.exists(sFilename_boundary):
    pDriver.DeleteDataSource(sFilename_boundary)

pSpatialRef = osr.SpatialReference()
pSpatialRef.ImportFromEPSG(4326)
#create a boundary for the mesh
pDataSource = pDriver.CreateDataSource(sFilename_boundary)
pLayer = pDataSource.CreateLayer('boundary', geom_type=ogr.wkbPolygon, srs=pSpatialRef)
pFeature = ogr.Feature(pLayer.GetLayerDefn())
pRing = ogr.Geometry(ogr.wkbLinearRing)
pRing.AddPoint(dLongitude_left_in, dLatitude_bot_in)
pRing.AddPoint(dLongitude_left_in + iNumber_cell_x*dResolution_x, dLatitude_bot_in)
pRing.AddPoint(dLongitude_left_in + iNumber_cell_x*dResolution_x, dLatitude_bot_in + iNumber_cell_y*dResolution_y)
pRing.AddPoint(dLongitude_left_in, dLatitude_bot_in + iNumber_cell_y*dResolution_y)
pRing.AddPoint(dLongitude_left_in, dLatitude_bot_in)
pPolygon = ogr.Geometry(ogr.wkbPolygon)
pPolygon.AddGeometry(pRing)
pFeature.SetGeometry(pPolygon)
pLayer.CreateFeature(pFeature)
pDataSource.Destroy()

sFilename_jigsaw_configuration = '/qfs/people/liao313/workspace/python/hrm/data/input/tin/pyflowline_jigsaw.json'
#check if it exists
if not os.path.exists(sFilename_jigsaw_configuration):
    print('The configuration file does not exist')
    exit()
else:
    print('The configuration file is: ', sFilename_jigsaw_configuration)
    #check it is a valid json file
with open(sFilename_jigsaw_configuration) as json_file:
    aConfig_jissaw = json.load(json_file)

sWorkspace_jigsaw = sWorkspace_output_tin

iFlag_use_mesh_dem = 0
iFlag_save_mesh = 1
#create a polygon based on
#read boundary
pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)

iFlag_global = 0
sFilename_mesh_netcdf= '/qfs/people/liao313/workspace/python/hrm/data/output/tin/tmp/mesh_triangles.nc'

aTin = create_tin_mesh( sFilename_mesh,
                       iFlag_global_in = 0,
                         iFlag_antarctic_in=0,
                         iFlag_arctic_in=0,
                         iFlag_run_jigsaw_in= 0,
                             sWorkspace_jigsaw_in = sWorkspace_jigsaw,
                         aConfig_jigsaw_in = aConfig_jissaw,
                         sFilename_mesh_netcdf_in = sFilename_mesh_netcdf,
                         pBoundary_in = pBoundary_wkt)
pass