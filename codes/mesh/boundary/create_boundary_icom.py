import sys, os, stat
import shutil
from shutil import copy2
import json
from osgeo import ogr, osr, gdal
from pyearth.system.define_global_variables import *
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh

sWorkspace_output = '/qfs/people/liao313/workspace/python/hrm/data/output'

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_output, 'boundary.geojson')

pDriver = ogr.GetDriverByName('GeoJSON')
if os.path.exists(sFilename_boundary):
    pDriver.DeleteDataSource(sFilename_boundary)

nrow = 9
ncolumn = 9
dResolution_x = 0.5
dResolution_y = 0.5
dLongitude_center = -77
dLatitude_center = 41

dLongitude_left_in = dLongitude_center - 0.5*ncolumn*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*nrow*dResolution_y

pSpatialRef = osr.SpatialReference()
pSpatialRef.ImportFromEPSG(4326)
#create a boundary for the mesh
pDataSource = pDriver.CreateDataSource(sFilename_boundary)
pLayer = pDataSource.CreateLayer('boundary', geom_type=ogr.wkbPolygon, srs=pSpatialRef)
pFeature = ogr.Feature(pLayer.GetLayerDefn())
pRing = ogr.Geometry(ogr.wkbLinearRing)
pRing.AddPoint(dLongitude_left_in, dLatitude_bot_in)
pRing.AddPoint(dLongitude_left_in + ncolumn*dResolution_x, dLatitude_bot_in)
pRing.AddPoint(dLongitude_left_in + ncolumn*dResolution_x, dLatitude_bot_in + nrow*dResolution_y)
pRing.AddPoint(dLongitude_left_in, dLatitude_bot_in + nrow*dResolution_y)
pRing.AddPoint(dLongitude_left_in, dLatitude_bot_in)
pPolygon = ogr.Geometry(ogr.wkbPolygon)
pPolygon.AddGeometry(pRing)
pFeature.SetGeometry(pPolygon)
pLayer.CreateFeature(pFeature)
pDataSource.Destroy()

print('done')