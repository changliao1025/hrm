import sys, os
from osgeo import osr, ogr
from pyearth.system.define_global_variables import *
from pyearth.gis.spatialref.get_utm_spatial_reference import get_utm_spatial_reference_wkt
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.gis.spatialref.conversion_between_degree_and_meter import degree_to_meter
from pyflowline.mesh.triangular.create_triangular_mesh import create_triangular_mesh



#design a simple method to create a latlon mesh on the sphere,
#this mesh will not cover the whole sphere, but only a small part of it

#the mesh is defined by the following parameters:
#location of the center point, which is the cell that will be highlighted
#number of cells in the x direction
#number of cells in the y direction
#size of the cell in the x direction
#size of the cell in the y direction
#for different context, we will hightlight different regions on earth

#use the delaware bay as an example

dLongitude_center = -77
dLatitude_center = 41

iNumber_cell_x = 9
iNumber_cell_y = 9
dResolution_x = 0.5
dResolution_y = 0.5

dLongitude_left_in = dLongitude_center - 0.5*iNumber_cell_x*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*iNumber_cell_y*dResolution_y

dResolution_degree_in = dResolution_x
dResolution = degree_to_meter(dLatitude_center, dResolution_x)

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output/triangular'
if not os.path.exists(sWorkspace_out):
    os.makedirs(sWorkspace_out)

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_out, 'boundary.geojson')

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

sFilename_output_in = os.path.join(sWorkspace_out, 'triangular.geojson')

pSpatialRef_source = osr.SpatialReference()
pSpatialRef_source.ImportFromEPSG(4326)

#find the best UTM zone for the center point
pProjection = get_utm_spatial_reference_wkt(dLongitude_center)
#conver the center point to UTM
uv_xcenter, uv_ycenter = reproject_coordinates(
        dLongitude_center, dLatitude_center, pSpatialRef_source.ExportToWkt(), pProjection_target_in = pProjection)

dX_left_in = uv_xcenter - 0.5*iNumber_cell_x*dResolution
dY_bot_in = uv_ycenter - 0.5*iNumber_cell_y*dResolution
dResolution_meter_in = dResolution
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y

pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)

create_triangular_mesh(dX_left_in, dY_bot_in,
                    dResolution_meter_in,
                    ncolumn_in, nrow_in,
                    sFilename_output_in,
                    pProjection  ,
                    pBoundary_in =  pBoundary_wkt    )