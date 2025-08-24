import sys, os
from osgeo import ogr, osr, gdal
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet
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

# lake michgan center
dLongitude_center = -87
dLatitude_center = 43.5

iNumber_cell_x = 9
iNumber_cell_y = 9
dResolution_x = 2
dResolution_y = 2

#1km resolution in degree
dResolution_1km = 30/3600.0 * 5
dResolution_x = dResolution_1km
dResolution_y = dResolution_1km
iNumber_cell_x = 9 * 10
iNumber_cell_y = 9 * 10

#convert the resolution from degree to km
dResolution = degree_to_meter(dResolution_x, dLatitude_center)

dLongitude_left_in = dLongitude_center - 0.5*iNumber_cell_x*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*iNumber_cell_y*dResolution_y
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = dResolution_x

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'

sWorkspace_output_latlon = os.path.join(sWorkspace_out, 'latlon')
if not os.path.exists(sWorkspace_output_latlon):
    os.makedirs(sWorkspace_output_latlon)

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_output_latlon, 'boundary_greatlakes_2.geojson')



ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = dResolution_x


sFilename_output_in = os.path.join(sWorkspace_output_latlon, 'latlon_globe.geojson')
pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)
dLongitude_left_in = -180
dLatitude_bot_in = -90
ncolumn_in = int(360 /10)
nrow_in = int(180/10)
dResolution_degree_in = 10.0
create_latlon_mesh(dLongitude_left_in,
                       dLatitude_bot_in,
                       dResolution_degree_in,
                       ncolumn_in, nrow_in,
                       sFilename_output_in)
#,
                 #      pBoundary_in=pBoundary_wkt)

sFilename_mesh_parquet = os.path.join(sWorkspace_output_latlon, 'latlon_greatlakes_1km.parquet')

convert_geojson_to_geoparquet(sFilename_output_in, sFilename_mesh_parquet)