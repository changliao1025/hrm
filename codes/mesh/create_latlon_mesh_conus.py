import sys, os
from osgeo import ogr, osr, gdal
from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.toolbox.data.geoparquet.convert_geojson_to_geoparquet import convert_geojson_to_geoparquet
from pyearth.gis.geometry.create_box_from_longitude_latitude import create_box_from_longitude_latitude
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



sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'

sWorkspace_output_latlon = os.path.join(sWorkspace_out, 'latlon')
if not os.path.exists(sWorkspace_output_latlon):
    os.makedirs(sWorkspace_output_latlon)

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'

pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)



#1km resolution in degree

dResolution_x = 1/8.0
dResolution_y = dResolution_x

dMinx, dMaxx, dMiny, dMaxy = aExtent
dLongitude_in = dMinx
dLatitude_in = dMiny
aBox_out, aCoordinates_out = create_box_from_longitude_latitude(dLongitude_in, dLatitude_in,
                                       dResolution_x, dResolution_y)


dLongitude_left_in =aBox_out[0]
dLatitude_bot_in = aBox_out[2]

dLongitude_in = dMaxx
dLatitude_in = dMaxy
aBox_out, aCoordinates_out = create_box_from_longitude_latitude(dLongitude_in, dLatitude_in,
                                       dResolution_x, dResolution_y)

dLongitude1 = aBox_out[0]
dLatitude1 = aBox_out[2]

ncolumn_in = int((dLongitude1 - dLongitude_left_in) / dResolution_x) + 1
nrow_in = int((dLatitude1 - dLatitude_bot_in) / dResolution_y) + 1

sFilename_output_in = os.path.join(sWorkspace_output_latlon, 'latlon_conus_8th.geojson')

dResolution_degree_in=dResolution_x
create_latlon_mesh(dLongitude_left_in,
                       dLatitude_bot_in,
                       dResolution_degree_in,
                       ncolumn_in, nrow_in,
                       sFilename_output_in,
                       pBoundary_in=pBoundary_wkt)

sFilename_mesh_parquet = os.path.join(sWorkspace_output_latlon, 'latlon_conus_8th.parquet')

convert_geojson_to_geoparquet(sFilename_output_in, sFilename_mesh_parquet)