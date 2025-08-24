import sys, os
sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from codes.shared.create_gradual_latlon_mesh import create_gradual_latlon_mesh

import numpy as np
from osgeo import ogr, osr, gdal
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
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

iNumber_cell_x = 13
iNumber_cell_y = 13
dResolution_x = 0.5
dResolution_y = 0.5

dLongitude_left_in = dLongitude_center - 0.5*iNumber_cell_x*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*iNumber_cell_y*dResolution_y
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = 0.2

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'

#create a sFilename_boundary using the bounding box of the mesh
#write a function to generate a NE30 mesh
#create a sFilename_boundary using the bounding box of the mesh
sWorkspace_output_latlon = os.path.join(sWorkspace_out, 'latlon')
if not os.path.exists(sWorkspace_output_latlon):
    os.makedirs(sWorkspace_output_latlon)

sFilename_output_in = os.path.join(sWorkspace_output_latlon, 'latlon_gradual.geojson')
sFilename_boundary = os.path.join(sWorkspace_output_latlon, 'boundary.geojson')




aLongitude_left_in = np.full(ncolumn_in + 1, dLongitude_center, dtype=float)

#create an array with different spacing in the x direction but it is symetric from left to right
center_index = ncolumn_in // 2

# Create an array with different spacing in the x direction that is symmetric from left to right
for i in range(1, center_index +1):
    dSpace = i * dResolution_degree_in
    aLongitude_left_in[center_index + i] = aLongitude_left_in[center_index + i - 1] + dSpace #right
    aLongitude_left_in[center_index - i] = aLongitude_left_in[center_index - i + 1] - dSpace

aLongitude_left_in[ncolumn_in] = aLongitude_left_in[ncolumn_in - 1] + dSpace

# Create an array with different spacing in the y direction that is symmetric from bottom to top
aLatitude_bottom_in = np.full(nrow_in + 1, dLatitude_center, dtype=float)
center_index = nrow_in // 2
for i in range(1, center_index+1 ):
    dSpace = i * dResolution_degree_in
    aLatitude_bottom_in[center_index + i] = aLatitude_bottom_in[center_index + i - 1] - dSpace #downward
    aLatitude_bottom_in[center_index - i] = aLatitude_bottom_in[center_index - i + 1] + dSpace #upward

aLatitude_bottom_in[nrow_in] = aLatitude_bottom_in[nrow_in - 1] - dSpace

print(aLongitude_left_in)
print(aLatitude_bottom_in)

pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)

create_gradual_latlon_mesh(dLongitude_left_in,
                       dLatitude_bot_in,
                       ncolumn_in, nrow_in,
                        aLongitude_left_in,
                       aLatitude_bottom_in,
                       sFilename_output_in,
                       pBoundary_in = pBoundary_wkt)