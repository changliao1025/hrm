import sys, os

from pyflowline.mesh.latlon.create_latlon_mesh import create_latlon_mesh

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
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = dResolution_x

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'
sFilename_output_in = os.path.join(sWorkspace_out, 'latlon.geojson')

create_latlon_mesh(dLongitude_left_in,
                       dLatitude_bot_in,
                       dResolution_degree_in,
                       ncolumn_in, nrow_in,
                       sFilename_output_in)