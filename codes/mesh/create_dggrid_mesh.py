import sys, os, stat
import shutil
from shutil import copy2
from osgeo import ogr, osr, gdal
from pyearth.system.define_global_variables import *
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_index_by_resolution
from pyflowline.mesh.dggrid.create_dggrid_mesh import create_dggrid_mesh



dLongitude_center = -77
dLatitude_center = 41

iNumber_cell_x = 9
iNumber_cell_y = 9
dResolution_x = 10
dResolution_y = 10
#convert the resolution from degree to km
dResolution = degree_to_meter(dResolution_x, dLatitude_center )

dLongitude_left_in = dLongitude_center - 0.5*iNumber_cell_x*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*iNumber_cell_y*dResolution_y
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = dResolution_x

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'

iFlag_save_mesh = 1
sWorkspace_output_dggrid = os.path.join(sWorkspace_out, 'dggrid')
if not os.path.exists(sWorkspace_output_dggrid):
    os.makedirs(sWorkspace_output_dggrid)

sFilename_mesh = os.path.join(sWorkspace_output_dggrid, 'dggrid.geojson')

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_output_dggrid, 'boundary.geojson')


sDggrid_type = 'ISEA3H'
iResolution_index = dggrid_find_index_by_resolution(sDggrid_type, dResolution)

sFilename_dggrid_bin = '/qfs/people/liao313/bin/dggrid'
#copy the binary to the target folder
sFilename_new = sWorkspace_output_dggrid + slash + 'dggrid'
copy2(sFilename_dggrid_bin, sFilename_new)
os.chmod(sFilename_new, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)
create_dggrid_mesh(1,
                   iFlag_save_mesh,
                   sFilename_mesh,
                   sWorkspace_output_dggrid,
                   iResolution_index_in= iResolution_index,
                   iFlag_antarctic_in=0,
                   iFlag_arctic_in=0,
                   sDggrid_type_in = sDggrid_type)
                   #sFilename_boundary_in = sFilename_boundary)