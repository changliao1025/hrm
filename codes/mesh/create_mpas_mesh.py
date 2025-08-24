import sys, os, stat
import shutil
from shutil import copy2
import json
from osgeo import ogr, osr, gdal
from pyearth.system.define_global_variables import *
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyflowline.mesh.mpas.create_mpas_mesh import create_mpas_mesh





dLongitude_center = -77
dLatitude_center = 41

iNumber_cell_x = 9
iNumber_cell_y = 9
dResolution_x = 0.5
dResolution_y = 0.5
#convert the resolution from degree to km
dResolution = degree_to_meter(dResolution_x, dLatitude_center )

dLongitude_left_in = dLongitude_center - 0.5*iNumber_cell_x*dResolution_x
dLatitude_bot_in = dLatitude_center - 0.5*iNumber_cell_y*dResolution_y
ncolumn_in = iNumber_cell_x
nrow_in = iNumber_cell_y
dResolution_degree_in = dResolution_x

sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'

iFlag_save_mesh = 1
sWorkspace_output_mpas = os.path.join(sWorkspace_out, 'mpas')
if not os.path.exists(sWorkspace_output_mpas):
    os.makedirs(sWorkspace_output_mpas)

sFilename_mesh = os.path.join(sWorkspace_output_mpas, 'mpas_tin.geojson')

#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_output_mpas, 'boundary.geojson')



sFilename_jigsaw_configuration = '/qfs/people/liao313/workspace/python/hrm/data/input/mpas/pyflowline_jigsaw_tin.json'
#check if it exists
if not os.path.exists(sFilename_jigsaw_configuration):
    print('The configuration file does not exist')
    exit()
else:
    print('The configuration file is: ', sFilename_jigsaw_configuration)
    #check it is a valid json file
with open(sFilename_jigsaw_configuration) as json_file:
    aConfig_jissaw = json.load(json_file)

sWorkspace_jigsaw = sWorkspace_output_mpas

iFlag_use_mesh_dem = 0
iFlag_save_mesh =1
#create a polygon based on
#read boundary
pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)

iFlag_global = 0
aMpas = create_mpas_mesh(sFilename_mesh,
                          iFlag_global_in = 0,
                         iFlag_run_jigsaw_in=1,
                         iFlag_antarctic_in=0,
                         iFlag_arctic_in=0,
                         sWorkspace_jigsaw_in = sWorkspace_jigsaw,
                         aConfig_jigsaw_in = aConfig_jissaw,
                         pBoundary_in = pBoundary_wkt)
pass