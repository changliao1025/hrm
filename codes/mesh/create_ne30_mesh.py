import sys, os, stat
import shutil
from shutil import copy2
import json
from osgeo import ogr, osr, gdal
from pyearth.system.define_global_variables import *
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary

from pyflowline.mesh.cubicsphere.create_cubicsphere_mesh import create_cubicsphere_mesh


dLongitude_center = -77
dLatitude_center = 0

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
sWorkspace_output_mpas = os.path.join(sWorkspace_out, 'cubicsphere')
if not os.path.exists(sWorkspace_output_mpas):
    os.makedirs(sWorkspace_output_mpas)

sFilename_mesh = os.path.join(sWorkspace_output_mpas, 'cubicsphere.geojson')

#create a sFilename_boundary using the bounding box of the mesh
#write a function to generate a NE30 mesh
#create a sFilename_boundary using the bounding box of the mesh
sFilename_boundary = os.path.join(sWorkspace_output_mpas, 'boundary.geojson')

pDriver = ogr.GetDriverByName('GeoJSON')


sWorkspace_output = '/qfs/people/liao313/workspace/python/hrm/data/output/cubicsphere'
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)

sFilename_output = os.path.join(sWorkspace_output, 'cubicsphere_globe.geojson')
pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)


aMpas = create_cubicsphere_mesh(1,
                            1,
                       dResolution,
                       sFilename_output,
                       sWorkspace_output)
                       #pBoundary_in=pBoundary_wkt)