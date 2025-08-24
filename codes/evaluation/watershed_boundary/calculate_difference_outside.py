import os, sys
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.visual.map.vector.map_multiple_vector_files import map_multiple_vector_files

from pyearth.toolbox.analysis.difference.difference_polygon_with_polygon_file import calculate_polygon_difference
from pyearth.toolbox.management.vector.merge_features import merge_features

from pyearth.toolbox.conversion.convert_polygon_to_polyline_file import convert_polygon_to_polyline_file

sWorkspace = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary'

sFilename_mpas  = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary/edge_cell.geojson'
sFilename_mpas = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501003/pyflowline/mpas.geojson'

sFilename_drt = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary/drt_16th_mesh_wo_holes.geojson'
sFilename_watershed_boundary = '/qfs/people/liao313/data/hexwatershed/mississippi/vector/mississippi_boundary.geojson'
aExtent_mississippi = gdal_get_vector_extent(sFilename_watershed_boundary)

#convert polygon to polyline
sFilename_watershed_boundary_new = os.path.join(sWorkspace, 'mississippi_boundary_polyline.geojson')

sFilename_drt_new = os.path.join(sWorkspace, 'drt_new.geojson')
merge_features(sFilename_drt, sFilename_drt_new)
sFilename_out = os.path.join(sWorkspace, 'drt_difference_outside.geojson')
calculate_polygon_difference( sFilename_watershed_boundary, sFilename_drt_new, sFilename_out  )

sFilename_mpas_new = os.path.join(sWorkspace, 'mpas_new.geojson')
merge_features(sFilename_mpas, sFilename_mpas_new)
sFilename_out = os.path.join(sWorkspace, 'mpas_difference_outside.geojson')
calculate_polygon_difference( sFilename_watershed_boundary, sFilename_mpas_new, sFilename_out  )

