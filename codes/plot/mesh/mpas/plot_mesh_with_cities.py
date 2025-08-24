import sys, os
sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from codes.shared.map_multiple_vector_files import map_multiple_vector_files

import numpy as np
from osgeo import ogr, osr, gdal
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.visual.map.pick_colormap import pick_colormap_hydrology

dLongitude_center = -77
dLatitude_center = 41
sColormap = pick_colormap_hydrology('nse')
iFiletype_in=1
sFilename_in='/qfs/people/liao313/workspace/python/hrm/data/output/mpas/mpas.geojson'
sFilename_boundary = '/qfs/people/liao313/workspace/python/hrm/data/output/mpas/boundary.geojson'
sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/mpas/mpas.png'


dBuffer = 0.75
#altlanta city -84.38798,33.74899
dLongitude_highlight = -84.3880
dLatitude_highlight = 33.7490
aExtent =  [-84.3880 - dBuffer, -84.3880 + dBuffer, 33.7490 - dBuffer, 33.7490 + dBuffer]

aFiletype_in = [3,  1]
sFilename_mesh= '/compyfs/liao313/04model/pyflowline/conus/pyflowline20241201019/mpas.geojson'

sFilename_city = '/compyfs/liao313/00raw/city/simplemaps_worldcities_basicv1.77/large_cities.geojson'
aFilename_in = [sFilename_mesh, sFilename_city]
sFilename_output_in = '/qfs/people/liao313/workspace/python/hrm/data/output/mpas/mpas_w_cites.png'

#use the McNary Dam location for highlight

map_multiple_vector_files(aFiletype_in,
                             aFilename_in,
                             iFlag_colorbar_in = None,
                             iFlag_title_in = 1,
                             iFlag_zebra_in=1,
                             iFlag_filter_in = 1,
                             dLongitude_highlight_in=dLongitude_highlight,
                             dLatitude_highlight_in=dLatitude_highlight,
                             aFlag_thickness_in = None,
                             aFlag_color_in = None,
                             aFlag_discrete_in = None,
                             aFlag_fill_in = None,
                             aVariable_in = None,
                             sFilename_output_in=sFilename_output_in,
                             iFlag_scientific_notation_colorbar_in = None,
                             iFlag_openstreetmap_in=1,
                             iBasemap_zoom_level_in=None,
                             iFont_size_in=None,
                             sColormap_in = None,
                             sTitle_in = 'Atlanta City',
                             iDPI_in = None,
                             iSize_x_in = None,
                             iSize_y_in = None,
                             aMissing_value_in = None,
                             aData_max_in = None,
                             aData_min_in = None,
                             sExtend_in =None,
                             sFont_in = None,
                             sUnit_in=None,
                             aLegend_in = None,
                             aExtent_in = aExtent,
                             pProjection_map_in=None,
                             pProjection_data_in=None,)
