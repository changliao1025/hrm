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


dBuffer = 1.0
#lake michigan -87.09209,43.73615
aExtent_lakemichigan = [-87.09209 - 1.2 *dBuffer, -87.09209 + 1.2 *dBuffer, 43.73615 - 2.0*dBuffer, 43.73615 +1.1* dBuffer]
aExtent = aExtent_lakemichigan
aFiletype_in = [3,  3]
sFilename_mesh= '/compyfs/liao313/04model/pyflowline/conus/pyflowline20241201019/mpas.geojson'
sFilename_river_network = '/compyfs/liao313/00raw/hydrology/hydroshed/hydrolake/global_lakes.geojson'
aFilename_in = [sFilename_mesh, sFilename_river_network]
sFilename_output_in = '/qfs/people/liao313/workspace/python/hrm/data/output/mpas/mpas_mesh_lakemichigan.png'

#use the McNary Dam location for highlight

map_multiple_vector_files(aFiletype_in,
                             aFilename_in,
                             iFlag_colorbar_in = None,
                             iFlag_title_in = 1,
                             iFlag_zebra_in=1,
                             iFlag_filter_in = 1,
                             aFlag_thickness_in = None,
                             aThickness_in = [0.5, 2],
                             aFlag_color_in = None,
                             aColor_in=['blue', 'red'],
                             aFlag_discrete_in = None,
                             aFlag_fill_in = None,
                             aVariable_in = None,
                             sFilename_output_in=sFilename_output_in,
                             iFlag_scientific_notation_colorbar_in = None,
                             iFlag_openstreetmap_in=0,
                             iFlag_esri_hydro_image_in=1,
                             iBasemap_zoom_level_in=None,
                             iFont_size_in=None,
                             sColormap_in = None,
                             sTitle_in = 'Lake Michigan',
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
