import sys, os
sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from codes.shared.map_vector_polygon_file import map_vector_polygon_file

import numpy as np
from osgeo import ogr, osr, gdal
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.visual.map.pick_colormap import pick_colormap_hydrology

dLongitude_center = -77
dLatitude_center = 41

iFiletype_in=1
sFilename_in='/qfs/people/liao313/workspace/python/hrm/data/output/tin/tin.geojson'
sFilename_boundary = '/qfs/people/liao313/workspace/python/hrm/data/output/tin/boundary.geojson'
sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/tin/tin.png'

aExtent = [-79.10236320495605, -74.35684242248536, 39.374137496948244, 42.9556583404541]
pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)
image_size = [1000, 1000]
dpi = 150
scale_denominator = calculate_scale_denominator(aExtent, image_size)
pSrc = osr.SpatialReference()
pSrc.ImportFromEPSG(3857) # mercator
pProjection = pSrc.ExportToWkt()
iFlag_openstreetmap_level = calculate_zoom_level(scale_denominator, pProjection, dpi=dpi)
print(iFlag_openstreetmap_level)

map_vector_polygon_file(iFiletype_in,
                          sFilename_in,
                        dLongitude_highligt_in=dLongitude_center,
                         dLatitude_highlight_in=dLatitude_center,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_terrain_image_in = 1,
                          iFlag_openstreetmap_level_in = iFlag_openstreetmap_level,
                          iFlag_scientific_notation_colorbar_in=None,
                          sTitle_in='Unstructured TIN Mesh',
                          iFlag_zebra_in = 1,
                          iFlag_fill_in = 0,
                          iDPI_in=None,
                          dMissing_value_in=None,
                          aLegend_in=['(c)'],
                          aExtent_in =aExtent )
