import sys, os
sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath

sys.path.append(sPath_project)
from codes.shared.map_vector_polygon_file import map_vector_polygon_file

import numpy as np
from osgeo import ogr, osr, gdal
import cartopy.crs as ccrs
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.visual.color.pick_colormap import pick_colormap_hydrology

dLongitude_center = -77
dLatitude_center = 41
sColormap = pick_colormap_hydrology('nse')
iFiletype_in=1
sFilename_in='/qfs/people/liao313/workspace/python/hrm/data/output/cubicsphere/cubicsphere_globe.geojson'
sFilename_boundary = '/qfs/people/liao313/workspace/python/hrm/data/output/cubicsphere/boundary.geojson'
sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/cubicsphere/cubicsphere.png'
sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/cubicsphere/cubicsphere_globe.png'


aExtent = [-79.10236320495605, -74.35684242248536, 39.374137496948244, 42.9556583404541]
aExtent_global = [-180,  180,-90, 90]
sWorkspace_out = '/qfs/people/liao313/workspace/python/hrm/data/output'
sWorkspace_output_mpas = os.path.join(sWorkspace_out, 'cubicsphere')
sFilename_boundary = os.path.join(sWorkspace_output_mpas, 'boundary.geojson')

pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)

image_size = [1000, 1000]
dpi = 150
scale_denominator = calculate_scale_denominator(aExtent, image_size)
pSrc = osr.SpatialReference()
pSrc.ImportFromEPSG(3857) # mercator
pProjection = pSrc.ExportToWkt()
iFlag_openstreetmap_level = calculate_zoom_level(scale_denominator, pProjection, dpi=dpi)
print(iFlag_openstreetmap_level)

pProjection_map = ccrs.Orthographic(central_longitude=0, central_latitude=0, globe=None)

map_vector_polygon_file(iFiletype_in,
                          sFilename_in,
                        #dLongitude_highligt_in=dLongitude_center,
                        # dLatitude_highlight_in=dLatitude_center,
                         iFlag_set_extent_in = None,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_coastline_in =1,
                          iFlag_terrain_image_in = 0,
                          iFlag_openstreetmap_level_in = None,
                          iFlag_scientific_notation_colorbar_in=None,
                          #sTitle_in='Structured Cubicsphere Mesh',
                          #iFlag_zebra_in = 1,
                          iFlag_fill_in = 0,
                          iDPI_in=None,
                          dMissing_value_in=None,
                          #aLegend_in=['(b)'],

                           pProjection_map_in=pProjection_map)
