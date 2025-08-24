import sys, os
import numpy as np
from osgeo import ogr, osr, gdal
import cartopy as cpl
import cartopy.crs as ccrs

from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from codes.shared.map_vector_polyline_file import map_vector_polyline_file


iFiletype_in=2
sFilename_in = '/compyfs/liao313/00raw/hydrology/hydroshed/hydroriver/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp'
#sFilename_in='/compyfs/liao313/04model/pyflowline/conus/flowline_hydroshed_simplified.geojson'
#sFilename_in= '/qfs/people/liao313/data/hexwatershed/conus/vector/river_networks_wo_lakes.geojson'

sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/original_river_networks.png'
#sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/river_networks_with_lakes.png'
#sFilename_output_in= '/qfs/people/liao313/workspace/python/hrm/data/output/river_networks_without_lakes.png'


aExtent = [-89.0, -85.0, 40.0, 44.0]

dLongitude_center = 0.5 *(aExtent[0] + aExtent[1])
dLatitude_center = 0.5 * (aExtent[2] + aExtent[3])

#create a wkt boundary from the extent

image_size = [1000, 1000]
dpi = 150
scale_denominator = calculate_scale_denominator(aExtent, image_size)
pSrc = osr.SpatialReference()
pSrc.ImportFromEPSG(3857) # mercator
pProjection = pSrc.ExportToWkt()
iFlag_openstreetmap_level = calculate_zoom_level(scale_denominator, pProjection, dpi=dpi)
print(iFlag_openstreetmap_level)

pProjection_map = cpl.crs.Orthographic(central_longitude =  dLongitude_center,
                                                central_latitude = dLatitude_center, globe=None)

sTitle = 'Original river networks with lakes'
#sTitle = 'Simplified river networks with lakes'
#sTitle = 'River networks without lakes'
map_vector_polyline_file(iFiletype_in,
                          sFilename_in,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_terrain_image_in = 1,
                          iFlag_openstreetmap_level_in = iFlag_openstreetmap_level,
                          iFlag_scientific_notation_colorbar_in=None,
                          sTitle_in=sTitle,
                          iFlag_zebra_in = 1,
                          iDPI_in=None,
                          dMissing_value_in=None,
                          aLegend_in=['(a)'],
                          pProjection_map_in = pProjection_map,
                          aExtent_in =aExtent )
