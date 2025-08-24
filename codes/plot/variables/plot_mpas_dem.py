import sys, os
import numpy as np
from osgeo import ogr, osr, gdal
import cartopy as cpl
import cartopy.crs as ccrs
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.visual.map.pick_colormap import pick_colormap_terrain

sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
#from codes.shared.map_vector_polygon_file import map_vector_polygon_file
from pyearth.visual.map.vector.map_vector_polygon_file import map_vector_polygon_file

sFilename_dem = '/qfs/people/liao313/workspace/python/hrm/data/output/mpas/dem_subset.geojson'
sFilename_output_in = '/qfs/people/liao313/workspace/python/hrm/data/output/mpas/dem_subset.png'
iFiletype_in=1
sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/mississippi/vector/mississippi_boundary.geojson'
pBoundary_wkt, aExtent = gdal_read_geojson_boundary(sFilename_boundary)
sColormap = pick_colormap_terrain('elevation')
map_vector_polygon_file(iFiletype_in,
                          sFilename_dem,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_terrain_image_in = 0,
                          iFlag_scientific_notation_colorbar_in=0,
                          sTitle_in='Elevation',
                          iFlag_zebra_in = 1,
                          iFlag_color_in=1,
                          iFlag_colorbar_in = 1,
                          iFlag_fill_in= 1,
                          iFlag_filter_in = 0,
                          iDPI_in=600,
                          dMissing_value_in=None,
                          sField_color_in='mean',
                          sColormap_in=sColormap,
                          sUnit_in='Unit: m',
                          aLegend_in=None,
                          aExtent_in =aExtent)