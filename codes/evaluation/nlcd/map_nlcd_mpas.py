import sys, os
import numpy as np
from osgeo import ogr, osr, gdal
import cartopy as cpl
import cartopy.crs as ccrs
from pyearth.gis.gdal.read.vector.gdal_read_geojson_boundary import gdal_read_geojson_boundary
from pyearth.visual.color.pick_colormap import pick_colormap_terrain
from pyearth.visual.formatter import OOMFormatter
sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
#from codes.shared.map_vector_polygon_file import map_vector_polygon_file
from pyearth.visual.map.vector.map_vector_polygon_file import map_vector_polygon_file
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
sFilename_dem = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/nlcd/nlcd_subset_mpas.geojson'

sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'

aExtent_conus = gdal_get_vector_extent(sFilename_boundary)
iFiletype_in=1
sFilename_output_in = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/nlcd/nlcd_subset_mpas.png'
sColormap = pick_colormap_terrain('lucc')

#create discrete labels for NLCD data
aDict_discrete_labels = {
    11: 'Open Water',
    12: 'Perennial Ice/Snow',
    21: 'Developed, Open Space',
    22: 'Developed, Low Intensity',
    23: 'Developed, Medium Intensity',
    24: 'Developed, High Intensity',
    31: 'Barren Land (Rock/Sand/Clay)',
    41: 'Deciduous Forest',
    42: 'Evergreen Forest',
    43: 'Mixed Forest',
    51: 'Dwarf Scrub',
    52: 'Shrub/Scrub',
    71: 'Grassland/Herbaceous',
    72: 'Sedge/Herbaceous',
    81: 'Pasture/Hay',
    82: 'Cultivated Crops',
    90: 'Woody Wetlands',
    95: 'Emergent Herbaceous Wetlands',
    0: 'No Data'
}
aColor_discrete = {
    11: '#466B9F',  # Open Water (NLCD RGB: 70, 107, 159)
    12: '#D1DEF8',  # Perennial Ice/Snow (NLCD RGB: 209, 222, 248)
    21: '#DEE4E7',  # Developed, Open Space (NLCD RGB: 222, 228, 231)
    22: '#E4CEAD',  # Developed, Low Intensity (NLCD RGB: 228, 206, 173)
    23: '#EA9E64',  # Developed, Medium Intensity (NLCD RGB: 234, 158, 100)
    24: '#C25D42',  # Developed, High Intensity (NLCD RGB: 194, 93, 66)
    31: '#EADCDC',  # Barren Land (Rock/Sand/Clay) (NLCD RGB: 234, 220, 220)
    41: '#A8D88E',  # Deciduous Forest (NLCD RGB: 168, 216, 142)
    42: '#38814E',  # Evergreen Forest (NLCD RGB: 56, 129, 78)
    43: '#DCCA8F',  # Mixed Forest (NLCD RGB: 220, 202, 143)
    51: '#F4EBBD',  # Dwarf Scrub (NLCD RGB: 244, 235, 189)
    52: '#B8D9E3',  # Shrub/Scrub (NLCD RGB: 184, 217, 227)
    71: '#AAACCA',  # Grassland/Herbaceous (NLCD RGB: 170, 172, 202)
    72: '#C9DD93',  # Sedge/Herbaceous (NLCD RGB: 201, 221, 147)
    73: '#B6D1B6',  # Lichens/Mosses (NLCD RGB: 182, 209, 182)
    74: '#C6DCD7',  # Emergent Herbaceous Wetlands (NLCD RGB: 198, 220, 215)
    81: '#76A8B4',  # Pasture/Hay (NLCD RGB: 118, 168, 180)
    82: '#C79D8E',  # Cultivated Crops (NLCD RGB: 199, 157, 142)
    90: '#7E967D',  # Woody Wetlands (NLCD RGB: 126, 150, 125)
    95: '#A2AD9C',  # Emergent Herbaceous Wetlands (NLCD RGB: 162, 173, 156) - Note: NLCD has two "Emergent Herbaceous Wetlands" categories (74 and 95) with slightly different colors.
    0:  '#FFFFFF'   # No Data / Background (White)
}
iFlag_plot_lucc = 0
if iFlag_plot_lucc == 1:
    map_vector_polygon_file(iFiletype_in,
                          sFilename_dem,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_scientific_notation_colorbar_in=0,
                          sTitle_in='Land cover (dominant type)',
                          iFlag_discrete_in=1,
                          iFlag_zebra_in = 1,
                          iFlag_color_in=1,
                          iFlag_colorbar_in = 1,
                          iFlag_fill_in= 1,
                          iFlag_filter_in = 0,
                          iDPI_in=600,
                          dMissing_value_in=None,
                          sField_color_in='type',
                          sColormap_in=sColormap,
                          sUnit_in='NLCD type',
                          aLegend_in=['(a)'],
                          aDict_discrete_labels_in = aDict_discrete_labels,
                          aDict_value_color_in = aColor_discrete,
                          aExtent_in =aExtent_conus)

sFilename_output_in = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/nlcd/nlcd_subset_ntype_mpas.png'

sColormap= 'Blues'
#format for integer
formatter = OOMFormatter(fformat = '%d')
map_vector_polygon_file(iFiletype_in,
                          sFilename_dem,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_terrain_image_in = 0,
                          iFlag_scientific_notation_colorbar_in=0,
                          sTitle_in='Number of land cover types',
                          iFlag_discrete_in=0,
                          iFlag_zebra_in = 1,
                          iFlag_color_in=1,
                          iFlag_colorbar_in = 1,
                          iFlag_fill_in= 1,
                          iFlag_filter_in = 0,
                          iDPI_in=600,
                          aLegend_in=['(c)', 'Unstructured MPAS mesh-based'],
                          dMissing_value_in=None,
                          sField_color_in='ntype',
                          sColormap_in=sColormap,
                          sFormat_colorbar_in=formatter,
                          sUnit_in='Count',
                          aExtent_in =aExtent_conus)