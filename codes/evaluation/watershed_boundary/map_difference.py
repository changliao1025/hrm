import os, sys
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.visual.map.vector.map_multiple_vector_files import map_multiple_vector_files


sWorkspace = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary'

sFilename_mpas  = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary/edge_cell.geojson'
sFilename_drt = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary/drt_16th_watershed_boundary.geojson'
sFilename_watershed_boundary = '/qfs/people/liao313/data/hexwatershed/mississippi/vector/mississippi_boundary.geojson'
aExtent_mississippi = gdal_get_vector_extent(sFilename_watershed_boundary)

aFlag_fill_in = [0,0,0]
aFiletype_in = [3,3,3]


aFilename_in = [ sFilename_drt, sFilename_mpas, sFilename_watershed_boundary]

aColor_in = ['red', 'blue', 'black']
aThickness_in = [1, 1, 1]

aExtent_outlet = [-92.45899893825303195,-89.79801124815780611, 29.5510526281927639,  31.39007657266963491]

sFilename_output_in = os.path.join(sWorkspace, 'watershed_boundary_difference_outlet.png')
aLegend_in = list()
aLegend_in.append('(b)')

map_multiple_vector_files(aFiletype_in,
                             aFilename_in,
                             aFlag_fill_in=aFlag_fill_in,
                             aColor_in=aColor_in,
                             aThickness_in = aThickness_in,
                             iFlag_esri_hydro_image_in = 1,
                             iFlag_zebra_in = 1,
                             iFlag_title_in =1,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in = 'Watershed boundary difference',
                             aLegend_in=aLegend_in,
                             aExtent_in=aExtent_outlet)