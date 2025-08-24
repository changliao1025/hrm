import os, sys
from pyearth.system.define_global_variables import *
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.visual.map.vector.map_vector_polygon_file import map_vector_polygon_file


sWorkspace = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation'

sFilename_mpas  = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation/dem_subset_mpas.geojson'
sFilename_drt = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/elevation/dem_subset_drt.geojson'
sFilename_watershed_boundary = '/qfs/people/liao313/data/hexwatershed/mississippi/vector/mississippi_boundary.geojson'
sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/conus/vector/conus_boundary.geojson'
aExtent_mississippi = gdal_get_vector_extent(sFilename_watershed_boundary)
aExtent_conus = gdal_get_vector_extent(sFilename_boundary)



aExtent_outlet = [-92.45899893825303195,-89.79801124815780611, 29.5510526281927639,  31.39007657266963491]
dBuffer = 1.0

#Sierra Nevada -119.85795,38.56839
aExtent_sierranevada = [-119.85795 - dBuffer, -119.85795 + dBuffer, 38.56839 - dBuffer, 38.56839 + dBuffer]
aExtent_rocky = [-108.84890749174736868,-101.40113177200181838, 36.41899616632803571,  40.89526136931795719]

aExtent = aExtent_conus

sFilename_output_in = os.path.join(sWorkspace, 'elevation_drt_conus.png')
#sFilename_output_in = os.path.join(sWorkspace, 'elevation_drt_sierranevada.png')
aLegend_in = list()
aLegend_in.append('(a)')
aLegend_in.append('Structured Latlon mesh-based')

dData_max= 4000
dData_min = 0
iFlag_drt = 1
if iFlag_drt == 1:
    map_vector_polygon_file(1, sFilename_drt,
                             iFlag_fill_in=1,
                             iFlag_zebra_in = 1,
                               iFlag_colorbar_in = 1,
                               iFlag_esri_hydro_image_in=0,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in = 'Surface elevation (mean)',
                             dData_max_in=dData_max,
                            dData_min_in=dData_min,
                            sColormap_in= 'terrain',
                            sField_color_in = 'mean',
                             sUnit_in='Unit: m',
                             aLegend_in=aLegend_in,
                             aExtent_in=aExtent)


aLegend_in = list()
sFilename_output_in = os.path.join(sWorkspace, 'elevation_mpas_conus.png')
#sFilename_output_in = os.path.join(sWorkspace, 'elevation_mpas_sierranevada.png')
aLegend_in.append('(b)')
aLegend_in.append('Unstructured MPAS mesh-based')
iFlag_mpas = 1
if iFlag_mpas ==1:
    map_vector_polygon_file(1, sFilename_mpas,
                             iFlag_fill_in=1,
                             iFlag_zebra_in = 1,
                               iFlag_colorbar_in = 1,
                               iFlag_esri_hydro_image_in=0,
                             dData_max_in=dData_max,
                            dData_min_in=dData_min,
                            sColormap_in= 'terrain',
                            sField_color_in = 'mean',
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in =  'Surface elevation (mean)',
                             sUnit_in='Unit: m',
                             aLegend_in=aLegend_in,
                             aExtent_in=aExtent)