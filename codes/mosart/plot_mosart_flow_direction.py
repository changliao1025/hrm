from pye3sm.mosart.map.structured.mosart_map_structured_flow_direction import mosart_map_flow_direction

sFilename_parameter_in = '/compyfs/inputdata/rof/mosart/MOSART_NLDAS_8th_20160426.nc'
sFilename_geojson_out = '/qfs/people/liao313/workspace/python/hrm/data/output/mosart/conus.geojson'
sFilename_png = '/qfs/people/liao313/workspace/python/hrm/figures/mosart/conus.png'

sFilename_parameter_in='/qfs/people/liao313/workspace/python/hrm/data/output/mosart/mosart_mississippi.nc'
sFilename_geojson_out='/qfs/people/liao313/workspace/python/hrm/data/output/mosart/mississippi.geojson'
sFilename_png='/qfs/people/liao313/workspace/python/hrm/figures/mosart/mississippi.png'

mosart_map_flow_direction(sFilename_parameter_in, sFilename_geojson_out, sFilename_png,
                           sRegion_in='Mississippi river basin', iFlag_esri_hydro_image_in =1)