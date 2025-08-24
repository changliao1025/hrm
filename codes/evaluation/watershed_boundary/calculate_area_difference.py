import os
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_file_area

sWorkspace = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary'
#drt
sFilename_out = os.path.join(sWorkspace, 'drt_difference_inside.geojson')
dArea_inside = calculate_polygon_file_area( sFilename_out  )
print('drt area:', dArea_inside)
sFilename_out = os.path.join(sWorkspace, 'drt_difference_outside.geojson')
dArea_outside = calculate_polygon_file_area( sFilename_out  )
print('drt area:', dArea_outside)

sFilename_out = os.path.join(sWorkspace, 'mpas_difference_inside.geojson')
dArea_inside  = calculate_polygon_file_area( sFilename_out  )
print('mpas area:', dArea_inside)
sFilename_out = os.path.join(sWorkspace, 'mpas_difference_outside.geojson')
dArea_outside = calculate_polygon_file_area( sFilename_out  )
print('mpas area:', dArea_outside)

#mpas