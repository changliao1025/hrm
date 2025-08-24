import os
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_file_area


sFilename_out = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501004/pyflowline/00000001/area_of_difference.geojson'
dArea_inside = calculate_polygon_file_area( sFilename_out  )
print('aof:', dArea_inside)


#mpas