from pyflowline.formats.read_flowline import read_flowline_geojson
sFilename = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501002/pyflowline/00000001/flowline_simplified.geojson'
sFilename = '/qfs/people/liao313/data/hexwatershed/mississippi/vector/flowline_hydroshed_simplified_1.0E4.geojson'
sFilename_out = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501004/pyflowline/flowline_filter.geojson'
aFlowline, pProjection_geojson = read_flowline_geojson(sFilename)

dLength = 0.0
for oFlowline in aFlowline:
    dLength += oFlowline.dLength

print('Total length of the river network: ', dLength)