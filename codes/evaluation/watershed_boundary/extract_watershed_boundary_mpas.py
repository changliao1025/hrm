import os
import json
from osgeo import gdal, ogr, osr
import numpy as np
from pyearth.system.define_global_variables import *
from pyflowline.classes.vertex import pyvertex
#can i use neighboring information to determind whethere a cell is at the edge of the watershed


sFilename_mesh_info = '/compyfs/liao313/04model/pyhexwatershed/mississippi/pyhexwatershed20250501003/pyflowline/mpas_mesh_info.json'

pDriver = ogr.GetDriverByName('GeoJSON')

sWorkspace_output = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary'
sFilename_edge = os.path.join(sWorkspace_output, 'edge_cell.geojson')
if os.path.exists(sFilename_edge):
    os.remove(sFilename_edge)

pSpatialRef = osr.SpatialReference()
pSpatialRef.ImportFromEPSG(4326)

pDataset = pDriver.CreateDataSource(sFilename_edge)
pLayer = pDataset.CreateLayer('edge_cell', pSpatialRef, ogr.wkbPolygon)
pLayer.CreateField(ogr.FieldDefn('ID', ogr.OFTInteger))

with open(sFilename_mesh_info) as json_file:
    data = json.load(json_file)
    ncell = len(data)
    lID = 1
    for i in range(ncell):
        pcell = data[i]
        nNeighbor = pcell['nNeighbor_land']
        nVertex = pcell['nVertex']
        lCellID = pcell['lCellID']

        if nNeighbor == nVertex:
            pass
        else:
            if nNeighbor < nVertex:
                #this is a cell at the edge of the watershed and we will save it

                #get all the vertex of the cell
                aVertex = pcell['aVertex']
                vVertex = list()
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pLinearRing = ogr.Geometry(ogr.wkbLinearRing)
                #do we assume that the point is CCW?
                #create a polygon from the vertex
                for j in range(nVertex):
                    point = dict()
                    point['dLongitude_degree'] = aVertex[j]['dLongitude_degree']
                    point['dLatitude_degree'] = aVertex[j]['dLatitude_degree']
                    pVertex = pyvertex(point)
                    vVertex.append(pVertex)
                    pLinearRing.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)
                    pass

                #add the first point to the end of the list
                pLinearRing.AddPoint(vVertex[0].dLongitude_degree, vVertex[0].dLatitude_degree)
                pPolygon.AddGeometry(pLinearRing)
                pFeature = ogr.Feature(pLayer.GetLayerDefn())
                pFeature.SetGeometry(pPolygon)
                pFeature.SetField('ID', lCellID)
                pLayer.CreateFeature(pFeature)
                pFeature.Destroy()
                lID += 1

                pass
            else:
                print("error")
                pass

pDataset = None
pDriver = None

print('Finished', lID)