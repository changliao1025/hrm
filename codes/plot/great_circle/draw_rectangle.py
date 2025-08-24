import os, sys
from osgeo import ogr, gdal, osr
import cartopy.crs as ccrs
#draw a rectangle on the sphere

#define several points

from pyearth.system.define_global_variables import *
from pyflowline.classes.vertex import pyvertex


sPath_project = '/qfs/people/liao313/workspace/python/hrm'
#add the project path of the pythonpath
sys.path.append(sPath_project)
from codes.shared.map_great_circle_path import map_great_circle_path

pt1= dict()
pt1['dLongitude_degree'] = -130
pt1['dLatitude_degree'] = 40

pt2= dict()
pt2['dLongitude_degree'] = -120
pt2['dLatitude_degree'] = 40

pt3= dict()
pt3['dLongitude_degree'] = -120
pt3['dLatitude_degree'] = 50

pt4= dict()
pt4['dLongitude_degree'] = -130
pt4['dLatitude_degree'] = 50

pVertex1 = pyvertex(pt1)
pVertex2 = pyvertex(pt2)
pVertex3 = pyvertex(pt3)
pVertex4 = pyvertex(pt4)

#save as a polygon using geojson format
pDriver = ogr.GetDriverByName('GeoJSON')
sWorkspace_output = '/qfs/people/liao313/workspace/python/hrm/data/output'
sFilename = os.path.join(sWorkspace_output, 'rectangle.geojson')
iFlag_create_cell = 1
if iFlag_create_cell == 1:


    if os.path.exists(sFilename):
        pDriver.DeleteDataSource(sFilename)

    #define a spatial reference
    pSpatialRef = osr.SpatialReference()
    pSpatialRef.ImportFromEPSG(4326)


    pDataSource = pDriver.CreateDataSource(sFilename)
    pLayer = pDataSource.CreateLayer('rectangle', geom_type=ogr.wkbPolygon, srs=pSpatialRef)
    pFeature = ogr.Feature(pLayer.GetLayerDefn())
    pRing = ogr.Geometry(ogr.wkbLinearRing)
    pRing.AddPoint(pVertex1.dLongitude_degree, pVertex1.dLatitude_degree)
    pRing.AddPoint(pVertex2.dLongitude_degree, pVertex2.dLatitude_degree)
    pRing.AddPoint(pVertex3.dLongitude_degree, pVertex3.dLatitude_degree)
    pRing.AddPoint(pVertex4.dLongitude_degree, pVertex4.dLatitude_degree)
    pRing.AddPoint(pVertex1.dLongitude_degree, pVertex1.dLatitude_degree)
    pPolygon = ogr.Geometry(ogr.wkbPolygon)
    pPolygon.AddGeometry(pRing)
    pFeature.SetGeometry(pPolygon)
    pLayer.CreateFeature(pFeature)
    pDataSource.Destroy()

iFiletype_in=1
aExtent = None #[-180, 180, -90, 90]
sFilename_output_in = os.path.join(sWorkspace_output, 'rectangle.png')

aVertex = list()
aVertex.append(pVertex1)
aVertex.append(pVertex2)
aVertex.append(pVertex3)
aVertex.append(pVertex4)
map_great_circle_path(aVertex,
                          sFilename_output_in= sFilename_output_in,
                          iFlag_terrain_image_in = None,
                          iFlag_openstreetmap_level_in = None,
                          iFlag_scientific_notation_colorbar_in=None,
                          sTitle_in='Mesh cell',
                          iFlag_zebra_in = 1,
                          iFlag_fill_in = False,
                          iThickness_in = 1,
                          iFlag_color_in = 0,
                          iDPI_in=600,
                          aLegend_in=None,
                          aExtent_in =aExtent)
