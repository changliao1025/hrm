import os
import datetime
import textwrap
import numpy as np
from osgeo import osr, gdal, ogr
from urllib.error import URLError
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from pyearth.system.define_global_variables import *
from pyearth.visual.map.zebra_frame import zebra_frame
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter
from pyearth.visual.map.map_servers import StadiaStamen
from pyearth.visual.map.map_servers import StadiaStamen, EsriTerrain, EsriHydro, Stadia_terrain_images, Esri_terrain_images, Esri_hydro_images
from pyproj import Geod
import shapely.geometry as sgeom
#osr.UseExceptions()
#use agg and backend
#mpl.use('agg')
#get the current year for openstreetmap copy right label
iYear_current = datetime.datetime.now().year
sYear = str(iYear_current)

def map_great_circle_path(aVertex_in,
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFlag_color_in = None,
                            iFlag_colorbar_in=None,
                            iFlag_zebra_in=None,
                            iFlag_fill_in=None,
                            iFont_size_in=None,
                            iFlag_openstreetmap_in = None,
                            iFlag_openstreetmap_level_in = None,
                            iFlag_terrain_image_in = None,
                            iThickness_in=None,
                            sTitle_in=None,
                            iDPI_in=None,
                            iSize_x_in=None,
                            iSize_y_in=None,
                            sExtend_in=None,
                            sFont_in=None,
                            aLegend_in=None,
                            aExtent_in=None,
                            pProjection_map_in=None,
                            pProjection_data_in = None):
    """
    plot vector data on a map
    currently only support geojson and shapefile
    by default, the program will plot all the polygons in the file
    in the furture, the program will support to plot only a subset of polygons

    Args:
        iFiletype_in (_type_): File format, geojson etc
        sFilename_in (_type_): _description_
        sFilename_output_in (_type_): _description_
        iFlag_scientific_notation_colorbar_in (_type_, optional): _description_. Defaults to None.
        sColormap_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
        iDPI_in (_type_, optional): _description_. Defaults to None.
        dMissing_value_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        sExtend_in (_type_, optional): _description_. Defaults to None.
        sUnit_in (_type_, optional): _description_. Defaults to None.
        aLegend_in (_type_, optional): _description_. Defaults to None.
    """
    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    pSRS_geodetic = ccrs.Geodetic()


    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 150

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 8

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 8

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if iThickness_in is not None:
        iThickness = iThickness_in
    else:
        iThickness = 0.25

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0

    if iFlag_fill_in is not None:
        iFlag_fill = iFlag_fill_in
    else:
        iFlag_fill = True





    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title = 1
    else:
        iFlag_title = 0
        sTitle = ''
    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend = 'max'

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    aValue = list()


    nPoint = len(aVertex_in)

    #for j in range(nFeature):
        #pFeature = pLayer.GetFeature(j)

    geod = Geod(ellps="WGS84")
    pGreat_circle_path = []

    for i in range(nPoint):
        pVertex = aVertex_in[i]
        dLon_max = float(np.max([dLon_max, pVertex.dLongitude_degree]))
        dLon_min = float(np.min([dLon_min, pVertex.dLongitude_degree]))
        dLat_max = float(np.max([dLat_max, pVertex.dLatitude_degree]))
        dLat_min = float(np.min([dLat_min, pVertex.dLatitude_degree]))

    for i in range(nPoint-1):
        pVertex1 = aVertex_in[i]
        pVertex2 = aVertex_in[i+1]
        lon1, lat1 = pVertex1.dLongitude_degree, pVertex1.dLatitude_degree
        lon2, lat2 = pVertex2.dLongitude_degree, pVertex2.dLatitude_degree
        points = geod.npts(lon1, lat1, lon2, lat2, 1000)  # 100 points along the path
        pGreat_circle_path.extend(points)

    #add the last point
    pVertex1 = aVertex_in[nPoint-1]
    pVertex2 = aVertex_in[0]
    lon1, lat1 = pVertex1.dLongitude_degree, pVertex1.dLatitude_degree
    lon2, lat2 = pVertex2.dLongitude_degree, pVertex2.dLatitude_degree
    points = geod.npts(lon1, lat1, lon2, lat2, 1000)  # 100 points along the path
    pGreat_circle_path.extend(points)

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        #pProjection_map = ccrs.Orthographic(central_longitude=0.50*(
        #    dLon_max+dLon_min),  central_latitude=0.50*(dLat_max+dLat_min), globe=None)
        pProjection_map = ccrs.Orthographic(central_longitude=0.50*(
            dLon_max+dLon_min),  central_latitude=0.0, globe=None)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84


    #pProjection_map._threshold /= 1.0E6
    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map)




    if aExtent_in is None:
        marginx = (dLon_max - dLon_min) / 20
        marginy = (dLat_max - dLat_min) / 20
        if (dLat_max + marginy)> 90:
            dLat_max = 90
        else:
            dLat_max = dLat_max + marginy
        if (dLat_min - marginy) < -90:
            dLat_min = -90
        else:
            dLat_min = dLat_min - marginy
        if (dLon_max + marginx) > 180:
            dLon_max = 180
        else:
            dLon_max = dLon_max + marginx
        if (dLon_min - marginx) < -180:
            dLon_min = -180
        else:
            dLon_min = dLon_min - marginx
        aExtent = [dLon_min, dLon_max, dLat_min, dLat_max]
    else:
        aExtent = aExtent_in

    print(aExtent)
    minx,  maxx, miny, maxy = aExtent
    ax.set_extent(aExtent, crs = pSRS_wgs84)



    ax.coastlines(linewidth=0.5, color='k', resolution='10m')
    try:
        dAlpha = 1.0
        if iFlag_openstreetmap_in is not None and iFlag_openstreetmap_in == 1:
            from cartopy.io.img_tiles import OSM
            if iFlag_openstreetmap_level_in is not None:
                iFlag_openstreetmap_level = iFlag_openstreetmap_level_in
            else:
                iFlag_openstreetmap_level = 9
                pass

            osm_tiles = OSM()
            #Add the OSM image to the map
            ax.add_image(osm_tiles, iFlag_openstreetmap_level) #, alpha=0.5
            sLicense_info = "© OpenStreetMap contributors "+ sYear + "." + " Distributed under the Open Data Commons Open Database License (ODbL) v1.0."
            sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=60))
            ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                    color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

            #we also need to set transparency for the image to be added
            dAlpha = 0.5
        else:
            pass

        if iFlag_terrain_image_in is not None and iFlag_terrain_image_in == 1:
            if iFlag_openstreetmap_level_in is not None:
                iFlag_zoom_level = iFlag_openstreetmap_level_in
            else:
                iFlag_zoom_level = 8
                pass

            esri_terrain = EsriTerrain()
            esri_hydro = EsriHydro()
            ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
            multi_poly = esri_hydro.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
            target_domain = multi_poly.geoms[0]
            #_, aExtent_stamen, _ = stamen_terrain.image_for_domain(target_domain, iFlag_zoom_level)
            _, aExtent_terrain, _ = esri_terrain.image_for_domain(target_domain, iFlag_zoom_level)
            _, aExtent_hydro, _ = esri_hydro.image_for_domain(target_domain, iFlag_zoom_level)
            #img_stadia_terrain = Stadia_terrain_images(aExtent, iFlag_zoom_level)
            img_esri_terrain  = Esri_terrain_images(aExtent, iFlag_zoom_level)
            img_eari_hydro  = Esri_hydro_images(aExtent, iFlag_zoom_level)
            #ax.imshow(img_stadia_terrain,  extent=aExtent_stamen, transform=stamen_terrain.crs, alpha=0.8)
            ax.imshow(img_esri_terrain,  extent=aExtent_terrain, transform=esri_terrain.crs)
            ax.imshow(img_eari_hydro,  extent=aExtent_hydro, transform=esri_hydro.crs)
            #add the license information
            sLicense_info = "© Stamen Design, under a Creative Commons Attribution (CC BY 3.0) license."
            sLicense_info = "© Esri Hydro Reference Overlay"
            sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=60))
            ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                    color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
            dAlpha = 0.8
        else:
            pass
    except URLError as e:
        print('No internet connection')
        dAlpha = 1.0

    lons, lats = zip(*pGreat_circle_path)
    ax.plot(lons, lats, transform=ccrs.Geodetic(), color='blue')
    #use the origina path to draw the path
    for i in range(nPoint-1):
        pVertex1 = aVertex_in[i]
        pVertex2 = aVertex_in[i+1]
        lons = [pVertex1.dLongitude_degree, pVertex2.dLongitude_degree]
        lats = [pVertex1.dLatitude_degree, pVertex2.dLatitude_degree]
        ax.plot(lons, lats, transform=pSRS_wgs84, color='red')
    #draw the last line
    pVertex1 = aVertex_in[nPoint-1]
    pVertex2 = aVertex_in[0]
    lons = [pVertex1.dLongitude_degree, pVertex2.dLongitude_degree]
    lats = [pVertex1.dLatitude_degree, pVertex2.dLatitude_degree]
    ax.plot(lons, lats, transform=pSRS_wgs84, color='red')


    ax.set_extent(aExtent, crs = pSRS_wgs84)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx+(maxx-minx)/9, (maxx-minx)/8),
                      ylocs=np.arange(miny, maxy+(maxy-miny)/9, (maxy-miny)/8))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}

    if iFlag_zebra == 1:
        ax.set_xticks(np.arange(minx, maxx+(maxx-minx)/11, (maxx-minx)/10))
        dummy = (maxy-miny)/10
        ax.set_yticks(np.arange(miny, maxy+(maxy-miny)/11, (maxy-miny)/10))
        ax.set_axis_off()

    if iFlag_title is None:
        ax.set_title( sTitle )
    else:
        if iFlag_title==1:
            ax.set_title( sTitle )
        else:
            pass
        ax.set_title(sTitle)

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(0.05, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=iFont_size + 2)


    if iFlag_zebra ==1:
        #ax.set_axis_off()
        ax.set_extent(aExtent, crs = pSRS_wgs84)
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    ax.set_extent(aExtent, crs = pSRS_wgs84)





    if sFilename_output_in is None:
        plt.show()
    else:
        if os.path.exists(sFilename_output_in):
            os.remove(sFilename_output_in)

        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == '.png':
            plt.savefig(sFilename_out, bbox_inches='tight')
        else:
            if sExtension == '.pdf':
                plt.savefig(sFilename_out, bbox_inches='tight')
            else:
                plt.savefig(sFilename_out, bbox_inches='tight', format='ps')

        plt.close('all')
        plt.clf()
