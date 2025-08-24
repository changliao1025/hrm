import os, sys
from pathlib import Path
from os.path import realpath
import importlib.util
from osgeo import osr, ogr
from shutil import copy2
from pyhexwatershed.configuration.read_configuration_file import pyhexwatershed_read_configuration_file
from pyhexwatershed.configuration.change_json_key_value import change_json_key_value
from pyearth.gis.gdal.read.vector.gdal_get_vector_extent import gdal_get_vector_extent
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area

#===================================
#set up workspace path
#===================================
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
sPath_data = realpath( sPath_parent +  '/data/mississippi' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/mississippi'
sWorkspace_output_backup = sWorkspace_output
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)

#===================================
#you need to update this file based on your own case study
#===================================
sFilename_configuration_in = realpath( sPath_parent +  '/data/mississippi/input/pyhexwatershed_mississippi_mpas_03.json' )
sFilename_basins_in = realpath( sWorkspace_input +  '/pyhexwatershed_mississippi_basins.json' )
sFilename_jigsaw_in = realpath( sWorkspace_input +  '/pyhexwatershed_jigsaw.json' )
if os.path.isfile(sFilename_configuration_in):
    pass
else:
    print('This configuration does not exist: ', sFilename_configuration_in )

#===================================
#setup case information
#===================================
iFlag_create_job = 0
iFlag_visualization = 1
iCase_index = 4
sMesh_type = 'mpas'
sDate='20250501'

#===================================
#setup output and HPC job
#===================================
sSlurm = 'short'
sSlurm = 'slurm'
sFilename = sWorkspace_output + '/' + sMesh_type + '.bash'
ofs = open(sFilename, 'w')
sLine  = '#!/bin/bash' + '\n'
ofs.write(sLine)

#===================================
#visualization spatial extent
#mississippi:-124.5604166666668675,24.0020833333330117 : -66.0062500000011312,49.9979166666659651
#===================================
dBuffer = 1.0

sFilename_mesh_boudary = '/qfs/people/liao313/data/hexwatershed/mississippi/vector/mississippi_boundary.geojson'
aExtent_mississippi = gdal_get_vector_extent(sFilename_mesh_boudary)


oPyhexwatershed = pyhexwatershed_read_configuration_file(sFilename_configuration_in, \
    iCase_index_in=iCase_index, sDate_in=sDate, sWorkspace_output_in=sWorkspace_output)

# Set the basin outlet coordinates -89.26249,29.10721
dLongitude_outlet_degree = -89.26249
dLatitude_outlet_degree = 29.10721

sWorkspace_output = oPyhexwatershed.sWorkspace_output

#we want to copy the example configuration file to the output directory
sFilename_configuration_copy= os.path.join( sWorkspace_output, 'pyhexwatershed_configuration_copy.json' )
copy2(sFilename_configuration_in, sFilename_configuration_copy)

#copy the basin configuration file to the output directory as well
sFilename_configuration_basins_copy = os.path.join( sWorkspace_output, 'pyhexwatershed_configuration_basins_copy.json' )
copy2(sFilename_basins_in, sFilename_configuration_basins_copy)

sFilename_jigsaw_configuration_copy = os.path.join( sWorkspace_output, 'jigsaw_configuration_copy.json' )
#copy2(sFilename_jigsaw_in, sFilename_jigsaw_configuration_copy)

#now switch to the copied configuration file for modification
sFilename_configuration = sFilename_configuration_copy
sFilename_basins = sFilename_configuration_basins_copy
sFilename_jigsaw = sFilename_jigsaw_configuration_copy

change_json_key_value(sFilename_configuration, 'sWorkspace_output', sWorkspace_output_backup) #output folder
change_json_key_value(sFilename_configuration, 'sFilename_basins', sFilename_basins) #basin configuration file
#no need to rerun the jigsaw
#change_json_key_value(sFilename_configuration, 'sFilename_jigsaw_configuration', sFilename_jigsaw) #basin configuration file


oPyhexwatershed = pyhexwatershed_read_configuration_file(sFilename_configuration,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate,
                    sMesh_type_in = sMesh_type)

if iFlag_create_job == 1:
    oPyhexwatershed._pyhexwatershed_create_hpc_job(sSlurm_in = sSlurm, hours_in = 10 )
    print(iCase_index)
    sLine  = 'cd ' + oPyhexwatershed.sWorkspace_output + '\n'
    ofs.write(sLine)
    sLine  = 'sbatch submit.job' + '\n'
    ofs.write(sLine)
else:
    #oPyhexwatershed.pyhexwatershed_export()
    pass

aExtent = None
aExtent_lake_michigan = [-89.8, -82.3, 40.1, 44.4]
aExtent_jeffson = [-85.6, -78.7, 34.8, 38.9]

aExtent_debug = [-93.8, -93.7, 32.5, 32.6]
aExtent_outlet = [-92.45899893825303195,-89.79801124815780611, 29.5510526281927639,  31.39007657266963491]
aLegend = list()

aLegend.append('Unstructured MPAS mesh-based')


if iFlag_visualization == 1:
    sFilename = os.path.join(  oPyhexwatershed.sWorkspace_output_hexwatershed, 'priority_flood.mp4' )
    sTitle = 'Stream burning seeded priority-flood'
    #oPyhexwatershed._animate(sFilename, iFlag_type_in =1,iFigwidth_in=5, iFigheight_in=3, sTitle_in = sTitle)
    sFilename = os.path.join(  oPyhexwatershed.sWorkspace_output_hexwatershed, 'priority_flood_track.gif' )
    #oPyhexwatershed._animate(sFilename, iFlag_type_in =2,iFigwidth_in=5, iFigheight_in=7)

    sFilename =  'filtered_flowline.png'
    #oPyhexwatershed.plot(sFilename, sVariable_in = 'flowline_filter', aExtent_in =aExtent_full  )
    sFilename = os.path.join( oPyhexwatershed.sWorkspace_output_hexwatershed, 'slope.png' )
    #oPyhexwatershed.plot(sFilename_output_in=sFilename, sVariable_in = 'slope',
    #                     iFlag_colorbar_in=1,
    #                     iFlag_title_in =1,
    #                     dData_min_in = 0.0, dData_max_in = 0.1,
    #                     sFilename_boundary_in=sFilename_mesh_boudary,
    #                     aExtent_in=aExtent_mississippi,
    #                     aLegend_in = aLegend)

    sFilename =  'mpas_mesh.png'
    sFilename = os.path.join(  oPyhexwatershed.sWorkspace_output, 'mpas_mesh.png' )
    sTitle = 'MPAS mesh'
    #oPyhexwatershed.plot(sFilename, sVariable_in = 'mesh', aExtent_in =aExtent, sFilename_output_in = sFilename,
    #                  iDPI_in=300, sTitle_in = sTitle,
    #                  iFlag_esri_hydro_image_in=0,
    #                   iFlag_terrain_image_in=0)

    sVariable_in = 'flow_direction'
    sFilename = os.path.join( oPyhexwatershed.sWorkspace_output_hexwatershed, 'flow_direction.png' )
    oPyhexwatershed.plot(sFilename_output_in=sFilename, iFlag_title_in =1,
                         iFlag_esri_hydro_image_in=0,
                         iFlag_filter_in=0,
                         iFlag_arrow_in=0,
                         sVariable_in = sVariable_in,
                         aExtent_in=aExtent_mississippi,
                         sFilename_boundary_in=sFilename_mesh_boudary,
                         aLegend_in = aLegend)

    sFilename = os.path.join( oPyhexwatershed.sWorkspace_output_hexwatershed, 'drainage_area.png' )
    #oPyhexwatershed.plot(sFilename_output_in=sFilename, sVariable_in = 'drainage_area',
    #                     iFlag_colorbar_in=1,
    #                     iFlag_title_in =1,
    #                     dData_min_in = 0.0, dData_max_in = 3.2E12,
    #                     sFilename_boundary_in=sFilename_mesh_boudary,
    #                     aExtent_in=aExtent_mississippi,
    #                     aLegend_in = aLegend)

    sFilename =  'meander.png'
    #oPyhexwatershed.plot(sFilename, iFlag_title=0, sVariable_in='overlap',    aExtent_in =aExtent_meander )

    sFilename =  'braided.png'
    #oPyhexwatershed.plot(sFilename, iFlag_title=0, sVariable_in='overlap',    aExtent_in =aExtent_braided )

    sFilename =  'confluence.png'
    #oPyhexwatershed.plot(sFilename, iFlag_title=0, sVariable_in='overlap',    aExtent_in =aExtent_confluence )

    sFilename =  'outlet.png'
    #oPyhexwatershed.plot(sFilename, iFlag_title=0, sVariable_in='overlap',    aExtent_in =aExtent_outlet )

    sFilename =  os.path.join( oPyhexwatershed.sWorkspace_output_pyflowline, 'area_of_difference.pdf' )

    aLegend_in = list()
    aLegend_in.append('(a)')
    #oPyhexwatershed.plot(sFilename_output_in = sFilename,
    #                     iFlag_title_in=1,
    #                     iFlag_fill_in=1,
    #                     iFlag_esri_hydro_image_in =0,
    #                     sColor_in = 'red', iDPI_in=300,
    #                     sVariable_in = 'area_of_difference',
    #                     aLegend_in = aLegend_in,
    #                     aExtent_in=aExtent_mississippi )
    pass

ofs.close()
print('Finished')