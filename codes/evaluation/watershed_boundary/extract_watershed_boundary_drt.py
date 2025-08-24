import os, sys
from osgeo import gdal, ogr, osr
import numpy as np
from scipy import ndimage  # For morphological operations like filling holes
import netCDF4 as nc
from pyearth.system.define_global_variables import *
from pyflowline.classes.vertex import pyvertex
#can i use neighboring information to determind whethere a cell is at the edge of the watershed
from pyearth.gis.geometry.create_box_from_longitude_latitude import create_box_from_longitude_latitude

def fill_holes_and_find_edge(mask):
    """
    Fills holes in a binary NumPy mask and then finds its edge.

    Args:
        mask (np.ndarray): A 2D binary NumPy array (dtype=bool or 0/1).

    Returns:
        np.ndarray: A boolean NumPy array of the same shape as the input mask,
                    where True indicates an edge pixel and False otherwise.
    """
    if mask.ndim != 2:
        raise ValueError("Input mask must be a 2D array.")

    # Convert to boolean if it's not already
    mask = mask.astype(bool)

    # Fill holes
    filled_mask = ndimage.binary_fill_holes(mask)

    # Pad the filled mask to handle edges correctly
    padded_mask = np.pad(filled_mask, pad_width=1, mode='constant', constant_values=False)

    # Create a result array initialized to False
    edge = np.zeros_like(filled_mask, dtype=bool)

    # Iterate through the original filled mask (excluding padding)
    for r in range(filled_mask.shape[0]):
        for c in range(filled_mask.shape[1]):
            # Check if the current pixel is part of the object (in the filled mask)
            if padded_mask[r + 1, c + 1]:
                # Check its neighbors in the padded mask
                neighbors = [
                    padded_mask[r, c + 1],     # Top
                    padded_mask[r + 2, c + 1],   # Bottom
                    padded_mask[r + 1, c],     # Left
                    padded_mask[r + 1, c + 2],   # Right
                    # Optional: Include diagonals if you define edges that way
                    padded_mask[r, c],         # Top-Left
                    padded_mask[r, c + 2],       # Top-Right
                    padded_mask[r + 2, c],       # Bottom-Left
                    padded_mask[r + 2, c + 2]   # Bottom-Right
                ]
                # If any neighbor is False (background), then the current pixel is an edge
                if not all(neighbors):
                    edge[r, c] = True

    return edge

def extract_drt_watershed_boundary(sFilename_parameter_in,
                                           sFilename_geojson_out,
                                           sLengend_in=None,
                                           iSize_x_in = None,
                                           iSize_y_in = None,
                                            dData_max_in = None,
                                            dData_min_in = None,
                                            aLegend_in=None,
                                           aExtent_in=None):

    if os.path.exists(sFilename_parameter_in):
        print("Yep, I can read that file!")
    else:
        print("Nope, the path doesn't reach your file. Go research filepath in python")
        print(sFilename_parameter_in)


    print(sFilename_parameter_in)

    aDatasets = nc.Dataset(sFilename_parameter_in)

    netcdf_format = aDatasets.file_format

    #output file
    dResolution_degree = 1.0/16
    dResolution_degree_half = dResolution_degree / 2.0
    # Copy variables
    for sKey, aValue in aDatasets.variables.items():
        #print(sKey, aValue)
        #print(aValue.datatype)
        #print( aValue.dimensions)
        if sKey == 'ID':
            aID =  (aValue[:]).data
        if sKey == 'dnID':
            aDnID =  (aValue[:]).data
        if sKey == 'fdir':
            aFdir =  (aValue[:]).data
        if sKey == 'latixy':
            aLatitude = (aValue[:]).data
        if sKey == 'longxy':
            aLongitude = (aValue[:]).data
        if sKey == 'areaTotal2':
            aAccu = (aValue[:]).data
            aAccu = aAccu / 1.0E6

    #get max and min latitude and longitude
    dLongtitue_min = np.min(aLongitude)
    dLongtitue_max = np.max(aLongitude)
    dLatitude_min = np.min(aLatitude)
    dLatitude_max = np.max(aLatitude)

    nRow = int((dLatitude_max - dLatitude_min) / dResolution_degree) + 1
    nColumn = int((dLongtitue_max - dLongtitue_min) / dResolution_degree) + 1

    #create a 2D array to store
    aMask = np.zeros((nRow, nColumn), dtype=int)
    aID_new = np.zeros((nRow, nColumn), dtype=int)

    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.CreateDataSource(sFilename_geojson_out)
    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326) # WGS84 lat/long
    pLayer = pDataset.CreateLayer('flowdir', pSrs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    nPoint = aID.size
    for i in np.arange(0, nPoint, 1):
        lID = int(aID[i])
        x_start = float(aLongitude[i])
        y_start = float(aLatitude[i])
        pBox_start, aCoordinates_out = create_box_from_longitude_latitude(x_start,    y_start,
                                                         dResolution_degree,
                                                           dResolution_degree)
        #calcuate the index of the cell
        left, right, bottom, top = aCoordinates_out
        iRow = int((dLatitude_max - y_start) / dResolution_degree)
        iColumn = int((x_start - dLongtitue_min) / dResolution_degree)
        aMask[iRow, iColumn] = 1
        aID_new[iRow, iColumn] = lID
        pass

    #find edge of the binary mask
    aMask_edge = fill_holes_and_find_edge(aMask)

    #now we need to
    for iRow in range(0, nRow, 1):
        for iColumn in range( 0, nColumn, 1):
            if aMask_edge[iRow, iColumn]== True:
                #now use this index to create a polygon for the left cell
                dLongitude_left = dLongtitue_min - dResolution_degree_half + iColumn * dResolution_degree
                dLongitude_right = dLongtitue_min - dResolution_degree_half + (iColumn + 1) * dResolution_degree
                dLatitude_bottom = dLatitude_max + dResolution_degree_half - (iRow + 1) * dResolution_degree
                dLatitude_top = dLatitude_max + dResolution_degree_half - (iRow ) * dResolution_degree
                #create a polygon
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pLinearRing = ogr.Geometry(ogr.wkbLinearRing)
                pLinearRing.AddPoint( float(dLongitude_left),  float(dLatitude_bottom))
                pLinearRing.AddPoint( float(dLongitude_left),  float(dLatitude_top))
                pLinearRing.AddPoint( float(dLongitude_right), float(dLatitude_top))
                pLinearRing.AddPoint( float(dLongitude_right), float(dLatitude_bottom))
                pLinearRing.AddPoint( float(dLongitude_left),  float(dLatitude_bottom))
                pPolygon.AddGeometry(pLinearRing)
                #add the polygon to the layer
                pFeature.SetGeometry(pPolygon)
                #set the ID
                pFeature.SetField('id', int(aID_new[iRow, iColumn]))
                pLayer.CreateFeature(pFeature)

    #Save and close everything
    pDataset = pLayer = pFeature  = None

if __name__ == '__main__':
    sFilename_parameter = '/compyfs/liao313/00raw/mosart/mosart_extract_16th.nc'
    sFilename_geojson = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/watershed_boundary/drt_16th_watershed_boundary.geojson'
    extract_drt_watershed_boundary(sFilename_parameter, sFilename_geojson)
    pass


print('Finished')