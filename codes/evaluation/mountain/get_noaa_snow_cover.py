import datetime
import os
import requests
import gzip
import shutil
from osgeo import gdal, ogr, osr

# rasterio and shapely are required for clipping.
# You may need to install them: pip install rasterio shapely
from pyearth.toolbox.analysis.extract.clip_raster_by_polygon_file import clip_raster_by_polygon_file


def download_and_clip_noaa_ims_snow_cover(
    target_date: datetime.date,
    bounding_box_latlon: list[float],  # [min_lon, min_lat, max_lon, max_lat] in WGS84
    output_directory: str,
    clipped_filename: str = "snow_cover_clipped.tif",
    keep_intermediate_files: bool = False
) -> str | None:
    """
    Downloads NOAA IMS Daily Northern Hemisphere Snow and Ice Analysis (1km, G02156, V4)
    for a specific date, unzips it, and clips it to the given bounding box.

    Data Source: NSIDC - https://n5eil01u.ecs.nsidc.org/MEASURES/G02156.004/

    Args:
        target_date (datetime.date): The date for the snow cover data.
        bounding_box_latlon (list[float]): A list [min_lon, min_lat, max_lon, max_lat]
                                           defining the area of interest in WGS84 coordinates.
        output_directory (str): Directory to save the final clipped file and intermediate files.
        clipped_filename (str): Filename for the final clipped GeoTIFF.
        keep_intermediate_files (bool): If True, raw downloaded .gz and unzipped .tif files
                                        will be kept. Defaults to False.

    Returns:
        str | None: Full path to the clipped GeoTIFF file if successful, otherwise None.
    """


    year = target_date.year
    month = target_date.month
    day = target_date.day
    day_of_year = target_date.timetuple().tm_yday

    base_url = "https://n5eil01u.ecs.nsidc.org/MEASURES/G02156.004"
    raw_filename_gz = f"IMS1KM.{year}{day_of_year:03}.NHEMI.TIF.gz"
    unzipped_filename_tif = f"IMS1KM.{year}{day_of_year:03}.NHEMI.TIF"
    url_path_date = f"{year}.{month:02}.{day:02}"
    full_url = f"{base_url}/{url_path_date}/{raw_filename_gz}"

    os.makedirs(output_directory, exist_ok=True)
    raw_filepath_gz = os.path.join(output_directory, raw_filename_gz)
    unzipped_filepath_tif = os.path.join(output_directory, unzipped_filename_tif)
    clipped_filepath_tif = os.path.join(output_directory, clipped_filename)

    # --- 1. Download the full .tif.gz file ---
    print(f"Attempting to download snow cover data for {target_date.strftime('%Y-%m-%d')} from:\n{full_url}")
    try:
        response = requests.get(full_url, stream=True, timeout=60) # Increased timeout
        response.raise_for_status()
        with open(raw_filepath_gz, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Successfully downloaded: {raw_filepath_gz}")
    except requests.exceptions.HTTPError as e:
        print(f"HTTP error occurred: {e}")
        if response.status_code == 404:
            print(f"Data not found for {target_date.strftime('%Y-%m-%d')}. It might not be available.")
        elif response.status_code in [401, 403]:
            print("Authorization error. Check if access requires login (e.g., Earthdata Login).")
        return None
    except requests.exceptions.RequestException as e:
        print(f"An error occurred during download: {e}")
        return None

    # --- 2. Unzip the .tif.gz file ---
    print(f"Unzipping {raw_filepath_gz} to {unzipped_filepath_tif}")
    try:
        with gzip.open(raw_filepath_gz, 'rb') as f_in:
            with open(unzipped_filepath_tif, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Successfully unzipped: {unzipped_filepath_tif}")
    except Exception as e:
        print(f"Error unzipping file {raw_filepath_gz}: {e}")
        if not keep_intermediate_files and os.path.exists(raw_filepath_gz):
            os.remove(raw_filepath_gz)
        return None

    # --- 3. Clip the unzipped .tif file ---
    print(f"Clipping {unzipped_filepath_tif} to bounding box: {bounding_box_latlon}")


    # Create a temporary polygon file for the bounding box
    min_lon,  max_lon, min_lat, max_lat = bounding_box_latlon

    sFilename_polygon_in = os.path.join(output_directory, "bounding_box.geojson")
    if os.path.exists(sFilename_polygon_in):
        os.remove(sFilename_polygon_in)
    driver = ogr.GetDriverByName('GeoJSON')
    data_source = driver.CreateDataSource(sFilename_polygon_in)
    layer = data_source.CreateLayer("bounding_box", geom_type=ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(min_lon, min_lat)
    ring.AddPoint(max_lon, min_lat)
    ring.AddPoint(max_lon, max_lat)
    ring.AddPoint(min_lon, max_lat)
    ring.AddPoint(min_lon, min_lat)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(polygon)
    layer.CreateFeature(feature)
    feature.Destroy()
    data_source.Destroy()

    clip_raster_by_polygon_file(unzipped_filepath_tif, sFilename_polygon_in,
                                clipped_filename )

    # --- 4. Clean up intermediate files ---
    if not keep_intermediate_files:
        print("Cleaning up intermediate files...")
        if os.path.exists(raw_filepath_gz):
            os.remove(raw_filepath_gz)
            print(f"Removed: {raw_filepath_gz}")
        if os.path.exists(unzipped_filepath_tif):
            os.remove(unzipped_filepath_tif)
            print(f"Removed: {unzipped_filepath_tif}")

    return clipped_filepath_tif


if __name__ == '__main__':
    # Example usage:
    # Note: Data availability might vary. Check the NSIDC website for coverage.
    # A date too recent might not have data available yet.
    sample_date = datetime.date(2024, 12, 10)  # Example: Feb 10, 2023

    # Bounding box for Sierra Nevada
    dBuffer = 1.0
    #Sierra Nevada -119.85795,38.56839
    aExtent_sierranevada = [-119.85795 - dBuffer, -119.85795 + dBuffer, 38.56839 - dBuffer, 38.56839 + dBuffer]

    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    data_output_dir = '/qfs/people/liao313/workspace/python/hrm/data/evaluation/mountain'

    print(f"Example: Downloading and clipping data for {sample_date.strftime('%Y-%m-%d')}")
    print(f"Bounding Box (Lat/Lon WGS84): {aExtent_sierranevada}")
    print(f"Output directory: {data_output_dir}")
    clipped_file = download_and_clip_noaa_ims_snow_cover(
        target_date=sample_date,
        bounding_box_latlon=aExtent_sierranevada,
        output_directory=data_output_dir,
        clipped_filename="snow_cover_clipped_denver.tif",
        keep_intermediate_files=False
    )



