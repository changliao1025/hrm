import numpy as np
from osgeo import gdal, osr
import math
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude

def get_raster_extent_and_crossing_length(sFilename_raster):
    """
    Calculate the line length needed to completely cross a raster file.

    Args:
        sFilename_raster (str): Path to the raster GeoTIFF file.
        center_point_wgs84 (tuple): Center point as (longitude, latitude) in WGS84.
        slope_degrees (float): Bearing in degrees (0° = North, 90° = East).

    Returns:
        tuple: (line_length_meters, raster_bounds_wgs84)
    """
    # Open the raster dataset
    pDataset_raster = gdal.Open(sFilename_raster, gdal.GA_ReadOnly)
    if pDataset_raster is None:
        raise FileNotFoundError(f"Could not open raster file: {sFilename_raster}")

    # Get raster georeferencing info
    gt = pDataset_raster.GetGeoTransform()
    proj = pDataset_raster.GetProjection()
    cols = pDataset_raster.RasterXSize
    rows = pDataset_raster.RasterYSize

    # Calculate raster bounds in raster CRS
    x_min = gt[0]
    x_max = gt[0] + cols * gt[1]
    y_max = gt[3]
    y_min = gt[3] + rows * gt[5]

    # Get raster CRS
    pSpatial_reference = osr.SpatialReference()
    pSpatial_reference.ImportFromWkt(proj)

    # Define WGS84 CRS
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.SetWellKnownGeogCS('WGS84')

    # Create coordinate transformation
    transform_to_wgs84 = osr.CoordinateTransformation(pSpatial_reference, wgs84_srs)

    # Transform raster corners to WGS84
    corners_raster = [
        (x_min, y_min),  # Bottom-left
        (x_max, y_min),  # Bottom-right
        (x_max, y_max),  # Top-right
        (x_min, y_max)   # Top-left
    ]

    corners_wgs84 = []
    for x, y in corners_raster:
        lon, lat, _ = transform_to_wgs84.TransformPoint(x, y)
        corners_wgs84.append((lon, lat))

    # Calculate raster bounds in WGS84
    lons = [corner[0] for corner in corners_wgs84]
    lats = [corner[1] for corner in corners_wgs84]

    raster_bounds_wgs84 = {
        'min_lon': min(lons),
        'max_lon': max(lons),
        'min_lat': min(lats),
        'max_lat': max(lats)
    }

    # Calculate diagonal distance of raster (maximum possible crossing distance)
    diagonal_distance = calculate_distance_based_on_longitude_latitude(
        raster_bounds_wgs84['min_lon'], raster_bounds_wgs84['min_lat'],
        raster_bounds_wgs84['max_lon'], raster_bounds_wgs84['max_lat']
    )

    # Add some buffer to ensure complete crossing (e.g., 20% extra)
    line_length_meters = diagonal_distance * 1.2

    pDataset_raster = None  # Close dataset

    return line_length_meters, raster_bounds_wgs84



def get_raster_cross_section_profile(sFilename_raster, center_point_wgs84, slope_degrees     ):
    """
    Get a complete cross-section profile that crosses the entire raster file.

    Args:
        sFilename_raster (str): Path to the raster GeoTIFF file.
        center_point_wgs84 (tuple): Center point as (longitude, latitude) in WGS84.
        slope_degrees (float): Bearing in degrees (0° = North, 90° = East).
        calculate_distance_based_on_longitude_latitude (callable): Your custom great circle distance function.

    Returns:
        tuple: Same as get_raster_pixels_intersecting_line_sphere function
    """

    # Calculate the line length needed to cross the entire raster
    line_length_meters, raster_bounds = get_raster_extent_and_crossing_length(
        sFilename_raster    )

    print(f"Raster bounds (WGS84): {raster_bounds}")
    print(f"Calculated line length to cross raster: {line_length_meters:.2f} meters")

    # Generate line coordinates
    line_coords_wgs84 = generate_line_from_center_and_slope(
        center_point_wgs84, slope_degrees, line_length_meters,
         nPoints=None
    )

    # Get the profile using the existing function
    return get_raster_pixels_intersecting_line_sphere(
        sFilename_raster, line_coords_wgs84   )


def clip_line_to_raster_bounds(line_coords_wgs84, raster_bounds_wgs84):
    """
    Clip the line coordinates to only include points within raster bounds.
    This is optional but can improve efficiency.

    Args:
        line_coords_wgs84 (list): List of (lon, lat) tuples
        raster_bounds_wgs84 (dict): Dictionary with min_lon, max_lon, min_lat, max_lat

    Returns:
        list: Clipped line coordinates
    """
    clipped_coords = []

    for lon, lat in line_coords_wgs84:
        if (raster_bounds_wgs84['min_lon'] <= lon <= raster_bounds_wgs84['max_lon'] and
            raster_bounds_wgs84['min_lat'] <= lat <= raster_bounds_wgs84['max_lat']):
            clipped_coords.append((lon, lat))

    return clipped_coords


# Modified version that automatically determines crossing length
def get_complete_raster_cross_section(sFilename_raster, center_point_wgs84, slope_degrees,
                                     calculate_distance_based_on_longitude_latitude, use_clipping=True):
    """
    Complete function that automatically calculates line length to cross entire raster.

    Args:
        sFilename_raster (str): Path to the raster GeoTIFF file.
        center_point_wgs84 (tuple): Center point as (longitude, latitude) in WGS84.
        slope_degrees (float): Bearing in degrees (0° = North, 90° = East).
        calculate_distance_based_on_longitude_latitude (callable): Your custom great circle distance function.
        use_clipping (bool): Whether to clip line to raster bounds for efficiency.

    Returns:
        tuple: (pixel_coords, geographic_coords, wgs84_coords, values, distances)
    """

    # Step 1: Calculate required line length
    line_length_meters, raster_bounds = get_raster_extent_and_crossing_length(
        sFilename_raster, center_point_wgs84, slope_degrees
    )

    print(f"Auto-calculated line length: {line_length_meters:.2f} meters")

    # Step 2: Generate line coordinates
    line_coords_wgs84 = generate_line_from_center_and_slope(
        center_point_wgs84, slope_degrees, line_length_meters,
        calculate_distance_based_on_longitude_latitude, nPoints=100  # Use enough points for good coverage
    )

    # Step 3: Optional clipping for efficiency
    if use_clipping:
        original_count = len(line_coords_wgs84)
        line_coords_wgs84 = clip_line_to_raster_bounds(line_coords_wgs84, raster_bounds)
        print(f"Clipped line from {original_count} to {len(line_coords_wgs84)} points")

    # Step 4: Get the profile
    return get_raster_pixels_intersecting_line_sphere(
        sFilename_raster, line_coords_wgs84, calculate_distance_based_on_longitude_latitude
    )

def generate_line_from_center_and_slope(center_point_wgs84, slope_degrees, line_length_meters,
                                        nPoints=None):
    """
    Generate line coordinates from center point and slope using great circle calculations.

    Args:
        center_point_wgs84 (tuple): Center point as (longitude, latitude) in WGS84.
        slope_degrees (float): Bearing in degrees (0° = North, 90° = East).
        line_length_meters (float): Total length of the line in meters.
        great_circle_distance_func (callable): Your custom great circle distance function.
        nPoints (int, optional): Number of points to generate along the line.

    Returns:
        list: List of (longitude, latitude) tuples defining the line.
    """

    if nPoints is None:
        # Default: one point per ~50 meters, minimum 10 points
        nPoints = max(10, int(line_length_meters / 50))

    center_lon, center_lat = center_point_wgs84
    half_length = line_length_meters / 2

    # Convert bearing to radians (geographic bearing: 0° = North)
    bearing_rad = math.radians(slope_degrees)

    # Calculate start and end points using spherical geometry
    start_point = calculate_destination_point_sphere(center_lon, center_lat,
                                                   bearing_rad + math.pi, half_length)  # Opposite direction
    end_point = calculate_destination_point_sphere(center_lon, center_lat,
                                                 bearing_rad, half_length)

    # Generate intermediate points along the great circle
    line_coords_wgs84 = []

    if nPoints == 1:
        line_coords_wgs84.append(center_point_wgs84)
    elif nPoints == 2:
        line_coords_wgs84.extend([start_point, end_point])
    else:
        # Use great circle interpolation for more than 2 points
        for i in range(nPoints):
            t = i / (nPoints - 1)  # Parameter from 0 to 1

            # Spherical linear interpolation (slerp) along great circle
            point = interpolate_great_circle(start_point, end_point, t)
            line_coords_wgs84.append(point)

    return line_coords_wgs84


def calculate_destination_point_sphere(lon_deg, lat_deg, bearing_rad, distance_meters):
    """
    Calculate destination point given start point, bearing and distance on a sphere.

    Args:
        lon_deg, lat_deg: Starting point in decimal degrees
        bearing_rad: Bearing in radians (0 = North, π/2 = East)
        distance_meters: Distance in meters

    Returns:
        tuple: (longitude, latitude) of destination point in decimal degrees
    """
    # Earth radius in meters (you might want to use the same value as your distance function)
    R = 6371000.0  # WGS84 mean radius

    # Convert to radians
    lat1_rad = math.radians(lat_deg)
    lon1_rad = math.radians(lon_deg)

    # Angular distance
    d_over_R = distance_meters / R

    # Calculate destination point
    lat2_rad = math.asin(
        math.sin(lat1_rad) * math.cos(d_over_R) +
        math.cos(lat1_rad) * math.sin(d_over_R) * math.cos(bearing_rad)
    )

    lon2_rad = lon1_rad + math.atan2(
        math.sin(bearing_rad) * math.sin(d_over_R) * math.cos(lat1_rad),
        math.cos(d_over_R) - math.sin(lat1_rad) * math.sin(lat2_rad)
    )

    # Convert back to degrees
    return (math.degrees(lon2_rad), math.degrees(lat2_rad))


def interpolate_great_circle(point1, point2, t):
    """
    Interpolate between two points along a great circle.

    Args:
        point1, point2: (longitude, latitude) tuples in decimal degrees
        t: Interpolation parameter (0.0 = point1, 1.0 = point2)

    Returns:
        tuple: Interpolated (longitude, latitude) in decimal degrees
    """
    lon1, lat1 = point1
    lon2, lat2 = point2

    # Convert to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Calculate angular distance between points
    d_lon = lon2_rad - lon1_rad

    a = (math.sin(lat1_rad) * math.sin(lat2_rad) +
         math.cos(lat1_rad) * math.cos(lat2_rad) * math.cos(d_lon))

    # Avoid numerical errors
    a = max(-1.0, min(1.0, a))
    d = math.acos(a)

    if abs(d) < 1e-10:  # Points are very close
        # Linear interpolation for very close points
        lat_interp = lat1 + t * (lat2 - lat1)
        lon_interp = lon1 + t * (lon2 - lon1)
        return (lon_interp, lat_interp)

    # Spherical interpolation
    A = math.sin((1 - t) * d) / math.sin(d)
    B = math.sin(t * d) / math.sin(d)

    x = A * math.cos(lat1_rad) * math.cos(lon1_rad) + B * math.cos(lat2_rad) * math.cos(lon2_rad)
    y = A * math.cos(lat1_rad) * math.sin(lon1_rad) + B * math.cos(lat2_rad) * math.sin(lon2_rad)
    z = A * math.sin(lat1_rad) + B * math.sin(lat2_rad)

    lat_interp_rad = math.atan2(z, math.sqrt(x*x + y*y))
    lon_interp_rad = math.atan2(y, x)

    return (math.degrees(lon_interp_rad), math.degrees(lat_interp_rad))


def get_raster_pixels_intersecting_line_sphere(sFilename_raster, line_coords_wgs84):
    """
    Extracts all raster pixel values that intersect with a given line on a sphere
    using custom great circle distance calculations.

    Args:
        sFilename_raster (str): Path to the raster GeoTIFF file.
        line_coords_wgs84 (list of tuples): A list of (longitude, latitude) tuples
                                           defining the vertices of the line in WGS84.
        great_circle_distance_func (callable): Your custom function that takes two points
                                              and returns distance. Should accept:
                                              great_circle_distance_func(lon1, lat1, lon2, lat2)
                                              and return distance in meters.

    Returns:
        tuple: A tuple containing:
               - pixel_coords (list of tuples): (row, col) pixel coordinates
               - geographic_coords (list of tuples): (x, y) coordinates in raster CRS
               - wgs84_coords (list of tuples): (lon, lat) coordinates in WGS84
               - values (np.array): Raster values at each intersecting pixel
               - distances (np.array): Cumulative great circle distances along the line (in meters)
    """

    # Open the raster dataset
    pDataset_raster = gdal.Open(sFilename_raster, gdal.GA_ReadOnly)
    if pDataset_raster is None:
        raise FileNotFoundError(f"Could not open raster file: {sFilename_raster}")

    # Get raster band and georeferencing info
    band = pDataset_raster.GetRasterBand(1)
    gt = pDataset_raster.GetGeoTransform()
    proj = pDataset_raster.GetProjection()

    # Get raster dimensions
    cols_raster = pDataset_raster.RasterXSize
    rows_raster = pDataset_raster.RasterYSize

    # Get raster CRS
    pSpatial_reference = osr.SpatialReference()
    pSpatial_reference.ImportFromWkt(proj)

    # Define WGS84 CRS for input coordinates
    wgs84_srs = osr.SpatialReference()
    wgs84_srs.SetWellKnownGeogCS('WGS84')

    # Create coordinate transformation objects
    transform_to_raster_crs = osr.CoordinateTransformation(wgs84_srs, pSpatial_reference)
    transform_to_wgs84_crs = osr.CoordinateTransformation(pSpatial_reference, wgs84_srs)

    # Read the entire raster array
    raster_array = band.ReadAsArray()
    nodata_value = band.GetNoDataValue()

    def world_to_pixel(x_world, y_world):
        """Convert world coordinates to pixel coordinates"""
        det = gt[1] * gt[5] - gt[2] * gt[4]
        pixel = (x_world * gt[5] - y_world * gt[2] - gt[0] * gt[5] + gt[3] * gt[2]) / det
        line = (y_world * gt[1] - x_world * gt[4] - gt[3] * gt[1] + gt[0] * gt[4]) / det
        return int(round(pixel)), int(round(line))

    def pixel_to_world(col, row):
        """Convert pixel coordinates to world coordinates (center of pixel)"""
        x_world = gt[0] + col * gt[1] + row * gt[2]
        y_world = gt[3] + col * gt[4] + row * gt[5]
        return x_world, y_world

    def bresenham_line(x0, y0, x1, y1):
        """Bresenham's line algorithm to get all pixels along a line"""
        pixels = []

        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx - dy

        x, y = x0, y0

        while True:
            pixels.append((x, y))

            if x == x1 and y == y1:
                break

            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                x += sx
            if e2 < dx:
                err += dx
                y += sy

        return pixels

    # Convert line coordinates to raster CRS and then to pixel coordinates
    all_pixel_coords = set()  # Use set to avoid duplicates

    for i in range(len(line_coords_wgs84) - 1):
        # Get start and end points of current segment
        lon1, lat1 = line_coords_wgs84[i]
        lon2, lat2 = line_coords_wgs84[i + 1]

        # Transform to raster CRS
        x1_raster, y1_raster, _ = transform_to_raster_crs.TransformPoint(lon1, lat1)
        x2_raster, y2_raster, _ = transform_to_raster_crs.TransformPoint(lon2, lat2)

        # Convert to pixel coordinates
        col1, row1 = world_to_pixel(x1_raster, y1_raster)
        col2, row2 = world_to_pixel(x2_raster, y2_raster)

        # Get all pixels along this line segment
        segment_pixels = bresenham_line(col1, row1, col2, row2)
        all_pixel_coords.update(segment_pixels)

    # Filter pixels that are within raster bounds
    valid_pixels = []
    for col, row in all_pixel_coords:
        if 0 <= row < rows_raster and 0 <= col < cols_raster:
            valid_pixels.append((col, row))

    # Sort pixels by their great circle distance from the first point
    if valid_pixels:
        first_lon, first_lat = line_coords_wgs84[0]

        def great_circle_distance_from_start(pixel):
            col, row = pixel
            x_world, y_world = pixel_to_world(col, row)
            # Transform pixel center back to WGS84 for distance calculation
            lon_wgs84, lat_wgs84, _ = transform_to_wgs84_crs.TransformPoint(x_world, y_world)
            return calculate_distance_based_on_longitude_latitude(first_lon, first_lat, lon_wgs84, lat_wgs84)

        valid_pixels.sort(key=great_circle_distance_from_start)

    # Extract data for each valid pixel
    pixel_coords = []
    geographic_coords = []
    wgs84_coords = []
    values = []
    distances = []

    cumulative_distance = 0.0
    prev_lon, prev_lat = None, None

    for col, row in valid_pixels:
        # Get pixel value
        pixel_value = raster_array[row, col]
        if nodata_value is not None and pixel_value == nodata_value:
            pixel_value = np.nan

        # Get world coordinates (center of pixel)
        x_world, y_world = pixel_to_world(col, row)

        # Transform back to WGS84
        lon_wgs84, lat_wgs84, _ = transform_to_wgs84_crs.TransformPoint(x_world, y_world)

        # Calculate cumulative great circle distance
        if prev_lon is not None:
            segment_dist = calculate_distance_based_on_longitude_latitude(prev_lon, prev_lat, lon_wgs84, lat_wgs84)
            cumulative_distance += segment_dist

        pixel_coords.append((row, col))
        geographic_coords.append((x_world, y_world))
        wgs84_coords.append((lon_wgs84, lat_wgs84))
        values.append(pixel_value)
        distances.append(cumulative_distance)

        prev_lon, prev_lat = lon_wgs84, lat_wgs84

    pDataset_raster = None  # Close the dataset

    return pixel_coords, geographic_coords, wgs84_coords, np.array(values), np.array(distances)


if __name__ == "__main__":
    from pyearth.gis.gdal.read.raster.gdal_get_raster_extent import gdal_get_raster_extent
    #test the function with a sample raster file
    sFilename_raster = '/compyfs/liao313/04model/subset/mississippi/elevation/river/000263140/dem_box.tif'
    aExtent = gdal_get_raster_extent(sFilename_raster)
    #get center
    center_lon = (aExtent[0] + aExtent[1]) / 2
    center_lat = (aExtent[2] + aExtent[3]) / 2
    center_point_wgs84 = (center_lon, center_lat)

    slope = 45  # degrees, northeast direction
    #call the function
    pixel_coords, geographic_coords, wgs84_coords, values, distances = get_raster_cross_section_profile(sFilename_raster, center_point_wgs84, slope)

    #plot the profile
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(distances, values, marker='o', linestyle='-', color='b')
    plt.title('Raster Cross-Section Profile')
    plt.xlabel('Cumulative Distance (meters)')
    plt.ylabel('Raster Value')
    plt.grid()
    #save as a png
    sFilename_png = '/qfs/people/liao313/workspace/python/hrm/figures/elevation/elevation_profile2.png'
    plt.savefig(sFilename_png, dpi=300)
