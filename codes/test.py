import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pyproj import Geod

# Define the polygon coordinates
polygon_coords = [
    (-130, 40),
    (-130, 50),
    (-120, 50),
    (-120, 40),
    (-130, 40)  # Closing the polygon
]

# Calculate the great circle path
geod = Geod(ellps="WGS84")
great_circle_path = []

for i in range(len(polygon_coords) - 1):
    lon1, lat1 = polygon_coords[i]
    lon2, lat2 = polygon_coords[i + 1]
    points = geod.npts(lon1, lat1, lon2, lat2, 100)  # 100 points along the path
    great_circle_path.extend(points)

# Plot the polygon
pProjection_map = ccrs.Robinson(central_longitude=0.0)
fig, ax = plt.subplots(subplot_kw={'projection': pProjection_map})
ax.set_extent([min(lon for lon, lat in polygon_coords), max(lon for lon, lat in polygon_coords),
               min(lat for lon, lat in polygon_coords), max(lat for lon, lat in polygon_coords)], crs=ccrs.PlateCarree())

lons, lats = zip(*great_circle_path)
ax.plot(lons, lats, transform=ccrs.Geodetic(), color='blue')
ax.plot(lons, lats, transform=ccrs.Geodetic(), color='blue')

ax.coastlines()
ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
# Save as a figure
sFilename = 'rectangle.png'
sWorkspace_output = '/qfs/people/liao313/workspace/python/hrm/data/output'

sFilename_out = os.path.join(sWorkspace_output, sFilename)
plt.savefig(sFilename_out, dpi=300)