import numpy as np
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import Polygon

from jmps_functions import Waypoint, get_distance, get_true_mag_course


def find_highest_point_in_polygon(dted_file, region_points):
    # Open the DTED file
    with rasterio.open(dted_file) as src:
        # Create a polygon from the region points
        polygon = Polygon(region_points)
        print(polygon)
        # Get the raster mask for the polygon
        mask = geometry_mask([polygon], out_shape=src.shape, transform=src.transform, invert=True)

        # Read the elevation data
        elevation = src.read(1, masked=True)
        print(elevation.shape)
        # Mask the elevation data
        masked_elevation = np.ma.masked_array(elevation, ~mask)
        print(masked_elevation.shape)
        # Find the highest point within the masked area
        highest_point = np.unravel_index(np.argmax(masked_elevation), masked_elevation.shape)
        print(highest_point)
        highest_elevation = masked_elevation[highest_point]
        print(highest_elevation)

        # Convert pixel coordinates to geographic coordinates
        lon, lat = src.xy(highest_point[0], highest_point[1])

    return lat, lon, highest_elevation


# Example usage:
dted_file = "C:\\Users\\drewh\\Documents\\Programming\\Python\\Palmetto Spark\\JMPS 2.0\\123086\\w078\\n34.dt2"
region_points = [
    (34.41178353127098, -78.95325998633476),
    (34.220878526601936, -78.09660954239314),
    (34.58636905551779, -77.97476195062806),
    (34.777274013206465, -78.83517633397904)
]
region_points = [
    (34.4, -78.7),
    (34.7, -78.6),
    (34.6, -78.2),
    (34.3, -78.3),
    (34.4, -78.7)
]


import matplotlib.pyplot as plt
plt.plot([i[1] for i in region_points], [i[0] for i in region_points])
plt.show()
# highest_point = find_highest_point_in_polygon(dted_file, region_points)
# print("Coordinates of the highest point:", highest_point)