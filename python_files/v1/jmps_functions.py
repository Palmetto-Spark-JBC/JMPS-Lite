"""---------------------------------------------------------------------------------------------------------------------
------------------------------------------------- JMPS v2.0 FUNCTIONS --------------------------------------------------
----------------------------------------------- Author: Drew Hollobaugh ------------------------------------------------
------------------------------------------------ Created: 6 April 2024 -------------------------------------------------
---------------------------------------------- Last Edited: 8 April 2024 -----------------------------------------------
---------------------------------------------------------------------------------------------------------------------"""


"""---------------------------------------------------------------------------------------------------------------------
---------------------------------------------------- MODULE IMPORTS ----------------------------------------------------
---------------------------------------------------------------------------------------------------------------------"""

from dted import LatLon, Tile
import geomag
import math
import matplotlib.pyplot as plt
import mgrs
import numpy as np
import openap as oap
from pathlib import Path
import rasterio
from rasterio.mask import mask
from shapely.geometry import Point, Polygon
import time


"""---------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------- CLASSES -------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------"""


class Waypoint(object):
    def __init__(self, name, _mgrs):
        self.name = name
        self.MGRS = _mgrs
        self.latlon = mgrs.MGRS().toLatLon(self.MGRS)
        # latitude ranges from 90S to 90N (y-axis)
        # longitude ranges from 180W to 180E (x-axis)
        self.closest_rad_dme = oap.nav.closest_fix(*self.latlon)

    def set_mgrs(self, _mgrs):
        self.MGRS = _mgrs
        self.latlon = mgrs.MGRS().toLatLon(self.MGRS)
        self.closest_rad_dme = oap.nav.closest_fix(*self.latlon)


"""---------------------------------------------------------------------------------------------------------------------
----------------------------------------------------- DTED FUNCTIONS ---------------------------------------------------
---------------------------------------------------------------------------------------------------------------------"""


def convert_meters_to_feet(meters):
    return meters * 3.2808416


def get_distance(p1, p2):
    """
    Get the distance in meters and nautical miles (nm) between two Waypoint Class objects or two tuples
    :param p1: Waypoint(object)
    :param p2: Waypoint(object)
    :return: (float) 42102.6, 22.72
    """
    if type(p1) == type(p2):
        if type(p1) == Waypoint:
            dist_meters = oap.aero.distance(*p1.latlon, *p2.latlon)  # in meters
        elif type(p1) == tuple:
            dist_meters = oap.aero.distance(*p1, *p2)  # in meters
        else:
            print(type(p1), type(p2))
            print(p1, p2)
            raise TypeError
    else:
        print(type(p1), type(p2))
        print(p1, p2)
        raise TypeError

    return round(dist_meters, 4), round(dist_meters * 0.0005396118, 4)


def get_true_mag_course(p1, p2):
    """
    Get the true course and mag course in between two Waypoint Class objects
    :param p1: Waypoint(object)
    :param p2: Waypoint(object)
    :return: 015.1, 022.2
    """
    if type(p1) == type(p2):
        if type(p1) == Waypoint:
            true_course = oap.aero.bearing(*p1.latlon, *p2.latlon)  # in true course
            mag_var_1 = geomag.declination(*p1.latlon)
            mag_var_2 = geomag.declination(*p2.latlon)
        elif type(p1) == tuple:
            true_course = oap.aero.bearing(*p1, *p2)  # in true course
            mag_var_1 = geomag.declination(*p1)
            mag_var_2 = geomag.declination(*p2)
        else:
            print(type(p1), type(p2))
            print(p1, p2)
            raise TypeError
    else:
        print(type(p1), type(p2))
        print(p1, p2)
        raise TypeError
    return round(true_course, 4), round(true_course - np.average([mag_var_1, mag_var_2]), 4)


def get_rectangle_points(p1, p2, corridor_size=5):
    tc, mc = get_true_mag_course(p1, p2)
    p1a = oap.aero.latlon(*p1.latlon, corridor_size * oap.aero.nm, (tc - 90 + 360) % 360)
    p1b = oap.aero.latlon(*p1.latlon, corridor_size * oap.aero.nm, (tc + 90 + 360) % 360)
    p2a = oap.aero.latlon(*p2.latlon, corridor_size * oap.aero.nm, (tc - 90 + 360) % 360)
    p2b = oap.aero.latlon(*p2.latlon, corridor_size * oap.aero.nm, (tc + 90 + 360) % 360)
    rect_pts = [(y, x) for x, y in [p1a, p1b, p2b, p2a, p1a]]
    rect_x = [j[0] for j in rect_pts]
    rect_y = [j[1] for j in rect_pts]
    plt.plot(rect_x, rect_y, lw=1)
    return rect_pts


def get_sub_points(p1, p2, increment=100):
    """
    Calculates the list of sub points on the path from p1 to p2, at the increment provided (in meters).
    :param p1: Way
    :param p2:
    :param increment:
    :return:
    """
    _distance_meters, _distance_nm = get_distance(p1, p2)
    print(_distance_meters, _distance_nm)
    tc, mc = get_true_mag_course(p1, p2)
    print(tc, mc)
    sub_points = [p1.latlon]
    for i in range(increment, int(_distance_meters-increment), increment):
        sub_points.append(oap.aero.latlon(*p1.latlon, i, tc))
    sub_points.append(p2.latlon)
    return sub_points


def is_point_inside_polygon(x, y, poly):
    """
    Determines if a point is inside a polygon.
    :param x: x-coordinate of the point
    :param y: y-coordinate of the point
    :param poly: list of (x, y) tuples representing the vertices of the polygon
    :return: (bool) True if the point is inside the polygon, False otherwise
    """
    x, y = float(x), float(y)
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def get_dted_filenames(_rect_points):
    """
    Calculates the folder and filenames needed for the given rectangle points
    :param _rect_points: (list) of (x, y) tuples. i.e: [(34.3, -78.5)]
    :return: set() of longitudes and set() of latitudes. i.e: {'w078'} {'n34'}
    latitude ranges from 90S to 90N (y-axis)
    longitude ranges from 180W to 180E (x-axis)
    """
    _latitudes = set()
    _longitudes = set()
    for point in _rect_points:
        x, y = point
        new_x = f"{'e' if x > 0 else 'w'}{abs(int(x-1)):0{3}}"
        new_y = f"{'n' if y > 0 else 's'}{abs(int(y)):0{2}}"
        _longitudes.add(new_x)
        _latitudes.add(new_y)
    return _longitudes, _latitudes


def find_highest_point(dted_file):
    """
    Finds the lat/lon coordinates and the elevation of the highest location within the DTED file.
    :param dted_file: i.e. n41.dt2
    :return: ((lat, lon), feet) i.e. ((34.98, -78.45), 88)
    latitude ranges from 90S to 90N (y-axis)
    longitude ranges from 180W to 180E (x-axis)
    """
    # Open the DTED file
    with rasterio.open(dted_file) as src:
        # Read the elevation data
        elevation = src.read(1)

        # Find the maximum elevation value
        max_elevation = elevation.max()

        # Find the index of the maximum elevation value
        max_index = np.unravel_index(np.argmax(elevation), elevation.shape)

        # Get the geographic coordinates of the highest point
        lon, lat = src.xy(max_index[0], max_index[1])

    return (lat, lon), max_elevation


def get_dted_corners(dted_file):
    """
    Returns the corners of the data contained in the dted file.
    :param dted_file: n35.dt2
    :return: ((-79.00, 35.00), (-77.99, 35.00), (-77.99, 33.99), (-79.00, 33.99), (-79.00, 35.00))
    """
    # Open the DTED file
    with rasterio.open(dted_file) as src:
        # Get the spatial extent (bounding box) of the DTED data
        bounds = src.bounds

        # Extract the coordinates of the four corners
        top_left = (bounds.left, bounds.top)
        top_right = (bounds.right, bounds.top)
        bottom_left = (bounds.left, bounds.bottom)
        bottom_right = (bounds.right, bounds.bottom)

    return top_left, top_right, bottom_right, bottom_left, top_left


def plot_terrain(x, y, z):
    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot surface
    surf = ax.plot_surface(x, y, z, cmap='viridis')

    # Add color bar
    fig.colorbar(surf)

    # Set labels and title
    # longitude ranges from 180W to 180E (x-axis)
    # latitude ranges from 90S to 90N (y-axis)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Elevation')
    ax.set_title('3D Plot with Elevation')

    # Show plot
    plt.show()


def find_highest_point_in_polygon(dted_file, rect_points):
    """
    Finds the lat/lon coordinates and the elevation of the highest location within the DTED file.
    :param dted_file: i.e. n41.dt2
    :param rect_points: <list> i.e. [(-78.6, 34.3), (-78.4, 34.2), (-78.3, 34.6), (-78.5, 34.7), (-78.6, 34.3)]
    :return: ((lat, lon), feet) i.e. ((34.4, -78.5) 59)
    latitude ranges from 90S to 90N (y-axis)
    longitude ranges from 180W to 180E (x-axis)
    """
    polygon_region = Polygon(rect_points)
    # Open the DTED file
    with rasterio.open(dted_file) as src:
        masked_data, masked_meta = mask(src, [polygon_region], crop=True)  # filter data by polygon region

        max_elevation = np.max(masked_data)  # get max elevation in filtered region
        _, height, width = masked_data.shape  # get shape of filtered data

        max_index = np.unravel_index(np.argmax(masked_data), (height, width))  # get index of max elevation point
        lon, lat = rasterio.transform.xy(masked_meta, max_index[0], max_index[1])  # Get the geographic coordinates of the highest point

    return (lat, lon), max_elevation


def run_leg_vertical_profile(p1, p2, _tfads_data, corridor_size=5):
    a = time.time()
    rect_points = get_rectangle_points(p1, p2, corridor_size=corridor_size)
    b = time.time()
    print(b-a)
    # Get Lat/Lon Coords and Elevation of the highest TFADS-O (obstacle) point:
    # TO DO: SPEED THIS UP, IT TAKES 4-5 SECONDS PER LEG WHICH IS WAYYYY TOO LONG
    max_obs_data = get_max_tfads_point_in_region(_tfads_data, rect_points)
    # print(max_obs_data[6], max_obs_data[5], max_obs_data[10])
    c = time.time()
    print(c - b)
    plt.plot(float(max_obs_data[6]), float(max_obs_data[5]), 'go')
    plt.text(float(max_obs_data[6]), float(max_obs_data[5]), max_obs_data[10])
    # plt.draw()

    # Get Lat/Lon Coords and Elevation of the highest DTED point:
    longitudes, latitudes = get_dted_filenames(rect_points)
    highest_elevation = -math.inf
    highest_point = None
    for folder in longitudes:
        for file in latitudes:
            dted_file = f"123086\\{folder}\\{file}.dt2"
            temp_highest_point, temp_highest_elevation = find_highest_point_in_polygon(dted_file, rect_points)
            if temp_highest_elevation > highest_elevation:
                highest_elevation = temp_highest_elevation
                highest_point = temp_highest_point
    highest_elevation = convert_meters_to_feet(highest_elevation)
    plt.plot(highest_point[1], highest_point[0], 'ro')
    plt.text(highest_point[1], highest_point[0], round(highest_elevation))
    # print(highest_point[1], highest_point[0], round(highest_elevation))
    plt.draw()
    return highest_point, highest_elevation


"""---------------------------------------------------------------------------------------------------------------------
--------------------------------------------------- TFADS-O FUNCTIONS --------------------------------------------------
---------------------------------------------------------------------------------------------------------------------"""


def get_max_tfads_point_in_region(_filtered_data, _region):
    sub_data = [i for i in _filtered_data if is_point_inside_polygon(i[6], i[5], _region)]  # takes 1s
    max_elevation_point = max(sub_data, key=lambda x: x[10])  # takes 0.0
    return max_elevation_point


def filter_tfads_data_by_country(_country='US'):
    tfads_file_name = f"TFADSP_LimDis-Filtered150+feet\\TFADSP_LimDis.txt"
    with open(tfads_file_name, "rb") as tfads_reader:
        tfads_data = tfads_reader.readlines()  # takes 0.45s
        sub_data = [i.decode("utf-8").split("\t") for i in tfads_data if i.decode("utf-8").split("\t")[2] == _country]  # takes 3.3s
        return sub_data

