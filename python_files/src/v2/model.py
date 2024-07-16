
import geomag
import geopandas as gpd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import matplotlib
import math
import mgrs
import numpy as np
import openap as oap
import pandas as pd
import rasterio
from rasterio.mask import mask
from shapely.geometry import Point, Polygon
import time
from tkinter import Label, Frame, Tk, StringVar, Entry, END, ttk
import tkinter as tk
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


def is_valid_mgrs(mgrs_coord):
    m = mgrs.MGRS()
    try:
        # Try to convert MGRS to lat/lon
        lat_lon = m.toLatLon(mgrs_coord)
        return True
    except (ValueError, mgrs.core.MGRSError):
        # If there is a ValueError, the MGRS coordinate is invalid
        return False


def convert_meters_to_nm(dist_meters):
    return round(dist_meters * 0.0005396118, 1)


def get_distance_nm(p1, p2):
    p1_latlong = mgrs.MGRS().toLatLon(p1)
    p2_latlong = mgrs.MGRS().toLatLon(p2)
    dist_meters = oap.aero.distance(*p1_latlong, *p2_latlong)  # in meters
    return convert_meters_to_nm(dist_meters)


def convert_meters_to_feet(meters):
    return meters * 3.2808416


def get_true_mag_course(p1, p2):
    if type(p1) == str:
        p1_latlong = mgrs.MGRS().toLatLon(p1)
    elif type(p1) == tuple:
        p1_latlong = p1
    else:
        raise TypeError(p1)

    if type(p2) == str:
        p2_latlong = mgrs.MGRS().toLatLon(p2)
    elif type(p2) == tuple:
        p2_latlong = p2
    else:
        raise TypeError(p2)

    try:
        true_course = oap.aero.bearing(*p1_latlong, *p2_latlong)  # in true course
        mag_var_1 = geomag.declination(*p1_latlong), 1
        mag_var_2 = geomag.declination(*p2_latlong), 1
        return round(true_course,1), round(true_course - np.average([mag_var_1, mag_var_2]),1)
    except mgrs.core.MGRSError as e:
        print(e)
        print(f"Unable to compute true/mag course for points:\np1: {p1}\np2: {p2}")
        return None, None


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


def get_max_tfads_point_in_region(_filtered_data, _region):
    sub_data = [i for i in _filtered_data if is_point_inside_polygon(i[6], i[5], _region)]  # takes 1s
    max_elevation_point = max(sub_data, key=lambda x: x[10])  # takes 0.0
    return max_elevation_point


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
    print(b - a)
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


def filter_tfads_data_by_country(_country='US'):
    tfads_file_name = f"TFADSP_LimDis-Filtered150+feet\\TFADSP_LimDis.txt"
    with open(tfads_file_name, "rb") as tfads_reader:
        tfads_data = tfads_reader.readlines()  # takes 0.45s
        sub_data = [i.decode("utf-8").split("\t") for i in tfads_data if i.decode("utf-8").split("\t")[2] == _country]  # takes 3.3s
        return sub_data


def get_closest_rad_dme(_lat, _lon):
    closest, distance_meters = oap.nav.closest_fix(_lat, _lon)
    true, mag = get_true_mag_course((closest[0], closest[1]), (_lat, _lon))
    return f"{closest[-1]}/R{mag:03.0f}/{convert_meters_to_nm(distance_meters):03.0f}"


class JMPS_GUI(Tk):
    def __init__(self):
        super().__init__()
        self._title = 'JMPS Lite'
        self.title(self._title)
        self.state('zoomed')

        self.point_headers1 = ["Point", "MGRS", "Latitude", "Mag Course", "Distance (NM)"]
        self.point_headers2 = ["", "Radial/DME", "Longitude", "True Course", ""]
        self.points_rows = []
        self.protocol("WM_DELETE_WINDOW", self.on_closing)  # Bind the close event
        self.initialize_page()

    def on_closing(self):
        self.destroy()  # Destroy the window, which ends the program

    def initialize_page(self):
        self.frame = ttk.Frame(self)
        self.frame.pack(fill=tk.BOTH, expand=1)

        self.left_frame = ttk.Frame(self.frame)
        self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        self.right_frame = ttk.Frame(self.frame)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1)

        self.initialize_menu()
        # self.initialize_toolbar()
        self.initialize_tabular_window()  # Initialize GUI input/display boxes

        # This defines the Python GUI backend to use for matplotlib
        matplotlib.use('TkAgg')

        # Initialize matplotlib figure for graphing purposes
        self.fig = plt.figure(1)

        # Special type of "canvas" to allow for matplotlib graphing
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.left_frame)
        self.canvas.draw()
        self.plot_widget = self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)
        # self.plot_widget.grid(row=0, column=0, rowspan=10)

        self.mainloop()

    def initialize_tabular_window(self):
        """Create the rows/labels/inputs on right side of app."""
        title = Label(self.right_frame, text=self._title, font=('Helvetica', 20), bd=10)
        title.grid(row=0, column=1, columnspan=6)

        for i, item in enumerate(self.point_headers1):
            label = Label(self.right_frame, text=item, bd=5)
            label.grid(row=2, column=2 + i, columnspan=1)

        for i, item in enumerate(self.point_headers2):
            label = Label(self.right_frame, text=item, bd=5)
            label.grid(row=3, column=2 + i, columnspan=1)

        self.row_start = 4
        self.add_point_row(self.right_frame, self.row_start)

    def add_point_row(self, frame, row, col=1):
        padx = 10
        pady = 2
        bd = 2
        point_number_label = Label(frame, text=f"{len(self.points_rows)+1}.")
        point_number_label.grid(row=row, column=col, padx=padx, pady=pady)

        name_input_var = StringVar()
        name_input_var.trace("w", lambda name, index, mode, var=name_input_var: self.point_callback(name_input_var))
        name_input = Entry(frame, textvariable=name_input_var, width=10, bd=bd)
        name_input.grid(row=row, column=col+1, padx=padx)

        mgrs_input_var = StringVar()
        mgrs_input_var.trace("w", lambda name, index, mode, var=mgrs_input_var: self.point_callback(mgrs_input_var))
        mgrs_input = Entry(frame, textvariable=mgrs_input_var, width=20, bd=bd)
        mgrs_input.grid(row=row, column=col+2, padx=padx)

        latitude_display = Entry(frame, width=20, bd=bd)
        latitude_display.grid(row=row, column=col + 3, padx=padx)
        latitude_display.config(state="readonly")

        mag_course_display = Entry(frame, width=7, bd=bd)
        mag_course_display.grid(row=row, column=col + 4, padx=padx)
        mag_course_display.config(state="readonly")

        distance_display = Entry(frame, width=7, bd=bd)
        distance_display.grid(row=row, column=col + 5, padx=padx)
        distance_display.config(state="readonly")

        rad_dme_display = Entry(frame, width=20, bd=bd)
        rad_dme_display.grid(row=row + 1, column=col + 2, padx=padx)
        rad_dme_display.config(state="readonly")

        longitude_display = Entry(frame, width=20, bd=bd)
        longitude_display.grid(row=row + 1, column=col + 3, padx=padx)
        longitude_display.config(state="readonly")

        true_course_display = Entry(frame, width=7, bd=bd)
        true_course_display.grid(row=row + 1, column=col + 4, padx=padx)
        true_course_display.config(state="readonly")

        # Add a horizontal separator
        separator = ttk.Separator(frame, orient='horizontal')
        separator.grid(row=row + 2, column=col, columnspan=6, sticky="ew", pady=10)

        self.points_rows.append([point_number_label, name_input, mgrs_input, latitude_display, mag_course_display, distance_display, rad_dme_display, longitude_display, true_course_display, separator])

    def initialize_menu(self):
        def donothing():
            print('clicked')
        menubar = tk.Menu(self.frame)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="New", command=donothing, accelerator="Ctrl+N")
        filemenu.add_command(label="Open", command=donothing, accelerator="Ctrl+O")
        filemenu.add_command(label="Save", command=donothing, accelerator="Ctrl+S")
        filemenu.add_command(label="Save as...", command=donothing, accelerator="Ctrl+Shift+S")
        filemenu.add_command(label="Close", command=donothing, accelerator="Ctrl+W")

        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.quit)
        menubar.add_cascade(label="File", menu=filemenu)
        editmenu = tk.Menu(menubar, tearoff=0)
        editmenu.add_command(label="Undo", command=donothing)
        editmenu.add_separator()
        editmenu.add_command(label="Cut", command=donothing)
        editmenu.add_command(label="Copy", command=donothing)
        editmenu.add_command(label="Paste", command=donothing)
        editmenu.add_command(label="Delete", command=donothing)
        editmenu.add_command(label="Select All", command=donothing)

        menubar.add_cascade(label="Edit", menu=editmenu)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Help Index", command=donothing)
        helpmenu.add_command(label="About...", command=donothing)
        menubar.add_cascade(label="Help", menu=helpmenu)

        self.config(menu=menubar)

    def initialize_toolbar(self):
        toolbar = self.right_frame
        button1 = tk.Button(toolbar, text='Button 1', font=('Helvetica', 20))
        button1.grid(row=5, column=2, columnspan=1)

        button2 = tk.Button(toolbar, text='Button 2', font=('Helvetica', 20))
        button2.grid(row=5, column=3, columnspan=1)

        button3 = tk.Button(toolbar, text='Button 3', font=('Helvetica', 20))
        button3.grid(row=5, column=4, columnspan=1)

        # toolbar.pack(side=tk.RIGHT, fill=tk.X)

    def delete_point_row(self):
        for i in self.points_rows[-1]:
            i.destroy()
        self.points_rows.pop()

    def update_plot(self):
        if all([is_valid_mgrs(i[2].get()) for i in self.points_rows[:-1]]):
            plt.clf()

            # update low level points
            p_longitudes = [mgrs.MGRS().toLatLon(i[2].get())[1] for i in self.points_rows[:-1] if is_valid_mgrs(i[2].get())]
            p_latitudes = [mgrs.MGRS().toLatLon(i[2].get())[0] for i in self.points_rows[:-1] if is_valid_mgrs(i[2].get())]
            p_labels = [i[1].get() for i in self.points_rows[:-1] if is_valid_mgrs(i[2].get())]
            if len(p_longitudes) > 0:
                geometry = [Point(xy) for xy in zip(p_longitudes, p_latitudes)]
                df = pd.DataFrame(p_longitudes, p_latitudes)
                crs = {'init': 'epsg:4326'}
                points_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
                points_df.plot(ax=self.fig.add_subplot(111), markersize=5, color='blue', marker='o', label='Neg', ls='-')
                plt.plot(p_longitudes, p_latitudes, lw=1, alpha=0.8, color='b')  # route centerline
                plt.plot(p_longitudes, p_latitudes, lw=20, alpha=0.4, color='m')  # route corridor

                [plt.text(long, lat, name, fontsize=10, ha='center') for long, lat, name in zip(p_longitudes, p_latitudes, p_labels)]

            self.fig.canvas.draw()

    def point_callback(self, var):
        if len(self.points_rows) > 1:
            # only look to delete the last row if > 1 row exist
            if self.points_rows[-1][1].get() == "" and self.points_rows[-1][2].get() == "":
                if self.points_rows[-2][1].get() == "" and self.points_rows[-2][2].get() == "":
                    self.delete_point_row()
                    print("last row deleted")

            else:
                if self.points_rows[-1][1].get() != "" or self.points_rows[-1][2].get() != "":
                    self.add_point_row(self.right_frame, len(self.points_rows) * 3 + 2 + self.row_start)
                    print("row added")

            for i in range(1, len(self.points_rows)):
                if is_valid_mgrs(self.points_rows[i - 1][2].get()):
                    lat, lon = mgrs.MGRS().toLatLon(self.points_rows[i - 1][2].get())

                    _closest = get_closest_rad_dme(lat, lon)
                    if _closest is not None:
                        self.points_rows[i - 1][6].config(state="normal")
                        self.points_rows[i - 1][6].delete(0, END)
                        self.points_rows[i - 1][6].insert(0, _closest)
                        self.points_rows[i - 1][6].config(state="readonly")

                    # print(_closest)

                    if lat is not None:
                        self.points_rows[i-1][3].config(state="normal")
                        self.points_rows[i-1][3].delete(0, END)
                        self.points_rows[i-1][3].insert(0, f"{lat}")
                        self.points_rows[i-1][3].config(state="readonly")
                    if lon is not None:
                        self.points_rows[i-1][7].config(state="normal")
                        self.points_rows[i-1][7].delete(0, END)
                        self.points_rows[i-1][7].insert(0, f"{lon}")
                        self.points_rows[i-1][7].config(state="readonly")

                    if is_valid_mgrs(self.points_rows[i][2].get()):
                        true_crs, mag_crs = get_true_mag_course(self.points_rows[i - 1][2].get(),
                                                                 self.points_rows[i][2].get())
                        if mag_crs is not None:
                            self.points_rows[i][4].config(state="normal")
                            self.points_rows[i][4].delete(0, END)
                            self.points_rows[i][4].insert(0, f"{mag_crs:03.0f}")
                            self.points_rows[i][4].config(state="readonly")

                            self.points_rows[i][8].config(state="normal")
                            self.points_rows[i][8].delete(0, END)
                            self.points_rows[i][8].insert(0, f"{true_crs:03.0f}")
                            self.points_rows[i][8].config(state="readonly")

                        distance = get_distance_nm(self.points_rows[i - 1][2].get(), self.points_rows[i][2].get())
                        if distance is not None:
                            self.points_rows[i][5].config(state="normal")
                            self.points_rows[i][5].delete(0, END)
                            self.points_rows[i][5].insert(0, f"{distance:,.1f}")
                            self.points_rows[i][5].config(state="readonly")

        else:
            # if len < 3, don't delete the last row
            if self.points_rows[-1][1].get() != "" or self.points_rows[-1][2].get() != "":
                self.add_point_row(self.right_frame, len(self.points_rows) * 3 + 2 + self.row_start)
        self.update_plot()
        # print(self.points_rows)


    def run_vertical_profile(self):
        tfads_data = filter_tfads_data_by_country()
        for i in range(0, len(self.points_rows) - 1):
            high_pt, high_elev = run_leg_vertical_profile(self.points_rows[i], self.points_rows[i + 1], tfads_data,
                                                          corridor_size=5)
            # print(high_pt, high_elev)
        plt.show()
