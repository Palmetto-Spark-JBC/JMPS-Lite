import math
import numpy as np
# import os
# import pickle
# import pyperclip
# import random
# import string
# import sys
import matplotlib.pyplot as plt
from tkinter import *
import geopandas as gpd
from jmps_functions import *

class JMPS_GUI(Tk):

    def __init__(self):
        super().__init__()
        self.points_rows = []
        self.point_headers = []
        self.top_frame = Frame(self)
        self.options = {}
        self.initializeUI()

    def initializeUI(self):
        self.title("JMPS 2.0")
        self.geometry(f"{self.winfo_screenwidth()//2-10}x{self.winfo_screenheight()-80}+{self.winfo_screenwidth()//2}+0")  # set screen to right half
        self.setupWindow()

    def add_point_row(self, frame, row):
        padx = 10
        pady = 2
        bd = 2
        point_label = Label(frame, text=f"{row}.")
        point_label.grid(row=row, column=0, padx=padx, pady=pady)

        point_input_var = StringVar()
        point_input_var.trace("w", lambda name, index, mode, var=point_input_var: self.point_callback(point_input_var))
        point_input = Entry(frame, textvariable=point_input_var, width=10, bd=bd)
        point_input.grid(row=row, column=1, padx=padx)

        mgrs_input_var = StringVar()
        mgrs_input_var.trace("w", lambda name, index, mode, var=mgrs_input_var: self.point_callback(mgrs_input_var))
        mgrs_input = Entry(frame, textvariable=mgrs_input_var, width=20, bd=bd)
        mgrs_input.grid(row=row, column=2, padx=padx)

        mag_course_display = Entry(frame, width=7, bd=bd)
        mag_course_display.grid(row=row, column=3, padx=padx)
        mag_course_display.config(state="readonly")

        distance_display = Entry(frame, width=7, bd=bd)
        distance_display.grid(row=row, column=4, padx=padx)
        distance_display.config(state="readonly")
        self.points_rows.append([point_label, point_input, mgrs_input, mag_course_display, distance_display])

    def delete_point_row(self):
        for i in self.points_rows[-1]:
            i.destroy()
        self.points_rows.pop()

    def setupWindow(self):
        """ Set up the widgets."""
        title = Label(self, text="JMPS 2.0", font=('Helvetica', 20), bd=10)
        title.grid(row=0, column=0, columnspan=6)

        self.top_frame = Frame(self)
        self.top_frame.grid(row=2, column=0, padx=10, pady=10)

        self.point_headers = ["Point", "MGRS", "Mag Course", "Distance (NM)"]
        self.points_rows = []
        for i, header in enumerate(self.point_headers):
            order_label = Label(self.top_frame, text=header, bd=10)
            order_label.grid(row=0, column=i + 1)

        self.add_point_row(self.top_frame, 1)

        self.initialize_plot()

    def point_callback(self, var):
        if len(self.points_rows) > 1:
            # only look to delete the last row if > 1 row exist
            if self.points_rows[-1][1].get() == "" and self.points_rows[-1][2].get() == "":
                if self.points_rows[-2][1].get() == "" and self.points_rows[-2][2].get() == "":
                    self.delete_point_row()

            else:
                if self.points_rows[-1][1].get() != "" or self.points_rows[-1][2].get() != "":
                    self.add_point_row(self.top_frame, len(self.points_rows) + 1)

            for i in range(1, len(self.points_rows)):
                if self.points_rows[i-1][2].get() != "" and self.points_rows[i][2].get() != "":
                    _, mag_course = get_true_mag_course(self.points_rows[i - 1][2].get(), self.points_rows[i][2].get())
                    self.points_rows[i][3].config(state="normal")
                    self.points_rows[i][3].delete(0, END)
                    self.points_rows[i][3].insert(0, f"{mag_course:03.0f}")
                    self.points_rows[i][3].config(state="readonly")

                    distance_meters,distance_nm = get_distance(self.points_rows[i-1][2].get(), self.points_rows[i][2].get())
                    self.points_rows[i][4].config(state="normal")
                    self.points_rows[i][4].delete(0, END)
                    self.points_rows[i][4].insert(0, f"{distance_nm:,.1f}")
                    self.points_rows[i][4].config(state="readonly")

        else:
            # if len < 3, don't delete the last row
            if self.points_rows[-1][1].get() != "" or self.points_rows[-1][2].get() != "":
                self.add_point_row(self.top_frame, len(self.points_rows) + 1)
        self.update_plot()

    def initialize_plot(self):
        px = 1/plt.rcParams['figure.dpi']
        w, h = self.winfo_screenwidth(), self.winfo_screenheight()
        self.fig, self.ax = plt.subplots(figsize=((w//2)*px, (h-80)*px))
        self.ax.plot()
        plt.draw()

    def update_plot(self):
        self.ax.clear()

        # update low level points
        p_longitudes = [mgrs.MGRS().toLatLon(i[2].get())[1] for i in self.points_rows[:-1] if i[2].get() != ""]
        p_latitudes = [mgrs.MGRS().toLatLon(i[2].get())[0] for i in self.points_rows[:-1] if i[2].get() != ""]
        if len(p_longitudes) > 0:
            geometry = [Point(xy) for xy in zip(p_longitudes, p_latitudes)]
            df = pd.DataFrame(p_longitudes, p_latitudes)
            crs = {'init': 'epsg:4326'}
            points_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)
            points_df.plot(ax=self.ax, markersize=5, color='blue', marker='o', label='Neg', ls='-')
            plt.plot(p_longitudes, p_latitudes, lw=1, alpha=0.8, color='b')  # route centerline
            plt.plot(p_longitudes, p_latitudes, lw=20, alpha=0.4, color='m')  # route corridor

        self.ax.plot()
        plt.show()


if __name__ == "__main__":
    app = JMPS_GUI()
    app.mainloop()

# 17squ2778500089
# 17squ5056535082
# 18std2599757311
# 18ste4601908063
# 18ste4979850611
# 18ste5818478635
# 18stf3435315366
# 17sqa3347240843
# 17sqa6051091243
# 18stg3710206584
