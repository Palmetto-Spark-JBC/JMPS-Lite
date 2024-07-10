# https://gis.stackexchange.com/questions/367496/plot-a-circle-with-a-given-radius-around-points-on-map-using-python
# plot threats like this

# https://towardsdatascience.com/geopandas-101-plot-any-data-with-a-latitude-and-longitude-on-a-map-98e01944b972


import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, Polygon
from jmps_functions import Waypoint, get_distance, get_true_mag_course, get_rectangle_points, get_sub_points

CEP = Waypoint('CEP', '17SQU2778500089')
A1 = Waypoint('A1', '17SQU3778640909')
B1 = Waypoint('B1', '18STD2599757311')
C1 = Waypoint('C1', '18STE4601908063')
D1 = Waypoint('D1', '18STE4979850611')
E1 = Waypoint('E1', '18STE5818478635')
F1 = Waypoint('F1', '18STF3435315366')
G1 = Waypoint('G1', '17SQA4148853803')
I04 = Waypoint('I04', '17SQA6051091243')

sub_points = get_sub_points(CEP, A1)
print("CEP", *CEP.latlon)
print(sub_points)
print(len(sub_points))
for i in range(len(sub_points)-1):
    print(sub_points[i], sub_points[i+1], get_distance(sub_points[i], sub_points[i+1]), get_true_mag_course(sub_points[i], sub_points[i+1]))
print("A1", *A1.latlon)
# PICK UP HERE - SPLICING A LEG INTO 100 METER CHUNKS TO THEN COMPARE DISTANCE TO EACH DTED POINT, BUT THE SUB POINTS LIST SEEMS TO GRADUALLY GET OFF THE HEADING very slightly. plot this to see significance
exit(123)
points = [CEP, A1, B1, C1, D1, E1, F1, G1, I04]
for i in points:
    print(i.name, i.MGRS, i.latlon)
p_longitudes = [i.latlon[1] for i in points]
p_latitudes = [i.latlon[0] for i in points]
geometry = [Point(xy) for xy in zip(p_longitudes, p_latitudes)]
df = pd.DataFrame(p_longitudes, p_latitudes)
crs = {'init': 'epsg:4326'}
points_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)


fig, ax = plt.subplots(figsize=(15, 15))
points_df.plot(ax=ax, markersize=5, color='blue', marker='o', label='Neg', ls='-')
[plt.text(pt.latlon[1], pt.latlon[0], f'{pt.name}', fontsize=10, ha='center') for pt in points]
plt.plot(p_longitudes, p_latitudes, lw=1, alpha=0.8, color='b')
plt.plot(p_longitudes, p_latitudes, lw=20, alpha=0.4, color='m')
plt.show()
