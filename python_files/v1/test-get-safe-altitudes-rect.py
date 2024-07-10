
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, Polygon
from jmps_functions import Waypoint, run_leg_vertical_profile, filter_tfads_data_by_country

CEP = Waypoint('CEP', '17SQU2778500089')
A1 = Waypoint('A1', '17SQU3778640909')
B1 = Waypoint('B1', '18STD2599757311')
C1 = Waypoint('C1', '18STE4601908063')
D1 = Waypoint('D1', '18STE4979850611')
E1 = Waypoint('E1', '18STE5818478635')
F1 = Waypoint('F1', '18STF3435315366')
G1 = Waypoint('G1', '17SQA4148853803')
I04 = Waypoint('I04', '17SQA6051091243')

points = [CEP, A1, B1, C1, D1, E1, F1, G1, I04]
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

tfads_data = filter_tfads_data_by_country()

for i in range(0, len(points)-1):
    high_pt, high_elev = run_leg_vertical_profile(points[i], points[i+1], tfads_data, corridor_size=5)
    # print(high_pt, high_elev)
plt.show()
