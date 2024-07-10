import matplotlib.pyplot as plt
from jmps_functions import is_point_inside_polygon

"""---------------------------------------------------------------------------------------------------------------------
---------------------------------------------------- TEST FUNCTIONS ----------------------------------------------------
---------------------------------------------------------------------------------------------------------------------"""


def convert_meters_to_nm(dist_meters):
    return round(dist_meters * 0.0005396118, 1)


def plot_region(points):
    fig, ax = plt.subplots(figsize=(15, 15))
    p_longitudes = [i[1] for i in points]
    p_latitudes = [i[0] for i in points]
    plt.plot(p_longitudes, p_latitudes, lw=1, alpha=0.8, color='b')
    plt.draw()


def plot_point(point, color='o'):
    plt.plot(point[1], point[0], color)
    plt.draw()


region_points = [(34.33008188232185, -78.58286403879423), (34.30404343551557, -78.46604946026831),
                         (34.66954400819279, -78.34581849511399), (34.69558245487992, -78.46314634675811),
                         (34.33008188232185, -78.58286403879423)]
plot_region(region_points)
test_pts = []
for x in range(3430, 3475, 5):
    for y in range(7830, 7865, 5):
        test_pts.append((x/100, -y/100))
print(test_pts)
inside_points = []
outside_points = []
for pt in test_pts:
    if is_point_inside_polygon(pt[0], pt[1], region_points):
        plot_point(pt, 'go')
        inside_points.append(pt)
    else:
        plot_point(pt, 'ro')
        outside_points.append(pt)
print(inside_points)
print(outside_points)
plt.show()


