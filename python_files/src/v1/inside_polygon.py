import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import time

def point_inside_region(x, y, region_points):
    # Initialize counter to count the number of intersections
    intersections = 0

    # Iterate through each pair of consecutive points in the region
    for i in range(len(region_points) - 1):
        x1, y1 = region_points[i]
        x2, y2 = region_points[i + 1]

        # Check if the point lies on the horizontal line passing through y
        if min(y1, y2) <= y <= max(y1, y2):
            # Check if the point lies to the right of the line segment
            if x < max(x1, x2):
                # Calculate the x-coordinate of the intersection point
                intersection_x = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                # If the point lies to the left of the intersection point, it's inside the region
                if x <= intersection_x:
                    intersections += 1

    # If the number of intersections is odd, the point is inside the region
    return intersections % 2 == 1


def point_inside_polygon(x, y, region_points):
    region_polygon = Polygon(region_points)
    _point = Point(x, y)
    return region_polygon.contains(_point)


def point_inside_polygon2(x, y, poly):
    """
    Determine if a point is inside a polygon.

    Arguments:
    x -- x-coordinate of the point
    y -- y-coordinate of the point
    poly -- list of (x, y) tuples representing the vertices of the polygon

    Returns:
    True if the point is inside the polygon, False otherwise
    """
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


# Example usage:
region_points = [
    (-78.7, 34.4),
    (-78.6, 34.7),
    (-78.2, 34.6),
    (-78.3, 34.3),
    (-78.7, 34.4)
]


# plt.plot([i[0] for i in region_points], [i[1] for i in region_points])
# plt.draw()

print(region_points)
start = time.time()
for i in range(1000):
    # Check if the test point is inside the region
    for x in range(-7900, -7800, 1):
        test_point = (x/100, 34.5)
        inside_region = point_inside_polygon2(test_point[0], test_point[1], region_points)
        color = 'go' if inside_region else 'ro'
        # plt.plot(test_point[0], test_point[1], color)
        # plt.draw()

    # Check if the test point is inside the region
    for y in range(3400, 3500, 1):
        test_point = (-78.5, y/100)
        # messes up at -78.5, 34.6
        inside_region = point_inside_polygon2(test_point[0], test_point[1], region_points)
        color = 'go' if inside_region else 'ro'
        # plt.plot(test_point[0], test_point[1], color)
        # plt.draw()

# plt.show()
end = time.time()
print(f"Func 3 took {end-start}")
# Func 1 took 0.33858776092529297 for 1000 iters
# Func 2 took 12.446444511413574 for 1000 iters
# Func 3 took 0.4298543930053711 for 1000 iters
# point_inside_polygon2 is the optimal function for determining if a point is inside the polygon
