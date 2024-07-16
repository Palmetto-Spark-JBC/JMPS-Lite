import unittest

from jmps_functions import Waypoint
from jmps_functions import get_distance, get_true_mag_course, get_rectangle_points, get_dted_corners, get_dted_filenames, is_point_inside_polygon
from shapely.geometry import Polygon


class TestGetDistance(unittest.TestCase):

    def test_distance(self):
        A1 = Waypoint('A1', '17SQU3778640909')
        B1 = Waypoint('B1', '18STD2599757311')
        self.assertTrue(get_distance(A1, B1) == (40399.9228, 21.8003))


class TestInsideRegion(unittest.TestCase):

    def test_inside_polygon(self):
        region_points = [(34.33008188232185, -78.58286403879423), (34.30404343551557, -78.46604946026831),
                         (34.66954400819279, -78.34581849511399), (34.69558245487992, -78.46314634675811),
                         (34.33008188232185, -78.58286403879423)]
        inside_points = [(34.35, -78.5), (34.35, -78.55), (34.4, -78.45), (34.4, -78.5), (34.4, -78.55),
                         (34.45, -78.45), (34.45, -78.5), (34.5, -78.45), (34.5, -78.5), (34.55, -78.4),
                         (34.55, -78.45), (34.55, -78.5), (34.6, -78.4), (34.6, -78.45), (34.65, -78.4),
                         (34.65, -78.45)]
        outside_points = [(34.3, -78.3), (34.3, -78.35), (34.3, -78.4), (34.3, -78.45), (34.3, -78.5), (34.3, -78.55),
                          (34.3, -78.6), (34.35, -78.3), (34.35, -78.35), (34.35, -78.4), (34.35, -78.45),
                          (34.35, -78.6), (34.4, -78.3), (34.4, -78.35), (34.4, -78.4), (34.4, -78.6), (34.45, -78.3),
                          (34.45, -78.35), (34.45, -78.4), (34.45, -78.55), (34.45, -78.6), (34.5, -78.3),
                          (34.5, -78.35), (34.5, -78.4), (34.5, -78.55), (34.5, -78.6), (34.55, -78.3),
                          (34.55, -78.35), (34.55, -78.55), (34.55, -78.6), (34.6, -78.3), (34.6, -78.35),
                          (34.6, -78.5), (34.6, -78.55), (34.6, -78.6), (34.65, -78.3), (34.65, -78.35),
                          (34.65, -78.5), (34.65, -78.55), (34.65, -78.6), (34.7, -78.3), (34.7, -78.35),
                          (34.7, -78.4), (34.7, -78.45), (34.7, -78.5), (34.7, -78.55), (34.7, -78.6)]
        for pt in inside_points:
            self.assertTrue(is_point_inside_polygon(pt[0], pt[1], region_points))
        for pt in outside_points:
            self.assertFalse(is_point_inside_polygon(pt[0], pt[1], region_points))


if __name__ == '__main__':
    unittest.main()

