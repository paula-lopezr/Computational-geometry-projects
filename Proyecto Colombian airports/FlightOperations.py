import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from itertools import combinations
from sknetwork.path import shortest_path
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


class FlightOperations:

    def __init__(self, border, airports):
        self.coords = np.loadtxt(airports, delimiter='   ', usecols=(0, 1, 2))
        self.names = np.loadtxt(airports, delimiter='   ',
                                dtype='str', usecols=(4, 5, 6))
        self.border = np.loadtxt(border, delimiter='  ')
        self.vor = Voronoi(np.fliplr(self.coords[:, :2]))
        self.ch = ConvexHull(self.vor.points)

    def plot_original(self, p_new=[]):
        """
        Plot the original map and Voronoi diagram.
        And plot and extra points if it desires.
        """

        min_ = np.min(self.vor.points)
        max_ = np.max(self.vor.points)
        norm = mpl.colors.Normalize(vmin=min_, vmax=max_, clip=True)
        fig = plt.figure(figsize=(8, 6))
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.Blues_r)
        voronoi_plot_2d(self.vor, show_vertices=False, line_colors='pink',
                        line_width=2, line_alpha=0.6, point_size=5)
        plt.plot(self.border[:, 1], self.border[:, 0],
                 alpha=1, zorder=10, c='black', lw=0.5)
        if len(p_new) > 0:
            try:
                plt.plot(p_new[:, 0], p_new[:, 1])
            except BaseException:
                plt.plot(p_new[0], p_new[1], '*')
        plt.show()

    def plotting_altitude(self, p_new=None):
        """
        Plots the location of the ariports and generates a scale color
        depending on the altitude of each airport.
        """

        box = [-100, -100], [-100, 100], [100, -100], [100, 100]
        box_e = [0, 0, 0, 0]
        points_vor = np.append(self.vor.points, box, axis=0)
        elevation = np.array(self.coords[:, 2])
        elevation = np.append(elevation, box_e, axis=0)
        min_ = np.min(elevation)
        max_ = np.max(elevation)
        vor2 = Voronoi(points_vor)
        norm = mpl.colors.Normalize(vmin=min_, vmax=max_, clip=True)
        fig = plt.figure(figsize=(8, 6))
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.Blues_r)
        voronoi_plot_2d(vor2, show_vertices=False, line_colors='pink',
                        line_width=2, line_alpha=0.6, point_size=5)
        for r in range(len(vor2.point_region)):
            region = vor2.regions[vor2.point_region[r]]
            if -1 not in region:
                polygon = [vor2.vertices[i] for i in region]
                plt.fill(*zip(*polygon), color=mapper.to_rgba(elevation[r]))
        plt.plot(self.border[:, 1], self.border[:, 0],
                 alpha=1, zorder=10, c='black', lw=0.5)
        plt.xlim(-82, -65)
        plt.ylim(-10, 20)
        if p_new:
            plt.plot(p_new[0], p_new[1], '*')
        plt.show()

    def _distancia(self, x1, x2, y1, y2):
        """
        Returns de Euclidean distance of a pair of points.
        """

        d = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        return d

    def _in_hull(self, hull, list_min, points, ver_vor):
        """
        Returns the maximum point in the convex hull.
        """

        polygon = Polygon(points[hull])
        point = Point(ver_vor[max(list_min)[1][1]])
        while not (polygon.contains(point)):
            list_min.remove(max(list_min))
            new_air = max(list_min)
            point = Point(ver_vor[max(list_min)[1][1]])
        return new_air

    def find_new_airport(self, plt_cond):
        """
        Returns the coordinates and the center of the new airport according to
        the largest circle centered in the convex hull and the border of Colombia.
        """

        points = self.vor.points
        region_vor = self.vor.regions
        ver_vor = self.vor.vertices
        point_reg = self.vor.point_region
        list_min = []
        hull = []
        for i in range(len(points)):
            lst = []
            index_region = point_reg[i]
            for j in region_vor[index_region]:
                if j != -1:
                    dist = self._distancia(
                        points[i][0], ver_vor[j][0], points[i][1], ver_vor[j][1])
                    lst.append(dist)
                else:
                    hull.append(i)
            if len(lst) != 0:
                list_min.append((min(lst), (i, j)))
        new_air = self._in_hull(hull, list_min, points, ver_vor)
        if plt_cond:
            p_new = ver_vor[new_air[1][1]][0], ver_vor[new_air[1][1]][1]
            self.plot_original(p_new)
        return "Center: " + \
            str(ver_vor[new_air[1][1]]) + " Radius: " + str(new_air[0])

    def farthest_airports(self):
        """
        Returns the farthest airports.
        """

        compair = 0
        vertex = self.ch.vertices
        for comb in combinations(vertex, 2):
            x1 = self.vor.points[comb[0]][0]
            x2 = self.vor.points[comb[1]][0]
            y1 = self.vor.points[comb[0]][1]
            y2 = self.vor.points[comb[1]][1]
            d = self._distancia(x1, x2, y1, y2)
            if compair < d:
                compair = d
                x = (comb, d)
        return "Los aeropuertos más lejanos son:" + \
            str(self.names[x[0][0]]) + str(self.names[x[0][1]])

    def nearest_airports(self):
        """
        Returns the nearest airports
        """

        edges = self.vor.ridge_points
        x0 = self.vor.points[edges[0]][0][0]
        x1 = self.vor.points[edges[0]][1][0]
        y0 = self.vor.points[edges[0]][0][1]
        y1 = self.vor.points[edges[0]][1][1]
        compair = self._distancia(x0, x1, y0, y1)
        for comb in edges:
            x1_ = self.vor.points[comb][0][0]
            x2_ = self.vor.points[comb][1][0]
            y1_ = self.vor.points[comb][0][1]
            y2_ = self.vor.points[comb][1][1]
            d = self._distancia(x1_, x2_, y1_, y2_)
            if compair >= d:
                compair = d
                x = (comb, d)
        return "Los aeropuertos más cercanos son:" + \
            str(self.names[x[0][0]]) + str(self.names[x[0][1]])

    def _path(self, threshold, origin, destination):
        """
        Return the index of the path from origin to destination.
        """

        matrix = np.zeros((len(self.vor.points), len(self.vor.points)))
        edg_neig = self.vor.ridge_points
        for i in edg_neig:
            altitude = self.coords[i][:, 2]
            if threshold < altitude[0] and threshold < altitude[1]:
                matrix[i[0], i[1]] = 1
                matrix[i[1], i[0]] = 1
        path = shortest_path(matrix, origin, destination)
        return path

    def path(self, threshold, origin, destination, plot=False):
        """
        Returns the names of the airport's path from origin to destination.
        And plots the path.
        """

        lst = self._path(threshold, origin, destination)
        if plot:
            self.plot_original(self.vor.points[lst])
        return self.names[lst]
