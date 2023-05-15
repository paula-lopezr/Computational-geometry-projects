import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree


class AutonFleet:

    def __init__(self, data):
        self.data = np.loadtxt(data)
        self.loc = self._pointVessels()
        self.kdt = KDTree(self.loc)
        self.new_data = self._new_data()
        self.kdt2 = KDTree(self.new_data)

    def _location(self, x1, x2, y1, y2):
        """
        Returns the mean point between two points.
        """

        p1 = np.mean([x1, x2])
        p2 = np.mean([y1, y2])
        return (p1, p2)

    def _distance(self, x1, x2, y1, y2):
        """
        Returns the distancce between two points.
        """

        d = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        return d

    def _pointVessels(self):
        """
        Returns the locations of the ships as a point.
        """

        loc = []
        for i in range(self.data.shape[0]):
            loc.append(self._location(self.data[i, 0], self.data[i, 1],
                                      self.data[i, 2], self.data[i, 3]))
        return np.array(loc)

    def _new_data(self):
        """
        Returns the new data transformed in terms of stern and bow to treat them
        as a segment into the KDTree.
        """

        new = self.data.copy()
        new[:, [1, 2]] = new[:, [2, 1]]
        new = new.reshape((new.shape[0] * 2, 2))
        return new

    def nearest_ship(self, s, b):
        """
        Returns the s closest ships according to the location of the ship b.
        Plots the ships as points inside a circle.
        """

        idx = self.kdt.query(self.loc[b], s+1)
        cir = plt.Circle(self.loc[b], idx[0][-1], fill=False,
                         color=[30 / 255, 51 / 255, 219 / 255])
        fig, ax = plt.subplots()
        ax.plot(self.loc[:, 0], self.loc[:, 1], 'o',
                alpha=1, zorder=0, c='black', lw=0.5)
        ax.plot(self.loc[b][0], self.loc[b][1], 'o', c='magenta', zorder=10)
        ax.plot(self.loc[idx[1]][:, 0], self.loc[idx[1]][:, 1], 'o',
                c=[1, 117 / 255, 20 // 255], zorder=0)
        ax.add_patch(cir)
        plt.show()
        return idx[1]

    def avoiding_collisions(self, r, b):
        """
        Plot the ships that are close to the ship b in a radius of r.
        """

        idx = self.kdt2.query_ball_point(self.loc[b], r)
        idx = sorted(idx)
        cir = plt.Circle(
            self.loc[b], r, fill=False, color=[
                0, 160 / 255, 227 / 255])
        fig, ax = plt.subplots()
        c = -1
        for i in range(len(self.new_data) - 1):
            if i in idx:
                c += 1
                if idx[c] % 2 == 0 and idx[c + 1] == idx[c] + 1:
                    p1 = idx[c]
                    p2 = idx[c + 1]
                    ax.plot((self.new_data[p1][0], self.new_data[p2][0]),
                            (self.new_data[p1][1], self.new_data[p2][1]),
                            c=[119 / 255, 0, 227 / 255])
                elif idx[c] % 2 == 0:
                    ax.plot((self.new_data[i][0],
                             self.new_data[i + 1][0]),
                            (self.new_data[i][1],
                             self.new_data[i + 1][1]),
                            c='black')
            elif i % 2 == 0:
                ax.plot((self.new_data[i][0],
                         self.new_data[i + 1][0]),
                        (self.new_data[i][1],
                         self.new_data[i + 1][1]),
                        c='black')
        ax.add_patch(cir)
        plt.show()

    def partner_lookup(self):
        """
        Plots all the closest-pair of locations of ships.
        """

        vic = []
        for i in range(len(self.loc)):
            vic.append((self.kdt.query(self.loc[i], 2))[1])
        plt.plot((self.data[:, 0], self.data[:, 1]), (self.data[:,
                 2], self.data[:, 3]), alpha=1, zorder=10, lw=0.5)
        for i in vic:
            plt.plot(self.loc[i][:, 0], self.loc[i][:, 1], '--', c='black')
        plt.show()

    def operation_radius(self):
        """
        Returns the approximate radius of
        the circle that contains all the ships.
        Plots the approximate circle.
        """

        mix = self.kdt2.mins
        may = self.kdt2.maxes
        radius = self._distance(mix[0], may[0], mix[1], may[1]) / 2
        center = self._location(mix[0], may[0], mix[1], may[1])
        circ = plt.Circle(
            center, radius, fill=False, color=[
                221 / 255, 69 / 255, 129 / 255])
        fig, ax = plt.subplots()
        for i in range(len(self.new_data)):
            if i % 2 == 0:
                ax.plot((self.new_data[i][0],
                         self.new_data[i + 1][0]),
                        (self.new_data[i][1],
                         self.new_data[i + 1][1]),
                        c='black')
        ax.add_patch(circ)
        plt.show()
        return radius
