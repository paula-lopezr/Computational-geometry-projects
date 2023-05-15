from scipy.spatial import Delaunay
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import networkx as nx


class TIN:
    def __init__(self, data):
        self.x = data[:, 0]
        self.y = data[:, 1]
        self.z = data[:, 2]
        self.data = np.array(data)
        self.tri = Delaunay(self.data)
        self.trian = Delaunay(self.data[:, :2])

    def __find_neighbors(self, pindex, trian):
        a, b = trian.vertex_neighbor_vertices
        return b[a[pindex]:a[pindex+1]]

    def plotting(self):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_trisurf(self.x, self.y, self.z, triangles=self.tri.simplices, cmap=cm.jet)
        surf = ax.plot_trisurf(self.x, self.y, self.z, cmap=cm.Pastel1, linewidth=0.1, vmin=0, vmax=1)
        plt.show()

    def elevation(self, pt, plt_cond=False):
        index = self.trian.find_simplex(pt)
        index_pts = self.trian.simplices[index]
        pt1 = self.data[index_pts[0]]
        pt2 = self.data[index_pts[1]]
        pt3 = self.data[index_pts[2]]
        pt_ref = np.cross((pt1 - pt2), (pt2-pt3))
        H = (((pt1[0]-pt[0])*pt_ref[0] + (pt1[1]-pt[1])*pt_ref[1])/pt_ref[2]) + pt1[2]
        if plt_cond:
            xp = [pt1[0], pt2[0], pt3[0]]
            yp = [pt1[1], pt2[1], pt3[1]]
            zp = [pt1[2], pt2[2], pt3[2]]
            self.__plot_pt(xp, yp, zp, pt[0], pt[1], H)
        return H

    def __maxi(self, index):
        if (self.data[index][2] >= max(self.data[self.__find_neighbors(index, self.trian)][:, 2])):
            return True
        else:
            return False

    def peak_point(self, pt):
        index = self.trian.find_simplex(pt)
        if (self.data[index][2] >= max(self.data[self.__find_neighbors(index, self.trian)][:, 2])):
            return True
        else:
            return False

    def relative_max(self):
        maxs, xmx, ymx, zmx = [], [], [], []
        for i in range(len(self.data[:, 2])):
            if self.__maxi(i):
                maxs.append(i)
        for i in maxs:
            xmx.append(self.data[i, 0])
            ymx.append(self.data[i, 1])
            zmx.append(self.data[i, 2])
        self.__plot_pt(self.x, self.y, self.z, xmx, ymx, zmx)

    def __minm(self, index):
        if (self.data[index][2] <= min(self.data[self.__find_neighbors(index, self.trian)][:, 2])):
            return True
        else:
            return False

    def anti_peak_point(self, pt):
        index = self.trian.find_simplex(pt)
        if (self.data[index][2] <= min(self.data[self.__find_neighbors(index, self.trian)][:, 2])):
            return True
        else:
            return False

    def __plot_pt(self, x, y, z, pt1, pt2, pt3):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_trisurf(x, y, z, triangles=self.tri.simplices, cmap=cm.jet, alpha=0.1)
        surf = ax.plot_trisurf(x, y, z, cmap=cm.Pastel1, linewidth=0.1, vmin=0, vmax=1, alpha=0.1)
        ax.scatter(pt1, pt2, pt3, marker='o')
        plt.show()

    def relative_min(self):
        mins, xmn, ymn, zmn = [], [], [], []
        for i in range(len(self.data[:, 2])):
            if self.__minm(i):
                mins.append(i)
        for i in mins:
            xmn.append(self.data[i, 0])
            ymn.append(self.data[i, 1])
            zmn.append(self.data[i, 2])
        self.__plot_pt(self.x, self.y, self.z, xmn, ymn, zmn)

    def __ady_matrix(self):
        ady = np.zeros((len(self.trian.simplices), len(self.trian.simplices)))
        c = 0
        y = 0
        for i in (self.trian.simplices):
            c += 1
            if c % len(self.trian.simplices) != 0:
                for j in self.trian.simplices:
                    y += 1
                    if y % len(self.trian.simplices) != 0:
                        if len(set(i).difference(set(j))) == 1:
                            ady[c-1, y-1] = 1
                        else:
                            ady[c-1, y-1] = 0
                    else:
                        y = 0
            else:
                c = 0
        return ady

    def graph(self):
        centroide = []
        for i in range(self.trian.simplices.shape[0]):
            x = (self.data[self.trian.simplices][i][0][0] + self.data[self.trian.simplices][i][1][0] + self.data[self.trian.simplices][i][2][0])/3
            y = (self.data[self.trian.simplices][i][0][1] + self.data[self.trian.simplices][i][1][1] + self.data[self.trian.simplices][i][2][1])/3
            centroide.append([x, y])
        centroide = np.array(centroide)
        dic = {i: centroide[i] for i in range(len(centroide))}
        rows, cols = np.where(self.__ady_matrix() == 1)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot()
        edges = zip(rows.tolist(), cols.tolist())
        gr = nx.Graph()
        gr.add_edges_from(edges)
        nx.draw(gr, pos=dic, node_size=5, ax=ax)
        ax.triplot(self.x, self.y, triangles=self.trian.simplices, color='pink')
        ax.scatter(centroide[:, 0], centroide[:, 1], s=5, color='black')
        plt.show()