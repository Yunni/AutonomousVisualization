import operator

import matplotlib.pyplot as plt
import matplotlib.delaunay as triang


def distance(x, y, i, j):
    '''The Euclidean distance of the i-th element and the j-th element

    '''
    return (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2


class UnionFind:
    '''A union-find structure to find clusters of a minimum spanning tree

    '''
    def __init__(self):
        self.elems = set()
        self.parent = {}

    def make_set(self, x):
        self.parent[x] = x
        self.elems.add(x)

    def find(self, x):
        if self.parent[x] == x:
            return x
        else:
            return self.find(self.parent[x])

    def union(self, x, y):
        if not x in self.parent:
            self.make_set(x)
        if not y in self.parent:
            self.make_set(y)
        x_root = self.find(x)
        y_root = self.find(y)
        self.parent[x_root] = y_root

    def second_largest_size(self):
        set_lengths = {}
        for elem in self.elems:
            root = self.find(elem)
            set_lengths[root] = set_lengths.get(root, 0) + 1
        set_lengths = sorted(set_lengths)
        if len(set_lengths) > 1:
            return set_lengths[1]
        return 0


class MinSpanTree:
    @classmethod
    def triangulation(cls, x, y):
        cens, edg, tri, neig = triang.delaunay(x, y)
        for t in tri:
            x_line = [x[t[0]], x[t[1]], x[t[2]], x[t[0]]]
            y_line = [y[t[0]], y[t[1]], y[t[2]], y[t[0]]]
            plt.plot(x_line, y_line)
        plt.plot(x, y, '.')
        plt.show()

    @classmethod
    def get_mst_score(cls, x, y):
        ''' Use the size of second largest connected component
        when building MST for Delaunay triangulation as score

        :param x: x-coordinates of all data
        :param y: y-coordinates of all data
        :return:
        '''
        # Do Delaunay triangulation, use length of the edges as weight
        cens, edg, tri, neig = triang.delaunay(x, y)
        edges = {}
        for t in tri:
            if not (t[0], t[1]) in edges.items():
                edges[(t[0], t[1])] = distance(x, y, t[0], t[1])
            if not (t[1], t[2]) in edges.items():
                edges[(t[1], t[2])] = distance(x, y, t[1], t[2])
            if not (t[2], t[0]) in edges.items():
                edges[(t[2], t[0])] = distance(x, y, t[2], t[0])

        # Use Kruskal to build a minimum spanning tree
        ordered_edges = sorted(edges.items(), key=operator.itemgetter(1))
        mst_edges = []
        lengths = []
        i = len(x)
        uf = UnionFind()
        for ((t1, t2), length) in ordered_edges:
            if MinSpanTree.is_connected(mst_edges, t1, t2):
                continue
            mst_edges.append((t1, t2))
            lengths.append(length)
            uf.union(t1, t2)
            i -= 1
            if length > lengths[0] * 4:
                break
        return uf.second_largest_size()

    @classmethod
    def is_connected(cls, edges, start, end):
        '''A depth first search of whether the nodes start and and is weak connected.

        :param edges: A list of tuples of edges (v1, v2)
        :param start: The start node
        :param end: The weak node
        :return:
        '''
        to_check = [start]
        checked = set()
        while to_check:
            node = to_check.pop()
            if node == end:
                return True
            checked.add(node)
            for t1, t2 in edges:
                if t1 == node and t2 not in checked:
                    to_check.append(t2)
                if t2 == node and t1 not in checked:
                    to_check.append(t1)
        return False