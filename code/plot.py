import matplotlib.pyplot as plt
from math import *

from hdt import DipTest

from mst import MinSpanTree


unaries = ["log", "sqrt", "exp"]
binaries = ["+", "-", "*", "/", "**"]


class ScatterPlot:
    '''A scatter plot class

    '''
    def __init__(self, x_expr, y_expr):
        '''Set the expressions

        :param x_expr:
        :param y_expr:
        :return:
        '''
        self.expr = (x_expr, y_expr)
        self.x = []
        self.y = []
        self.score = float("-inf")

    def set_points(self, features, data):
        '''Set points on the graph

        :param features:
        :param data:
        :return:
        '''
        for row in data:
            v = {}
            for i in range(len(features)):
                v[features[i]] = eval(row[i])
            try:
                x = eval(self.expr[0])
                y = eval(self.expr[1])
            except OverflowError:
                return False
            self.x.append(x)
            self.y.append(y)
        try:
            x_min = min(self.x)
            y_min = min(self.y)
        except TypeError:
            return False
        self.x = [record - x_min + 1 for record in self.x]
        self.y = [record - y_min + 1 for record in self.y]
        return True

    def plot_data(self):
        '''Show the scatter plot

        :return:
        '''
        MinSpanTree.get_mst(self.x, self.y)
        plt.xlabel(self.expr[0])
        plt.ylabel(self.expr[1])
        plt.plot(self.x, self.y, '.')
        plt.show()

    def calculate_score1(self):
        '''Caculate the score using Hartigen's Dip test

        :return:
        '''
        self.score = (DipTest(self.x) + DipTest(self.y))
        return self.score

    def calculate_score2(self):
        '''Caculate the score using Minimum Spanning Tree

        :return:
        '''
        self.score = MinSpanTree.get_mst_score(self.x, self.y)
        return self.score