import numpy as np
import matplotlib.pyplot as plt
import matplotlib.delaunay as triang
import csv
from math import log, exp, sqrt
from hdt import DipTest
from MST import MST

unaries = ["log", "sqrt", "exp"]
binaries = ["+", "-", "*", "/", "**"]

class Plot:
    def __init__(self, xexpr, yexpr):
        self.expr = (xexpr, yexpr)
        self.x = []
        self.y = []
        #self.complexity = 0
        #for opr in unaries + binaries:
            #self.complexity += xexpr.count(opr)
            #self.complexity += yexpr.count(opr)
        self.score = float('-inf')
        
    def setData(self, features, data):
        for row in data:
            v = {}
            for i in range(len(features)):
                v[features[i]] = eval(row[i])
            x = eval(self.expr[0])
            y = eval(self.expr[1])
            self.x.append(x)
            self.y.append(y)
            
        xmin = min(self.x)
        ymin = min(self.y)
        self.x = list(map(lambda d: d - xmin + 1, self.x))
        self.y = list(map(lambda d: d - ymin + 1, self.y))
            
    def plotData(self):
        MST.getMST(self.x,self.y)
        plt.xlabel(self.expr[0])
        plt.ylabel(self.expr[1])
        plt.plot(self.x, self.y, '.')
        plt.show()
        
    def calculateScore1(self):
        self.score = (DipTest(self.x) + DipTest(self.y)) #* (512 - 2 ** self.complexity)#        self.score = score
        return self.score

    def calculateScore2(self):
        self.score = MST.getMST(self.x,self.y)
        return self.score