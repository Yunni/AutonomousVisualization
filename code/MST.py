import matplotlib.pyplot as plt
import matplotlib.delaunay as triang
import operator

class UnionFind:
    def __init__(self):
        self.elems = set()
        self.parent = {}
    
    def makeSet(self, x):
        self.parent[x] = x
        self.elems.add(x)
    
    def find(self, x):
        if self.parent[x] == x:
            return x
        else:
            return self.find(self.parent[x])
    
    def union(self, x, y):
        if not x in self.parent:
            self.makeSet(x)
        if not y in self.parent:
            self.makeSet(y)
        xRoot = self.find(x)
        yRoot = self.find(y)
        self.parent[xRoot] = yRoot
        
    def secondLargestSet(self):
        L = {}
        for elem in self.elems:
            s = self.find(elem)
            if not s in L:
                L[s] = 0
            L[s] += 1
        L = sorted(L)
        result = 0
        if len(L) > 1:
            result = L[1]
        return result

class MST:
    def triangulation(x, y):
        cens,edg,tri,neig = triang.delaunay(x,y)
        for t in tri:
            xT = [x[t[0]], x[t[1]], x[t[2]], x[t[0]]]
            yT = [y[t[0]], y[t[1]], y[t[2]], y[t[0]]]
            plt.plot(xT, yT)
        plt.plot(x, y, '.')
        plt.show()
            
    def getMST(x, y):
        cens,edg,tri,neig = triang.delaunay(x,y)
        edges = {}
        for t in tri:
            if not (t[0],t[1]) in edges.items():
                edges[(t[0],t[1])] = (x[t[0]] - x[t[1]])**2 + (y[t[0]] - y[t[1]])**2
            if not (t[1],t[2]) in edges.items():
                edges[(t[1],t[2])] = (x[t[1]] - x[t[2]])**2 + (y[t[1]] - y[t[2]])**2
            if not (t[2],t[0]) in edges.items():
                edges[(t[2],t[0])] = (x[t[2]] - x[t[0]])**2 + (y[t[2]] - y[t[0]])**2
        orderedEdges = sorted(edges.items(), key=operator.itemgetter(1))
        E = []
        L = []
        i = len(x)
        uf = UnionFind()
        for ((t1, t2), length) in orderedEdges:
            if MST.BFS(E, t1, t2):
                continue
            E.append((t1, t2))
            L.append(length)
            uf.union(t1, t2)
            i = i - 1
            if length > L[0]*4: break
        #print(E)
        #for (t1, t2) in E:
            #xT = [x[t1], x[t2]]
            #yT = [y[t1], y[t2]]
            #plt.plot(xT, yT)
        #plt.plot(x, y, '.')
        #plt.show()
        return uf.secondLargestSet()
        
    def BFS(graph, start, end):
        todo = [(start, [start])]
        while 0 < len(todo):
            (node, path) = todo.pop(0)
            next_nodes = []
            for (t1, t2) in graph:
                if t1 == node:
                    next_nodes.append(t2)
                if t2 == node:
                    next_nodes.append(t1)
            for next_node in next_nodes:
                if next_node in path:
                    continue
                elif next_node == end:
                    return True
                else:
                    todo.append([next_node, path + [next_node]])
        return False