# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

il = lambda l, offset=0: zip(np.arange(len(l))+offset,l)

class Node:
    def __init__(self,postring,i):
        x, y = postring.strip().split()
        self.x, self.y = [float(x), float(y)]
        self.i = i

class Element:
    def __init__(self,nos,i):
        args = nos.strip().split()
        self.nodes = [int(a) for a in args[0:3]]
        self.material = args[3]

class Mesh:
    def __init__(self, filepath):
        with open(filepath) as f:
            ls = f.readlines()
            l = ls[0].strip().split()
            ns = int(l[0])
            es = int(l[1])
            cs = int(l[2])
            self.nodes = [Node(a[1],a[0]) for a in il(ls[1:ns+1],1)]
            self.elements = [Element(a[1],a[0]) for a in il(ls[ns+1:ns+es+1],1)]
    
    def getNode(self, index):
        nodes = filter(lambda n: n.i==index, self.nodes)
        if len(nodes) == 1:
            return nodes[0]
        else:
            raise
    
    def plotMesh(self):
        cbase = colors.cnames.keys()
        fig, ax = plt.subplots()
        for e in self.elements:
            x = []
            y = []
            for ni in e.nodes:
                n = self.getNode(ni)
                x.append(n.x)
                y.append(n.y)
            x.append(x[0])
            y.append(y[0])
            ax.plot(x, y, color=cbase[int(e.material)], marker='o')
        ax.axis('equal')
        plt.show()
        pass
        
m = Mesh("""/home/leo/Documents/Master/Eletromagnetismo Computacional 2/vanti/dados.malha""")
m.plotMesh()