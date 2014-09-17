# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Polygon

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
        self.material = int(args[3])
        self.source = float(args[4])

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
        cbase = 'green white green blue coral violet navy brown orange beige'.split()
        fig, ax = plt.subplots()
        for e in self.elements:
            x = []
            y = []
            for ni in e.nodes:
                n = self.getNode(ni)
                x.append(n.x)
                y.append(n.y)
            c = cbase[e.material] if e.source == 0.0 else 'grey'  # Cinza se for fonte.
            p = Polygon(zip(x,y), closed=True, fc=c, ec='k', alpha=.3, linewidth=.4)
            ax.add_patch(p)
        ax.axis('equal')
        plt.show()

m = Mesh("""/home/leo/Documents/Master/Eletromagnetismo Computacional 2/vanti/dados.malha""")
m.plotMesh()