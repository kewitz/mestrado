#!/usr/
# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Thu May 15 11:50:45 2014
"""
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.mlab import griddata
from matplotlib.patches import Polygon

il = lambda l, offset=0: zip(arange(len(l))+offset,l)

class Node:
    def __init__(self,postring,i):
        x, y = postring.strip().split()
        self.x, self.y = [float(x), float(y)]
        self.i = i

class Element:
    dN = matrix('-1 1 0;-1 0 1', dtype=double)
    w = double(.5)
    
    def __init__(self,nos,i,mesh):
        params = nos.strip().split()
        self.nodes = [mesh.getNode(int(a)) for a in params[0:3]]
        self.material = int(params[3])
        self.source = float(params[4])
    
    def calc(self,sm):
        n1, n2, n3 = self.nodes
        # Calcular jacobiana.
        J = matrix([[n2.x-n1.x, n2.y-n1.y],[n3.x-n1.x, n3.y-n1.y]])
        det = linalg.det(J)
        Jinv = matrix([[n3.y-n1.y, n1.y-n2.y],[n1.x-n3.x, n2.x-n1.x]])/det
        self.M = (Jinv*self.dN)
        self.M = self.M.T * self.M * self.w * det
        # Preenche a matriz de rigidez.
        for li in range(3):
            gi = self.nodes[li].i - 1
            for lj in range(3):
                gj = self.nodes[lj].i - 1
                sm[gi,gj] += self.M[li,lj]

class Mesh:
    def __init__(self, filepath):
        with open(filepath) as f:
            ls = f.readlines()
            l = ls[0].strip().split()
            ns = int(l[0])
            es = int(l[1])
            cs = int(l[2])
            self.nodes = [Node(a[1],a[0]) for a in il(ls[1:ns+1],1)]
            self.elements = [Element(a[1],a[0],self) for a in il(ls[ns+1:ns+es+1],1)]
            self.bounds = map(lambda b: (int(b.strip().split()[0]),double(b.strip().split()[1])), ls[ns+es+1:ns+es+cs+1])
    
    def getNode(self, index):
        nodes = filter(lambda n: n.i==index, self.nodes)
        if len(nodes) == 1:
            return nodes[0]
        else:
            raise
    
    def solve(self):
        ns = len(self.nodes)
        self.stiff = zeros((ns,ns))
        self.values = zeros(ns)
        # Calcula jacobiana de cada elemento e preenche a matriz de rigidez.
        map(lambda e: e.calc(self.stiff), self.elements)
        # Aplica condição de contorno.
        for bound in self.bounds:
            i, v = bound
            i -= 1
            line = zeros(ns)
            line[i] = 1
            self.stiff[i,:] = line
            self.values[i] = v
        self.result = linalg.solve(self.stiff,self.values.T)
        print "Done!"
    
    def plotMesh(self):
        cbase = 'green white green blue coral violet navy brown orange beige'.split()
        fig, ax = plt.subplots()
        for e in self.elements:
            x = []
            y = []
            for n in e.nodes:
                x.append(n.x)
                y.append(n.y)
            c = cbase[e.material] if e.source == 0.0 else 'grey'  # Cinza se for fonte.
            p = Polygon(zip(x, y), closed=True, fc=c, ec='k', alpha=.3, linewidth=.4)
            ax.add_patch(p)
        ax.axis('equal')
        self.plotResult()
        plt.show()
    
    def plotResult(self):
        x = []
        y = []
        z = []
        for n in self.nodes:
            x.append(n.x)
            y.append(n.y)
            z.append(self.result[n.i-1])
        tri.Triangulation(x, y)
        plt.tricontour(x, y, z, 15, linewidths=0.5, colors='k')
        plt.tricontourf(x, y, z, 15, cmap=plt.cm.rainbow)
        plt.plot(x, y, 'ko', ms=3)
        plt.colorbar()
        