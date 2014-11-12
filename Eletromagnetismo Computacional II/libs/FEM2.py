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
from multiprocessing import Pool, Array

il = lambda l, offset=0: zip(arange(len(l))+offset,l)
vmod = lambda p1, p2: sqrt(sum((p1 - p2)**2))
materiais = {
    1: {'mur':1.0, 'sigma': 0.0, 'color':'white'},
    2: {'mur':1000.0, 'sigma': 0.0, 'color':'blue'},
    5: {'mur':1000.0, 'sigma': 1E7, 'color':'red'}
}
mu = 4E-7*pi

class Node:
    def __init__(self,postring,i):
        x, y = postring.strip().split()
        self.x, self.y = [float(x), float(y)]
        self.i = i
        self.p = array([self.x,self.y])

class Element:
    dN = matrix('-1 1 0;-1 0 1', dtype=double)
    w = double(.5)
    
    def __init__(self,nos,i,mesh):
        params = nos.strip().split()
        self.dT = dT = mesh.dT
        self.nodes = [mesh.getNode(int(a)) for a in params[0:3]]
        self.material = materiais[int(params[3])]
        self.source = float(params[4])
        self.input = self.source > 0
        n1, n2, n3 = self.nodes
        a, b, c = vmod(n1.p,n2.p), vmod(n1.p,n3.p), vmod(n3.p,n2.p)
        s = (a+b+c)/2
        self.D = D = sqrt(s*(s-a)*(s-b)*(s-c))*2
        q1, q2, q3 = n2.y - n3.y, n3.y - n1.y, n1.y - n2.y
        r1, r2, r3 = n3.x - n2.x, n1.x - n3.x, n2.x - n1.x
        self.qr = matrix([[q1*q1+r1*r1, q1*q2 + r1*r2, q1*q3+r1*r3],
                          [q1*q2 + r1*r2, q2*q2 + r2*r2, q2*q3+r2*r3],
                          [q1*q3+r1*r3, q3*q2 + r3*r2, q3*q3+r3*r3]]) * -1.0/(mu*self.material['mur']*2*D)
        self.qr2 = (-self.material['sigma']*D/(12*dT))*matrix([[1, 0.5, 0.5],[0.5, 1, 0.5], [0.5, 0.5, 1]])

    def makeSS(self, ns):
        SS = zeros((ns,ns))
        # Preenche a matriz de rigidez.
        for li in range(3):
            gi = self.nodes[li].i - 1
            for lj in range(3):
                gj = self.nodes[lj].i - 1
                SS[gi,gj] += self.qr[li,lj] + self.qr2[li,lj]
        return SS
    
    def makeA(self, ns, **k):
        assert 'A' in k and 'JE' in k, "Need to specify current vector A and JE"
        A = zeros(ns)
        A1, A2, A3 = [k['A'][n.i-1] for n in self.nodes]
        D, dT = self.D, self.dT
        Q = (self.material['sigma']*D/(12*dT))*(matrix([[1, 0.5, 0.5],[0.5, 1, 0.5], [0.5, 0.5, 1]]).dot(array([A1, A2, A3])))
        if self.input:
            Q += k['JE']*D/6
        N1, N2, N3 = [n.i-1 for n in self.nodes]
        A[N1] = Q[0,0]
        A[N2] = Q[0,1]
        A[N3] = Q[0,2]
        return A
        

class Mesh:
    def __init__(self, filepath, **kwargs):
        self.dT = kwargs['dT'] if 'dT' in kwargs else 1
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
    
    def makeStiff(self,iterations=1,**kargs):
        ns = len(self.nodes)
        self.stiff = zeros((ns,ns))
        for s in map(lambda e: e.makeSS(ns), self.elements):
            self.stiff += s
        # Aplica condição de Contorno:
        for i,v in self.bounds:
            i -= 1
            self.stiff[i,:] = zeros(ns)
            self.stiff[i,i] = 1.0
        return self.stiff
    
    def plotMesh(self,**k):
        cbase = 'green white green blue coral violet navy brown orange beige'.split()
        fig, ax = plt.subplots()
        for e in self.elements:
            x = []
            y = []
            for n in e.nodes:
                x.append(n.x)
                y.append(n.y)
            c = e.material['color'] if e.source == 0.0 else 'grey'  # Cinza se for fonte.
            p = Polygon(zip(x, y), closed=True, fc=c, ec='k', alpha=.1, linewidth=.4)
            ax.add_patch(p)
        ax.axis('equal')
        if 'r' in k:
            self.plotResult(r=k['r'])
        plt.show()
    
    def plotResult(self, **k):
        x = []
        y = []
        z = []
        for n in self.nodes:
            x.append(n.x)
            y.append(n.y)
            z.append(k['r'][n.i-1])
        t = tri.Triangulation(x, y)
        plt.tricontour(t, z, 15, linewidths=0.5)
        #plt.tricontourf(t, z, 15, cmap=plt.cm.rainbow)
        plt.colorbar()
        #plt.plot(x, y, 'ko', ms=3)
        