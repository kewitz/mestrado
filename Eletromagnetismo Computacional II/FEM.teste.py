# -*- coding: utf-8 -*-
from numpy import *
from libs import FEM2 as FEM

m = FEM.Mesh("""/home/leo/Documents/Master/Eletromagnetismo Computacional 2/vanti/rele.malha""", dT = .5E-2)
ss = m.makeStiff()
ns = len(m.nodes)

time = arange(0,4E-1,m.dT)
A = zeros((ns,len(time)))
JE = [t*200E3 if t < .01 else 2E3 for t in time]
for i,t in enumerate(time[1:],1):
    Q = zeros(ns)
    for e in m.elements:
        Ae = e.makeA(ns, A=A[:,i-1], JE=JE[i])
        Q += Ae
    for ib, bv in m.bounds:
        Q[ib-1] = bv
    A[:,i] = linalg.solve(ss,Q)
#m.plotResult(r=A[:,38])
m.plotMesh(r=A[:,50])