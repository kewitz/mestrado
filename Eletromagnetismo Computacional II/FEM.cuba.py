# -*- coding: utf-8 -*-
from numpy import *
from libs import FEMgmsh as FEM

m = FEM.Mesh(file="""/home/leo/Documents/Master/FEMRafa/cuba.msh""",
             verbose=True)
inCirc = list(set([n.i for el in m.elements if el.dim == 1
                   and 2 in el.tags for n in el.nodes]))
outCirc = list(set([n.i for el in m.elements if el.dim == 1
                    and 1 in el.tags for n in el.nodes]))


def elementCalc(e):
    dN = matrix('-1 1 0;-1 0 1', dtype=double)
    w = double(.5)
    n1, n2, n3 = e.nodes
    # Calcular jacobiana.
    J = matrix([[n2.x-n1.x, n2.y-n1.y], [n3.x-n1.x, n3.y-n1.y]])
    det = linalg.det(J)
    Jinv = matrix([[n3.y-n1.y, n1.y-n2.y], [n1.x-n3.x, n2.x-n1.x]])/det
    M = (Jinv*dN)
    e.M = M.T * M * w * det

ss = m.getStiff(elementCalc)

a = zeros(ss.shape[0])
r = zeros(ss.shape[0])

for i in inCirc:
    i -= 1
    ss[i, :] = zeros(ss.shape[0])
    ss[i, i] = 1.0
    a[i] = 100.0
for i in outCirc:
    i -= 1
    ss[i, :] = zeros(ss.shape[0])
    ss[i, i] = 1.0

r[1:] = linalg.solve(ss[1:, 1:], a[1:])
m.plotMesh(result=r)
