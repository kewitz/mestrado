# -*- coding: utf-8 -*-
from numpy import *
from libs import FEMgmsh as FEM

m = FEM.Mesh(file="""/home/leo/Documents/Master/FEMRafa/cuba.msh""",
             verbose=True)
inCirc = m.nodesOnLine([2, 7], True)
outCirc = m.nodesOnLine(1, True)

# Parse dos dados de excitação.
with open("""/home/leo/Documents/Master/FEMRafa/input.csv""") as f:
    lines = f.readlines()
    sike = [l.strip().split(',') for l in lines]
    sike = array([[float(t)*1E-6, float(v)] for t, v in sike])

# Define função auxiliar e constantes.
vmod = lambda p1, p2: sqrt(sum((p1 - p2)**2))
sigma, eps, dt = 1.0/300.0, 10*8.8541878176E-12, 1E-8


def elementCalc(e):
    n1, n2, n3 = e.nodes
    a, b, c = vmod(n1.p, n2.p), vmod(n1.p, n3.p), vmod(n3.p, n2.p)
    s = (a+b+c)/2
    D = sqrt(s*(s-a)*(s-b)*(s-c))*2
    q1, q2, q3 = n2.y - n3.y, n3.y - n1.y, n1.y - n2.y
    r1, r2, r3 = n3.x - n2.x, n1.x - n3.x, n2.x - n1.x
    X = matrix([[q1*q1 + r1*r1, q1*q2 + r1*r2, q1*q3 + r1*r3],
                [q1*q2 + r1*r2, q2*q2 + r2*r2, q2*q3 + r2*r3],
                [q1*q3 + r1*r3, q3*q2 + r3*r2, q3*q3 + r3*r3]])
    e.D = (eps/(2*D*dt)) * X
    e.M = (sigma/(2*D))*X + e.D


def gradV(e, V):
    indexes = [n.i-1 for n in e.nodes]
    n1, n2, n3 = e.nodes
    v1, v2, v3 = v.take(indexes)
    q1, q2, q3 = n2.y - n3.y, n3.y - n1.y, n1.y - n2.y
    r1, r2, r3 = n3.x - n2.x, n1.x - n3.x, n2.x - n1.x
    DET = n2.x*n3.y + n1.x*n2.y + n3.x*n1.y\
        - n1.x*n3.y - n3.x*n2.y - n2.x*n1.y
    EX = -(q1*v1 + q2*v2 + q3*v3)/DET
    EY = -(r1*v1 + r2*v2 + r3*v3)/DET
    return (EX, EY)

#%% Cria Matriz de rigidez e domínio V.
X = m.getStiff(elementCalc)
V = zeros((X.shape[0], len(sike)+1))
ns = X.shape[0]
nin = len(inCirc)
non = len(outCirc)

for bi in inCirc+outCirc:
    X[bi, :] = zeros(ns)
    X[bi, bi] = 1.0

for i, Vi in enumerate(sike[:, 1]):
    va = zeros(ns)
    for e in [e for e in m.elements if e.dim is 2]:
        indexes = [n.i-1 for n in e.nodes]
        ve = e.D.dot(V[:, i].take(indexes))
        va.put(indexes, ve.A[0]+va.take(indexes))
    # Seta condições de contorno.
    va.put(inCirc, zeros(nin)+Vi)
    va.put(outCirc, zeros(non))
    # Solve.
    V[1:, i+1] = linalg.solve(X[1:, 1:], va[1:])
    if i % 10 == 0:
        print i
m.plotResult(result=V[:, -1])
save("""/home/leo/Documents/Master/FEMRafa/resultados.npy""", V)
