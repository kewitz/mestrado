# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:09:54 2014

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# Defs
irange = lambda a: zip(range(len(a)),a)
pi = np.pi

# Parâmetros
L =     1.0 # Comprimento do fio em m.
a =     1E-3 # Raio do condutor.
Delta = .02 # Delta.
V0 = 1.0
Eps = 8.854E-12

# Funções
l = lambda m, n: Delta/np.abs(y[m-1]-y[n-1]) if m != n else 2*np.log(Delta/a)
g = lambda m:    4*pi*Eps*V0

# Domínios e constantes
y = np.arange(0,L,Delta)
ns = y.size
ms = ns
L = np.matrix(np.zeros((ms,ns)))
G = np.matrix(np.zeros((ns)))

# Processamento
for i in range(ms):
    m = float(i+1.0)
    G[0,i] = g(m)  # Monta o vetor de tensões
    for j in range(ns):
        n = float(j+1.0)
        L[i,j] = l(m,n)  # Monta a matriz de ?

rho = np.linalg.solve(L,G.T)  # Obtem o vetor de Permeabilidades
Q = rho.sum() * Delta
print "Carga total no condutor %.3epC" % Q

# Plots
p1.plot(y,rho)