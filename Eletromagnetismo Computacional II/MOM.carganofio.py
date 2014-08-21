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
L =     np.float64(1.0) # Comprimento do fio em m.
a =     np.float64(1E-3) # Raio do condutor.
Delta = np.float64(.01) # Delta.
V0 = np.float64(1.0)
Eps = np.float64(8.854E-12)

# Funções
l = lambda m, n: Delta/np.abs(y[m-1]-y[n-1]) if m != n else 2*np.log(Delta/a)
g = lambda m:    4*pi*Eps*V0

# Domínios e constantes
y = np.arange(0,L,Delta, dtype=np.float64)
ns = y.size
ms = ns
L = np.matrix(np.zeros((ms,ns)))
G = np.matrix(np.zeros((ns)))

# Processamento
for i in range(ms):
    m = np.float64(i+1.0)
    G[0,i] = g(m)  # Monta o vetor de tensões
    for j in range(ns):
        n = np.float64(j+1.0)
        L[i,j] = l(m,n)  # Monta a matriz de ?

rho = np.linalg.solve(L,G.T)  # Obtem o vetor de Permeabilidades
Q = rho.sum() * Delta
print "Carga total no condutor %.3eC" % Q

# Plots
plt.plot(y,rho)