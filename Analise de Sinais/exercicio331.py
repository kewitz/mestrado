# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:07:28 2014

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
exp = np.exp
pi = np.pi


def discplot(x,y, color='k', marker='o'):
    plt.plot(x, y, ''+color+marker, markersize=4)    
    plt.vlines(x, [0], y, ''+color, lw=2)    
    plt.grid(True)
				
def regen(w,aks,ns):
	ks = len(aks)
	fn = lambda n: np.array([aks[k]*exp(1j*k*n*w) for k in range(ks)]).sum()
	return [fn(n) for n in ns]
	
def makeAk(funcn, periodo, w):
	weights = [funcn(n) for n in range(periodo)]
	somatorio = lambda k: np.array([weights[n]*exp(-1j*k*n*w) for n in range(periodo)]).sum()
	return [(1.0/periodo)*(somatorio(k)) for k in range(periodo)]

fsize = (8,4)
picfolder = '/home/leo/Documents/Sinais e Sistemas Lineares/exercicios/'
save = True

# 3.31
ns = np.arange(0,30)
# a)
xn = lambda n: 1 if n%10<8 else 0
gn = lambda n: xn(n) - xn(n-1)
plt.figure(figsize=fsize)
discplot(ns, map(gn,ns))
plt.xlim(ns.min()-1,ns.max()+1)
plt.ylim(-1.5,1.5)
plt.xlabel('$n$')
plt.ylabel('$g[n]$')
if save:
	plt.savefig(picfolder+'3-31a.eps')
	plt.close()

# b)
w = 2*pi/10
aks = makeAk(gn,10,w)
plt.figure(figsize=fsize)
discplot(range(len(aks)),aks)
plt.xlim(-1,10)
plt.xlabel('$k$')
plt.ylabel('$real(Ak)$')
if save:
	plt.savefig(picfolder+'3-31b.eps')
	plt.close()

# c)