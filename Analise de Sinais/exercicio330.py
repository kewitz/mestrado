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

# 3.30
w = 2*pi/6

# a)
akx = makeAk(lambda n: 1 + np.cos(2*(pi/6)*n), 6, w)
# b)
aky = makeAk(lambda n: np.sin(2*n*(pi/6) + pi/4), 6, w)
# d)
akz = makeAk(lambda n: np.sin(2*(pi/6)*n + (pi/4))*(1 + np.cos(2*(pi/6)*n)), 6, w)
# c)
akz2 = np.convolve(akx,aky)

# Plots
fsize = (8,4)
picfolder = '/home/leo/Documents/Sinais e Sistemas Lineares/exercicios/'
for tp in [('a',akx),('b',aky),('d',akz),('c',akz2)]:
	fname = '3-30'+tp[0]+'.eps'
	print 'Saving', fname
	plt.figure(figsize=fsize)
	ks = np.arange(len(tp[1]))
	discplot(ks, tp[1])
	plt.xlim(ks.min()-1,ks.max()+1)
	plt.xlabel('$k$')
	plt.ylabel('$real(Ak)$')
	plt.savefig(picfolder+fname)
	plt.close()

# Comparação entre c) e d)
plt.figure(figsize=fsize)
x = range(20)
z1 = regen(w,akz,range(20))
z2 = regen(w,akz2,range(20))
plt.plot(x,z1,'b+',x,z2,'rx', markersize=10)
plt.xlim(min(x)-1, max(x)+1)
plt.grid()
plt.xlabel('$n$')
plt.legend(('z[n] (d)','z[n] (c)'))
fname = '3-30cd.eps'
plt.savefig(picfolder+fname)
plt.close()