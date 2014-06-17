# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:07:28 2014

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
exp = np.exp
pi = np.pi

def discplot(x,y, color='k'):
    plt.plot(x, y, 'o'+color)    
    plt.vlines(x, [0], y, ''+color, lw=2)    
    plt.grid(True)
				
def regen(ak,n,w,ks):
    return np.matrix([ak(k*w) * exp(1j*k*w*n) for k in range(ks)]).sum()
				
# 3.30
plt.subplot(3, 1, 1)
w = 2*pi/6
akx = lambda kw: (1.0/6) * (2 + (3.0/2)*exp(-1j*kw) + (1.0/2)*exp(-2j*kw) + (1.0/2)*exp(-4j*kw) + (3.0/2)*exp(-5j*kw) )
xn = [regen(akx,n,w,7) for n in range(20)]
discplot(range(len(xn)), xn, 'r')

plt.subplot(3, 1, 2)
aky = lambda kw: (1.0/6) * (0.707 + 0.707*exp(-1j*kw) - 0.707*exp(-2j*kw) - 0.707*exp(-3j*kw) + 0.707*exp(-4j*kw) + 0.707*exp(-5j*kw) )
yn = [regen(aky,n,w,7) for n in range(20)]
discplot(range(len(yn)), yn, 'r')

plt.subplot(3, 1, 3)
akz = lambda kw: (1.0/6) * (1.4142 + 1.06*exp(-1j*kw) - 0.353*exp(-2j*kw) + 0.353*exp(-4j*kw) + 1.06*exp(-5j*kw) )
zn = [regen(akz,n,w,7) for n in range(20)]
discplot(range(len(zn)), zn, 'r')