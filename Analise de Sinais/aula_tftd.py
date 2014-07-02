# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 15:03:33 2014

@author: leo
"""
import numpy as np
pi = np.pi
exp = np.exp


def hjw (w):
	if abs(w) < pi/3: return 2
	elif (abs(w) <= (2*pi/3)): return 1.0 - abs(w) + pi/3
	else: return 0

xjw = lambda w: -exp(-3j*w)+2*exp(-1j*w)-1 

wd = np.arange(-pi,pi,.001)

resultado = np.array(map(hjw,wd)) * xjw(wd)

plot(wd,np.abs(resultado))
plt.xticks( [-2*pi/3,-pi/3, 0, pi/3, 2*pi/3], ['$-2\pi / 3$','$-\pi / 3$','$0$','$\pi/3$','$2\pi/3$'])