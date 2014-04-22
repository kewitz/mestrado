# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 20:18:07 2014

@author: leo
"""
import matplotlib.pyplot as plt
import numpy as np


#Define time
t = arange(-50,50,dtype=float64)

#Define functions
ut = lambda t: 1 if t >= 0 else 0
xt = lambda t: 2*t*sin(3*t - 2)*ut(t)

#Define signals
x1 = [xt(a) for a in t]
ev = [(xt(a)+xt(-a))/2 for a in t]
od = [(xt(a)-xt(-a))/2 for a in t]

#Plot x(t)
lx1, = plt.plot(t, x1)
#Plot Even and Odd components
lev, = plt.plot(t, ev)
lod, = plt.plot(t, od)

#Add legend and stuff...
plt.legend( (lx1, lev, lod), (r'$x(t)$',r'$Ev\{x(t)\}$', r'$Od\{x(t)\}$'), 'best', shadow=True)
plt.title(r"$x(t) = 2 t \sin(3t-2) u(t)$")
plt.xlabel(r"$t$")
plt.grid(True)
plt.show()