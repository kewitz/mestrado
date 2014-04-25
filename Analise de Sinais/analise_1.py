# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 20:18:07 2014

@author: leo
"""
import matplotlib.pyplot as plt
import numpy as np


#Define time
t = np.arange(-50,50,.1,dtype=np.float64)

#Define functions
ut = lambda t: 1 if t >= 0 else 0
xt = lambda t: 2*t*np.sin(3*t - 2)*ut(t)

#Define signals
x1 = [xt(a) for a in t]
ev = [(xt(a)+xt(-a))/2 for a in t]
od = [(xt(a)-xt(-a))/2 for a in t]

#Setup figure.
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#Plot x(t)

lx1, = ax1.plot(t, x1)

#Plot Even and Odd components
lev, = ax2.plot(t, ev, 'r-', alpha=0.6)
lod, = ax3.plot(t, od, 'g-', alpha=0.6)

#Add legend and stuff...
ax1.set_title(r"$x(t) = 2 t \sin(3t-2) u(t)$")
ax2.set_title(r"$Ev\{x(t)\}$")
ax3.set_title(r"$Od\{x(t)\}$")
plt.xlabel(r"$t$")

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
plt.show()