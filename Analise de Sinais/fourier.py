# -*- coding: utf-8 -*-
"""
Created on Tue May 27 16:08:43 2014

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Snippets and Constants
indexed = lambda l, offset=0: zip(np.arange(len(l))+offset,l)
f = np.float64  # Default format.
pi = f(np.pi)  #Double precision pi.
e = np.exp(f(1))  # Exponential as unit.

a0 = 1.0/3.0;
w = (2*pi)/6.0
ak = lambda k: (1/(3*k*w))*np.sin(k*w) + (2/(3j*k*w))*(-exp(-1j*k*w)+2*exp(-
    2j*k*w) - exp(3j*k*w))

dom = np.arange(0,21)
#%%

# Config and plot first frame.
fig, ax = plt.subplots()
line, = ax.plot(dom, np.ones(dom.shape))
plt.ylim(-3,3)
# Animation function
def animate(k):
    xt = zeros(dom.shape) + a0
    for i, x in indexed(dom):
        soma = map(lambda t: ak(k)*exp(1j*k*w*x) + ak(-k)*exp(-1j*k*w*x), range(1,k))
        xt[i] += np.sum(soma)
    line.set_ydata(xt)
    return line,
# Plot animation
ani = animation.FuncAnimation(fig, animate, range(1,10), interval=1000)
plt.show()