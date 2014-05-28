# -*- coding: utf-8 -*-
"""
Simulação do erro de dispersão do FDTD quando não considerado o passo mágico.
Created on Wed May 28 11:11:30 2014
@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from libs.FDTD import FDTD

a = FDTD()
a.setFreq(2.4E6)
a.dt = a.dt/(1+2e-1)
a.run()

b = FDTD()
b.setFreq(2.4E6)
b.run()

#%% Plot
# Config and plot first frame.
fig, ax = plt.subplots()
line, = ax.plot(a.xd, np.ones(a.xd.shape),'r')
line2, = ax.plot(b.xd, np.ones(b.xd.shape), 'b')
plt.ylim(-1.5,1.5)
plt.xlim(a.xd.min(),a.xd.max())
plt.xlabel(r"$x (m)$")
plt.legend(("Not Magic", "Magic Step"))
# Animation function
def animate(t):
    line.set_ydata(a.st[t,:])
    line2.set_ydata(b.st[t,:])
    return line,
# Plot animation
ani = animation.FuncAnimation(fig, animate, range(len(a.td)), interval=20)
plt.show()