# -*- coding: utf-8 -*-
"""
Simulação do erro de dispersão do FDTD quando não considerado o passo mágico em
um impulso gaussiano modulado em seno.

Created on Wed May 28 11:11:30 2014
@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

from libs.FDTD import Yee

def gaussFactory(fdtd):
    """Gera função excitação `f(x)` de um impulso gaussiano modulado em sin."""
    width = (2*np.power(fdtd.tal,2))
    omega = 6*np.pi*fdtd.fop
    return lambda t: np.sin(omega*t)\
            * np.exp(-np.power(t-2*fdtd.tc,2) / width)

# Simulação sem passo mágico.
a = Yee()
a.setFreq(2.4E6)
a.run(gaussFactory(a))

#%% Plot
fig = plt.figure()
ax = p3.Axes3D(fig)

line1, = ax.plot(a.sd, a.et[50,:], zs=0, zdir='z', label='Ex')
line2, = ax.plot(a.sd, a.ht[50,:], zs=0, zdir='y', label='Hy')
ax.legend()

ax.set_xlabel('Z')
ax.set_ylabel('Ex')
ax.set_zlabel('Hy')

ax.set_ylim3d(a.et.min(), a.et.max())
ax.set_zlim3d(a.ht.min(), a.ht.max())

def animate(t):
    line1.set_data(a.sd, a.et[t,:])
    line1.set_3d_properties(0)
    line2.set_data(a.sd, a.ht[t,:])
    line2.set_3d_properties(a.ht[t,:])
    return True
    
# Plota animação
ani = animation.FuncAnimation(fig, animate, range(len(a.td)), blit=False,
                              interval=33)
plt.show()
