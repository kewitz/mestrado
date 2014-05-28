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
from libs.FDTD import FDTD

def gaussFactory(fdtd):
    """Gera função excitação `f(x)` de um impulso gaussiano modulado em sin."""
    width = (2*np.power(fdtd.tal,2))
    omega = 6*np.pi*fdtd.fop
    return lambda t: np.sin(omega*t)\
            * np.exp(-np.power(t-2*fdtd.tc,2) / width)

# Simulação sem passo mágico.
a = FDTD()
a.setFreq(2.4E6)
a.dt = a.dt/(1+2e-1)
a.run(gaussFactory(a))

# Simulação com passo mágico.
b = FDTD()
b.setFreq(2.4E6)
b.run(gaussFactory(b))

# Cria plot e as primeiras linhas.
fig, ax = plt.subplots()
line, = ax.plot(a.xd, np.ones(a.xd.shape),'r')
line2, = ax.plot(b.xd, np.ones(b.xd.shape), 'b')

#Configura legendas e limites.
plt.ylim(-1.5,1.5)
plt.xlim(a.xd.min(),a.xd.max())
plt.xlabel(r"$x (m)$")
plt.legend(("Not Magic", "Magic Step"))

# Animação das linhas.
def animate(t):
    line.set_ydata(a.st[t,:])
    line2.set_ydata(b.st[t,:])
    return line,
    
# Plota animação
ani = animation.FuncAnimation(fig, animate, range(len(a.td)), interval=20)
plt.show()