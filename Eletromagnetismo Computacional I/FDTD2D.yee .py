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

from libs.FDTD2D import Yee

def gauss(k, t, fdtd):
    width = (2*np.power(fdtd.tal,2))
    omega = 6*np.pi*fdtd.fop
    func = lambda t: np.sin(omega*t) * np.exp(-np.power(t-2*fdtd.t0,2) / width)
    fdtd.Ez[k,-1,:] = func(t)

# Simulação sem passo mágico.
a = Yee()
a.setFreq(2.4E6)
a.bound['Ez'][0,:] = 0
a.bound['Ez'][-1,:] = 0
a.bound['Ez'][20:50+1,40:60+1] = 0

a.bound['Hx'][0,:] = 0
a.bound['Hx'][-1,:] = 0
a.bound['Hx'][20,40:60+1] = 0
a.bound['Hx'][50,40:60+1] = 0

a.bound['Hy'][:,0] = 0
a.bound['Hy'][:,-1] = 0
a.bound['Hy'][20:50+1,40] = 0
a.bound['Hy'][20:50+1,60] = 0
a.run(gauss)

plt.plot(a.Ez[139,:,40])