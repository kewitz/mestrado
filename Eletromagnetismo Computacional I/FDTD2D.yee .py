# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Simulação do erro de dispersão do FDTD quando não considerado o passo mágico em
um impulso gaussiano modulado em seno.

Created on Wed May 28 11:11:30 2014
@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

from src.Yee2D.fastYee import FYee

def gauss(k, t, fdtd):
	width = (2*np.power(fdtd.tal,2))
	omega = 6*np.pi*fdtd.fop
	func = lambda t: np.exp(-np.power(t-2*fdtd.t0,2) / width)
	fdtd.Ez[k,1,:] = func(t)

a = FYee()
a.setFreq(2.4E9)

a.bound['Ez'][0,:] = 0
a.bound['Ez'][-1,:] = 0
a.bound['Ez'][20:50+1,40:60+1] = 0

a.bound['Hx'][0,:] = 0
a.bound['Hx'][-1,:] = 0
a.bound['Hx'][:,0] = 0
a.bound['Hx'][:,-1] = 0
a.bound['Hx'][20,40:60+1] = 0
a.bound['Hx'][50,40:60+1] = 0

a.bound['Hy'][20:50+1,40] = 0
a.bound['Hy'][20:50+1,60] = 0

a.run(gauss,t=3000)

#%%Plot
fig = plt.figure()
ims = []


for k in a.Ez[::15,:,:]:
    im = plt.imshow(k)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=30, blit=True, repeat_delay=0)

plt.show()

#%% Save Plot
#Writer = animation.writers['mencoder_file']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
#ani.save('img.mp4', writer=writer)