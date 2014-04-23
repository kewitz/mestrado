# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 15:00:06 2014
@author: leo
"""
import matplotlib.pyplot as plt
import numpy as np


class MDF:
    """Uma classe de cálculo do Método das Diferenças Finitas."""
    alpha = 1
    def __init__(self,x,y):
        self._x = y
        self._y = x
        self.k = 0
        self.delta = 0.01   #Max value delta for convergence
        self.X, self.Y = np.meshgrid(x, y)
        self.space = np.zeros(self.X.shape)
        self.bounds = {'top':y*0, 'bot':y*0, 'left':x*0, 'right':x*0}   # Use 'n' if Newmann
    
    def run(self,iteractions=None):
        if iteractions:           
            for k in range(iteractions):
                self.iterate()
        else:
            this = self.space.copy()
            while True:
                self.iterate()
                d = max(abs((this-self.space)).flatten())
                if d<self.delta:                    
                    break
                this = self.space.copy()
            print "Converged in %d iteractions" % self.k
            
    def bound(self,i,b):
        r = self.bounds[b][i]
        if r == 'n':
            if b == 'top':
                return self.space[i,1]
            elif b == 'bot':
                return self.space[i,self.space.shape[1]-1]
            elif b == 'left':
                return self.space[1,i]
            elif b == 'right':
                return self.space[self.space.shape[0]-1,i]
        return r
        
    def iterate(self):
        sx = range(self.space.shape[0])
        sy = range(self.space.shape[1])
        for ix in sx:
            for iy in sy:
                vt = self.space[ix,iy-1] if iy > min(sy) else self.bound(ix,
                    'top')
                vb = self.space[ix,iy+1] if iy < max(sy) else self.bound(ix,
                    'bot')
                vl = self.space[ix-1,iy] if ix > min(sx) else self.bound(iy,
                    'right')
                vr = self.space[ix+1,iy] if ix < max(sx) else self.bound(iy,
                    'left')
                #v = v1 + (self.alpha/4.0) * ( ((vt+vb+vl+vr)/4.0) - 4.0*v1 )
                v = (vt+vb+vl+vr)/4.0   # Calculate Vk+1
                self.space[ix,iy] = v
        self.k += 1
                
    def plot(self):
        #h = plt.contour(self.X,self.Y,np.rot90(self.space))
        h = plt.imshow(np.rot90(self.space, 3), interpolation='nearest', 
                       cmap=cm.afmhot)
        cbar = plt.colorbar(h)
        #plt.clabel(h, inline=1, fontsize=10)
        
plt.subplot(2,2,1)
ex1 = MDF(np.arange(0,10,.5), np.arange(0,10,.5))
ex1.bounds['left'] = ['n' for a in ex1.bounds['left']]
ex1.bounds['right'] = ['n' for a in ex1.bounds['right']]
ex1.bounds['top'] += 100
ex1.run()
ex1.plot()
plt.title("a)")

plt.subplot(2,2,2)
ex2 = MDF(np.arange(0,10,.5), np.arange(0,10,.5))
ex2.bounds['top'] += 100
ex2.bounds['top'][0:1] = 0
ex2.bounds['top'][19:20] = 0
ex2.run()
ex2.plot()
plt.title("b)")

plt.subplot(2,2,3)
ex3 = MDF(np.arange(0,10,.5), np.arange(0,10,.5))
ex3.bounds['top'][2:18] = 100
ex3.bounds['left'][2:18] = 60
ex3.bounds['bot'][2:18] = 40
ex3.bounds['right'][2:18] = 20
ex3.run()
ex3.plot()
plt.title("c)")