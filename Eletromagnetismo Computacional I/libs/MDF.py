# -*- coding: utf-8 -*-
"""
The MIT License (MIT)
Copyright (c) 2014 Leonardo Kewitz

Created on Wed Apr 16 15:00:06 2014
"""
import matplotlib.pyplot as plt
import numpy as np

class MDF:
    """Uma classe de cálculo do Método das Diferenças Finitas."""
    def __init__(self,x,y):
        self._x = x
        self._y = y
        self.k = 0
        self.verb = False
        self.alpha = 1.0
        self.delta = 1E-8   #Max value delta for convergence
        self.X, self.Y = np.meshgrid(x, y)
        self.space = np.zeros(self.X.shape, dtype=np.float64)
        self.bounds = {'top':y*0, 'bot':y*0, 'left':x*0, 'right':x*0}   # Use 'n' if Newmann
    
    def run(self,maxiteractions=None):
        this = self.space.copy()
        while True:
            self.iterate()
            d = max(abs((this-self.space)).flatten())
            if d <= self.delta or (maxiteractions != None and self.k >= maxiteractions):                    
                break
            this = self.space.copy()
        if self.verb:
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
                v1 = self.space[ix,iy]
                vt = self.space[ix,iy-1] if iy > min(sy) else self.bound(ix,
                    'top')
                vb = self.space[ix,iy+1] if iy < max(sy) else self.bound(ix,
                    'bot')
                vl = self.space[ix-1,iy] if ix > min(sx) else self.bound(iy,
                    'left')
                vr = self.space[ix+1,iy] if ix < max(sx) else self.bound(iy,
                    'right')
                v = (1-self.alpha)*v1 + self.alpha*((vt+vb+vl+vr)/4.0)
                #v = (vt+vb+vl+vr)/4.0   # Calculate Vk+1
                self.space[ix,iy] = v
        self.k += 1
                
    def plot(self):        
        h = plt.contour(self.X, self.Y, np.rot90(self.space,1))
        plt.clabel(h, inline=1, fontsize=10)