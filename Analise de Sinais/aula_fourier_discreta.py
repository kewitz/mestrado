# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:07:28 2014

@author: leo
"""

w = 2*pi/9
ak = lambda k: (1.0/9) * (3 + exp(-1j*k) + 2*exp(-3j*k) - exp(-4j*k) -2*exp(-6j*k))
harms = lambda k: arange(-k,k+1)

def discplot(x,y, color='k'):
    plt.plot(x, y, 'o'+color)    
    plt.vlines(x, [0], y, ''+color, lw=2)    
    plt.grid(True)

def xn(n):
    return np.matrix([ak(k*w) * exp(1j*k*w*n) for k in range(9)]).sum()

plt.figure()
discplot(range(20), [xn(n) for n in range(20)], 'r')

#%%
ak1 = lambda k: (2.0/15) * sin(2*pi*k*5.0/30)/sin(pi*k/15)
ak2 = lambda k: (4.0/15) * sin(2*pi*k*5.0/30)/sin(pi*k/15) * exp(-5j*k*2*pi/15)
ak3 = lambda k: (-3.0/15) * sin(2*pi*k*3.0/30)/sin(pi*k/15) * exp(-9j*k*2*pi/15)

ak = lambda k: 1.4 if k == 0 else ak1(k)+ak2(k)+ak3(k)

xn1 = lambda n: np.matrix([ak(k) * exp(1j*k*(2.0*pi/15)*n) for k in range(15)]).sum()
plt.figure()
discplot(range(20), [xn1(n) for n in range(20)], 'r')