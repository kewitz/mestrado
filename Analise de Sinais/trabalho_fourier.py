# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 15:55:18 2014
@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt

# Macros
pi = np.pi; exp = np.exp; arange = np.arange; zeros = np.zeros
indexed = lambda l, offset=0: zip(np.arange(len(l))+offset,l)

# Constantes
w = 2.0*pi*0.25
a0 = 6.0/4.0
# Funções
ak = lambda k: a0 if k == 0 else (1.0/(4*pow(1j*k*w,2))) * (2.0*exp(-3j*k*w)*(3j*k*w + 1) - 2.0*exp(-2j*k*w)*(2j*k*w + 1) - 2.0*exp(-1j*k*w)*(1j*k*w + 1) + 2) \
            + (1.0/(4j*k*w))* (-6.0*exp(-3j*k*w) + 2.0*exp(-2j*k*w) + 4.0*exp(-1j*k*w))
harms = lambda k: arange(-k,k+1)

# Domínio do Tempo
td = arange(0,6,.05)
# X(t) com n harmônicas.
def xt(k):
    xt = zeros(td.shape)
    for n, t in indexed(td):
        xt[n] = np.matrix([ak(h)*exp(1j*h*w*t) for h in harms(k)]).sum()
    return np.abs(xt)


# Plota
def ploth(harm):   
    plt.figure()
    plt.plot(td, xt(harm))
    plt.grid(True)
    plt.title('$x(t)$')
    plt.xlabel('$t$')

ploth(5)
ploth(10)
ploth(20)

plt.figure()
plt.vlines(harms(10), [0], [abs(ak(k)) for k in harms(10)], 'k', lw=2)
plt.plot(harms(10), [abs(ak(k)) for k in harms(10)], 'ko')
plt.xlim(-11,11)

plt.grid(True)
plt.title('$a_k$')
plt.legend()