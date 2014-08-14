# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:09:54 2014

@author: leo
"""
import numpy as np
import matplotlib.pyplot as plt

from libs.MOM import MOM

# Constantes
x = np.arange(0,1,.01)  # Domínio de x

#%% Exemplo 1:
lmn = "(m*n) / (m+n+1)"
gm = "( m*(m+1) ) / ( (m+2)*(m+3) )"
E1 = MOM(lmn,gm,"x*(1-x**n)")
print E1.getAlphas(2).alphas

#%% Exemplo 2:
lmn = "( (m/n)**(n-1) )*n*(n+1)"
gm = "4 - 12*( (m / N)**2 )"

E2 = MOM(lmn,gm,"x*(1-x**n)")
p1 = E2.getAlphas(1).solveFor(x)
p2 = E2.getAlphas(2).solveFor(x)
p3 = E2.getAlphas(3).solveFor(x)
plt.figure(1)
plt.plot(x,p1,x,p2,x,p3)

#%% Exemplo 3:
lmn = "(m*n*(3*m*n + 7*m + 7*n + 15))/( 3*(m+3)*(n+3)*(m+n+3) )"
gm = "(-8*m*(2*m+1))/(15*(m+3)*(m+5))"
E3 = MOM(lmn,gm,"x*(1-x**n)")
p2 = E3.getAlphas(2).solveFor(x)
p3 = E3.getAlphas(3).solveFor(x)
p4 = E3.getAlphas(4).solveFor(x)
plt.figure(2)
plt.plot(x,p2,x,p3,x,p4)