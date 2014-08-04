# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:09:54 2014

@author: leo
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

irange = lambda a: zip(range(len(a)),a)

def alphas(ns,lmn,gm):
	"""
	Retorna os valores de alpha dado a quantidade de `ns` e as funções 
	lambdas: `lmn` e `gm`.
	
	Parâmetros
	----------
	ns : int
		Número de N.
	lmn, gm : string
		Função utilizada para definir `l[m,n]` e `g[m]` respectivamente.
		
	Exemplo
	-------
	>>> lmn = "( (m/n)**(n-1) )*n*(n+1)"
	>>> gm = "4 - 12*( (m / N)**2 )"
	>>> alphas(2,lmn,gm)
	matrix([[ 5.],
		[-3.]])
	"""
	x, m, n, N = sp.symbols("x m n N")
	ms = ns
		
	lmn = sp.sympify(lmn.replace("N","ns")) if type(lmn) == str else lmn
	gm = sp.sympify(gm.replace("N","ns")) if type(gm) == str else gm
	
	assert dir(lmn).count('evalf') == 1, "lmn must be string or sympy.core.mul.Mul"
	assert dir(gm).count('evalf') == 1, "gm must be string or sympy.core.mul.Mul"
	
	L = np.matrix(np.zeros((ms,ns)))
	G = np.matrix(np.zeros((ns)))
	for i in range(ms):
		m = float(i+1.0)	
		G[0,i] = gm.evalf(subs={'m':m,'ns':ns})
		for j in range(ns):
			n = float(j+1.0)			
			L[i,j] = lmn.evalf(subs={'m':m,'ns':ns,'n':n})
	return np.linalg.solve(L,G.T)

def solve(alphas,fx,domain):
	"""
	Resolve a função fx para o domínio especificado.
	
	Example
	-------
	
	>>> lmn = "( (m/n)**(n-1) )*n*(n+1)"
	>>> gm = "4 - 12*( (m / N)**2 )"
	>>> x = np.arange(0,1,.1)
	>>> p1 = solve(alphas(1,lmn,gm),"x*(1-x**n)",x)
	array([ 0.  , -0.36, -0.64, -0.84, -0.96, -1.  , -0.96, -0.84, -0.64, -0.36])
	"""
	ns = len(alphas)
	func = sp.sympify(fx.replace("N","ns"))
	result = np.zeros(domain.shape)
	for ix,x in irange(domain):
		result[ix] = float()
		for ia, a in irange(alphas):
			result[ix] += a*func.evalf(subs={'ns':ns,'n':ia+1,'x':x})
	return result

	
# Exemplo 1:
#lmn = "(m*n) / (m+n+1)"
#gm = "( m*(m+1) ) / ( (m+2)*(m+3) )"
#print alphas(2,lmn,gm)

#% Exemplo 2:
lmn = "( (m/n)**(n-1) )*n*(n+1)"
gm = "4 - 12*( (m / N)**2 )"
x = np.arange(0,1,.01)
p1 = solve(alphas(1,lmn,gm),"x*(1-x**n)",x)
p2 = solve(alphas(2,lmn,gm),"x*(1-x**n)",x)
p3 = solve(alphas(3,lmn,gm),"x*(1-x**n)",x)
plt.plot(x,p1,x,p2,x,p3)