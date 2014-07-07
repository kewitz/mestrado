# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 15:03:33 2014

@author: leo
"""
import numpy as np
pi = np.pi
exp = np.exp

# Exercício em sala TFTD 1
def hjw (w):
	if abs(w) < pi/3: return 2
	elif (abs(w) <= (2*pi/3)): return 1.0 - abs(w) + pi/3
	else: return 0
xjw = lambda w: -exp(-3j*w)+2*exp(-1j*w)-1 
wd = np.arange(-pi,pi,.001)
resultado = np.array(map(hjw,wd)) * xjw(wd)

# Exercício em sala TFTD 2
hjw2 = lambda w: 16.0/(8 - (6*exp(-1j*w)) + exp(-2j*w))
wd2 = np.arange(-pi,pi,.001)
wd2phase = np.arange(-2*pi,2*pi,.001)
r2 = hjw2(wd2)
r3 = np.angle(hjw2(wd2phase))

# Plot
fsize = (8,4)
picfolder = '/home/leo/Documents/Sinais e Sistemas Lineares/exercicios/'
save = True
ex1 = False
ex2 = True

if ex1:
	plt.figure(figsize=fsize)
	plot(wd,np.abs(resultado))
	plot(wd,np.abs(xjw(wd)), '-.')
	plot(wd,np.abs(np.array(map(hjw,wd))),'--')
	plt.xticks( [-2*pi/3,-pi/3, 0, pi/3, 2*pi/3], ['$-2\pi / 3$','$-\pi / 3$','$0$','$\pi/3$','$2\pi/3$'])
	plt.xlabel('$w$')
	plt.legend(('a','b','c'))
	plt.grid()
	if save:
		plt.savefig(picfolder+'tftd-1.eps')
		plt.close()
if ex2:
	plt.figure(figsize=fsize)	
	plot(wd2,np.abs(r2))
	plt.xticks( [-2*pi,-pi, 0, pi, 2*pi], ['$-2\pi$','$-\pi$','$0$','$\pi$','$2\pi$'])
	if save:
		plt.savefig(picfolder+'tftd-21.eps')
		plt.close()
	
	plt.figure(figsize=fsize)
	plot(wd2phase,r3)
	plt.xticks( [-2*pi,-pi, 0, pi, 2*pi], ['$-2\pi$','$-\pi$','$0$','$\pi$','$2\pi$'])
	plt.xlabel('$w$')
	plt.grid()
	if save:
		plt.savefig(picfolder+'tftd-22.eps')
		plt.close()