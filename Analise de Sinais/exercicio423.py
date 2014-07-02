# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 15:41:46 2014

@author: leo
"""
import numpy as np

td = np.arange(0,50,.05)
Xt = map(lambda x: x*exp(-x) if x <= 1 and x >= 0 else 0, td)
#Xjw = fftshift(fft(Xt))
Xjw = map(lambda jw: (-exp(-1.0-jw)/(1.0+jw)) 
				- (-exp(-1.0-jw)/np.power(1.0+jw,2.0))
				+ 1.0/np.power(1+jw,2),
		np.arange(-10*pi,10*pi,.2)*1j)

fsize = (8,4)
picfolder = '/home/leo/Documents/Sinais e Sistemas Lineares/exercicios/'
save = True

figure(figsize=fsize)
plot(td[0:100],Xt[0:100])
if save:
	plt.savefig(picfolder+'4-23da.eps')
	plt.close()

figure(figsize=fsize)
plot(Xjw)
if save:
	plt.savefig(picfolder+'4-23db.eps')
	plt.close()