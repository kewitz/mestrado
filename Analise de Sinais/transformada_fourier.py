# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 15:33:47 2014

@author: leo
"""

from multiprocessing import Pool
import matplotlib.pyplot as plt
from numpy import *

def a(x):
	return -exp(-1j*x)/(-2-1j*x)

def b(x):
	return exp(-1j*x)/(2-1j*x) - exp(-1j*x)/(-2-1j*x)	

if __name__ == '__main__':
	pool = Pool(processes=4)              # start 4 worker processes
	fa = pool.map(a, np.linspace(20*pi,-20*pi,101).tolist())
	fb = pool.map(b, np.linspace(20*pi,-20*pi,101).tolist())
	plt.plot(absolute(fa), label='a)')
	plt.plot(absolute(fb), label='b)')
	plt.legend()