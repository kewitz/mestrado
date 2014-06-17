# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 11:08:26 2014

@author: leo
"""

from multiprocessing import Pool, Array

def f(x):
	i = x[0]
	return x[1]

if __name__ == '__main__':
	pool = Pool(processes=4)              # start 4 worker processes
	ar = np.ndarray((3,3,3))
	r = pool.map(f, zip(range(10),[ar]*10))
