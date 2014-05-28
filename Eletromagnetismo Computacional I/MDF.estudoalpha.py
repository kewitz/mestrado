# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 17:09:45 2014

@author: leo
"""

import matplotlib.pyplot as plt
import numpy as np
import libs.MDF as m
import multiprocessing

threads = multiprocessing.cpu_count()

def f(alpha):
    test = m.MDF(np.arange(0,10), np.arange(0,10))
    test.delta = 1E-8
    test.alpha = np.float64(alpha)
    test.bounds['left'] = ['n' for a in test.bounds['left']]
    test.bounds['right'] = ['n' for a in test.bounds['right']]
    test.bounds['top'] += 100
    test.run()
    print "Alpha %.2f -> %i iterations." % (alpha, test.k)
    return test.k
    

if __name__ == '__main__':
    multiprocessing.freeze_support()
    alphas = np.arange(1,2,.01)
    pool = multiprocessing.Pool(processes=threads)    # start 4 worker processes
    x = pool.map(f, alphas)          # prints "[0, 1, 4,..., 81]"
    bestA, bestK = (alphas[x.index(min(x))],min(x))
    plt.plot(alphas, x)
    plt.xlabel(r'$\alpha$')
    plt.ylabel('Iterations')
    plt.title(r'$\alpha$ impact on performance.')    
    plt.annotate(r'$\alpha=$'+str(bestA), xy=(bestA, bestK),  xycoords='data',
                xytext=(30, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3,rad=.2")
                )