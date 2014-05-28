# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 16:17:46 2014

@author: leo
"""
import matplotlib.pyplot as plt
import numpy as np
import libs.MDF as m


plt.subplot(2,2,1)
ex1 = m.MDF(np.arange(0,10,.5), np.arange(0,10,.5))
ex1.bounds['left'] = ['n' for a in ex1.bounds['left']]
ex1.bounds['right'] = ['n' for a in ex1.bounds['right']]
ex1.bounds['top'] += 100
ex1.run()
ex1.plot()
plt.title("a)")

plt.subplot(2,2,2)
ex2 = m.MDF(np.arange(0,10,.5), np.arange(0,10,.5))
ex2.alpha = 1.1
ex2.bounds['top'] += 100
ex2.bounds['top'][0:1] = 0
ex2.bounds['top'][19:20] = 0
ex2.run()
ex2.plot()
plt.title("b)")

plt.subplot(2,2,3)
ex3 = m.MDF(np.arange(0,10,.3), np.arange(0,10,.3))
ex3.bounds['top'] = ['n']+[100 for a in range(len(ex3._x))]+['n']
ex3.bounds['left'] = ['n']+[60 for a in range(len(ex3._y))]+['n']
ex3.bounds['bot'] = ['n']+[40 for a in range(len(ex3._x))]+['n']
ex3.bounds['right'] = ['n']+[20 for a in range(len(ex3._y))]+['n']
ex3.run()
ex3.plot()
plt.title("c)")