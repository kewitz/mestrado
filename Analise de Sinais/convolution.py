# -*- coding: utf-8 -*-
"""
Created on Tue May 13 16:59:49 2014

@author: leo
"""

import matplotlib.pyplot as plt
import numpy as np

#Constants
e = np.exp(1)
#Define Time
t = np.arange(-10,10,.05,dtype=np.float64)

#Define Function
def yt(t):
    if t < 1:
        return -e**(2*t-4) + (e**(2*t-10))/2 + (e**(2*t))/2
    elif 1 <= t <= 3:
        return (-e**(2*t-4))/2 + (e**(2*t-10))/2 + (e**2)/2 - (e**(2*t-4))/2
    elif 3 < t <= 6:
        return -(e**2)/2 + (e**(2*t-10))/2
    elif t > 6:
        return 0

y = [yt(a) for a in t]
plt.plot(t, y, 'b')