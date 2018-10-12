# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 16:07:10 2018

@author: GA28573
"""

import numpy as np
from numpy import sin, cos, pi
import pylab as plt



def fit(X, A, Y, W, T):
    value = T + (A*X)*(Y*cos(Y*T)*sin(W*T) - W*cos(W*T)*sin(Y*T))/(W**2 - Y**2)
    return value


X = np.linspace(0, 8000, 100)
A = 0.1
Y = 2*pi*1.5 + 1
W = 2*pi*1.5
T = 60

plt.plot(X, fit(X, A, Y, W, T))
plt.grid()
plt.show()