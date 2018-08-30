# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 14:56:03 2018

@author: GA28573
"""

import numpy as np
import pylab as plt
from scipy.fftpack import dct, idct

def dct_stats(start, end, num_points, freq,  nCounts):
    time = np.linspace(start, end, num_points)
    x = np.cos(2*np.pi*time)
    null_hypothesis = np.mean(x)/nCounts
    input_array = (x - nCounts*null_hypothesis)/np.sqrt(nCounts*null_hypothesis * (1 - null_hypothesis))
    modes = dct(input_array, norm='ortho')
    total_power_plotted = np.sum(modes)
    peak_power = np.max(modes)
    return peak_power, total_power_plotted

def make_binary_data(start, stop, points, freq, nCounts):
    time = np.linspace(start, end, points)
    probs = 0.5*np.cos(2*np.pi*time) + 0.5
    binary_points = np.random.binomial(nCounts, probs)
    return binary_points

time = np.linspace(0, 4, 200)
nCounts = 1
plt.plot(time, make_binary_data(0,4,200,0.25, 5))
plt.show()

#x = 0.5*np.cos(2*np.pi*time) + 0.5
x = np.cos(2*np.pi*time)
null_hypothesis = np.mean(x)/nCounts
input_array = (x - nCounts*null_hypothesis)/np.sqrt(nCounts*null_hypothesis * (1 - null_hypothesis))




start = 0
end = 4
freq = 0.25
peak_power = []
total_power = []
fractions = []
nrange = range(1,600)
for N in nrange:
    peak, total = dct_stats(start, end, N, freq, 1)
    peak_power.append(peak)
    total_power.append(total)
    fractions.append(peak/total)


