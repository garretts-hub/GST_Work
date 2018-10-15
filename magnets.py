# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:11:06 2018

@author: GA28573
"""
import numpy as np
import pylab as plt
from scipy.optimize import curve_fit


u = 4*np.pi*1e-7 #kg*m/s2
R = 0.30 #m
N = 10 #loops
I = 1 #amps
top_loop = 0.08 #m
bottom_loop = -1*top_loop #m

def B_z(abs_z, loop_loc):
    r = R**2 + (abs_z - loop_loc)**2
    return (u*I*N*R**2)/(2*np.pi*r**(3/2))

locs = np.linspace(8*bottom_loop, 8*top_loop, 300)
vals = B_z(locs, bottom_loop) + B_z(locs, top_loop)

def gaussian(x, a, b, c):
    return a*np.exp(-(x-b)**2/(2*c**2))
#p, o = curve_fit(gaussian, locs, vals)
#theory_vals = gaussian(locs, p[0], p[1], p[2])

plt.plot(locs, vals)
#plt.plot(locs, theory_vals)
plt.axvline(x=top_loop, color='orange', ls='dashed')
plt.axvline(x=bottom_loop, color='orange', ls='dashed')
plt.xlabel("Location between the loops, m")
plt.ylabel("z-component of B-field, T")
plt.title("Magnetic field along central axis between two conducting loops")
plt.show()