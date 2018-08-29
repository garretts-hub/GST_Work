# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:18:37 2018

@author: GA28573
"""


import sys
#sys.path.append("../../transmon-sim/physical_sim")
sys.path.append("C:/Users/GA28573/AppData/Local/Continuum/anaconda32/Lib/site-packages/qutip-4.3.1-py3.6-win-amd64.egg/qutip")


import warnings
import pygsti
from pygsti.objects import dataset as _ds
import Tomography as _tm
import NoiseSignal as _ns
import qutip as _qt
import numpy as np
import pylab as plt
from overrotation_probabilities import plot_rotation_error


rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
time_per_trial = 0.012 #s
num_trials = 10
time = np.arange(0, num_trials*time_per_trial, time_per_trial)
probs = []
angles = []

xerr = (0.2, 0.8)
yerr = (0.1, 1e8)
zerr = (0.1, 1e8)


for t in time:
    angle = np.pi/2 + xerr[0]*np.cos(2*np.pi*xerr[1]*t)
    angles.append(angle)
    rho = _qt.to_super(_qt.rx(angle)) * rho0
    p1 = (rho.dag()*rho1).norm()
    probs.append(p1)
    #print("Time {:.4f}, Prob {:.4f}, Angle {:.4f}".format(t, p1, angle))
    
plt.scatter(time, probs)
plt.grid()
plt.title("Probabilities")
plt.ylim(0,1)
plt.show()

'''plt.scatter(time, angles)
plt.grid()
plt.title("Angles")
plt.ylim(0, np.pi)
plt.show()'''


angle = np.pi/2 #+ xerr[0]*np.cos(2*np.pi*xerr[1]*t)
gx = _qt.to_super(_qt.rx(angle))
rho = gx * gx * rho0
p1 = (rho.dag()*rho1).norm()
print(p1)