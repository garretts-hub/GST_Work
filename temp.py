# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:25:26 2018

@author: GA28573
"""
import pygsti
from pygsti.construction import std1Q_XYI
from pygsti.objects import DataSet
from pygsti.objects import GateString
from create_data_garrett import generate_model_data, generate_model_data_rotation_error, drifted_rotation

import matplotlib.pyplot as plt
import numpy as np
import NoiseSignal as _ns
import qutip as qt


time = np.linspace(0, 1e-6, 100000) #secoinds
xerr = (0.001, 8e8) #amplitude in radians, frequency in Hz
yerr = (0.001, 8e8)
zerr = (0.001, 8e8) 
gate = 'Gx'


def calculate_drift_probabilities(t, xerr, yerr, zerr, gate_of_interest):
    standard_gates = dict()
    standard_gates['Gi'] = qt.to_super(qt.qeye(2))
    standard_gates['Gx'] = qt.to_super(qt.rx(np.pi/2))# + drifted_rotation(t, xerr)))
    standard_gates['Gy'] = qt.to_super(qt.ry(np.pi/2 + drifted_rotation(t, yerr)))
    standard_gates['Gz'] = qt.to_super(qt.rz(np.pi/2 + drifted_rotation(t, yerr)))

    rho0 = qt.operator_to_vector(qt.ket2dm(qt.basis(2,0)))
    rho1 = qt.operator_to_vector(qt.ket2dm(qt.basis(2,1)))
    rho = standard_gates[gate_of_interest] * rho0
    p0 = (rho.dag()*rho0).norm()
    p1 = (rho.dag()*rho1).norm()
    
    return (p0, p1) #zero probability, one probability (in tuple form)


calculate_drift_probabilities(time,xerr,yerr,zerr,gate)

#plt.plot(time, calculate_drift_probabilities(time,xerr,yerr,zerr,gate)[1] )
#plt.show()

