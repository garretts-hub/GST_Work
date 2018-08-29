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

def plot_rotation_error(time_array, xerr, yerr, zerr, selected_gate_str):
    time = time_array #seconds
    gate = selected_gate_str
    rho0 = qt.operator_to_vector(qt.ket2dm(qt.basis(2,0)))
    rho1 = qt.operator_to_vector(qt.ket2dm(qt.basis(2,1)))
    rho = 0
    rho_plain = 0
    zero_probs = []
    one_probs = []
    one_probs_plain = []
    angles = []
    amp = 0
    freq = 0
    for t in time:
        if gate == 'Gx':
            angle = np.pi/2 + drifted_rotation(t, xerr)
            angles.append(angle)
            rho = qt.to_super(qt.rx(angle)) * rho0
            rho_plain = qt.to_super(qt.rx(np.pi/2)) * rho0
            amp = xerr[0]
            freq = xerr[1]
        elif gate == 'Gy':
            angle = np.pi/2 + drifted_rotation(t, yerr)
            angles.append(angle)
            rho = qt.to_super(qt.ry(angle)) * rho0
            rho_plain = qt.to_super(qt.ry(np.pi/2)) * rho0
            amp = yerr[0]
            freq = yerr[1]
        elif gate == 'Gz':
            angle = np.pi/2 + drifted_rotation(t, zerr)
            angles.append(angle)
            rho = qt.to_super(qt.rz(angle)) * rho0
            rho_plain = qt.to_super(qt.rz(np.pi/2)) * rho0
            amp = zerr[0]
            freq = zerr[1]
        p0 = (rho.dag()*rho0).norm()
        p1 = (rho.dag()*rho1).norm()
        p1_plain = (rho_plain.dag()*rho1).norm()
        zero_probs.append(p0)
        one_probs.append(p1)
        one_probs_plain.append(p1_plain)
    
    plot_time = time #plots time in units of seconds -->change the string below accordingly
    plot_time_label = "time, seconds"
    

    plt.plot(plot_time, one_probs, label="With Drift")
    plt.plot(plot_time, one_probs_plain, linestyle='dashed', label="Without drift")
    plt.xlabel(plot_time_label)
    plt.ylabel("probability of being in 1-state")
    plt.title("Time-Dependent Overrotation\nby {} rads at {:.1e} Hz for {}".format(amp, freq, gate))
    plt.legend(loc='lower right')
    plt.grid()
    plt.ylim(0,1)
    plt.show()
    
    '''
    plt.plot(plot_time, angles, label="Overrotated angle")
    plt.plot(plot_time, np.asarray([np.pi/2]*len(plot_time)), linestyle='dashed', label="Expected angle")
    plt.xlabel(plot_time_label)
    plt.ylabel("Rotation angle from gate, radians")
    plt.title("Time-Dependent Overrotation\nby {} rads at {:.1e} Hz for {}".format(amp, freq, gate))
    plt.grid()
    plt.legend(loc="lower right")
    plt.ylim(np.pi/2 - np.pi/3, np.pi/2 + np.pi/3)    
    plt.show()'''
    
    #for i in range(len(time)):
        #print("Time {:.3f} s, Probability {:.3f}".format(time[i], one_probs[i]))
    

if __name__=='__main__':
    time = np.linspace(0, 3, 1000) #seconds
    xerr = (0.2, 8e6) #amplitude in radians, frequency in Hz
    yerr = (0.001, 8e8)
    zerr = (0.001, 8e8) 
    gate = 'Gx'
    plot_rotation_error(time, xerr, yerr, zerr, gate)