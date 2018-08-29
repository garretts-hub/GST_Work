# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 15:24:07 2018

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


def gate_string_to_list(gate_string):
    #returns gate sequence as a list
    #ex, 'GxGxGx' --> ['Gx', 'Gx', 'Gx', 'Gx']
    #"Gx(Gy)^2Gz" --> ['Gx', 'Gy', 'Gy', 'Gz']
    #do not raise anything to the 1 exponent, only 2 or more; always use parenthesis
    #do not put more than one gate in parentheses for right now
    gate_list = []
    for i in range(len(gate_string)):
        char = gate_string[i]
        if char == "G":
            gate_list.append('G'+ gate_string[i+1])
        elif char.isnumeric():
            for count in range(int(char)-1):
                gate_list.append('G' + gate_string[i-3])
    return gate_list

def gate_list_to_string(gate_list):
    gate_string = ''
    for gate in gate_list:
        gate_string = gate_string + gate
    return gate_string


def create_data(time_per_shot, num_samples, num_counts, num_shots, gate_list, xerr, yerr, zerr):
    #time_per_shot: time in seconds for a single (prep-gate-measure+delay)
    #num_samples: how many timestamps and strings of counts you want to have
    #num_counts: how many data points to create (how many ones and zeros) per sample (i.e. per timestamp) --> affects both time of a sample and precision
    #num_shots: how many times (shots) you apply (prep-gate_list-meas) to get one count (a single 0 or 1) --> determines the time of one count, but won't affect precision
    #gate_list: gates you want to do for your operation, entered as a list of strings
    #xerr,yerr,zerr: 2D tuples with overrotation amplitude in radians and frequency in Hz
    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
    zero_counts = []
    one_counts = []
    time_per_count = num_shots*time_per_shot #the time required to get a single count (a single 0 or 1) within a given timestamped sample
    timestep = num_counts*time_per_count #the time to get a full bitstring of zeros and ones for one sample, i.e. one timestamp
    timestamps = np.arange(0, num_samples*timestep, timestep) #array of 1*timestep, 2*timestep,....(num_samples)*timestep
    plot_times = []
    probs = []
    angles = []

    for sample in range(num_samples):
        #print("Starting timestamp {:.3f} s".format(time))
        sample_time = sample*timestep
        print("sample time is {}".format(sample_time))
        rho = rho0
        sample_bitlist = []
        
        for count in range(num_counts):
            count_time = count*time_per_count
            count_bitlist = []
            print("Count time is {}".format(count_time))
            
            for shot in range(num_shots):
                shot_time = shot*time_per_shot
                current_time = sample_time + count_time + shot_time
                print("Current shot time is {}".format(current_time))
                plot_times.append(current_time)
                for gate in gate_list:
                    if gate == 'Gx':
                        angle = np.pi/2 + xerr[0]*np.cos(2*np.pi*xerr[1]*current_time)
                        angles.append(angle)
                        rho = _qt.to_super(_qt.rx(angle)) * rho
                    elif gate == 'Gy':
                        angle = np.pi/2 + yerr[0]*np.cos(2*np.pi*yerr[1]*current_time)
                        angles.append(angle)
                        rho = _qt.to_super(_qt.ry(angle)) * rho
                    elif gate == 'Gz':
                        angle = np.pi/2 + zerr[0]*np.cos(2*np.pi*zerr[1]*current_time)
                        angles.append(angle)
                        rho = _qt.to_super(_qt.rz(angle)) * rho
                #print(rho)
                #calculate probabilities of being in 1 after the shot has been applied
                p1 = (rho.dag()*rho1).norm()
                if p1 >= 1:
                    p1 = 1
                elif p1<0:
                    p1 = 0
                #print("*****Time {:.3f}, prob {:.3f}".format(time, p1))
                #print(p1)
                probs.append(p1)
                print("Shot gave probability {:.3f}".format(p1))
                shot_value = np.random.binomial(1, p1) #simulates a summation of the number of 1-counts you get in one bitstring sample
                count_bitlist.append(shot_value)
            
            #determine whether the majority of shots in this count were 0 or 1, and add that bit to the sample's bitstring:
            if count_bitlist.count(1) > count_bitlist.count(0):
                sample_bitlist.append(1)
            else:
                sample_bitlist.append(0)
        
        #go back to the sample bitstring, and count the 1s and 0s occurences at this timestamp
        one_counts.append(sample_bitlist.count(1))
        zero_counts.append(sample_bitlist.count(0))
        
    plt.plot(plot_times, probs)
    plt.ylim(0,1)
    plt.xlabel("Time, seconds")
    plt.ylabel("Probability of Measuring State {1}")
    plt.title("Simulated {} with Drifting Rotation Angle\n $\Theta = \pi/2 \pm {}cos(2\pi *{:.2f}*t)$".format(gate_list_to_string(gate_list), xerr[0], xerr[1]))
    plt.grid()
    plt.show()
        
    return (np.asarray(one_counts), np.asarray(zero_counts), np.asarray(timestamps), probs)



if __name__=='__main__':
    time = 0.016 #time for a single experiment/shot, in seconds
    num_samples = 10   #number of timestamped bitstrings of 0s and 1s
    num_counts = 1       #number of 0s and 1s per sample/timestamp
    num_shots = 2     #number of times you do the experiment to get a single count (a single 0 or 1)
    gate_list = ['Gx']
    xerr = (0.2, 0.8)
    yerr = (0.1, 1e8)
    zerr = (0.1, 1e8)
    counts, zero_counts, times, probs = create_data(time, num_samples, num_counts, num_shots, gate_list, xerr, yerr, zerr)
    #plt.scatter(times, counts)
    #plt.show()
    #print(counts)
    #print(times)
    #print(gate_string_to_list("Gx(GyGi)^2Gz"))
    #plot_probabilities(gate_list, 0, 3, xerr, yerr, zerr, 500)
    #time = np.linspace(0, 1.4, num)
    #plot_rotation_error(time, xerr, yerr, zerr, gate_string)

    
    
    
    