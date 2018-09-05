# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:07:03 2018

@author: GA28573
"""

import sys
#sys.path.append("../../transmon-sim/physical_sim")
sys.path.append("C:/Users/GA28573/AppData/Local/Continuum/anaconda32/Lib/site-packages/qutip-4.3.1-py3.6-win-amd64.egg/qutip")
import NoiseSignal2 as _ns
import qutip as _qt
import numpy as np
import pylab as plt
from helpers import *

def create_data(time_per_count, num_samples, num_counts, gate_list, time_unit, noise_type=None, walking_amp=None, telegraph_amp=None, \
                res=None, freq_list=None, amp_list=None, phase_list=None, start_f=None, stop_f=None, fluctuators=None, plot_noise=False, \
                add_noise=False, noise_object=None, dc_angle_offset=0, constant_linear_drift=0):
    #time_per_shot: time in seconds for a single (prep-gate-measure+delay)
    #num_samples: how many timestamps and strings of counts you want to have
    #num_counts: how many data points to create (how many ones and zeros) per sample (i.e. per timestamp) --> affects both time of a sample and precision
    #num_shots: how many times (shots) you apply (prep-gate_list-meas) to get one count (a single 0 or 1) --> determines the time of one count, but won't affect precision
    #gate_list: gates you want to do for your operation, entered as a list of strings
    #xerr,yerr,zerr: 2D tuples with overrotation amplitude in radians and frequency in Hz
    #constant linear drift: enter in rads/second
    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
    zero_counts = []
    one_counts = []
    timestep = num_counts*time_per_count #the time to get a full bitstring of zeros and ones for one sample, i.e. one timestamp
    timestamps = np.arange(timestep, num_samples*timestep, timestep) #array of 1*timestep, 2*timestep,....(num_samples)*timestep
    probs = []
    total_time = (time_per_count*num_counts*num_samples)/time_unit
    
    sig = 0
    if noise_type == "Sine":
        #this function returns the noise object so you can enter it back in as a parameter
        # in the event that you call the function repeatedly for an identical set of parameters
        if noise_object != None:
            sig = noise_object
            #print("REUSING NOISE OBJECT")
            while total_time > sig.times[-1] + timestep:
                sig.next_interval()
                #print("Doing a NEXT INTERVAL")
        else:
            #print("INITIALIZING NEW NOISE")
            sig = _ns.NoiseSignalSine(time_unit=time_unit)
            sig.configure_noise(resolution_factor=res, freq_list=freq_list, amp_list=amp_list, phase_list=phase_list, total_time=total_time)
            sig.init()
            if add_noise != None:  
                sig.add_random_noise(add_noise) #add normal noise with specified std deviation if requested
    elif noise_type == "Random Walk":
        sig = _ns.NoiseSignalRandomWalk(initial_seed=1234, time_unit=time_unit)
        sig.configure_noise(walking_amp, res, total_time)
        sig.init() 
    elif noise_type == "Telegraph":
        sig = _ns.NoiseSignalTelegraph(initial_seed=1234, time_unit=time_unit)
        sig.configure_noise(exponent=1, amplitude=telegraph_amp, total_fluctuators=fluctuators, start_freq=start_f, stop_freq=stop_f, total_time=total_time)
        sig.init()
        sig.interpolation_settings(do_interpolation=True, resolution_factor=res)
    
    if plot_noise==True:
        sig.plot_noise_signal()
        
    angle_list = []
    expected_angle_list = []
        
    compressed_gate_list = compress_gate_list(gate_list)
    for time in timestamps:
        noise_at_time = 0
        if noise_type != None:
            #print(time/time_unit)
            noise_at_time = sig[time/time_unit]
        rho = rho0
        total_angle = 0
        total_ideal_angle = 0
        
        for gate in compressed_gate_list:
            gate_name = gate[0]
            gate_repetitions = gate[1]
            '''
            Next step: calculate change in rotation error within each shot. Currently takes the time at the start of the experiment shot
            and applies that to all gates in one shot. Depending on the timescale of the error and time per shot, this simplification may need
            to be addressed so that each gate, say each Gx in (Gx)^11, has an error associated with its specific time, not the same error for
            all 11 Gx gates.
            '''
            
            if gate_name == 'Gx':
                angle = (np.pi/2 + noise_at_time)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                ideal_angle = (np.pi/2)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                rho = (_qt.to_super(_qt.rx(angle))) * rho
                #this section just keeps the angle between 0 and pi
                angle = angle % (2*np.pi)
                ideal_angle = ideal_angle % (2*np.pi)
                if angle > np.pi:
                    angle = 2*np.pi - angle
                if angle < 0:
                    angle = 0 + abs(angle)
                if ideal_angle > np.pi:
                    ideal_angle = 2*np.pi - ideal_angle
                if ideal_angle < 0:
                    ideal_angle = 0 + abs(ideal_angle)
            elif gate_name == 'Gy':
                angle = (np.pi/2 + noise_at_time)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                ideal_angle = (np.pi/2)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                rho = (_qt.to_super(_qt.rx(angle))) * rho
                #this section just keeps the angle between 0 and pi
                angle = angle % (2*np.pi)
                ideal_angle = ideal_angle % (2*np.pi)
                if angle > np.pi:
                    angle = 2*np.pi - angle
                if angle < 0:
                    angle = 0 + abs(angle)
                if ideal_angle > np.pi:
                    ideal_angle = 2*np.pi - ideal_angle
                if ideal_angle < 0:
                    ideal_angle = 0 + abs(ideal_angle)
            elif gate_name == 'Gz':
                angle = (np.pi/2 + noise_at_time)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                ideal_angle = (np.pi/2)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                rho = (_qt.to_super(_qt.rx(angle))) * rho
                #this section just keeps the angle between 0 and pi
                angle = angle % (2*np.pi)
                ideal_angle = ideal_angle % (2*np.pi)
                if angle > np.pi:
                    angle = 2*np.pi - angle
                if angle < 0:
                    angle = 0 + abs(angle)
                if ideal_angle > np.pi:
                    ideal_angle = 2*np.pi - ideal_angle
                if ideal_angle < 0:
                    ideal_angle = 0 + abs(ideal_angle)
            elif gate_name == "Gi":
                #apply only the oscillating drift angle, don't add it to pi/2
                angle = (noise_at_time)*gate_repetitions + dc_angle_offset + constant_linear_drift*time
                ideal_angle = 0 + dc_angle_offset + constant_linear_drift*time
                rho = (_qt.to_super(_qt.rx(angle))) * rho
                angle = angle % (2*np.pi)
                ideal_angle = ideal_angle % (2*np.pi)
                if angle > np.pi:
                    angle = 2*np.pi - angle
                if angle < 0:
                    angle = 0 + abs(angle)
                if ideal_angle > np.pi:
                    ideal_angle = 2*np.pi - ideal_angle
                if ideal_angle < 0:
                    ideal_angle = 0 + abs(ideal_angle)
            total_angle += angle
            total_ideal_angle += ideal_angle
        
        #append the total rotation angle after timestamp 'time'
        angle_list.append(total_angle)
        expected_angle_list.append(total_ideal_angle)        
        #calculate probabilities of being in 1 after the experiment has been applied
        p1 = (rho.dag()*rho1).norm()
        #fix the p1 if it exceeds 1 due to rounding error
        if p1 > 1:
            #print("p1 exceeds 1 for time {}".format(time))
            #print("prob is {}".format(p1))
            #print("Resetting to {}".format(2-p1))
            p1 = 2 - p1
        probs.append(p1)
        one_count = np.random.binomial(num_counts, p1) #simulates a summation of the number of 1-counts you get in one bitstring sample
        zero_count = num_counts - one_count #simulates summation of the number of 0-counts in one bitstring sample
        one_counts.append(one_count)
        zero_counts.append(zero_count)
            
    
    if plot_noise == True:
        plt.plot(timestamps, np.asarray(expected_angle_list),label="Ideal Angle",ls='dashed',color='orange')
        plt.plot(timestamps, np.asarray(angle_list), label="Drifting Angle")
        plt.legend()
        plt.xlabel("Time, seconds")
        plt.yticks(np.linspace(0, np.pi, 5), ['0', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'])
        plt.ylabel("Angle, radians\n(Displayed between 0 and $\pi$)")
        plt.title("Time-Dependent Rotation Angle after each {}".format(gate_list_to_string(gate_list)))
        plt.grid()
        plt.show()        
        
        plt.plot(timestamps, probs)
        plt.ylim(0,1)
        plt.xlabel("Time, seconds")
        plt.ylabel("Probability of Measuring State {1}")
        plt.title("Simulated {} with {} Noise".format(gate_list_to_string(gate_list), noise_type))
        plt.grid()
        plt.show()
        
    return (np.asarray(one_counts), np.asarray(zero_counts), np.asarray(timestamps), probs,np.asarray(expected_angle_list), np.asarray(angle_list),  sig)

def ramsey_experiment(left_gates, right_gates, L, detuning, nSamples, nCounts, time_per_gate, switching_time, time_units, noise_type, freq_list,\
                      amp_list, phase_list, start_f=0, stop_f=0, fluctuators=0,plot_noise=False,add_noise=False,noise_object=None):



if __name__=='__main__':

    gate_string = "GxGxGxGx"                                
    #print("Start with string {}".format(gate_string))
    #gate_list = gate_string_to_list(gate_string)
    #print("Input list of length {}".format(len(gate_list)))
    #print("Compressed form: {}".format(compress_gate_list(gate_list)))
    #print("Converted back to string: {}".format(gate_list_to_string(gate_list)))
    nSamples = 4  #total samples to take for each measurement
    nCounts = 1      #total shots to take at one; =nSamples: same noise, probabilities for all repeats; =1, new experiment & noise for each count
    time_per_count = 1/60 #seconds
    time_units = 1e-3 #seconds
    amp = 0.05
    noise_type='Sine' #Sine, Random Walk, Telegraph
    plot_noise=True
    walking_amp = 0.001
    telegraph_amp = 0.02
    res = 4
    low_frequency_noise = [np.round(0.0015*i, 4) for i in range(5)]
    low_frequency_noise = [] #uncomment this to eliminate low_frequency noise
    low_frequency_amps = [0.005*i for i in range(len(low_frequency_noise))]
    low_frequency_phase = [0]*len(low_frequency_noise)
    freq_list=[0.5]#1.2, 6, 8.4, 9.6] + low_frequency_noise
    amp_list=[0.00]#.002, 0.002, 0.0015, 0.0015] + low_frequency_amps
    phase_list=[0]#,0,0,0] + low_frequency_phase
    dc_angle_offset = 0
    constant_linear_drift = 0
    add_noise=False#0.005
    start_f = 0.1
    stop_f = 2
    fluctuators= 40

    freq_list=tuple(freq_list)
    amp_list = tuple(amp_list)
    phase_list= tuple(phase_list)
    #create_1f_data(time_per_count, nSamples, nCounts, gate_list, amp, fluctuators, start_f, stop_f, time_units)
    ones, zeros, times, probs, expected_angles, angles, sig = create_data(time_per_count, nSamples, nCounts, gate_list, time_units, noise_type, walking_amp, telegraph_amp, \
                res, freq_list, amp_list, phase_list, start_f, stop_f, fluctuators,plot_noise,add_noise,noise_object=None,dc_angle_offset=dc_angle_offset, constant_linear_drift=constant_linear_drift)
    
    
    #print("times have {} points".format(len(times)))
    #print("ones list has {} points".format(len(ones)))
    #print("there were {} input points".format(nSamples))
    print("Max probability is {}".format(max(probs)))
    print("Min probability is {}".format(min(probs)))
    print("Peak-to-Trough Amplitude of probability oscillation is {}".format(max(probs) - min(probs) ))
    
    plt.figure(figsize=(10,4))
    plt.plot(times, ones/nCounts,marker='.')
    plt.ylim(0,1)
    plt.title("1-state Probability Using Simulated Data\nDrift Freq {} Hz at {} rads, averaging {}-counts per point".format(freq_list,amp_list,nCounts))
    plt.ylabel("Probability\n(Total 1-counts/total counts per sample)")
    plt.show()
