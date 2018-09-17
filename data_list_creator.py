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
                rho = (_qt.to_super(_qt.ry(angle))) * rho
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
                rho = (_qt.to_super(_qt.rz(angle))) * rho
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
                rho = (_qt.to_super(_qt.rz(angle))) * rho
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
        
        plt.figure(figsize=(15,3))
        plt.plot(timestamps, probs)
        plt.ylim(0,1)
        plt.xlim(timestamps[0], timestamps[-1])
        plt.xlabel("Time, seconds")
        plt.ylabel("Probability of Measuring State {1}")
        plt.title("Simulated {} with {} Noise".format(gate_list_to_string(gate_list), noise_type))
        plt.grid()
        plt.show()
        
    return (np.asarray(one_counts), np.asarray(zero_counts), np.asarray(timestamps), probs,np.asarray(expected_angle_list), np.asarray(angle_list),  sig)

def ramsey_experiment(left_gate_list, right_gate_list, L, field_f, transition_f, nCounts, time_per_gate, gate_switching_time, \
                      experiment_sample_time, time_units=1e-6, noise_type=None, freq_list=0,amp_list=0, phase_list=0,\
                      start_f=0, stop_f=0, fluctuators=0, plot_noise = False):
    #left and right gate_lists: lists of the gates sandwiching the Gi gates
    #L: the number of Gi gates. If you want to vary this, make it a list
    #field_f: frequency of external field in Hz. To vary this, make it a list. Either L or field_f must vary.
    #transition_f: frequency of the qubit transition in Hz.
    #time_per_gate: total time to complete one Gx, Gy, Gi, or Gz
    #switching_time: additional time before you start any new gate or change gates
    #experiment_sample_time: time at which you trigger a single expeimrent (single count), usually 60 Hz triggered
    #nCounts: number of times to repeat one experiment (a single set of parameters)
    #time_units: baseline units (in seconds) for any additional drift noise signal you may want to add. Defaults to ms.
    
    #check that the input has something to be varied
    if type(L) == type(field_f):
        print("Error: Either L or field_f must vary over a range!")
        return None
    #this list will contain the varied parameter, either detuning in Hz or delay time in s
    varied_param = []

    if type(L) == list:
        total_experiments = len(L)
        experiment_list = []
        for l in L:
            experiment_list.append(left_gate_list + ['Gi']*l + right_gate_list)
        field_f_list = [field_f]*total_experiments
        varied_param = [(time_per_gate*l + gate_switching_time) for l in L]
        
    else:
        total_experiments = len(field_f)
        experiment_list = [left_gate_list + ['Gi']*L + right_gate_list]*total_experiments
        field_f_list = field_f
        varied_param = delta_list
        
    #create a noise object:
    if noise_type != None and noise_type == "Sine":
        total_time = total_experiments*nCounts*experiment_sample_time + 1#assumes that every count is taken every 0.016 seconds; no experiment takes longer
        sig = _ns.NoiseSignalSine(time_unit=time_units)
        sig.configure_noise(resolution_factor=1, freq_list=freq_list, amp_list=amp_list, phase_list=phase_list, total_time=total_time/time_units)
        sig.init()
        
        if plot_noise==True:
            sig.plot_noise_signal()
        
    probs = [] #a list of total_expeirment* elements. Has the probability for each experiment set (assumes constant p for all nCounts)
    time_per_experiment_set = [] #has the time at the start of each experiment set
    ones_counts = [] #number of 1s counted for each experiment set
    angles = [] #total theta rotation for each experiment set
    ideal_angles = [] #has the ideal theta rotation for each experiment set
    transition_f_list = [] #list of the transition frequency at the start of each experiment set
    detuning_list = [] #list of the detuning at the start of each experiment set
    absolute_time = 0 #seconds
    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
    
    for experiment_index in range(total_experiments):
        experiment = experiment_list[experiment_index]
        compressed_experiment = compress_gate_list(experiment)
        #the following line assumes that each count for each experiment set is triggered at 0.0167 seconds (no single gate sequence can be longer than 1/60)
        absolute_time += experiment_sample_time*nCounts
        #print("Starting experiment set {} at {} s".format(experiment_index, absolute_time))
        rho = rho0
        total_angle = 0
        total_ideal_angle = 0
        modified_transition_f = transition_f
        detuning = field_f_list[experiment_index] - modified_transition_f
        for gate_name, repetitions in compressed_experiment:
            if noise_type != None:
                if absolute_time >= total_time:
                    print("Abs time: {}, total_time: {}".format(absolute_time, total_time))
                detuning_noise = sig[absolute_time/time_units]
            else:
                detuning_noise = 0
            if gate_name == 'Gx':
                angle = (np.pi/2)*repetitions
                ideal_angle = (np.pi/2)*repetitions
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
                angle = (np.pi/2)*repetitions
                ideal_angle = (np.pi/2)*repetitions
                rho = (_qt.to_super(_qt.ry(angle))) * rho
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
                #print("starting {} with z-rotation 2pi*{:.2f}".format(gate_list_to_string(experiment), detuning*(time_per_gate*repetitions + gate_switching_time)))
                #make the transition frequency oscillate between a fraction of its nominal value
                modified_transition_f = transition_f*(1+detuning_noise)
                detuning = field_f_list[experiment_index] - modified_transition_f
                angle = 2*np.pi*detuning*(time_per_gate*repetitions + gate_switching_time)
                rho = (_qt.to_super(_qt.rz(angle))) * rho
            total_angle += angle
            total_ideal_angle += ideal_angle
           
        #calculate probabilities of being in 1 after the all the gates in the experiment have been applied
        p1 = (rho.dag()*rho1).norm()
        #fix the p1 if it exceeds 1 due to rounding error
        if p1 > 1:
            p1 = 2 - p1
        #get nCounts of data (bits) for the experiment
        one_counts = np.random.binomial(nCounts,p1)      
        angles.append(total_angle)
        ideal_angles.append(total_ideal_angle)
        probs.append(p1)
        ones_counts.append(one_counts)
        time_per_experiment_set.append(absolute_time)
        transition_f_list.append(modified_transition_f)
        detuning_list.append(detuning)

    return np.asarray(ones_counts), np.asarray(time_per_experiment_set), np.asarray(probs), np.asarray(varied_param)

if __name__=='__main__':

    gate_string = "(Gx)^49"                                
    #print("Start with string {}".format(gate_string))
    gate_list = gate_string_to_list(gate_string)
    #print("Input list of length {}".format(len(gate_list)))
    #print("Compressed form: {}".format(compress_gate_list(gate_list)))
    #print("Converted back to string: {}".format(gate_list_to_string(gate_list)))
    nSamples = 8000  #total samples to take for each measurement
    nCounts = 1      #total shots to take at one; =nSamples: same noise, probabilities for all repeats; =1, new experiment & noise for each count
    time_per_count = 1/60 #seconds
    time_units = 1e-3 #seconds
    noise_type='Sine' #Sine, Random Walk, Telegraph
    plot_noise=True
    walking_amp = 0.001
    telegraph_amp = 0.02
    res = 4
    low_frequency_noise = [np.round(0.0015*i, 4) for i in range(5)]
    low_frequency_noise = [] #uncomment this to eliminate low_frequency noise
    low_frequency_amps = [0.005*i for i in range(len(low_frequency_noise))]
    low_frequency_phase = [0]*len(low_frequency_noise)
    freq_list=[1.2]#1.2, 6, 8.4, 9.6] + low_frequency_noise
    amp_list=[0.002]#.002, 0.002, 0.0015, 0.0015] + low_frequency_amps
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
    #ones, zeros, times, probs, expected_angles, angles, sig = create_data(time_per_count, nSamples, nCounts, gate_list, time_units, noise_type, walking_amp, telegraph_amp, \
                #res, freq_list, amp_list, phase_list, start_f, stop_f, fluctuators,plot_noise,add_noise,noise_object=None,dc_angle_offset=dc_angle_offset, constant_linear_drift=constant_linear_drift)
    
    
    #print("times have {} points".format(len(times)))
    #print("ones list has {} points".format(len(ones)))
    #print("there were {} input points".format(nSamples))

    p_amplitude = (max(probs) - min(probs))/2
    print(p_amplitude)

    
    left_gate_list = ['Gx']
    right_gate_list = left_gate_list
    L = [i for i in range(1,300,3)]
    transition_f = 700e12
    field_f = transition_f + 5e4
    nCounts = 50
    time_per_gate = 6e-6
    gate_switching_time = 3e-6
    experiment_sample_time = 1/60 #s
    ###parameters for drifting noise object
    noise = 'Sine'
    plotted_noise = True
    freqs = [1.2]
    amps = [1e-10] #percent by which you want to oscillate the transition frequency
    phases = [0]
        
    rams_ones, rams_times, rams_probs, varied_params = ramsey_experiment(left_gate_list, right_gate_list, L, field_f, transition_f, \
                      nCounts, time_per_gate, gate_switching_time, experiment_sample_time, time_units=1e-3, noise_type=noise, \
                      freq_list=freqs,amp_list=amps, phase_list=phases, plot_noise=plotted_noise)
    
    plt.figure(figsize=(10,4))
    plt.plot(varied_params, rams_ones/nCounts, marker='.', label="Averaged over {}-Counts".format(nCounts))
    plt.xlabel("Varied Parameter (either delay time in s, or detuning in Hz)")
    plt.ylabel("1-state Probability")
    plt.title("Ramsey Fringes WITH Oscillatory Detuning")
    plt.legend(loc="lower right")
    plt.grid()
    plt.show()
    
    
    rams_ones, rams_times, rams_probs, varied_params = ramsey_experiment(left_gate_list, right_gate_list, L, field_f, transition_f, \
                      nCounts, time_per_gate, gate_switching_time, experiment_sample_time, time_units=1e-3, noise_type=None, \
                      freq_list=freqs,amp_list=amps, phase_list=phases, plot_noise=plotted_noise)
    
    plt.figure(figsize=(10,4))
    plt.plot(varied_params, rams_ones/nCounts, marker='.', label="Averaged over {}-Counts".format(nCounts))
    plt.xlabel("Varied Parameter (either delay time in s, or detuning in Hz)")
    plt.ylabel("1-state Probability")
    plt.title("Ramsey Fringes without Oscillatory Detuning")
    plt.legend(loc="lower right")
    plt.grid()
    plt.show()
