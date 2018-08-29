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



def gate_string_to_list(gate_string):
    #returns gate sequence as a list
    #ex, 'GxGxGx' --> ['Gx', 'Gx', 'Gx', 'Gx']
    #"Gx(Gy)^2Gz" --> ['Gx', 'Gy', 'Gy', 'Gz']
    #do not raise anything to the 1 exponent, only 2 or more; always use parenthesis
    #currently can only handle exponents up to 99
    #do not put more than one gate in parentheses for right now
    gate_list = []
    for i in range(len(gate_string)):
        char = gate_string[i]
        if char == "G":
            gate_list.append('G'+ gate_string[i+1])
        elif char.isnumeric():
            if i != (len(gate_string)-1): 
                next_char = gate_string[i+1]
            else:
                next_char = ''
            prev_char = gate_string[i-1]
            if prev_char.isnumeric():
                pass
            if next_char.isnumeric():
                char = int(gate_string[i] + gate_string[i+1])
            for count in range(int(char)-1):
                gate_list.append('G' + gate_string[i-3])
    return gate_list

def gate_list_to_string(gate_list):
    gate_string = ''
    for gate in gate_list:
        gate_string = gate_string + gate
    return gate_string

def compress_gate_list(gate_list):
    #this will compress the gate_list down to a list of 2-element lists, where the first element is the gate
    # and the second element is the number of consecutive appearances. The for loop below runs faster
    # when it can read the gates this way.
    num_gates = len(gate_list)
    compressed_list = []
    for i in range(num_gates):
        if i==0:
            compressed_list.append([gate_list[i], 1])
        else:
            if gate_list[i]==gate_list[i-1]:
                compressed_list[-1][1] += 1
            else:
                compressed_list.append([gate_list[i], 1])
    return compressed_list

def SNR(frequencies, powers, signal_frequencies, signal_bandwidth):
    #frequencies & powers from the FT; assumes evenly-spaced frequencies
    #signal frequencies are the frequencies you're interested in
    #bandwidth is the amount of frequency plus/minus your signal frequencies that 
    #   you're willing to count as part of your signal power
    frequencies = np.asarray(frequencies)
    powers = np.asarray(powers)
    mean_power = np.mean(powers)
    signals = []
    signal_power_list = [] #includes one element for the total power for each signal frequency of interest
    for f in signal_frequencies:
        lower_bound = f - signal_bandwidth
        upper_bound = f + signal_bandwidth
        signals.append((lower_bound, upper_bound))
    for signal in signals:
        signal_power = 0
        for i in range(len(frequencies)):
            f = frequencies[i]
            p = powers[i]
            if f >= signal[1]:
                break
            if f >= signal[0] and f <= signal[1]:
                signal_power += p
        signal_power_list.append(signal_power)
    SNR = sum(signal_power_list)/mean_power
    return SNR
    
    

def create_data(time_per_count, num_samples, num_counts, gate_list, time_unit, noise_type=None, walking_amp=None, telegraph_amp=None, \
                res=None, freq_list=None, amp_list=None, phase_list=None, start_f=None, stop_f=None, fluctuators=None, plot_noise=False, add_noise=False, noise_object=None):
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
    
    compressed_gate_list = compress_gate_list(gate_list)
    
    for time in timestamps:
        noise_at_time = 0
        if noise_type != None:
            #print(time/time_unit)
            noise_at_time = sig[time/time_unit]
        rho = rho0
        
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
                angle = np.pi/2 + noise_at_time
                rho = (_qt.to_super(_qt.rx(angle)))**gate_repetitions * rho
            elif gate_name == 'Gy':
                angle = np.pi/2 + noise_at_time
                rho = (_qt.to_super(_qt.ry(angle)))**gate_repetitions * rho
            elif gate_name == 'Gz':
                angle = np.pi/2 + noise_at_time
                rho = (_qt.to_super(_qt.rz(angle)))**gate_repetitions * rho
        #calculate probabilities of being in 1 after the experiment has been applied
        p1 = (rho.dag()*rho1).norm()
        probs.append(p1)
        one_count = np.random.binomial(num_counts, p1) #simulates a summation of the number of 1-counts you get in one bitstring sample
        zero_count = num_counts - one_count #simulates summation of the number of 0-counts in one bitstring sample
        one_counts.append(one_count)
        zero_counts.append(zero_count)
    
    if plot_noise == True:
        plt.plot(timestamps, probs)
        plt.ylim(0,1)
        plt.xlabel("Time, seconds")
        plt.ylabel("Probability of Measuring State {1}")
        plt.title("Simulated {} with {} Noise".format(gate_list_to_string(gate_list), noise_type))
        plt.grid()
        plt.show()
        
    return (np.asarray(one_counts), np.asarray(zero_counts), np.asarray(timestamps), probs, sig)



if __name__=='__main__':

    gate_string = "(Gx)^5"
    gate_list = gate_string_to_list(gate_string)
    #print(compress_gate_list(gate_list))
    nSamples = 1000  #total samples to take for each measurement
    nCounts = 1      #total shots to take at one; =nSamples: same noise, probabilities for all repeats; =1, new experiment & noise for each count
    time_per_count = 0.016 #seconds
    time_units = 1e-3 #seconds
    amp = 0.05
    noise_type='Sine' #Sine, Random Walk, Telegraph
    plot_noise=True
    walking_amp = 0.001
    telegraph_amp = 0.02
    res = 1
    freq_list=(0.8,)
    amp_list=(0.04,)
    phase_list=(0,)
    add_noise=0.01
    start_f = 0.1
    stop_f = 2
    fluctuators= 40
    
    
    
    #create_1f_data(time_per_count, nSamples, nCounts, gate_list, amp, fluctuators, start_f, stop_f, time_units)
    #create_data(time_per_count, nSamples, nCounts, gate_list, time_units, noise_type, walking_amp, telegraph_amp, \
    #            res, freq_list, amp_list, phase_list, start_f, stop_f, fluctuators,plot_noise,add_noise)
    
    
    
    