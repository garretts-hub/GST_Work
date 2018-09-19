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
from pygsti.extras import drift
import pylab as plt

def gate_string_to_list(gate_string):
    #returns gate sequence as a list
    #ex, 'GxGxGx' --> ['Gx', 'Gx', 'Gx', 'Gx']
    #"Gx(Gy)^2Gz" --> ['Gx', 'Gy', 'Gy', 'Gz']
    #do not raise anything to the 1 exponent, only 2 or more; always use parenthesis
    #currently can only handle exponents up to 99
    #do not put more than one gate in parentheses for right now
    gate_list = []
                
    string_list = list(gate_string)
    list_len = len(string_list)
    for char_index in range(list_len):
        char = string_list[char_index]
        #print("Index {}: {}".format(char_index, char))
        if char == ")":
            num = ''
            start_num = False
            #do a loop to find the start and end of the exponent number
            for sub_index in range(char_index, list_len):
                sub_char = string_list[sub_index]
                if sub_char.isnumeric():
                    num = num + sub_char
                    start_num = True
                elif not(sub_char.isnumeric()) and start_num == True:
                    break
            #print("Found number " + num)
            last_gate = gate_list[-1]
            gate_list += [last_gate]*(int(num)-1)
        if char == "G":
            gate = char + string_list[char_index + 1]
            gate_list.append(gate)
    return gate_list

def gate_list_to_string(gate_list):
    gate_string = ''
    #check if gate_list is in compressed form
    if type(gate_list[0]) == list:
        for gate in gate_list:
            num = gate[1]
            name = gate[0]
            if num > 1:
                gate_string = gate_string + '(' + name + ')^' + str(num)
            if num == 1:
                gate_string = gate_string + name
    #treat the gate_list as a list with each element as a single gate in string form
    if type(gate_list[0]) == str:
        total_gates = len(gate_list)
        count = 1
        for gate_index in range(total_gates):
            gate = gate_list[gate_index]
            #print("Gate {}: {}".format(gate_index, gate))
            if gate_index != (total_gates-1) and gate_list[gate_index + 1] == gate:
                count += 1
            elif count > 1:
                gate_string = gate_string + "(" + gate + ")^" + str(count)
                count = 1
            else:
                gate_string = gate_string + gate
                count = 1
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

def create_sorted_tuples(x, y, x_limits=None):
    #pass in two length L lists, x and y
    #returns a length-L list of 2 element tuples, (x_i, y_i),
    # (or returns only a range of sorted tuples selected from the x_limits specified)
    #sorted from greatest to lowest y value
    low_index = 0
    high_index = len(x)
    if x_limits != None:
        if len(x_limits) != 2:
            print("Error: x_limits must be an iterable of length 2 (min and max frequency to scan through!)")
            return None
        for i in range(len(x)):
            if x[i] > x_limits[0] and low_index == 0:
                low_index = i
            if x[i] > x_limits[1]:
                high_index = i
                break
    
    grouped = []
    if len(x) != len(y):
        print("Error: both input lists must be of the same length!")
        return None
    for i in range(low_index, high_index):
        grouped.append((x[i], y[i]))
    sorted_groups = sorted(grouped, key=lambda tup: tup[1], reverse=True)
    return sorted_groups

def select_tuples(sorted_list_of_tuples, x_values_list):
    #picks out N-tuples from the N-length values list
    #picks out the tuples whose x-value is in x_values_list
    tuples_of_interest = []
    for x in x_values_list:
        difference_list = [] #a list that will be filled with abs(x - tup[0]) for all tuples
        for tup in sorted_list_of_tuples:
            difference_list.append(abs(x - tup[0]))
        index_of_closest_x = difference_list.index(min(difference_list))
        tuples_of_interest.append(sorted_list_of_tuples[index_of_closest_x])
    return tuples_of_interest

def RMS(powers_list):
    #calculates RMS average for a list of powers
    # this function is used in other functions
    N = len(powers_list)
    squares_list = [val**2 for val in powers_list]
    sum_of_squares  = sum(squares_list)
    RMS = np.sqrt(sum_of_squares/N)
    return RMS

def SNR(frequencies, powers, central_freq, freq_band):
    #frequencies and powers are lists.
    # central freq is the frequency you're interested in, and
    # freq_band is how many frequencies plus or minus the central freq
    # that you're willing to consider. This is important for frequency
    # step size in the fourier transform.
    # Divides the total power within your frequency band by RMS mean.
    signal_power = 0
    lower_f = central_freq - freq_band
    upper_f = central_freq + freq_band
    for index in range(len(frequencies)):
        f = frequencies[index]
        if f > upper_f:
            break
        if f > lower_f:
            signal_power += powers[index]
    rms = RMS(powers)
    return signal_power/rms

def find_closest_frequency(frequencies, freq_of_interest):
    #finds the frequency and its index from a list closest to that of interest
    frequencies = np.asarray(frequencies)
    differences = abs(frequencies - freq_of_interest)
    minimum_difference = min(differences)
    minimum_index = list(differences).index(minimum_difference)
    closest_freq = frequencies[minimum_index]
    return closest_freq, minimum_index
                

def find_max_power(frequencies, powers, low_bound, hi_bound):
    #Finds the highest power within band of frequencies
    low_index = 0
    high_index = 0
    for f_i in range(len(frequencies)):
        freq = frequencies[f_i]
        if freq > low_bound and low_index == 0:
            low_index = f_i
        if freq > hi_bound:
            high_index = f_i
            break
    #print(frequencies[low_index:high_index])
    max_power = max(powers[low_index:high_index])
    return max_power

def find_multi_max_power(frequencies, powers, low_bound, high_bound, num_powers):
    #sums together the num_powers-highest powers in the band of frequencies
    sorted_tuples = create_sorted_tuples(frequencies, powers, (low_bound, high_bound))
    #print(sorted_tuples)
    power_sum = 0
    for i in range(num_powers):
        power_sum += sorted_tuples[i][1]
        #print("Adding {} from x-value {}".format(sorted_tuples[i][1], sorted_tuples[i][0]))
    return power_sum

def find_band_power(frequencies, powers, low_bound, high_bound):
    #sums the powers within a band of frequencies
    #returns the sum of those powers, and the range of frequencies that it spans
    power = 0
    low_freq = 0
    high_freq = 0
    for fi in range(len(frequencies)):
        freq = frequencies[fi]
        if freq > low_bound and freq < high_bound:
            power += powers[fi]
            if low_freq == 0:
                #print("Setting low freq at {}".format(freq))
                low_freq = freq
        elif freq > high_bound:
            high_freq = frequencies[fi - 1]
            #print("Setting high freq at {}".format(frequencies[fi - 1]))
            break
    freq_range = high_freq - low_freq
    return (power, freq_range)

def above_and_below(central_freq, frequencies, powers):
    '''
    finds the powers for the frequency just below and just above central_freq.
    Returns: 
        -the summed powers
        -lower frequency
        -upper frequency
        -the index of that closest frequency
    '''
    closest_index = 0
    summed_power = 0
    for fi in range(len(frequencies)):
        freq = frequencies[fi]
        if freq > central_freq:
            upper = freq
            lower = frequencies[fi - 1]
            summed_power = abs(powers[fi]) + abs(powers[fi - 1])
            #for the case that the upper frequency is closer
            if upper-central_freq < central_freq-lower:
                closest_index = fi
            #for all others, i.e. the lower freq is closer
            else:
                closest_index = fi - 1
            break
    return summed_power, lower, upper, closest_index
            

def single_frequency_reconstruction(drifted, low_freq, high_freq, \
                                    print_info=False, plot_original=False, plot_results=False): 
    #requires a Drift object as the input
    '''
    Takes the frequency with the maximum power within a band of specified frequencies, 
    and does the IDCT of just that frequency and its power. Returns a plot showing the
    probability as a function of time.
    '''
    #copy relevant info from the drift object
    all_frequencies = list(drifted.frequencies)
    all_modes = drifted.pspepo_modes[0,0,1,:]
    all_powers = list(all_modes**2)
    nSamples = len(drifted.data[0,0,1,:])
    nCounts = drifted.number_of_counts
    
    if plot_original:
        plt.figure(figsize=(15,3))
        plt.plot(all_frequencies, np.asarray(all_powers)/nSamples)
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Normalized Power, a.u.")
        plt.title("ORIGINAL Normalized Power Spectrum")
        plt.show()
        
        plt.figure(figsize=(15,3))
        times = np.linspace(0, drifted.timestep*nSamples, nSamples)
        plt.plot(times, drifted.pspepo_reconstruction[0,0,1,:])
        plt.xlabel("Time, seconds")
        plt.ylabel("1-State Probability")
        plt.ylim(0,1)
        plt.grid()
        plt.xlim(0, times[-1])
        plt.title("ORIGINAL pyGSTi Probability Reconstruction\n(using the entire power spectrum)")
        plt.show()
        
        
    #find the max power within the frequency range and its associated index and frequency
    max_power = find_max_power(all_frequencies, all_powers, low_freq, high_freq)
    max_power_index = all_powers.index(max_power)
    max_frequency = all_frequencies[max_power_index]
    if print_info:
        print("Max normalized power is {:.2f} at a frequency of {:.4f} Hz.".format(max_power/nSamples, max_frequency))
    
    #all modes except for the maximum one in the specified range are set to zero
    modified_modes = [0]*len(all_modes)
    modified_modes[max_power_index] = all_modes[max_power_index]
    
    #do the reconstruction here
    the_null_hypothesis = np.mean(drifted.data[0,0,1,:])*np.ones(nSamples,float)/nCounts
    my_reconstruction = drift.IDCT(modified_modes, null_hypothesis=the_null_hypothesis, counts=nCounts)/nCounts
    my_reconstruction = drift.renormalizer(my_reconstruction, method='logistic')
    
    if print_info:
        plt.figure(figsize=(15,3))
        plt.plot(all_frequencies, (np.asarray(modified_modes)**2)/nSamples)
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Normalized Power, a.u.")
        plt.title("Normalized Power Spectrum going into the IDCT\n(Peak at {:.4f} Hz)".format(max_frequency))
        plt.show()
    
    if plot_results:
        plt.figure(figsize=(15,3))
        plt.title("Reconstructed Probability Plot using IDCT\n(Using {:.4f} Hz)".format(max_frequency))
        plt.ylabel("1-State Probability")
        plt.xlabel("Time, seconds")
        plt.ylim(0,1)
        times = np.linspace(0, drifted.timestep*nSamples, nSamples)
        plt.grid()
        plt.xlim(0, times[-1])
        plt.plot(times, my_reconstruction, label='Reconstructed Probability: {:.3f}sin(2$\pi$ft) + {:.3f}'.format(max(my_reconstruction) - the_null_hypothesis[0],\
                 the_null_hypothesis[0]))
        plt.plot(times, the_null_hypothesis, label='Null Hypothesis: {:.3f}'.format(the_null_hypothesis[0]))
        plt.legend(loc="lower right")
        plt.show()
        
    amplitude = max(my_reconstruction) - the_null_hypothesis[0]
        
    return my_reconstruction, amplitude






def multi_frequency_reconstruction(drifted, central_freq, tolerance_band,\
                                    print_info=False, plot_original=False, plot_results=False): 
        #requires a Drift object as the input
    '''
    Looks within the range of (tolerance band - central freq) to (central freq + tolerance band).
    Takes the two highest peaks in that band, and associates their power with the peak closest to the
      central_freq power.
    Sets the power of all other frequencies to zero, then does the IDCT of the modified power spectrum.
    Returns a plot showing the probability as a function of time and a calculated amplitude of the
    probability oscillation about the mean.
    '''
    #copy relevant info from the drift object
    all_frequencies = list(drifted.frequencies)
    all_modes = drifted.pspepo_modes[0,0,1,:]
    all_powers = list(all_modes**2)
    nSamples = len(drifted.data[0,0,1,:])
    nCounts = drifted.number_of_counts
    
    if plot_original:
        plt.figure(figsize=(15,3))
        plt.plot(all_frequencies, np.asarray(all_powers)/nSamples)
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Normalized Power, a.u.")
        plt.title("ORIGINAL Normalized Power Spectrum")
        plt.show()
        
        plt.figure(figsize=(15,3))
        times = np.linspace(0, drifted.timestep*nSamples, nSamples)
        plt.plot(times, drifted.pspepo_reconstruction[0,0,1,:])
        plt.xlabel("Time, seconds")
        plt.ylabel("1-State Probability")
        plt.ylim(0,1)
        plt.grid()
        plt.xlim(0, times[-1])
        plt.title("ORIGINAL pyGSTi Probability Reconstruction\n(using the entire power spectrum)")
        plt.show()
        
        
    #Find the two highest frequencies within the frequency band of interest, and associate their combined
    # power with the frequency closest to the central frequency
    info_tuples = []
    for fi in range(len(all_frequencies)):
        freq = all_frequencies[fi]
        if (central_freq - tolerance_band) < freq < (central_freq + tolerance_band):
            info_tuples.append((fi, freq, all_modes[fi], all_powers[fi]))
        elif freq > (central_freq + tolerance_band):
            break
    #sort them descending by their power, the 4th element in the tuple
    sorted_tuples = sorted(info_tuples, key=lambda tup: tup[3], reverse=True)
    summed_mode = abs(sorted_tuples[0][2]) + abs(sorted_tuples[1][2])
    summed_power = summed_mode**2
    input_modes = [0]*len(all_modes)
    closest_freq, closest_index = find_closest_frequency(all_frequencies, central_freq)
    input_modes[closest_index] = summed_mode
    if print_info:
        print("Using {:.4f} and {:.4f} Hz, with summed power {}".format(sorted_tuples[0][1], sorted_tuples[1][1], summed_power))
    
   
    #do the reconstruction here
    the_null_hypothesis = np.mean(drifted.data[0,0,1,:])*np.ones(nSamples,float)/nCounts
    my_reconstruction = drift.IDCT(input_modes, null_hypothesis=the_null_hypothesis, counts=nCounts)/nCounts
    my_reconstruction = drift.renormalizer(my_reconstruction, method='logistic')
    
    if print_info:
        plt.figure(figsize=(15,3))
        plt.plot(all_frequencies, (np.asarray(input_modes)**2)/nSamples)
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Normalized Power, a.u.")
        plt.title("Normalized Power Spectrum going into the IDCT\n(Merging power at {:.4f} and {:.4f} Hz)".format(sorted_tuples[0][1], sorted_tuples[1][1]))
        plt.show()
    
    if plot_results:
        plt.figure(figsize=(15,3))
        plt.title("Reconstructed Probability Plot using IDCT\n(Merging power at {:.4f} and {:.4f} Hz)".format(sorted_tuples[0][1], sorted_tuples[1][1]))
        plt.ylabel("1-State Probability")
        plt.xlabel("Time, seconds")
        plt.ylim(0,1)
        times = np.linspace(0, drifted.timestep*nSamples, nSamples)
        plt.grid()
        plt.xlim(0, times[-1])
        plt.plot(times, my_reconstruction, label='Reconstructed Probability: {:.3f}sin(2$\pi$ft) + {:.3f}'.format(max(my_reconstruction) - the_null_hypothesis[0],\
                 the_null_hypothesis[0]))
        plt.plot(times, the_null_hypothesis, label='Null Hypothesis: {:.3f}'.format(the_null_hypothesis[0]))
        plt.legend(loc="lower right")
        plt.show()
        
    amplitude = max(my_reconstruction) - the_null_hypothesis[0]
    
    return my_reconstruction, amplitude
        
    
if __name__=='__main__':
    f = np.arange(0, 20)
    p = -(f - 10)**2 + 100
    print(f)
    print(p)
    plt.plot(f, p)
    
    maxpow = find_multi_max_power(f, p, 6, 15, 2)
    print("Summed max power is {}".format(maxpow))
    
    
    
    
    
    
    
    
    
    
    