# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 15:25:23 2018

@author: GA28573
"""

import numpy as np
from scipy.fftpack import dct, idct
import pylab as plt
import math


class Drift(object):
    def __init__(self, timestamps_array, data_array, nCounts=None):
        '''Input an array of of the num of 1 counts per timestamp'''
        self.data = data_array
        self.times = timestamps_array
        self.samples = len(data_array)
        self.counts = nCounts
        
        summed_timestep = 0
        for i in range(1, self.samples):
            summed_timestep += self.times[i] - self.times[i-1]
        self.avg_timestep = summed_timestep/self.samples
        
        #check if input is valid
        for point in self.data:
            if nCounts != None:
                if point > self.counts:
                    raise ValueError("Cannot have a sample with greater than nCounts!")
                elif point < 0:
                    raise ValueError("Cannot have a sample with negative counts!")
                    
        if len(self.data) != len(self.times):
            raise ValueError("Must have the same number of data points and timestamps!")
        
        #initialize variables
        self.timesteps = [0]*(len(self.times)-1) #will list the time spacing between all 
        for i in range(0, len(self.times)-1 ):
            self.timesteps[i] = self.times[i+1] - self.times[i]
        self.avg_timestep = np.mean(self.timesteps)
        self.frequencies = np.zeros(self.samples)
        self.modes = np.zeros(self.samples)
        self.powers = np.zeros(self.samples)
        
    def __repr__(self):
        return "Drift(Len={}, start_time={}s, end_time={}s, counts_per_sample={}, average_timestep={:.3}s".format(self.samples, self.times[0], self.times[-1], self.counts, self.avg_timestep)
    
    def _dct(self, print_details=False, null_hypothesis=None):
        """    
        x : array; Data string, on which the normalization and discrete cosine transformation is performed. If
            counts is not specified, this must be a bit string.
        null_hypothesis : array, optional
            If not None, an array to use in the normalization before the DCT. If None, it is
            taken to be an array in which every element is the mean of x.
        """
        delta_freq = 1/(2*self.avg_timestep*self.samples)
        self.frequencies = np.arange(0, self.samples*delta_freq, delta_freq)
        if print_details:
            print("Frequency spacing is {:.4f} Hz".format(delta_freq))
        x = self.data
        nCounts = self.counts
        x_mean = np.mean(x)
        N = len(x) 
        # If the null hypothesis is not specified, we take our null hypothesis to be a constant bias
        # coin, with the bias given by the mean of the data / number of counts.
        if null_hypothesis is None:    
            null_hypothesis = x_mean/nCounts
            if null_hypothesis <= 0 or null_hypothesis >= 1:
                return np.zeros(N)
        input_array= (x - nCounts*null_hypothesis)/np.sqrt(nCounts*null_hypothesis * (1 - null_hypothesis))
        self.modes = dct(input_array, norm='ortho')
        self.powers = self.modes**2
    
    def _manual_ndft(self, print_details=False):
        times = self.times
        vals = self.data
        nCounts = self.counts
        N = self.samples
        null_hypothesis = np.mean(vals)/nCounts
        input_array= (vals - nCounts*null_hypothesis)/np.sqrt(nCounts*null_hypothesis * (1 - null_hypothesis))
        T = times[-1]
        if print_details: print("Will return {}-frequencies (N/2), spaced at {} Hz".format(math.ceil(N/2), (1/T)))
        frequencies = np.arange(N)*(1/T)
        modes = np.zeros(N)
        for m in range(len(frequencies)):
            summed = 0
            for n in range(N):
                summed += input_array[n]*np.exp(-2j*np.pi*m/T*times[n])
            modes[m] = summed.real
        #since the ndft covers -N/2 to N/2 (not 0 -> N) integer multiples of 1/T for frequency,
        #this is because of the Nyquist theorem. Although the timestep may not be
        # perfectly defined, and therefore the sample rate varies, calculating an
        #average sample rate will only give you half the frequencies.
        #The frequencies should be symmetric about the halfway frequency, so add their
        #powers together.
        halfN = math.ceil(N/2)
        new_frequencies = []
        new_powers = []
        for n in range(halfN):
            new_frequencies.append(frequencies[n])
            power = (modes[n])**2 + (modes[len(modes) - 1 - n])**2 #sum the first and last powers, then the second & second to last..etc.
            #power = modes[n]**2 #this is what it would be if you didn't add the reflection, and just truncated at N/2
            new_powers.append(power)
        self.frequencies = new_frequencies
        self.powers = new_powers
        if print_details: 
            print("Average timestep is {:.4f} s".format(self.avg_timestep))
            print("Average sample rate is {:.3f} Hz".format(1/self.avg_timestep))
            print("Max frequency (under Nyquist limit) is {:.3f} Hz".format(self.frequencies[-1]))
        return self.frequencies, self.powers
    
    def plot_input(self):
        plt.plot(self.times, self.data,marker='.')
        plt.grid()
        plt.xlabel("Time, s")
        plt.ylabel("Measured Value, a.u.")
        plt.title("Data Set to be Transformed")
        plt.show()
    
        
    def plot_power_spectrum(self):
        plt.plot(self.frequencies, self.powers, marker='.')
        plt.grid()
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Power, a.u.")
        plt.title("Power Spectrum")
        plt.show()
    
    
    