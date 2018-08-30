# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 20:24:35 2018

@author: Garrett
"""

import numpy as np
import pylab as plt
import nfft

def make_time(start, stop, num, sigma):
    t = np.linspace(start, stop, num)
    for i in range(len(t)):
        t[i] += np.random.normal(0,sigma)
        while t[i] < start or t[i] > stop:
            t[i] += np.random.normal(0,sigma)
    return t

def get_powers(times, vals):
    N = len(times)
    T = times[-1]
    print("Will return {}-frequencies, spaced at {} Hz".format(N, (1/T)))
    freqs = np.arange(N)*(1/T)
    powers = np.zeros(N)
    for m in range(len(powers)):
        summed = 0
        for n in range(N):
            summed += vals[n]*np.exp(-2j*np.pi*m/T*times[n])
        powers[m] = summed
    return freqs, (powers)**2


N = 50
freq = 4.5
start_t = 0
stop_t = 3
num_points = 40
sigma = (stop_t/num_points)*0.5
x = make_time(start_t, stop_t, num_points, sigma)
x0 = np.linspace(start_t, stop_t,800)
f = np.sin(2*np.pi*freq*x)
f0 = np.sin(2*np.pi*freq*x0)

plt.plot(x,f, label="Irregular",marker='.',markersize=8,ls="None")
plt.plot(x0,f0, label="Regular")
plt.grid()
plt.legend(loc='lower right')
plt.xlabel("Time, seconds")
plt.title("Data Points")
plt.show()

f_try, p_try = get_powers(x, f)
plt.plot(f_try, p_try, marker='.')
plt.grid()
plt.xlabel("Time, seconds")
plt.title("Attempted Manual ndft")
plt.show()

k = np.arange(-N//2, N//2)
powers = np.absolute(nfft.nfft_adjoint(x,f,N))
plt.plot(k,powers.real,label="Real",marker='.')
#plt.plot(k, powers.imag, label="Imaginary")
plt.title("Non-Equispaced Fast Fourier Transform")
plt.legend()
plt.grid()
plt.xlim(0, k[-1])
plt.show()