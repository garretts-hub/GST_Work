# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 12:34:44 2018

@author: Garrett
"""

from scipy.fftpack import dct, idct, fft, ifft, rfft, irfft
import numpy as np
import pylab as plt

def sine(x, f, a):
    return 0.5 + a*np.sin(2*np.pi*f*x)

time = np.linspace(0, 10, 2000)
spacing = time[1] - time[0]
vals = sine(time, 5, 0.30)
normalized_vals = (vals - np.average(vals))
plt.plot(time, vals)
plt.xlabel("Time, s")
plt.ylabel("Arbitary Units")
plt.title("Sine Test Function")
plt.grid()
plt.ylim(0,1)
plt.show()

dct_modes = dct(normalized_vals, norm='ortho')
frequencies = [t/(2*len(time)*spacing) for t in range(len(time))]

plt.plot(dct_modes.real, label="Real")
plt.plot(dct_modes.imag, label='Imaginary')
plt.title("DCT Modes with Normalized Input")
plt.xlabel("Mode Index")
plt.ylabel("Mode value")
plt.legend()
plt.grid()
plt.show()


plt.plot(frequencies, np.absolute(dct_modes)**2/len(time))
plt.title("DCT Modes Squared, Divided by num points")
plt.xlabel("Mode Index")
plt.ylabel("Mode value")
plt.grid()
plt.show()

fft_modes = fft(normalized_vals)
rfft_modes = rfft(normalized_vals)
plt.plot(rfft_modes)
plt.grid()
plt.title("Using Real FFT with Normalized Input")
plt.xlabel("Mode Index")
plt.ylabel("Mode Value")
plt.show()

plt.plot(rfft_modes**2/len(time))
plt.title("rFFT modes squared divided by num samples")
plt.xlabel("Mode Index")
plt.ylabel("Mode Value")
plt.grid()
plt.show()

plt.plot(fft_modes.real, marker='.', label="Real")
plt.plot(fft_modes.imag, marker='.', label="Imaginary")
plt.title("FFT Modes with Normalized Input")
plt.xlabel("Mode Index")
plt.ylabel("Mode Value")
plt.legend()
#plt.xlim(4,6)
plt.grid()
plt.show()

plt.plot(np.absolute(fft_modes)**2/len(time))
plt.title("FFT Modes, Norm Squared, Divided by Num Samples")
plt.xlabel("Mode Index")
plt.ylabel("Mode Value")
plt.grid()
plt.show()


print("Sum of DCT Modes squared and divided by num points is {}".format(np.sum(np.absolute(dct_modes)**2)/len(time)))
print("Sum of rFFT modes squared divided by num points is {}".format(np.sum(rfft_modes**2)/len(time)))
print("Sum of FFT Modes squared divided by num points is {}".format(np.sum(np.absolute(fft_modes)**2)/len(time)))
print("Sum of input values squared is {}".format(np.sum(normalized_vals**2)))

amplitudes = np.arange(0.01, 0.3, 0.01)
peak_powers_fft = [] #modes squared, then divided by numpoints
peak_powers_dct = []
f = 5 #hz
num_points = 2000
total_time = 10 #seconds
for amp in amplitudes:
    vals = sine(np.linspace(0, total_time, num_points), f, amp)
    shifted_vals = vals - np.average(vals)
    fft_modes = fft(shifted_vals)
    peak_power_fft = np.max(np.absolute(fft_modes)**2/num_points)
    peak_powers_fft.append(peak_power_fft)
    dct_modes = dct(shifted_vals, norm='ortho')
    peak_power_dct = np.max(np.absolute(dct_modes)**2/num_points)
    peak_powers_dct.append(peak_power_dct)
    
plt.plot(amplitudes, peak_powers_fft)
plt.xlabel("Signal Amplitude, a.u.")
plt.ylabel("Peak Power divided by num points")
plt.title("FFT Scaling of Power to Sine Amplitude")
plt.grid()
plt.show()

plt.plot(amplitudes, peak_powers_dct)
plt.xlabel("Signal Amplitude, a.u.")
plt.ylabel("Peak Power divided by num points")
plt.title("DCT Scaling of Power to Sine Amplitude")
plt.grid()
plt.show()




