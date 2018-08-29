# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 14:45:18 2018

@author: GA28573
"""
import sys
sys.path.append("C:/Users/GA28573/AppData/Local/Continuum/anaconda32/Lib/site-packages/qutip-4.3.1-py3.6-win-amd64.egg/qutip")
import pygsti
from pygsti.construction import std1Q_XYI
from create_data import generate_model_data, generate_model_data_rotation_error
import matplotlib.pyplot as plt
import NoiseSignal2 as _ns

maxLengths = [1,2,4]
gs_target = std1Q_XYI.gs_target
prep_fiducials = std1Q_XYI.prepStrs[0:3]
meas_fiducials = std1Q_XYI.effectStrs[0:3]
germ_strs = ['Gi', 'Gx', 'Gy']
germs = pygsti.construction.gatestring_list(germ_strs)
listOfExperiments = pygsti.construction.make_lsgst_experiment_list(gs_target, prep_fiducials, meas_fiducials, germs, maxLengths)
gx3 = pygsti.objects.GateString(('Gx','Gx','Gx'))
gxgi8gx = pygsti.objects.GateString(None, "Gx(Gi)^8Gx")
listOfExperiments = [gx3]
print(listOfExperiments)

sequence = pygsti.objects.GateString(('Gx','Gx','Gx'))
nSamples = 50  #total samples to take for each measurement
shots = 1      #total shots to take at one; =nSamples: same noise, probabilities for all repeats; =1, new experiment & noise for each count
nRepeats = nSamples   #num times to repeat a sequence before moving to the next
L = 32 #max length of germ sequences
error_type = 'amplitude'
amplitude_type = 'non_markovian_1f'#'fixed' #or 'non_markovian_1f'
error_amp = 0.10 #amplitude factor for generated noise
start_f = 0.1
stop_f = 100
fluctuators = 40
time_units = 1e-3 #ms
time_per_exp = 0.016 #time in seconds (will be converted to ns in the function) for one prep-sequence-meas
time_reprogram = 0.1 #time in seconds (will be converted to ns) to reprogram a new sequence
seeds=[1234]

dataset, ones, zeros, timestamps, probs = \
    generate_model_data_rotation_error([sequence], nSamples, shots, nRepeats, L, error_type, amplitude_type, error_amp, \
                                       start_f, stop_f, time_units, fluctuators, time_per_exp*1e9, time_reprogram*1e9, seeds, do_recalibration=True)

