# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 10:16:22 2018

@author: GA28573
"""
import pygsti
from pygsti.objects import DataSet
from pygsti.objects import GateString
from create_data_garrett import generate_model_data, generate_model_data_rotation_error
import matplotlib.pyplot as plt

listOfExperiments = [GateString(('Gx',)), GateString(('Gy',))]#, GateString(('Gi',))]
print(listOfExperiments)

samples = 9
shots = 1
repeats = 3#shots
maxLengths = [1]#,2,4,8,16,32]
time_per_sequence = 5e-4 #nanoseconds
time_to_reset = 0.1 #nanoseconds, time to reprogram a new sequence
error_amplitude = 0
seeds = [1234]
xerr = (0.004,7e8)
yerr = (0.0026,1e9)
zerr=(0,0)
output = generate_model_data_rotation_error(False, listOfExperiments, samples, shots, repeats, max(maxLengths), 'phase_damping', 'fixed', error_amplitude, time_per_sequence, time_to_reset, seeds)#, xerr, yerr, zerr)
ds = output[0] #this is the conventional DataSet object
timestamped = output[1] #this is the dictionary of timestamped counts
print("\n")
print(timestamped)
print("\n")

tdds = pygsti.objects.DataSet(outcomeLabels=['0','1'])

for gate in timestamped.keys():
    sequence_list = []
    count_list = [] #list of count dictionaries of the form {'0':10, '1':90} for a tdds object
    time_list = [] #list of timestamps; one per dictionary
    item = timestamped[gate]
    gate_name = item[0]
    for i in range(1, len(item)):
        
        zero_counts = item[i][0]
        one_counts = item[i][1]
        timestamp = item[i][2]

        count_list.append({'0':zero_counts, '1':one_counts})
        time_list.append(timestamp)  
    tdds.add_series_data((gate_name,), count_list, time_list)
tdds.done_adding_data()
print(tdds)

