# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:51:04 2018

@author: GA28573
"""

import numpy as np
import pylab as plt
from drift_file_io import drift_input, merge_lines


file_loc = 'N:/Programs/Ions Share/Gate-Set Tomography/DriftAnalysis/2018_08_27 Gi Data/2018_08_27_1627_08_DRIFT.txt'
raw_bit_array, ones_count_array, zeros_count_array, timestamp_array = drift_input(file_loc)
nCounts = ones_count_array[0] + zeros_count_array[0]
nSamples = len(ones_count_array)

plt.figure(figsize=(10,4))
plt.title("Experiment Raw Data")
plt.ylabel("1-State Probability,\naveraged over the entire row")
plt.xlabel("Absolute time, s (Timestamp of Row)")
plt.plot(timestamp_array - timestamp_array[0], ones_count_array/nCounts)
plt.show()