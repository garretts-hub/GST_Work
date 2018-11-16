# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 09:20:17 2018

@author: GA28573
"""

import numpy as np

def drift_input(file_loc):
    '''
    Returns a 4-D tuple with an array of the raw bit strings,
    an array of the total 1 counts per timestamp, an array of the
    total 0 counts per timestamp, and an array of the timestamp for
    each total count/bitstring
    '''
    raw_bit_list = []
    ones_counts_per_stamp = []
    zeros_counts_per_stamp = []
    timestamps = []
    with open(file_loc, 'r') as input:
        for row in input:
            row = row.split("\t")
            raw_bits = row[0].strip()
            time = row[1][11:].strip()
            hour = int(time[:2])
            minutes = int(time[2:4])
            seconds = int(time[5:7])
            msecs = int(time[9:])
            time_in_seconds = hour*60*60 + minutes*60 + seconds + msecs/1000
            #print("{} {} {} {}".format(hour, minutes, seconds, msecs))
            #print("Time in seconds {}".format(time_in_seconds))
            raw_bit_list.append(raw_bits)
            timestamps.append(time_in_seconds)
            ones = raw_bits.count('1')
            zeros = raw_bits.count('0')
            ones_counts_per_stamp.append(ones)
            zeros_counts_per_stamp.append(zeros)
    return (np.asarray(raw_bit_list), np.asarray(ones_counts_per_stamp), np.asarray(zeros_counts_per_stamp), np.array(timestamps))  

def merge_lines(file_loc, timestep, num_rows=None):
    '''
    Assumes equal spacing between each point, regardless of how many counts per line, how many lines, or the timestamps
    '''
    ones_counts_per_stamp = []
    zeros_counts_per_stamp = []
    timestamps = []
    with open(file_loc, 'r') as input:
        row_index = 1
        for row in input:
            row = row.split('\t')
            row = row[0].strip()
            for char in row:
                char = int(char)
                if char == 0:
                    zeros_counts_per_stamp.append(1)
                    ones_counts_per_stamp.append(0)
                elif char == 1:
                    zeros_counts_per_stamp.append(0)
                    ones_counts_per_stamp.append(1)
            row_index += 1
            if num_rows != None:
                if row_index > num_rows:
                    break
    num_samples = len(ones_counts_per_stamp)
    timestamps = [i*timestep for i in range(1, num_samples)]
    return (np.asarray(ones_counts_per_stamp), np.asarray(zeros_counts_per_stamp), np.array(timestamps))  

def sixteen_bit_lines_backwards(file_loc, timestep, num_experiments):
    '''
    Assumes 16 bits per line. Merges all lines into one bit string of 
    equidistant points. If more than one experiment is merged into the line, the long bitstring
    will be split up into num_experiments of equal length.
    '''
    datasets = []
    bitstring = []
    with open(file_loc, 'r') as file:
        for row in file:
            row = [int(bit) for bit in row[0:16]]
            #print("Forward row",row)
            row = row[::-1]#due to a LabView error, the bits are read in backwards
            #print("Reversed row", row_list)
            #print("Row list is",row_list)
            bitstring.extend(row)
    
    total_points = len(bitstring)
    points_per_experiment = total_points/num_experiments
    
    #print("total points is",total_points)
    #print("points per experiment is ",points_per_experiment)
    if points_per_experiment != round(points_per_experiment, 0):
        raise ValueError("Points per experiment is a decimal: Either the input number of experiments or the data file is incorrect.")
    else:
        points_per_experiment = int(points_per_experiment)
    for index in range(num_experiments):
        data = np.asarray(bitstring[index*points_per_experiment:(index + 1)*points_per_experiment])
        timestamp_array = np.asarray([i*timestep for i in range(1, points_per_experiment + 1)])
        set_tuple = (data, 1 - data, timestamp_array)
        datasets.append(set_tuple)
    return datasets #returns a num_experiment-length list with tuples containing (ones_count_array, zeros_count_array, timestamp_array)
    
    

def experiment_per_line(file_loc, timestep):
    '''
    Reads in a file with N-lines of C-bits each; 
    treats each line as an independent experiment;
    assumes all C-bits in all N-lines have equal spacing (timestep)
    '''
    datasets = [] #each element will be a tuple consisting of three lists: ones counts, zeros counts, and timesteps
    with open(file_loc, 'r') as file:
        for row in file:
            zeros_counts_per_stamp = []
            ones_counts_per_stamp = []
            row = row.split('\t')
            row = row[0].strip()
            for char in row:
                char = int(char)
                if char == 0:
                    zeros_counts_per_stamp.append(1)
                    ones_counts_per_stamp.append(0)
                elif char == 1:
                    zeros_counts_per_stamp.append(0)
                    ones_counts_per_stamp.append(1)
            num_samples = len(ones_counts_per_stamp)
            timestamps = [i*timestep for i in range(1, num_samples)]
            row_tup = (np.asarray(ones_counts_per_stamp), np.asarray(zeros_counts_per_stamp), np.asarray(timestamps))
            datasets.append(row_tup)
    
    return datasets #returns a list of tuples, with each tuple containing the ones array, zeros array, and timestamps for an experiment

'''def experiment_per_multiple_lines(file_loc, timestep, rows_per_experiment):
    #need documentation here
    #find how many rows are in the file
    with open(file_loc, 'r') as file:
        for i, l in enumerate(file):
            pass
        rows_in_file = i
        
    #find how many experiments you'll be able to get, given the number of rows in the file and the requested rows per experiment
    total_experiments = rows_in_file % rows_per_experiment
    
    datasets = [] #each element is a tuple containing three lists: ones counts, zeros counts, and timesteps
    with open(file_loc, 'r') as file:
        for exper in range(total_experiments):'''
            
          
    
        
    

def calculate_average_timestep(timestamp_array):
    differences = []
    for index in range(len(timestamp_array) - 1):         
        diff = float(timestamp_array[index + 1]) - float(timestamp_array[index])
        differences.append(diff)
    average = sum(differences)/len(differences)
    return average
            
if __name__=='__main__':
    raw, ones, zeros, times = drift_input("N:/Programs/Ions Share/Gate-Set Tomography/DriftAnalysis/PrelimData/2018_08_14/2018_08_14_1030_13_DRIFT_MOD.txt")
    #ones, zeros, times = manual_read("N:/Programs/Ions Share/Gate-Set Tomography/DriftAnalysis/PrelimData/2018_08_10_1322_27_DRIFT.txt", 1/60)

    #for i in range(len(zeros)):
        #print("Time {:.3f}, zero: {}, one: {}".format(times[i], zeros[i], ones[i]))
    #print(calculate_average_timestep(times))
    file_loc = 'N:/Programs/Ions Share/Gate-Set Tomography/DriftAnalysis/2018_11_16_64000/2018_11_16_1134_47_DRIFT.txt'
    data_list = sixteen_bit_lines_backwards(file_loc, 1/60, 3)
    print(len(data_list[0][0]))
    print(len(data_list))
