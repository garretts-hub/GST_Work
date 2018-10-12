# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 11:25:52 2018

@author: GA28573
"""

import numpy as np
import pylab as plt
from scipy import signal

def p1(t, nominal):
    #probability of being in the 1-state at time t
    return (0.5 + nominal[1]*np.sin(2*np.pi*nominal[0]*t + nominal[2]))

def p0(t, nominal):
    #probability of being on 0-state at time t
    return 1 - p1(t, nominal)

def p1_sine(times, nominal):
    return 0.5 + nominal[1]*np.sin(2*np.pi*nominal[0]*times + nominal[2])
def p0_sine(times, nominal):
    return 0.5 - nominal[1]*np.sin(2*np.pi*nominal[0]*times + nominal[2])

def p1_square(times, nominal, tpp=None):
    #make the time_per_pulse a third of the period by default
    if tpp == None:
        tpp = 1/(3*nominal[0])
    return 0.5 + nominal[1]*signal.square(2*np.pi*nominal[0]*times + nominal[2], tpp)
def p0_square(times, nominal, tpp=None):
    #make the time_per_pulse a third of the period by default
    if tpp == None:
        tpp = 1/(3*nominal[0])
    return 0.5 - nominal[1]*signal.square(2*np.pi*nominal[0]*times + nominal[2], tpp)

def p1_saw(times, nominal):
    #to make this more of a triangle, increase tpp to 0.5
    tpp = 0
    return 0.5 + nominal[1]*signal.sawtooth(2*np.pi*nominal[0]*times + nominal[2], tpp)
def p0_saw(times, nominal):
    #to make this more of a triangle, increase tpp to 0.5
    tpp = 0
    return 0.5 - nominal[1]*signal.sawtooth(2*np.pi*nominal[0]*times + nominal[2], tpp)

def loss(time_data, y_data, nominal, form='sine', tpp=None):
    if form == 'sine':
        p1 = p1_sine(time_data, nominal)
        p0 = p0_sine(time_data, nominal)
    elif form == 'square':
        p1 = p1_square(time_data, nominal, tpp)
        p0 = p0_square(time_data, nominal, tpp)
    elif form == 'saw':
        p1 = p1_saw(time_data, nominal)
        p0 = p0_saw(time_data, nominal)
    sum_term = 0
    for i in range(len(time_data)):
        sum_term += np.log( (1-y_data[i])*p0[i] + y_data[i]*p1[i])
    return -1*sum_term

def bit_flip(input_bit, flip_probability):
    out_bit = np.random.choice([input_bit, not(input_bit)], 1, p=[1-flip_probability, flip_probability])
    return out_bit[0]

def variable_loss(time_data, y_data, variable_name, variable_array, nominal_variables, form='sine', tpp = None):
    losses=[]
    for variable in variable_array:
        f = nominal_variables[0]
        amp = nominal_variables[1]
        phase = nominal_variables[2]
        tpp = tpp
        
        if variable_name == 'frequency':
            f = variable
        elif variable_name == 'amplitude':
            amp = variable
        elif variable_name == 'phase':
            phase = variable
        elif variable_name == 'tpp':
            tpp = variable
        
        l = loss(time_data, y_data, (f, amp, phase), form=form, tpp=tpp)
        losses.append(l)
    
    variable_array = list(variable_array)
    for i in range(len(losses) - 1, -1, -1):
        if np.isnan(losses[i]):
            del losses[i]
            del variable_array[i]
    minimum_index = losses.index(min(losses))
    return variable_array, losses, minimum_index

def MLE(times, vals, nominal_params, f_range, a_range, p_range, form, tpp=None, input_f=None, input_a=None, input_p=None, plot_range=None):
    if input_f == None:
        variable_array = f_range
        variable = 'frequency'
        x, y, index = variable_loss(times, vals, variable, variable_array, nominal_params, form, tpp)
        plt.plot(x, y, marker='.')
        plt.grid()
        plt.xlabel(variable.capitalize())
        plt.title("Max Likelihood Estimation for Sinusoidal Binomial Probability\n--Frequency--")
        plt.show()
        input_f = x[index]
    
    if input_a == None:
        variable_array = a_range
        variable = 'amplitude'
        x, y, index = variable_loss(times, vals, variable, variable_array, nominal_params, form, tpp)
        plt.plot(x, y, marker='.')
        plt.grid()
        plt.xlabel(variable.capitalize())
        plt.title("Max Likelihood Estimation for Sinusoidal Binomial Probability\n--Amplitude--")
        plt.show()
        input_a = x[index]
    
    if input_p == None:
        variable_array = p_range
        variable = 'phase'
        x, y, index = variable_loss(times, vals, variable, variable_array, nominal_params, form, tpp)
        plt.plot(x, y, marker='.')
        plt.grid()
        plt.xlabel(variable.capitalize())
        plt.title("Max Likelihood Estimation for Sinusoidal Binomial Probability\n--Phase--")
        plt.show()
        input_p = x[index]
    
    optimal_params = (input_f, input_a, input_p)
    
    if form == 'sine':
        reconst = p1_sine(times, optimal_params)
    
    if form == "square":
        x, y, index = variable_loss(times, vals, 'tpp', np.linspace(0, 1/input_f, 50), optimal_params, form)
        plt.plot(x, y, marker='.')
        plt.grid()
        plt.title("Max Likelihood Estimation for Sinusoidal Binomial Probability\n--Time per Pulse--")
        plt.show()
        tpp = x[index]
        reconst = p1_square(times, optimal_params, tpp)
        
    if form == 'saw':
        reconst = p1_saw(times, optimal_params)
        
    plt.plot(times, vals, marker='.', ls='None', label="Data Points")
    plt.plot(times, reconst, label="Reconstruction")
    if plot_range != None:
        plt.xlim(plot_range[0], plot_range[1])
    else:
        plt.xlim(0, 3)
    plt.grid()
    plt.xlabel("Time, seconds")
    plt.legend(loc="lower right")
    plt.title("Reconstructed 1-State Probability")
    plt.show()
    
    print("Frequency: {:.3f} Hz\nAmplitude: {:.3f}\nPhase: {:.3} radians\n".format(input_f, input_a, input_p))
    if form == 'square': 
        print("Time per pulse: {:.4f} seconds".format(tpp))
    
    return times, reconst

def two_dimensional_optimization(times, vals, f, a_range, p_range, form, verbose=True):
    num_a = len(a_range)
    num_p = len(p_range)
    losses = np.ndarray(shape=(num_a, num_p))
    for a_i in range(num_a):
        a = a_range[a_i]
        for p_i in range(num_p):
            p = p_range[p_i]
            l = loss(times, vals, (f, a, p), form=form, tpp=None)
            losses[a_i, p_i] = l
            
    minimum = np.min(losses)
    min_loc = np.where(losses == minimum)
    min_loc = np.asarray(min_loc).T
    a_index = min_loc[0, 0]
    p_index = min_loc[0, 1]
    optimized_tuple = (f, a_range[a_index], p_range[p_index])
    
    if form == 'sine':
        prob = p1_sine(times, (f, a_range[a_index], p_range[p_index]))
    if form == 'saw':
        prob = p1_saw(times, (f, a_range[a_index], p_range[p_index]))
    if form == 'square':
        prob = p1_square(times, (f, a_range[a_index], p_range[p_index]), tpp=None)
    
    return losses, prob, optimized_tuple
                

def three_dimensional_optimization(times, vals, f_range, a_range, p_range, form, verbose=True):
    num_f = len(f_range)
    num_a = len(a_range)
    num_p = len(p_range)
    losses = np.ndarray(shape=(num_f, num_a, num_p))
    for f_i in range(num_f):
        f = f_range[f_i]
        print("Scanning amp & phase with F = {:.3f} Hz".format(f))
        for a_i in range(num_a):
            a = a_range[a_i]
            for p_i in range(num_p):
                p = p_range[p_i]
                l = loss(times, vals, (f, a, p), form=form, tpp=None)
                losses[f_i, a_i, p_i] = l
            
    minimum = np.min(losses)
    min_loc = np.where(losses == minimum)
    min_loc = np.asarray(min_loc).T
    f_index = min_loc[0, 0]
    a_index = min_loc[0, 1]
    p_index = min_loc[0, 2]
    optimized_tuple = (f_range[f_index], a_range[a_index], p_range[p_index])
    
    if form == 'sine':
        prob = p1_sine(times, (f_range[f_index], a_range[a_index], p_range[p_index]))
    if form == 'saw':
        prob = p1_saw(times, (f_range[f_index], a_range[a_index], p_range[p_index]))
    if form == 'square':
        prob = p1_square(times, (f_range[f_index], a_range[a_index], p_range[p_index]), tpp=None)
    
    return losses, prob, optimized_tuple
                

import tensorflow as tf
def tensorflow_optimization(times, vals, init_f, init_a, init_p, tpp=None, nepochs=2, neval_period=10, 
                     learning_rate=0.001, optimizer="gd",
                     mini_batch_size=512, verbose=True, do_plot=True):
  '''
  Tensorflow model of sinusoid evolution, with J as the latent variable.
  
  Runs maximum likelihood fitting, returns J.
  
  data = dict with data
  nepochs = number of "epochs" of fitting, where each epoch goes through all the data
  
  This version uses minibatching, to allow randomness to encourage jumps out of
  shallow local minima.
  '''
  print("Model training for %s epochs, with evaluation every %s steps" % (nepochs,neval_period))
  
  data = {'samples': vals, 'time': times}
  
  samples = np.array(data['samples']).astype(np.float32)  
  tpts = np.array(data['time']).astype(np.float32)  
  batch_size = len(tpts)
  
  xy = np.stack([tpts, samples], axis=1)
  print("Input data xy shape=%s" % str(xy.shape))
  n_mini_batches = int(xy.shape[0] / mini_batch_size)
  print("Each epoch has %d mini_batches of size %s" % (n_mini_batches, mini_batch_size))  
  print("Input data shapes samples=%s, tpts=%s" % (samples.shape, tpts.shape))
  
  # graph input
  X = tf.placeholder(tf.float32, name="X_time_points")
  Y = tf.placeholder(tf.float32, name="Y_samples")        # y = 0 or 1
  print("Input data placeholders X=%s, Y=%s" % (X, Y))
  
  # model variables
  F = tf.Variable(init_f, name="f", dtype=tf.float32)
  A = tf.Variable(init_a, name='a', dtype=tf.float32)
  P = tf.Variable(init_p, name='p', dtype=tf.float32)
  
  # model to predict sample probability (prob(k=0))
  Y_sp = A*tf.sin(2*np.pi*F*X + P) + 0.5
  
  # loss function (binary cross-entropy)
  Ybin = Y  # 0 or 1
  loss = tf.log( Ybin * Y_sp + (1-Ybin) * (1-Y_sp) + 1.0e-8)   # log likelihood
  loss = tf.reduce_mean(loss)  # take mean over batch
  
  # optimizer
  if optimizer=="gd":
    gdo = tf.train.GradientDescentOptimizer(learning_rate=learning_rate)
    optimizer_op = gdo.minimize(loss)
    op_name = "GD"
  elif optimizer=="adagrad":
    optimizer_op = tf.train.AdagradOptimizer(learning_rate).minimize(loss)
    op_name = "Adagrad"
  else:
    optimizer_op = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)
    op_name = "Adam"
  print("Using %s optimmizer, learning rate=%s" % (op_name, learning_rate))
  
  losses = []    # record losses at each epoch step
  steps = []
  
  # run MLE
  if verbose:
    print("Running MLE over %s datapoints with %s epochs" % (batch_size, nepochs))
  with tf.Session() as sess:
    sess.run(tf.global_variables_initializer())
    for k in range(nepochs):
      np.random.shuffle(xy)     # random in-place permutation of first dimension
      for n in range(n_mini_batches):
        n0 = n*mini_batch_size
        sess.run(optimizer_op, feed_dict={X:xy[n0:n0+mini_batch_size, 0], 
                                          Y:xy[n0:n0+mini_batch_size, 1]})      

      if not (k % neval_period):
        results = sess.run([loss, F, A, P], feed_dict={X:tpts, Y:samples})
        print("    Epoch %s: loss=%s, F=%s, A=%s, P=%s" % tuple([k] + results))
        losses.append(results[0])
        steps.append(k)
        if np.isnan(results[0]):
          raise Exception("loss is NaN, quitting!")
  
    results = sess.run([loss, F, A, P], feed_dict={X:tpts, Y:samples})
    
  m_loss, m_F, m_A, m_P = results
  if verbose:
    print("Results from ML regression: loss=%s, F=%s, A=%s, P=%s" % (m_loss, m_F, m_A, m_P))

  
  if do_plot:
    plt.plot(steps, losses, 'go')
    plt.plot(steps, losses)
    plt.grid(True)
    plt.xlabel("Epoch step number")
    plt.ylabel("Loss (negative log likelihood)")
    plt.suptitle("Tensorflow MLE on dataset with %s samples using %s optimizer" % 
              (samples.size, op_name))
  
  return {'losses': losses,
         'steps': steps,
         'results': results}
  


if __name__=='__main__':
    
    #### Create data with the lines below
    nominal = (1.21, 0.19, np.pi/2)
    times = np.arange(0, 30, 1/60)
    form = 'sine'
    tpp = None
    
    if form == 'sine':
        vals = np.random.binomial(1, p1_sine(times, nominal))
        prob = p1_sine(times, nominal)
    if form == 'saw':
        vals = np.random.binomial(1, p1_saw(times, nominal))
        prob = p1_saw(times, nominal)
    if form == 'square':
        vals = np.random.binomial(1, p1_square(times, nominal, tpp))
        prob = p1_square(times, nominal, tpp)
    if True:
        #hard code True or False if you want bit flip noise added in
        flip_prob = 0.05
        vals = [bit_flip(val, flip_prob) for val in vals]
        
    plt.plot(times, vals, ls="None", marker='.', label="Data Points")
    plt.plot(times, prob, label="Probability")
    plt.title("Input Data Set\nF = {} Hz, A = {}, P = {:.3f} Radians".format(nominal[0], nominal[1], nominal[2]))
    plt.xlim(0, 5)
    plt.xlabel("Time, s")
    plt.show()
        
    if False:
        resolution = 20
        a_range=np.linspace(0.14, 0.22, resolution)
        p_range=np.linspace(0, 2*np.pi, resolution)
        freq = 1.21
        losses, opt_prob, opt_params = two_dimensional_optimization(times, vals, freq, a_range, p_range, form=form)
        
        plt.figure(figsize=(6,6))
        plt.pcolormesh(a_range,p_range,losses)
        plt.xlabel("Amplitude")
        plt.ylabel("Phase")
        plt.colorbar()
        plt.show()
        
        plt.plot(times, prob, ls='dashed', label="Input")
        plt.plot(times, opt_prob, label="Optimized Fit")
        plt.xlabel("Time, s")
        plt.title("2-D MLE with Amplitude and Phase Variation\nInput: F = {:.3f} Hz, A = {:.3f}, P = {:.3f} Radians\nOutput: F = {:.3f} Hz, A = {:.3f}, P = {:.3f} Radians".format(nominal[0], nominal[1], nominal[2], opt_params[0], opt_params[1], opt_params[2]))
        plt.xlim(0, 6)
        plt.ylim(0, 1)
        plt.grid()
        plt.legend(loc='lower right')
        plt.show()
        
    if False:
        f_range=np.linspace(1.18, 1.22, 14)
        a_range=np.linspace(0.14, 0.22, 20)
        p_range=np.linspace(0, 2*np.pi, 20)
        losses, opt_prob, opt_params = three_dimensional_optimization(times, vals, f_range, a_range, p_range, form=form)
        
        plt.plot(times, prob, ls='dashed', label="Input")
        plt.plot(times, opt_prob, label="Optimized Fit")
        plt.xlabel("Time, s")
        plt.title("3-D MLE\nInput: F = {:.3f} Hz, A = {:.3f}, P = {:.3f} Radians\nOutput: F = {:.3f} Hz, A = {:.3f}, P = {:.3f} Radians".format(nominal[0], nominal[1], nominal[2], opt_params[0], opt_params[1], opt_params[2]))
        plt.xlim(0, 6)
        plt.ylim(0, 1)
        plt.grid()
        plt.legend(loc='lower right')
        plt.show()
        
    
    if False:    
        ##### Analyze data using the lines below     
        opt_times, opt_prob = MLE(times, vals, nominal, np.linspace(0.1, 2, 200), np.linspace(0, 0.3, 200), np.linspace(0, 1, 200),\
            form, tpp=None, input_f=None, input_a=None, input_p=None)
        
        #plot the actual data and optimized results below
        plt.plot(times, prob, label="Input Function")
        plt.plot(times, opt_prob, label="Reconstructed")
        plt.legend(loc='lower right')
        plt.grid()
        plt.xlabel("Time, s")
        plt.ylabel("1-State Probability")
        plt.legend(loc='lower right')
        plt.title("Combined Plot of Actual Input and Optimized Output")
        plt.xlim(0, 5)
        plt.ylim(0,1)
        plt.show()
    
    if False:
        init_f = 1.3
        init_a = 0.1
        init_p = 0.1
        tensorflow_optimization(times, vals, init_f, init_a, init_p, nepochs=4000, neval_period=10,\
                                learning_rate=0.0001, optimizer="adam", mini_batch_size=512, verbose=True, do_plot=True)
    
