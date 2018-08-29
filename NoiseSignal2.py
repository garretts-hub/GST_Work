import numpy as np
import matplotlib.pyplot as plt
import bisect

# Base class for all noisy signals
# Not to be used directly. Use the derived classes instead
# 
class NoiseSignal(object):
    """
    Generic noisy signal
    """

    def __init__(self, initial_seed=1234, resolution=1.0, time_unit=1e-9):
        self.resolution_factor = resolution
        self.amplitude = 1.0
        self.random_seed = initial_seed
        self.is_init = False
        self.interpolation_settings()
        self.time_unit=time_unit
        self.time_units_per_sec = 1/time_unit

        # This base random seed is used to create the individual seeds used to generate each instance of the noise
        self.seed_generator = np.random.RandomState(self.random_seed)

    def __repr__(self):
        return "NoiseSignal(resolution_factor={0}, amplitude={1}, seed={2}".format(self.resolution_factor, self.amplitude, self.random_seed)
    def __str__(self):
        return "NoiseSignal(resolution_factor={0}, amplitude={1}, seed={2}".format(self.resolution_factor, self.amplitude, self.random_seed)
    def __getitem__(self, key):
        if self.do_interpolation:
            return self.eval_interpolated(key)
        else:
            return self.eval(key)

    def reset(self):
        self.seed_generator = np.random.RandomState(self.random_seed)

    def get_new_seed(self):
        return self.seed_generator.randint(np.iinfo(np.int32).max)

    def interpolation_settings(self, do_interpolation=False, resolution_factor=1.0, itime=10.0):
        # enable flag
        self.do_interpolation = do_interpolation
        # resolution of samples from the signal
        self.interpolation_resolution = resolution_factor
        # window to interpolate over
        self.interpolation_time = itime

    def eval_interpolated(self, t):
        if self.do_interpolation:
            time_steps, data_samples = self.get_sample_data(t)
            ival = np.interp(t, time_steps, data_samples)*self.amplitude
            return ival
        else:
            return self.eval(t)
        
class NoiseSignalRandomWalk(NoiseSignal):
    '''
    Creates a random walk function, starting at 0, over a specified time interval.
    Create the object, configure it, then init. Do more inits or next_intervals as necessary
    '''
    def __str__(self):
        return "NoiseSignalRandomWalk(resolution_factor={0}, amplitude={1}, seed={2}".format(self.resolution_factor, self.amplitude, self.random_seed)

    def configure_noise(self, amplitude, resolution_factor, total_time):
        '''total_time will remain constant after configuration, but the times and points will get data appended to
        each time you call init or next_interval.
        amplitude: the step size of each walk; the value at a given time is the last time's value plus or minus this amplitude
        resolution: how many values per time_unit
        '''
        self.amplitude = amplitude
        self.resolution_factor = resolution_factor
        self.total_time = total_time #number of time_units per .init() call
        self.points = []
        self.times = []
        self.num_points = int(self.total_time*self.resolution_factor)
        self.next_interval_start_time = 0
        
    def init(self, random_seed=None):
        '''
        Starts a random walk from zero, up to total_time with num_points in between.
        '''
        if random_seed==None:
            random_seed = self.get_new_seed()
        rand_switch = np.random.RandomState(random_seed)
        num_points = self.num_points
        binaries = [0] + [self.amplitude if rand_switch.randint(2)==1 else -1*self.amplitude for i in range(num_points-1)]
        points = [0]*num_points
        points[0] = 0
        for i in range(1, num_points):
            points[i] = points[i-1] + binaries[i]
        self.points += points
        step = (self.total_time - 0)/self.num_points
        if len(self.times) == 0:
            self.times += [i*step for i in range(num_points)]  #array of times from 0 to total_time, units of time_unit
            self.next_interval_start_time = num_points*step
        else:
            self.times += [self.next_interval_start_time+i*step for i in range(num_points)]
            self.next_interval_start_time = self.next_interval_start_time + num_points*step
        self.is_init = True
        return random_seed
    
    def eval(self, t):
        '''Returns noise value at time t (specified in time_unit) in units of time_unit'''
        if t>self.total_time:
            raise ValueError("NoiseSignalNormal(): Requested time greater than total simulation time")
        if t<0:
            raise ValueError("NoiseSignalNormal(): Requested time < 0")

        tidx = int(t*self.resolution_factor)
        return self.points[tidx]
    
    def get_sample_data(self):
        '''returns a tuple with a time list and associated noise values at each time point'''
        times = self.times #returns all the calculated times, including any .next_intervals() you may do
        vals = self.points
        return times, vals
    
    def next_interval(self,random_seed=None):
        '''Does self.init again, but starts where the last init ended, not from zero.
        Will will add total_time amount of time to self.time (a total of num_points added)
        '''
        if random_seed==None:
            random_seed = self.get_new_seed()
        rand_switch = np.random.RandomState(random_seed)
        num_points = int(self.total_time * self.resolution_factor)
        binaries = [self.amplitude if rand_switch.randint(2)==1 else -1*self.amplitude for i in range(num_points)]
        points = [0]*num_points
        points[0] = binaries[0] + self.points[len(self.points) -1]
        for i in range(1, num_points):
            points[i] = points[i-1] + binaries[i]
        self.points += points
        step = (self.total_time - 0)/self.num_points
        self.times += [self.next_interval_start_time+i*step for i in range(num_points)]
        self.next_interval_start_time = self.next_interval_start_time + num_points*step
        return random_seed
    
    def plot_noise_signal(self):
        if self.is_init == False:
            raise ValueError("Must be initialized first!")
        marker = '.'
        linewidth = 1
        markersize = 1
        
        t, n = self.get_sample_data()
        plt.plot(t,n, marker=marker, linewidth=linewidth, markersize=markersize)
        plt.grid()
        plt.title("Random Walk Noise\n$\pm${} per step, {:.1f}(x{:.1e}s) per step".format(self.amplitude, self.total_time/self.num_points, self.time_unit) )
        plt.ylabel("Arbitrary Units")
        plt.xlabel("Time, in units of {:.1e} seconds".format(self.time_unit))
        plt.show()
        
class NoiseSignalSine(NoiseSignal):
    '''
    Creates a noise signal composed of a summation of sine waves over a specified time interval.
    create the object, configure it, then do an init and you're good to go. Do a next_interval or another init if necessary
    '''
    def __str__(self):
        return "NoiseSignalSine(resolution_factor={0}, seed={1}".format(self.resolution_factor, self.random_seed)

    def configure_noise(self, resolution_factor, freq_list, amp_list, phase_list, total_time):
        '''
        resolution factor: how many points per time_unit
        freq_list: a list or tuple of the all the desired frequencies, in Hz
        amp_list: a list or tuple of the amplitudes for all frequencies
        phase_list: a list or tuple of the phases of all the waves (all lists should be the same length)
        total_time: total # of time_units you want to cover per .init call
        '''
        self.resolution_factor = resolution_factor
        self.total_time = total_time #number of time_units per .init() call
        self.list_of_point_lists = []
        self.summed_points = []
        self.times = []
        self.num_points = int(self.total_time*self.resolution_factor)
        self.step = self.total_time/self.num_points #the time interval between time points, in units of time_units
        self.next_interval_start_time = 0
        self.freq_list = freq_list
        self.amp_list = amp_list
        self.phase_list = phase_list
        
    def init(self):
        times = [i*self.step for i in range(0,self.num_points+1)] #step = totaltime/numpoints, numpoits = int(totaltime*resfactor)
        #print(times)
        list_of_point_lists = [0]*len(self.freq_list)
        for i in range(len(self.freq_list)):
            '''find values for each frequency, then su mthem up later'''
            freq = self.freq_list[i]
            amp = self.amp_list[i]
            phase = self.phase_list[i]
            value_list = []
            for t in times:
                val = amp*np.sin(2*np.pi*freq*t*self.time_unit + phase)
                value_list.append(val)
            list_of_point_lists[i] = value_list
        self.list_of_point_lists = list_of_point_lists
        #each frequency's list is done, now add them all together in one list
        summed_vals = [0]*len(times)
        for t_index in range(len(times)):
            for wave_index in range(len(self.list_of_point_lists)):
                summed_vals[t_index] += self.list_of_point_lists[wave_index][t_index]
        #print(summed_vals)
        self.summed_points += summed_vals
        #add these lines at the end of the init call
        self.times += times
        self.next_interval_start_time = self.times[-1] + self.step
                
    def get_sample_data(self):
        '''returns a tuple with a time list and associated noise values at each time point'''
        times = self.times #returns all the calculated times, including any .next_intervals() you may do
        vals = self.summed_points #returns the summed sine wave for each point in self.times
        return times, vals
    
    def add_random_noise(self, sigma):
        '''add normally distributed noise to the overall summed_values
        Specify the std deviation, sigma, of the noise, which centers on zero'''
        for i in range(len(self.summed_points)):
            self.summed_points[i] += np.random.normal(0,sigma)
    
    def next_interval(self):
        times = [self.next_interval_start_time+i*self.step for i in range(self.num_points)] 
        #print(times)
        for i in range(len(self.freq_list)):
            '''find values for each frequency, then su mthem up later'''
            freq = self.freq_list[i]
            amp = self.amp_list[i]
            phase = self.phase_list[i]
            value_list = []
            for t in times:
                val = amp*np.sin(2*np.pi*freq*t*self.time_unit + phase)
                value_list.append(val)
            self.list_of_point_lists[i].extend(value_list)

        #each frequency's list is done, now add them all together in one list
        summed_vals = [0]*len(times)
        for t_index in range(len(times)):
            for wave_index in range(len(self.list_of_point_lists)):
                summed_vals[t_index] += self.list_of_point_lists[wave_index][t_index]
        #print(summed_vals)
        self.summed_points += summed_vals
        #add these lines at the end of the init call
        self.times += times
        self.next_interval_start_time = self.times[-1] + self.step
        
    def plot_noise_signal(self, show_components=True, seconds=False):
        '''Shows a plot of the individual waves, and a plot of the summed wave, if components = True'''
        if seconds: 
            times = self.times*self.time_unit
            time_label = "s"
        else:
            times = self.times
            time_label = "time unit"
        times = self.times
        summed = self.summed_points
        components = self.list_of_point_lists
        marker ='.'
        linewidth=1
        markersize = 1
        if show_components == True:
            for wave_index in range(len(components)):
                wave= components[wave_index]
                freq = self.freq_list[wave_index]
                plt.plot(times, wave, label="{:.2f} Hz".format(freq), marker=marker, linewidth=linewidth, markersize=markersize)
            plt.grid()
            plt.xlabel("Time in {}".format(time_label))
            plt.ylabel("Amplitude, a.u.")
            plt.legend()
            plt.title("Component Waves")
            plt.show()
        
        plt.plot(times, summed, marker=marker, linewidth=linewidth,markersize=markersize)
        plt.grid()
        plt.xlabel("Time in {}".format(time_label))
        plt.ylabel("Amplitude, a.u.")
        plt.title("Sinusoidal Noise Signal")
        plt.show()
    
    def eval(self, t):
        if t>self.total_time:
            raise ValueError("NoiseSignalSine(): Requested time greater than total simulation time")
        if t<0:
            raise ValueError("NoiseSignalSine(): Requested time < 0")

        tidx = int(t*self.resolution_factor)
        #print("time t is {}".format(t))
        #print("tindex is {}".format(tidx))
        #print("total time is {}".format(self.total_time))
        return self.summed_points[tidx]
        
            
    

class NoiseSignalNormal(NoiseSignal):

    def __str__(self):
        return "NoiseSignalNormal(resolution_factor={0}, amplitude={1}, seed={2}, mean{3}, variance{4})".format(self.resolution_factor, self.amplitude, self.random_seed, self.mean, self.variance)
    def __repr__(self):
        if self.is_init:
            str_val = "normal_deviants = [ "
            for ii in range(len(self.points)):
                str_val += "{0} ".format(self.points[ii])
            str_val += "]"

            return str_val
        else:
            return "normal_deviants = []"

    def plot_amplitudes(self):
        times = [ti/self.resolution_factor for ti in range(int(self.total_time*self.resolution_factor))]
        plt.plot(times, [sv*self.amplitude for sv in self.points])
        plt.show()

    def configure_noise(self, amplitude, mean, variance, resolution_factor, total_time):
        self.mean = mean
        self.variance = variance
        self.amplitude = amplitude
        self.resolution_factor = resolution_factor
        self.std_dev = np.sqrt(variance)
        self.total_time = total_time

    def init(self, random_seed=None):
        if random_seed==None:
            random_seed = self.get_new_seed()
        rand_switch = np.random.RandomState(random_seed)
        num_points = int(self.total_time * self.resolution_factor)

        self.points = [self.std_dev*rand_switch.randn() + self.mean for ii in range(num_points)]
        self.is_init = True

        return random_seed

    def eval(self, t):
        if t>self.total_time:
            raise ValueError("NoiseSignalNormal(): Requested time greater than total simulation time")
        if t<0:
            raise ValueError("NoiseSignalNormal(): Requested time < 0")

        tidx = int(t*self.resolution_factor)
        return self.points[tidx]*self.amplitude

    def get_sample_data(self, t):
        lidx = int(np.floor(t*self.resolution_factor))
        ltime = lidx/self.resolution_factor
        uidx = int(np.ceil(t*self.resolution_factor))
        utime = uidx/self.resolution_factor

        sample_times = [ltime]
        sample_vals = [self.points[lidx]]
        if uidx<len(self.points):
            sample_times += [utime]
            sample_vals += [self.points[uidx]]

        while (lidx>0 or (uidx+1)<len(self.points)) and ((len(sample_times)<2) or (utime-ltime)<self.interpolation_time):
            if lidx>0:
                lidx -= 1
                ltime = lidx/self.resolution_factor
                sample_vals = [self.points[lidx]] + sample_vals
                sample_times = [ltime] + sample_times
            if (uidx+1)<len(self.points):
                uidx += 1
                utime = uidx/self.resolution_factor
                sample_vals += [self.points[uidx]]
                sample_times += [utime]

        if len(sample_times)<2:
            raise ValueError("Cannot find enough sample data")

        return sample_times, sample_vals

class NoiseSignalTelegraph(NoiseSignal):

    def __str__(self):
        return "NoiseSignalTelegraph(resolution_factor={0}, amplitude={1}, seed={2}, exponent={3}, total_fluctuators={4}, start_freq={5}, stop_freq={6}, total_time={7})".format(self.resolution_factor, self.amplitude, self.random_seed, self.exponent, self.total_fluctuators, self.start_freq, self.stop_freq, self.total_time)
    def __repr__(self):
        if self.is_init:
            str_val = "flux_switch_times = [ "
            for ii in range(len(self.switching_times)):
                str_val += "{0}:{1} ".format(self.switching_times[ii],self.amplitude_sums[ii])
            str_val += "]"
        else:
            return "flux_switch_times = []"

        return str_val

    def plot_amplitudes(self):
        prev_time = -1
        switch_times = []
        switch_amps = []
        for ii in range(len(self.switching_times)):
            sw_time = self.switching_times[ii]
            sw_sum = self.amplitude*self.amplitude_sums[ii]

            if prev_time >= 0:
                switch_times += [sw_time]
                switch_amps += [prev_amp]
            switch_times += [sw_time]
            switch_amps += [sw_sum]
            prev_time = sw_time
            prev_amp = sw_sum

        plt.plot(switch_times, switch_amps)
        plt.show()
        
    def plot_noise_signal(self):
        if self.is_init != True:
            raise ValueError("Noise signal must first be initialized!")
        num_points = int(self.total_time/self.resolution_factor)
        step = self.total_time/num_points
        times = [i*step for i in range(num_points)] #step = totaltime/numpoints, numpoits = int(totaltime*resfactor)
        vals = []
        for t in times:
            vals.append(self[t])
        plt.plot(times, vals,marker='.',linewidth=1)
        plt.xlabel("Time in time_unit")
        plt.ylabel("Noise, arb. units")
        plt.title("Telegraph Noise\nFrom {} Hz to {} Hz".format(self.start_freq, self.stop_freq))
        plt.grid()
        plt.show()

    def configure_noise(self, exponent, amplitude, total_fluctuators, start_freq, stop_freq, total_time):
        self.exponent = float(exponent)
        self.amplitude = float(amplitude / (2.0 * np.pi))
        self.total_fluctuators = total_fluctuators
        self.start_freq = float(start_freq)
        self.stop_freq = float(stop_freq)
        self.total_time = float(total_time)
        self.create_fluctuators()
        
    # create structures needed for noise generator
    def create_fluctuators(self):

        freq_area = np.log(self.stop_freq/self.start_freq)/((self.total_fluctuators-1))
        freq_fact = np.exp(freq_area)

        if ( self.exponent<0.0 or self.exponent>2.0 ):
            raise ValueError("Invalid exponent for Telegraph noise, should be > 0.0 and <= 2.0")

        self.freq_points = [self.start_freq*np.power(freq_fact,ii) for ii in range(self.total_fluctuators+1)]
        if self.exponent==1.0:
            self.freq_amps = [np.sqrt(2.0 * (np.log(self.freq_points[ii+1])-np.log(self.freq_points[ii]))) for ii in range(self.total_fluctuators)]
        else:
            self.freq_amps = [np.sqrt(1.0/(1.0-self.exponent) * (np.power(self.freq_points[ii+1],1.0-self.exponent) - np.power(self.freq_points[ii],1.0-self.exponent))) for ii in range(self.total_fluctuators)]
        del self.freq_points[self.total_fluctuators]

    # initialize a new set of fluctuators
    def init(self, random_seed=None):
        if random_seed==None:
            random_seed = self.get_new_seed()
        self.rand_switch = np.random.RandomState(random_seed)

        # initialize the state of each fluctuator (+-1.0)
        self.f_state = [1.0 if self.rand_switch.randint(2)==0 else -1.0 for ii in range(self.total_fluctuators)]

        # initialize the switching time for each fluctuator
        # next switching time is exponentially distributed, so we can
        # generate a uniform rv and transform it to an exponential rv.
        # also convert from seconds to nanoseconds (1 s = 10^9 ns) and start at a random point in the interval [0,1]
        self.f_next_switch = [(self.time_units_per_sec)*self.rand_switch.uniform()*(-np.log(1-self.rand_switch.uniform())/self.freq_points[ii]) for ii in range(self.total_fluctuators)]
         
        # create an array of switching times (for any fluctuator) and the sum of all the fluctuators at that time
        #
        self.amplitude_sums = [sum(self.f_state[ii]*self.freq_amps[ii] for ii in range(self.total_fluctuators))]
        self.switching_times = [0.0]

        tval = 0.0;
        while ( tval < self.total_time ):
            # find the fluctuator that is to switch next
            next_switch_idx = self.f_next_switch.index(min(self.f_next_switch))
            tval = self.f_next_switch[next_switch_idx]

            # update the switching time for this fluctuator
            self.f_next_switch[next_switch_idx] += self.time_units_per_sec*(-np.log(1-self.rand_switch.uniform())/self.freq_points[next_switch_idx])

            # update the state of this fluctuator
            self.f_state[next_switch_idx] *= -1.0

            # add the switching time and new sum to the lists
            self.switching_times += [tval]
            self.amplitude_sums += [self.amplitude_sums[len(self.amplitude_sums)-1] + self.f_state[next_switch_idx]*2.0*self.freq_amps[next_switch_idx]]

        self.is_init = True
        return random_seed

    # move to a new interval of time, i.e., continue the 1/f signal without re-initializing
    def next_interval(self):

        self.f_next_switch = [next_switch-self.total_time for next_switch in self.f_next_switch]

        self.amplitude_sums = [sum(self.f_state[ii]*self.freq_amps[ii] for ii in range(self.total_fluctuators))]
        self.switching_times = [0.0]

        tval = 0.0;
        while ( tval < self.total_time):
            # find the fluctuator that is to switch next
            next_switch_idx = self.f_next_switch.index(min(self.f_next_switch))
            tval = self.f_next_switch[next_switch_idx]

            # update the switching time for this fluctuator
            self.f_next_switch[next_switch_idx] += self.time_units_per_sec*(-np.log(1-self.rand_switch.uniform())/self.freq_points[next_switch_idx])

            # update the state of this fluctuator
            self.f_state[next_switch_idx] *= -1.0

            # add the switching time and new sum to the lists
            self.switching_times += [tval]
            self.amplitude_sums += [self.amplitude_sums[len(self.amplitude_sums)-1] + self.f_state[next_switch_idx]*2.0*self.freq_amps[next_switch_idx]]

    # return the value of the noise at the indicated time
    def eval(self, t):
        tdx = bisect.bisect_left(self.switching_times, t)
        if ( tdx >= len(self.switching_times) ):
             tdx -= 1

        if tdx==0:
            e_val = self.amplitude_sums[0]*self.amplitude
        else:
            e_val = self.amplitude_sums[tdx-1]*self.amplitude
        return e_val


    def get_sample_data(self, t):
        uidx = 0
        while uidx<len(self.switching_times) and self.switching_times[uidx]<t:
            uidx += 1

        if uidx>0:
            lidx = uidx - 1
        else:
            lidx = 0

        ltime = self.switching_times[lidx]
        utime = self.switching_times[uidx]
        sample_times = [ltime]
        sample_vals = [self.amplitude_sums[lidx]]
        if not ltime == utime:
            sample_times += [utime]
            sample_vals += [self.amplitude_sums[uidx]]

        while (lidx>0 or (uidx+1)<len(self.switching_times)) and ((len(sample_times)<2) or (utime-ltime)<self.interpolation_time):
            if lidx>0:
                lidx -= 1
                ltime = self.switching_times[lidx]
                sample_vals = [self.amplitude_sums[lidx]] + sample_vals
                sample_times = [ltime] + sample_times
            if (uidx+1)<len(self.switching_times):
                uidx += 1
                utime = self.switching_times[uidx]
                sample_vals += [self.amplitude_sums[uidx]]
                sample_times += [utime]

        if len(sample_times)<2:
            raise ValueError("Cannot find enough sample data")

        return sample_times, sample_vals

if __name__=='__main__':
    total_time = 1000
    time_unit = 1e-3

    sig = NoiseSignalTelegraph(initial_seed=1234, time_unit=time_unit, resolution=1)
    sig.configure_noise(exponent=1, amplitude=1, total_fluctuators=40, start_freq=0.1, stop_freq=5, total_time=total_time)
    sig.init()
    sig.plot_noise_signal()
    
    amp = 1
    res = 1
    freq_list = (2,5, 8) #Hz
    amp_list = (0.1, 0.06, 0.04) #arb. units
    phase_list = (0,0,0)
    sig = NoiseSignalRandomWalk(initial_seed = 1, time_unit=time_unit)
    sig.configure_noise(amp, res, total_time)
    sig.init()
    sig.plot_noise_signal()
    
    sine = NoiseSignalSine(time_unit = 1e-3)
    sine.configure_noise(resolution_factor=0.1, freq_list=freq_list, amp_list=amp_list, phase_list=phase_list, total_time=1000)
    sine.init()
    sine.next_interval()
    sine.plot_noise_signal(True)



