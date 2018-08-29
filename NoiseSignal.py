import sys
import numpy as np
import matplotlib.pyplot as plt

# Base class for all noisy signals
# Not to be used directly. Use the derived classes instead
# 
class NoiseSignal(object):
    """
    Generic noisy signal
    """

    def __init__(self, initial_seed=1234, resolution=1.0):
        self.resolution_factor = resolution
        self.amplitude = 1.0
        self.random_seed = initial_seed
        self.is_init = False
        self.interpolation_settings()

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
        #print("Freq area is: {}".format(freq_area))
        #print("Freq fact is: {}".format(freq_fact))

        if ( self.exponent<0.0 or self.exponent>2.0 ):
            raise ValueError("Invalid exponent for Telegraph noise, should be > 0.0 and <= 2.0")

        self.freq_points = [self.start_freq*np.power(freq_fact,ii) for ii in range(self.total_fluctuators+1)]
        #print(self.freq_points)
        if self.exponent==1.0:
            self.freq_amps = [np.sqrt(2.0 * (np.log(self.freq_points[ii+1])-np.log(self.freq_points[ii]))) for ii in range(self.total_fluctuators)]
        else:
            self.freq_amps = [np.sqrt(1.0/(1.0-self.exponent) * (np.power(self.freq_points[ii+1],1.0-self.exponent) - np.power(self.freq_points[ii],1.0-self.exponent))) for ii in range(self.total_fluctuators)]
        del self.freq_points[self.total_fluctuators]
        #print(self.freq_amps)

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
        self.f_next_switch = [1e9*self.rand_switch.uniform()*(-np.log(1-self.rand_switch.uniform())/self.freq_points[ii]) for ii in range(self.total_fluctuators)]
         
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
            self.f_next_switch[next_switch_idx] += 1e9*(-np.log(1-self.rand_switch.uniform())/self.freq_points[next_switch_idx])

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
            self.f_next_switch[next_switch_idx] += 1e9*(-np.log(1-self.rand_switch.uniform())/self.freq_points[next_switch_idx])

            # update the state of this fluctuator
            self.f_state[next_switch_idx] *= -1.0

            # add the switching time and new sum to the lists
            self.switching_times += [tval]
            self.amplitude_sums += [self.amplitude_sums[len(self.amplitude_sums)-1] + self.f_state[next_switch_idx]*2.0*self.freq_amps[next_switch_idx]]

    # return the value of the noise at the indicated time
    def eval(self, t):
        tdx = 0
        while tdx<len(self.switching_times) and self.switching_times[tdx]<t:
            tdx += 1

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
    max_time = 5000 #you'll get data points from 0 to this value at time intervals specified by the decimal precision below, seconds
    decimals = 3
    signal = NoiseSignalTelegraph(1)
    exponent = 1.0
    amplitude = 0.01
    total_fluctuators = 100
    start_freq = 1
    stop_freq = 100e3
    signal.configure_noise(exponent, amplitude, total_fluctuators, start_freq, stop_freq, total_time)
    
    #signal = NoiseSignalNormal(1)
    #signal.configure_noise(1e-5, 1.0, 1e-2, 1.0, 200)
    total_time = max_time/10
    total_time10 = total_time*10**decimals
    
    
    sample_val = []
    for ii in range(1):
        signal.init()
        signal.interpolation_settings(True)
        tvals = [t/10**decimals for t in range(int(total_time10))] #gives you time from 0 to one order of magnitude less than total_time
        
        for t in tvals:
            sample_val += [signal.eval_interpolated(t)]
            #print("Time is {:.4}".format(t))
            #print("Signal_val is {}".format([signal.eval_interpolated(t)]))
        tt = tvals
        #print(tt)
            
    for ii in range(1,10):
        signal.next_interval()
        
        tvals = [t/10**decimals for t in range(int(total_time10))] #gives you time from 0 to one order of magnitude less than total_time
        for t in tvals:
            sample_val += [signal.eval(t)]
        tt += [t+ii*total_time for t in tvals]
    
    plt.plot(tt, sample_val)
    plt.grid()
    plt.xlabel("Time, seconds")
    plt.ylabel("Noise, arbitrary units")
    plt.title("$1/f$ Noise with $f$: {} to {} Hz".format(start_freq, stop_freq))
    plt.show()
    print("Length of values {}".format(len(sample_val)))






