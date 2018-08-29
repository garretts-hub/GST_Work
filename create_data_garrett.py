import sys
#sys.path.append("../../transmon-sim/physical_sim")
sys.path.append("C:/Users/GA28573/AppData/Local/Continuum/anaconda32/Lib/site-packages/qutip-4.3.1-py3.6-win-amd64.egg/qutip")


import warnings
import pygsti
from pygsti.objects import dataset as _ds
import Tomography as _tm
import NoiseSignal as _ns
import qutip as _qt
import numpy as np

# Return coherent noise that is meant to be added after a perfect gate
# The error_type can be one of:
#     'amplitude': add a rotation about the axis specified by the gateName
#     'coherent':  add a rotaion about all XYZ axes
#
# The error_generator object returns a random deviant that is used
# to calculate the error to generate
#
def coherent_noise(error_type, gateName, noise_signals, tm, error_amp):
    epsilons = [noise_signals[ii][tm]*error_amp for ii in range(len(noise_signals))]

    if error_type=='amplitude':
        if gateName=='Gx':
            SO = _qt.to_super(_qt.rx(np.pi*epsilons[0]))
        elif gateName=='Gy':
            SO = _qt.to_super(_qt.ry(np.pi*epsilons[0]))
        else:
            SO = _qt.to_super(_qt.qeye(2))
    elif error_type=='coherent':
        if len(epsilons)<3:
            raise Exception("ERROR: >= 3 noise signals required for XYZ coherent errors")
        S1 = _qt.to_super(_qt.rx(np.pi*epsilons[0]))
        S2 = _qt.to_super(_qt.ry(np.pi*epsilons[1]))
        S3 = _qt.to_super(_qt.rz(np.pi*epsilons[2]))
        SO = S1*S2*S3
    else:
        warnings.warn("Unsupported error_type specified")
        SO = _qt.to_super(_qt.qeye(2))

    return SO

# return Superoperator (in standard basis) for single qubit amplitude damping
#
def amplitude_damping(gamma):
    E0 = _qt.ket2dm(_qt.basis(2,0)) + np.sqrt(1-gamma)*_qt.ket2dm(_qt.basis(2,1))
    E1 = np.sqrt(gamma)*_qt.sigmap()
    
    SE = _qt.to_super(E0) + _qt.to_super(E1)

    return SE

# return Superoperator (in standard basis) for single qubit phase damping
#
def phase_damping(gamma):
    E0 = _qt.ket2dm(_qt.basis(2,0)) + np.sqrt(1-gamma)*_qt.ket2dm(_qt.basis(2,1))
    E1 = np.sqrt(gamma)*_qt.ket2dm(_qt.basis(2,1))
    
    SE = _qt.to_super(E0) + _qt.to_super(E1)

    return SE


def drifted_rotation(time, err):
    e = err[0]
    w = err[1]
    return e*np.cos(2*np.pi*w*time)

# Generate GST data under a specified model
# Create noise generators for each independent random variable based on the model
#
# Parameters:
#      gatestring_list: list of GateString() objects where each element specifies a different measurement to take
#      nSamples: total samples to take for each measurement
#      shots:    total shots to take at once (i.e. using the same output state, shots=nSamples => each experiment is run once and we use the probabilities to determine 0/1 counts,
#                                             shots=1 => a new experiment (with potentially different noise) is run to obtain each new count)
#      nRepeats:  Total number of times to repeat a sequence before moving on to the next
#      L: maximum length of germ sequences
#      error_type: (amplitude, coherent)
#      amplitude_type: type of amplitude error (fixed, non_markovian_1f)
#      error_amp:  amplitude factor (strength) of generated noise
#      time_per_exp: time for each complete experiment (i.e. one sequence of prep+fiducial+germs+fiducial+measurement) (in nanoseconds)
#      time_reprogram: time to reprogram a new sequence
#      seeds: list of seeds, one for each independent random variable
#
def generate_model_data_rotation_error(give_updates,gatestring_list, nSamples, shots, nRepeats, L, error_type, amplitude_type, error_amp, time_per_exp, time_reprogram, seeds, do_recalibration=True, xerr=(0,0), yerr=(0,0), zerr=(0,0)):
    time_samples = len(gatestring_list)*(L+6)*nSamples
    #print("total time_samples is",time_samples)
    if error_type=='amplitude':
        noise_sig = _ns.NoiseSignalNormal(initial_seed=seeds[0], resolution=1.0)
        noise_sig.configure_noise(1.0, 1.0, 1.0, 1.0, time_samples)
        noise_sig.init()
        noise_signals = [noise_sig]
    elif error_type=='coherent':
        noise_signals = []
        for ii in len(seeds):
            noise_sig = _ns.NoiseSignalNormal(initial_seed=seeds[ii], resolution=1.0)
            noise_sig.configure_noise(1.0, 1.0, 1.0, 1.0, time_samples)
            noise_sig.init()
            noise_signals += [noise_sig]
    elif error_type=='amplitude_damping':
        SError = amplitude_damping(error_amp)
    elif error_type=='phase_damping':
        SError = phase_damping(error_amp)

    if amplitude_type=='non_markovian_1f':
        calibration_interval = (time_per_exp*nRepeats + time_reprogram)*len(gatestring_list)
        error_amp_signal = _ns.NoiseSignalTelegraph(initial_seed=seeds[-1], resolution=1e-3)
        error_amp_signal.configure_noise(1.0, error_amp, 40, 1, 1e6, calibration_interval)

    if shots>nRepeats:
        shots = nRepeats
    if give_updates: print("{} Shots will be taken".format(shots))
    
    dataset = _ds.DataSet( collisionAction='aggregate' )
    labels = [('0',), ('1',)]
    ol0,ol1 = labels[0], labels[1]
    standard_gates = dict()
    standard_gates['Gi'] = _qt.to_super(_qt.qeye(2))
    standard_gates['Gx'] = _qt.to_super(_qt.rx(np.pi/2))
    standard_gates['Gy'] = _qt.to_super(_qt.ry(np.pi/2))
    standard_gates['Gz'] = _qt.to_super(_qt.rz(np.pi/2))

    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))

    nidx = 0
    samples_per_iteration = [nRepeats for ii in range(nSamples // nRepeats)]
    if sum(samples_per_iteration) < nSamples:
        samples_per_iteration += [nSamples - sum(samples_per_iteration)]
    if give_updates: print("there are {} samples per iteration".format(samples_per_iteration))
    
    
    first_loop = True
    gs_counts = {}
    gs_bitstring = {}
    gs_timestamped_counts = {} #counts the number of 1's and the timestep at which that was taken per gateset
    for gs in gatestring_list:
        gs_timestamped_counts[gs] = [str(gs)]
    
    tm = 0
    T=1
    for num_repeats in samples_per_iteration:
        if give_updates: print("************Beginning Tth iteration of each sequence, T = {}**************".format(T))
        # each of the outer loops involves a new calibration. model as a reconfiguration of the 1/f noise amplitude
        if amplitude_type=='non_markovian_1f':
            if first_loop or do_recalibration:
                error_amp_signal.init()
            else:
                error_amp_signal.next_interval()

        
        for gs in gatestring_list:
            #print("looking at gate",gs)
            numS = 0
            if first_loop:
                gs_counts[gs] = (0,0)
                gs_bitstring[gs] = ""
            while numS<num_repeats:
                c_shots = shots
                if numS+c_shots > num_repeats:
                    c_shots = num_repeats - numS

                rho = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))            
                if amplitude_type=='non_markovian_1f':
                    err_amp = error_amp_signal[tm]   # assume the noise is not faster than a gate
                else:
                    err_amp = error_amp
                for s in gs:
                    standard_gates = dict()
                    standard_gates['Gi'] = _qt.to_super(_qt.qeye(2))
                    standard_gates['Gx'] = _qt.to_super(_qt.rx(np.pi/2 + drifted_rotation(tm, xerr)))
                    standard_gates['Gy'] = _qt.to_super(_qt.ry(np.pi/2 + drifted_rotation(tm, yerr)))
                    standard_gates['Gz'] = _qt.to_super(_qt.rz(np.pi/2 + drifted_rotation(tm, zerr)))
                    G = standard_gates[s]
                    if error_type=='amplitude_damping' or error_type=='phase_damping':
                        Ge = SError
                    else:
                        Ge = coherent_noise(error_type, s, noise_signals, nidx, err_amp)
                    rho = Ge * G * rho
                    nidx += 1    
                p0 = (rho.dag()*rho0).norm()
                p1 = (rho.dag()*rho1).norm()
                c0 = np.random.binomial(c_shots, p0)
                c1 = c_shots - c0
                numS += c_shots
                # update running time
                tm += time_per_exp*c_shots
                if give_updates: print("at {:.3f} ns, {} shot(s) of {} gave {} zeros and {} ones.".format(tm, c_shots, str(gs), c0, c1))
                gs_tuple = (gs_counts[gs][0]+c0, gs_counts[gs][1]+c1)
                gs_counts[gs] = gs_tuple
                gs_bitstring[gs] = gs_bitstring[gs] + str(c1)
                gs_timestamped_counts[gs].append((c0, c1, tm*1e-9))
                
                #print(gs,gs_counts[gs],gs_bitstring[gs], gs_timestamped_counts[gs])
            tm += time_reprogram

        first_loop = False
        T += 1

    for gs in gatestring_list:
        counts = {}
        counts[ol0] = gs_counts[gs][0]
        counts[ol1] = gs_counts[gs][1]
        if give_updates: print(gs, counts)
        dataset.add_count_dict(gs, counts)

    dataset.done_adding_data()
    return (dataset,gs_timestamped_counts)

def generate_model_data(gateset, gates, gatestring_list, nSamples):
    num_qubits = 1
    num_basis = (2**num_qubits)*(2**num_qubits)
    standard_gates = dict()

    S = _qt.to_super(_qt.qeye(2**num_qubits))

    # pull off the gates from the gateset and convert to the standard basis
    for lbl,obj in gateset.iter_objs():
        if lbl in gates:
            T = _tm.Tomography(num_qubits, 2, 'pauli')
            data = obj.to_vector().reshape(num_basis, num_basis)
            SO = _qt.Qobj(inpt=data, shape=S.shape, dims=S.dims)
            T.set_super(SO)
            T.change_basis(new_basis='std')
            standard_gates[lbl] = T.get_super()

    dataset = _ds.DataSet( collisionAction='aggregate' )
    labels = [('0',), ('1',)]
    ol0,ol1 = labels[0], labels[1]

    rho0 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
    rho1 = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,1)))
    for gs in gatestring_list:

        rho = _qt.operator_to_vector(_qt.ket2dm(_qt.basis(2,0)))
        for s in gs:
            G = standard_gates[s]
            rho = G * rho

        p0 = (rho.dag()*rho0).norm()
        p1 = (rho.dag()*rho1).norm()
        counts = {}
        counts[ol0] = p0*nSamples
        counts[ol1] = p1*nSamples
        dataset.add_count_dict(gs, counts)
        #print("rhoOut =",rho, "probs (0,1) = ({0},{1})".format((rho.dag()*rho0).norm(), (rho.dag()*rho1).norm()))

    dataset.done_adding_data()
    return dataset


    
