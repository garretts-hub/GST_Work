from qutip import *
from copy import copy
import sys
sys.path.append("C:/Users/GA28573/AppData/Local/Continuum/anaconda32/Lib/site-packages/qutip-4.3.1-py3.6-win-amd64.egg/qutip")
from pygsti import *
import numpy as np
import scipy.linalg as spl
from pygsti.tools.gatetools import fidelity

# Class to manage the tomography of simulations
# 
# To use:
#     1) T = Tomography(num_qubits=num_qubits, levels=levels)   # Initialize a Tomography object
#     2) tomography_states = T.tomography_states()              # Generate standard tomography states for the size and number of levels in system
#     3) for tomograpy_state in tomography_states:
#     4)      rho = (sample output state from simulation or experiment for process applied to tomography_state)
#     5)      T.accumulate_result(rho, tomography state #)
#     6) T.calculate_superoperator()   # calculates superoperator for given samples and tomography states
#     7) T.get_super()    # returns superoperator, other formats also supported: (get_chi(), get_choi(), get_krause())
#     Multiple samples can be added for each tomography state
#    
class Tomography(object):

    def __init__(self, num_qubits=1, levels=2, basis='std', use_entangled_states=False):
        self.num_qubits = num_qubits
        self.levels = levels

        # initialize the result of each tomography run
        self.num_states = 1 << 2*self.num_qubits
        num_basis_ml = pow(self.levels, self.num_qubits)
        psi0 = Qobj(shape=(num_basis_ml,1))
        psi0.dims = [[self.levels for ii in range(self.num_qubits)], 
                    [1 for ii in range(self.num_qubits)]]
        rho0 = ket2dm(psi0)
        self.tomo_results = [rho0 for ii in range(1 << 2*self.num_qubits)]
        self.valid_super = False
        self.basis = basis
        self.num_tomo_samples = [0 for ii in range(pow(levels, 2*num_qubits))]
        self.use_entangled_states = use_entangled_states

    # change the basis of a number, i.e. take the digits of the number interpreted in the in_base and convert to numbers interpreted in the output base
    # e.g. the number (2 in decimal) interpreted in binary would be 10, and the same number interpreted in base 3 would be 3 (i.e. 3*1 + 1*0)
    def change_base(self, num, in_base, out_base):
        in_num = num
        out_num = 0
        out_pos = 1
        while in_num>0:
            digit = in_num % in_base
            out_num += digit*out_pos
            in_num = in_num // in_base
            out_pos *= out_base
            
        return out_num

    # Generate tomography states (using entangled or single-qubit rotations)
    def tomography_states(self, return_dm=False, return_two_level=False):
        if self.use_entangled_states:
            return self.tomography_states_entangled(return_dm, return_two_level)
        else:
            return self.tomography_states_expt(return_dm)

    # generate the standard set of tomography states
    def tomography_states_entangled(self, return_dm=False, return_two_level=False):
        num_basis = 1 << self.num_qubits
        num_states = 1 << 2*self.num_qubits

        tomo_states = []
        basisp = [0,0]
        basism = [0,0]

        for snum in range(num_states):
            #dims_v = [[self.levels for li in range(self.num_qubits)], [1 for li in range(self.num_qubits)]]
            dims_v = [[2 for li in range(self.num_qubits)], [1 for li in range(self.num_qubits)]]
            shape_v = (num_basis, 1)
            if snum < num_basis:
                # states with a single basis state
                psi = Qobj(inpt=basis(num_basis, snum), shape=shape_v, dims=dims_v)
            elif snum < (num_basis + (num_basis*(num_basis-1))/2):
                basisp[1] += 1
                if ( basisp[1] == num_basis ):
                    basisp[0] += 1
                    basisp[1] = basisp[0] + 1
                psi = Qobj(inpt=(basis(num_basis, basisp[0]) + basis(num_basis, basisp[1]))/np.sqrt(2), shape=shape_v, dims=dims_v)
            else:
                basism[1] += 1
                if ( basism[1] == num_basis ):
                    basism[0] += 1
                    basism[1] = basism[0] + 1
                psi = Qobj(inpt=(basis(num_basis, basism[0]) + 1j*basis(num_basis, basism[1]))/np.sqrt(2), shape=shape_v, dims=dims_v)
            if return_dm:
                tomo_states += [ket2dm(psi)]
            else:
                tomo_states += [psi]

        if self.levels==2 or return_two_level:
            self.I = qeye(num_basis)
            return tomo_states
        else:
            # return the qubit basis states in kets of the proper shape
            num_basis_ml = pow(self.levels, self.num_qubits)
            # create transformation matrix
            # old code
            #identity_list = [[1 if ii==jj else 0 for ii in range(num_basis)] for jj in range(num_basis_ml)]
            # new code
            identity_list = [[0 for ii in range(num_basis)] for jj in range(num_basis_ml)]
            for ii in range(num_basis):
                jj = self.change_base(ii, 2, self.levels)
                identity_list[jj][ii] = 1
            dims_v = [[self.levels for ii in range(self.num_qubits)],[2 for ii in range(self.num_qubits)]]
            I = Qobj(inpt=identity_list, shape=(num_basis_ml, num_basis), dims=dims_v, type='oper')
            multi_tomo_states = []
            for tomo_state in tomo_states:
                psi = I*tomo_state
                psi.dims = [[self.levels for ii in range(self.num_qubits)], 
                    [1 for ii in range(self.num_qubits)]]
                if return_dm:
                    multi_tomo_states += [ket2dm(psi)]
                else:
                    multi_tomo_states += [psi]
            return multi_tomo_states

    def tomography_states_expt(self, return_dm=False):
        """generate precursor states that are fully unentangled"""

        if self.levels != 2:
            raise ValueError('This method only supports two-level systems')

        tomo_states = []
        num_states = 1 << 2*self.num_qubits
        unitaries = [qeye(2),  # prepare 0
                     rotation(sigmax(), np.pi),  # prepare 1
                     rotation(sigmay(), np.pi/2),  # prepare X
                     rotation(sigmax(), -np.pi/2)]  # prepare Y

        # Make the initial state in a way that keeps quTiP happy
        gs = basis(2, 0)
        init_state = tensor([gs]*self.num_qubits)

        for snum in range(num_states):
            # Get the indices for which gates on which qubits
            indices = self.get_indices(snum)
            # build the gate
            gate_list = [unitaries[i] for i in indices]
            gate = tensor(gate_list)

            # Operate on the starting state
            fin_state = gate * init_state
            if return_dm:
                tomo_states.append(ket2dm(fin_state))
            else:
                tomo_states.append(fin_state)

        return tomo_states

    # add to the running tomography results
    def accumulate_result(self, state, snum):
        if state.type=='ket':
            rho = ket2dm(state)
        elif state.type=='oper':
            rho = state
        else:
            raise TypeError("tomography result must be ket or operator")

        self.tomo_results[snum] += rho
        self.num_tomo_samples[snum] += 1

    # collapse a multi_level operator to a qubit operator
    def collapse_operator_multi_level(self, G):
        if self.levels==2:
            return G

        num_basis = 1 << self.num_qubits
        num_basis_ml = pow(self.levels, self.num_qubits)

        GdataIn = G.full().reshape(num_basis_ml*num_basis_ml,1)
        GdataOut = np.array([0+0j for ii in range(num_basis) for jj in range(num_basis)])
        
        for di in range(num_basis*num_basis):
            rdx = self.change_base(di // num_basis, 2, self.levels)
            cdx = self.change_base(di % num_basis, 2, self.levels)
            GdataOut[di] = GdataIn[rdx*num_basis_ml+cdx][0]

        Gout = Qobj(inpt=GdataOut.reshape(num_basis,num_basis),type='oper')
        return Gout

    # convert a Unitary to a multi-level Hilbert space
    def convert_operator_multi_level(self, G, annihilate=True):
        if self.levels==2:
            return G

        num_basis = 1 << self.num_qubits
        num_basis_ml = pow(self.levels, self.num_qubits)

        # data from the input operator
        GdataIn = G.full().reshape(num_basis*num_basis,1)

        # data for output operator
        # not sure if we should annihilate the states that are not in the Unitary or do nothing
        if annihilate:
            # zeros on diagonal for basis states including a higher level
            GdataOut = np.array([0+0j for ii in range(num_basis_ml) for jj in range(num_basis_ml)])
        else:
            # ones on diagonal for basis states including a higher level
            GdataOut = np.array([1+0j if ii==jj else 0 for ii in range(num_basis_ml) for jj in range(num_basis_ml)])

        for di in range(num_basis*num_basis):
            rdx = self.change_base(di // num_basis, 2, self.levels)
            cdx = self.change_base(di % num_basis, 2, self.levels)
            GdataOut[rdx*num_basis_ml+cdx] = GdataIn[di][0]

        Gout = Qobj(inpt=GdataOut.reshape(num_basis_ml,num_basis_ml),type='oper')
        Gout.dims = [[self.levels for ii in range(self.num_qubits)], 
                    [self.levels for ii in range(self.num_qubits)]]
        return Gout

    def get_tomo_offset(self, ii, jj):
        num_basis = 1 << self.num_qubits
        idxs = [1,0]
        offset = 0
        while not (idxs[0]==jj and idxs[1]==ii):
            offset += 1
            idxs[0] += 1
            if idxs[0]==num_basis:
                idxs[1] += 1
                idxs[0] = idxs[1] + 1
        return offset

    def calculate_superoperator(self, convert_to_qubit_basis=False):
        if self.use_entangled_states:
            self.calculate_superoperator_entangled(convert_to_qubit_basis)
        else:
            self.calculate_superoperator_expt(convert_to_qubit_basis)

    def calculate_superoperator_expt(self, convert_to_qubit_basis=False):
        """Calculate superoperator based on input states made from single-qubit
        operations.  """

        # First check if two level
        if self.levels != 2:
            raise ValueError("Method only supports two-level systems")

        minv = min(self.num_tomo_samples[:(1<<2*self.num_qubits)-1])
        maxv = max(self.num_tomo_samples[:(1<<2*self.num_qubits)-1])
        if not minv==maxv:
            raise Exception("Number of samples for each tomography state should be the same")
        if not minv>=1:
            raise Exception("Number of tomography samples should be greater than one")
        nmc = minv

        num_states = 1 << 2*self.num_qubits

        # Generate the input state to density-matrix-element conversion matrix
        input_states = self.tomography_states_expt(return_dm=True)
        coeff_list = self.find_linear_combos()

        dm_rhos = []

        # Initialize a blank density matrix
        psi0 = Qobj(shape=[2,1])
        ket = tensor([psi0]*self.num_qubits)

        rho_blank = ket2dm(ket)

        # Match shape of rho_blank to the input tomography data
        rho_blank.dims = self.tomo_results[0].dims

        # Make the appropriate linear combinations of tomography data
        for ii in range(num_states):
            # Outer loop: one linear combo per each density matrix element
            dmat = copy(rho_blank)
            for jj in range(num_states):
                # Inner loop: make appropriate linear combo of tomo results
                dmat += self.tomo_results[jj]/nmc * coeff_list[ii][jj][0]
            dm_rhos.append(dmat)

            vec = operator_to_vector(dmat)
            selector = basis(num_states, ii)
            selector.dims = vec.dims

            if ii == 0:
                SS = vec * selector.dag()
            else:
                SS += vec * selector.dag()

        # Do some reformatting to make qutip happy
        I = tensor([qeye(2)] * self.num_qubits)
        Sop = to_super(I)
        Sop.data = SS.data

        self.superoperator = Sop
        self.valid_super = True

    # normalize the results and generate the superoperator
    def calculate_superoperator_entangled(self, convert_to_qubit_basis=False):

        minv = min(self.num_tomo_samples[:(1<<2*self.num_qubits)-1])
        maxv = max(self.num_tomo_samples[:(1<<2*self.num_qubits)-1])
        
        if not minv==maxv:
            raise Exception("Number of samples for each tomography state should be the same")
        if not minv>=1:
            raise Exception("Number of tomography samples should be greater than one")
        nmc = minv

        num_basis = 1 << self.num_qubits
        num_states = 1 << 2*self.num_qubits
        minus_offset = (num_basis >> 1)*(num_basis-1)
        dm_rhos = []

        # normalize out the number of MC runs, does not necessarily make trace one
        for dmop in self.tomo_results:
            dmop /= nmc
            dm_rhos += [dmop]

        if convert_to_qubit_basis and self.levels>2:
            for ii in range(len(dm_rhos)):
                dm_rhos[ii] = self.collapse_operator_multi_level(dm_rhos[ii])

        for ii in range(num_basis):
            for jj in range(num_basis):
                if ii==jj:
                    dm_rho = dm_rhos[ii]
                elif ii < jj:
                    # construct non-physical DM using tomography states
                    pm_offset = self.get_tomo_offset(ii, jj)
                    dmp = dm_rhos[num_basis + pm_offset]
                    dmm = dm_rhos[num_basis + minus_offset + pm_offset]
                    dm_rho = dmp - 1j*dmm - ((1 - 1j)/2)*(dm_rhos[ii] + dm_rhos[jj])
                else:
                    # same as above, but need to swap indicies and use conjugate version
                    pm_offset = self.get_tomo_offset(jj, ii)
                    dmp = dm_rhos[num_basis + pm_offset]
                    dmm = dm_rhos[num_basis + minus_offset + pm_offset]
                    dm_rho = dmp + 1j*dmm - ((1 + 1j)/2)*(dm_rhos[ii] + dm_rhos[jj])

                dm_vec = operator_to_vector(dm_rho)
                sel_vec = basis(num_states, ii*num_basis + jj).dag()

                if ii==0 and jj==0:
                    SS = dm_vec * sel_vec
                else:
                    SS += dm_vec * sel_vec

        # identity operator with one's only for qubit basis states
        #
        Ilist = []
        if convert_to_qubit_basis:
            num_levels = 2
        else:
            num_levels = self.levels

        num_states_ml = pow(num_levels, 2*self.num_qubits)
        num_basis_ml = pow(num_levels, self.num_qubits)

        for ii in range(num_states):
            idx = self.change_base(ii, 2, num_levels)
            pl = [1 if idx==ii else 0 for ii in range(num_states_ml)]
            Ilist += [pl]

        I = Qobj(inpt=Ilist, shape=(num_states, num_states_ml), type='oper')
        SopData = SS*I

        SopDataObj = Qobj(inpt=SopData, shape=(num_states_ml, num_states_ml))
        Id_list = [qeye(num_levels) for idx in range(self.num_qubits)]
        Sop = to_super(tensor(Id_list))
        IdI = Qobj(I.dag()*I, shape=Sop.shape, dims=Sop.dims)
        Sop2 = Sop-IdI
        Sop.data = SopDataObj.data
        #Sop += Sop2

        self.superoperator = Sop
        self.valid_super = True

    def set_super(self, Sop):
        num_basis_ml = pow(self.levels, self.num_qubits)
        if not Sop.shape == (num_basis_ml*num_basis_ml, num_basis_ml*num_basis_ml):
            raise Exception("Invalid shape for superopertor")

        self.superoperator = Sop
        self.valid_super = True

    def get_choi(self):
        if not self.valid_super:
            raise Exception("No valid superoperator, must call calculate_superoperator() first")

        return to_choi(self.superoperator)

    def get_super(self):
        if not self.valid_super:
            raise Exception("No valid superoperator, must call calculate_superoperator() first")

        return self.superoperator

    def get_kraus(self):
        if not self.valid_super:
            raise Exception("No valid superoperator, must call calculate_superoperator() first")
        #TODO: ensure Kraus operators have correct Qobj.dims
        return to_kraus(self.superoperator)

    def get_chi(self):
        if not self.valid_super:
            raise Exception("No valid superoperator, must call calculate_superoperator() first")

        if self.levels==2:
            norm_factor = pow(2,self.num_qubits)*pow(2,self.num_qubits)
            return to_chi(self.superoperator)/norm_factor
        else:
            # convert the multi-level superoperator to a 2-level one and then return the chi matrix for this
            # we only keep states that are valid 2-level and do not normalize
            norm_factor = pow(2,self.num_qubits)*pow(2,self.num_qubits)
            Id_list = [qeye(2) for ii in range(self.num_qubits)]
            Sop = to_super(tensor(Id_list))
            SopML = self.get_super()
            for rdx in range(Sop.shape[0]):
                for cdx in range(Sop.shape[1]):
                    rdx_ml = self.change_base(rdx, 2, self.levels)
                    cdx_ml = self.change_base(cdx, 2, self.levels)
                    val = SopML.data[(rdx_ml,cdx_ml)]
                    if np.abs(val) > 0:
                        Sop.data[(rdx,cdx)] = val
            chi = to_chi(Sop)/norm_factor
            return to_chi(chi)

    # Pauli transform matrix to convert to/from the Pauli basis
    def pauli_transform_matrix(self, S=None):
        if S==None:
            S = self.get_super()

        # create conversion operators
        I = ket2dm(basis(2,0)) + ket2dm(basis(2,1))
        X = basis(2,0)*basis(2,1).dag() + basis(2,1)*basis(2,0).dag()
        Y = (-basis(2,0)*basis(2,1).dag() + basis(2,1)*basis(2,0).dag())*1j
        Z = ket2dm(basis(2,0)) - ket2dm(basis(2,1))

        bstates = [0 for ii in range(self.num_qubits)]
        num_basis = pow(4,self.num_qubits)
        OPs = [I,X,Y,Z]
        Odata = []
        for ii in range(pow(4,self.num_qubits)):
                    
            O = OPs[bstates[-1]]
            for bi in range(self.num_qubits-2,-1,-1):
                O = tensor(O,OPs[bstates[bi]])
            F = operator_to_vector(O).full()
            Odata += [[F[fi][0] for fi in range(F.shape[0])]]

            inc = False
            vi = 0
            while (vi < self.num_qubits) and (not inc):
                bstates[vi] = (bstates[vi] + 1) % 4
                if not bstates[vi]==0:
                    inc = True
                vi += 1

        OD = Qobj(inpt=Odata, dims=S.dims, shape=S.shape, superrep=S.superrep)
        ODt = Qobj(inpt=OD.trans(), superrep=S.superrep)
        return ODt

    # convert superoperator from standard basis to pauli basis
    def to_pauli_basis(self, S):
        PT = self.pauli_transform_matrix()
        PT = Qobj(inpt=PT, dims=PT.dims, shape=PT.shape, superrep=S.superrep)
        #SP = PT.dag()*S*PT/(2**self.num_qubits)
        SP = PT.dag()*S*PT/(1 << (2*self.num_qubits))
        return SP

    # convert from pauli basis to standard basis
    def to_standard_basis(self, S):
        PT = self.pauli_transform_matrix()
        PT = Qobj(inpt=PT, dims=PT.dims, shape=PT.shape, superrep=S.superrep)
        #SP = PT*S*PT.dag()/(2**self.num_qubits)
        SP = PT*S*PT.dag()/(1 << (2*self.num_qubits))
        return SP

    # change the basis, current basis's supported are the standard basis and the Pauli basis
    def change_basis(self, new_basis='std'):
        if new_basis=='pauli' or new_basis=='Pauli' or new_basis=='p':
            new_basis = 'pauli'
        else:
            new_basis = 'std'

        if not self.valid_super:
            print('Cannot change basis, superoperator is not valid')
            return
        if not self.levels==2:
            print('Cannot change basis, only two-level operators are supported')

        if new_basis==self.basis:
            return
        else:
            S = self.get_super()
            PT = self.pauli_transform_matrix(S)
            if new_basis=='pauli' and self.basis=='std':
                # change from standard basis to Pauli one
                SP = PT.dag()*S*PT/(2**self.num_qubits)
            elif new_basis=='std' and self.basis=='pauli':
                # change from Pauli basis to standard basis
                SP = PT*S*PT.dag()/(2**self.num_qubits)
            else:
                raise Exception("Invalid conversion basis")

            self.basis = new_basis
            self.set_super(SP)

    # average gate fidelity of a superoperator
    # the superoperator is assumes to be for the identity operator (gate must be removed first)
    # returns the average fidelity over all tomography states as well as the average leakage
    def average_gate_fidelity_old(self, S=None):
        if S==None:
            S = self.get_super()
        if not S.type=='super':
            raise Exception("Input to average gate fidelity should be a superoperator")
        fid_sum = 0
        num_ts = 0
        trace_sum = 0
        num_basis = 1 << self.num_qubits
        Sop = self.get_super()

        for ts in self.tomography_states():
            vs = operator_to_vector(ket2dm(ts))
            vs2 = Sop*vs
            fid_sum += (vs2.dag()*S*vs).tr() 
            num_ts += 1

            if self.levels>2:
                rhop = self.collapse_operator_multi_level(vector_to_operator(S*vs))
                trace_sum += rhop.tr()
            else:
                trace_sum += 1

        return (fid_sum/num_ts, 1-trace_sum/num_ts)

    # Calculate the process fidelity, which is defined as:
    # F(A,B) = Tr( sqrt( sqrt(A) * B * sqrt(A) ) ) **2
    # Where the process matricies A,B are Choi matricies in the Pauli basis
    #
    def process_fidelity(self, target):
        if not self.valid_super:
            raise Exception("No valid superoperator, must call calculate_superoperator() first")

        if target.type=='oper':
            TS = sprepost(target,target)
        elif target.type=='super':
            TS = target
        else:
            raise TypeError("Invalid operator type in fidelity calculation")

        # Choi matrix of objects SO in pauli basis
        A = to_choi(self.to_pauli_basis(self.get_super()))
        # Choi matrix of target SO in pauli basis
        B = to_choi(self.to_pauli_basis(TS))
        A = self.to_pauli_basis(to_choi(self.get_super()))
        B = self.to_pauli_basis(to_choi(TS))

        # uses fidelity function from pyGSTi
        # this function can return the fidelity value in objects of different shapes
        # depending on how it is calculated, so flatten to a single value
        mf = fidelity(A.full(), B.full())
        indicies = tuple([0 for ii in range(len(shape(mf)))])
        mf_val = mf[indicies]
        d = int(pow(self.levels, self.num_qubits))

        return mf_val
        
    # Calculate the average gate fidelity between the superoperator stored in the object 
    # and a target gate, which is defined as:
    # (d * fid_process + 1)/(1 + d)
    # Where 'd' is the dimension of the Hibert space, and fid_process is the process fidelity
    #
    def average_gate_fidelity(self, target):
        process_fid = self.process_fidelity(target)
        d = int(pow(self.levels, self.num_qubits))
        avg_gate_fid = (d*process_fid + 1)/(1+d)

        return(avg_gate_fid)

    def average_gate_infidelity(self, target):
        return 1.0-self.average_gate_fidelity(target)

    # Calculate the normalized diamond distance, defined as:
    # diamond_norm(S - Target)/2
    #
    def diamond_distance(self, target):
        if not self.valid_super:
            raise Exception("No valid superoperator, must call calculate_superoperator() first")

        if target.type=='oper':
            TS = sprepost(target,target)
        elif target.type=='super':
            TS = target
        else:
            raise TypeError("Invalid operator type in fidelity calculation")
        S = self.get_super()

        dd = ((TS-S).dnorm()/2.0)
        return(dd)

    def get_indices(self, snum):
        """Helper function to generate which-gate indices for tomo states

        Parameters:
        ---------
        snum : int
            State preparation index

        Returns indices, a list of length self.num_qubits that gives the
        correct state preparation pulse for each qubit
        """
        loop = self.num_qubits-1
        indices = []
        while loop >= 0:
            indices.append(snum // 4**loop)
            snum -= indices[-1] * 4**loop
            loop -= 1

        return indices

    def find_linear_combos(self):
        """Find linear combinations of prepared qubit states that
        generate a linearly independent set rho_j that spans the n-qubit
        density matrix space.

        Returns solns, a list of 4**n elements where each element is a
        4**n complex-valued np.array whose elements are the coefficients
        for the i-th n-qubit basis states.
        """
        n_states = 1 << 2 * self.num_qubits

        unitaries = [qeye(2),  # prepare 0
                     rotation(sigmax(), np.pi),  # prepare 1
                     rotation(sigmay(), np.pi/2),  # prepare X
                     rotation(sigmax(), -np.pi/2)]  # prepare Y

        state = basis(2,0)
        self.basis_states = [U * state for U in unitaries]

        # Generate the qubit basis
        tomo_states = self.tomography_states_expt(return_dm=True)

        # Make the coefficient matrix representing the density matrix elements
        coeffs = np.zeros([n_states, n_states], dtype='complex')
        for k in range(n_states):
            coeffs[:, k] = np.squeeze(operator_to_vector(tomo_states[k]).full())

        # Find the linear combination coefficients that allow us to generate
        # density matrix elements.
        solns = []
        for k in range(n_states):
            # initialize all elements at zero
            target = np.zeros([n_states, 1])
            # mark only the density matrix element of interest
            target[k] = 1
            # solve for the linear combination of input states
            b = np.linalg.solve(coeffs, target)
            solns.append(b)

        return solns
