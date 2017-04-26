import pdb
import sys
sys.path.append('../functions/.')
from CEMA import *

import pyjacob

# import py_eval_jacobian 
import cantera as ct
import numpy as np
import timeit

import scipy.linalg as LA 




#create gas from original mechanism file gri30.cti
gas = ct.Solution('Skeletal29_N.cti')
#reorder the gas to match pyJac
n2_ind = gas.species_index('N2')
# print(n2_ind)
specs = gas.species()[:]
# reorder to have the last entry as N2 -> when NOx not considered?
gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
        species=specs[:n2_ind] + specs[n2_ind + 1:] + [specs[n2_ind]],
        reactions=gas.reactions())


#set the gas state
P = 101325 # Pa
phi = 0.10 # Most reactive mixture fraction
T = 800  # K

gas.set_equivalence_ratio(phi, 'C2H4', 'O2:1.0, N2:3.76')
gas.TP = T,P

D,l,r = EigPbJac(gas)



# print D


real_lambda = D.real
max_eig_idx = np.argmax(real_lambda)
print max_eig_idx, ' maximum eigenvalue position'

# print max eigenvector
# print l[max_eig_idx,:]

pdb.set_trace()

ei = EI(D,l,r,max_eig_idx)

important_spec_idx = np.where(ei>1e-3) 
print important_spec_idx[0], " important species index"

real_idx = important_spec_idx[0] - 1
print real_idx, " real index of species (-1 for T)"

important_spec_ei = ei[real_idx]
# print important_spec_idx
# print important_spec_ei
important_spec_idx = np.asarray(important_spec_idx)


for i in range(len(real_idx)):
    
    print gas.species_name(real_idx[i])


## display species order in gas object
# specs = gas.species()[:]

# for s in specs:
# 	print s
































