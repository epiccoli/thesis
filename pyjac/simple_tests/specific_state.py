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



QSS_species = ['N2','H','H2','O','OH','O2','H2O2','H2O','HO2','CO','CH3','CH2O','CO2','CH4','C2H2','C2H4','CH2CO','C2H6']

#create gas from original mechanism file gri30.cti
gas = ct.Solution('Skeletal29_N.cti')
#reorder the gas to match pyJac
N2_ind = gas.species_index('N2')
# print(N2_ind)
specs = gas.species()[:]
# reorder to have the last entry as N2 -> when NOx not considered?
gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
        species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
        reactions=gas.reactions())


#set the gas state
P = 101325 # Pa
phi = 0.10 # Most reactive mixture fraction
T = 800  # K

gas.set_equivalence_ratio(phi, 'C2H4', 'O2:1.0, N2:3.76')
gas.TP = T,P

y = np.zeros(gas.n_species)
y_massfr = np.concatenate([gas.Y[0:N2_ind], gas.Y[N2_ind+1:]])

# TEST EIG without QSS species

y_massfr_counter = 0
for i in range(len(QSS_species)):
	if i < N2_ind:
		# check if species is in QSS, if not, strip element from y_massfr
		if gas.species_name(i) not in QSS_species:
			print gas.species_name(i), " is not in QSS"
			y_massfr = np.delete(y_massfr,y_massfr_counter)

		y_massfr_counter += 1

	if i > N2_ind:
		# check if species is in QSS, if not, strip element from y_massfr
		if gas.species_name(i) not in QSS_species:
			print gas.species_name(i), " is not in QSS"
			y_massfr = np.delete(y_massfr,y_massfr_counter)
		y_massfr_counter += 1

pdb.set_trace()
y[0] = T
y[1:] = y_massfr

jac = create_jacobian(T,P,y)

D, vl, vr = LA.eig(jac, left = True)

print np.amax(D)

pdb.set_trace()



# print D,


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
































