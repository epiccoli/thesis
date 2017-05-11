import pdb
import sys
sys.path.append('../functions/')
from CEMA import *

import pyjacob

import cantera as ct
import numpy as np
import timeit
import matplotlib.pyplot as plt
import scipy.linalg as LA 

#create gas from original mechanism file gri30.cti
gas = ct.Solution('Skeletal29_N.cti')

#set the gas state
P = 101325 # Pa
phi = 1.0 # Most reactive mixture fraction
T = 1400  # K
fuel_species = 'C2H4'

gas.set_equivalence_ratio(phi, fuel_species, 'O2:1.0, N2:3.76')
gas.TP = T,P

# 1D flame simulation
width = 0.03
initial_grid = np.linspace(0,width, 7) # coarse first
f = ct.FreeFlame(gas, grid=initial_grid)
f.set_refine_criteria(ratio=2, slope=0.05, curve=0.1)

f.solve(loglevel=0, auto=True)
# print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))

# Chemical explosive mode analysis for 1D flame
eigenvalues, global_expl_indices, track_species= solve_eig_flame(f,gas)

# print eigenvalues[2,:]
# pdb.set_trace()

spec_dictionary = get_species_names(track_species,gas)
print spec_dictionary



manual_spec = ['CH', 'OH', 'H', 'HCO', 'H2O2']
# manual_spec = ['CO', 'HCO', 'O', 'O2', 'C']

# TO SEE WHICH SPECIES ARE IN DICTIONARY OF MOST EXPLOSIVE SPECIES
print spec_dictionary

# {'C': 18, 'CH': 19, 'CO': 9, 'H2O2': 6, 'OH': 4, 'H2': 2, 'H': 1, 'O': 3, 'HCO': 20, 'C2H5O': 28, 'CH3': 10, 'TXCH2': 21, 'SXCH2': 22, 'O2': 5}

## EIG FIGURE ## 

# TEST one eigenvalue at a time
# plt.figure()
# plt.plot(f.grid,np.log10(eigenvalues[0,:]+1),'--')
# PLOT all eigenvalues stored
for i in range(len(eigenvalues[:,0])):
	plt.plot(f.grid,eigenvalues[i,:]/1e6,'--',label=str(i))

# plt.plot(f.grid,np.log10(manual_eig[:]+1),'--',label='manual')

plt.legend()
plt.xlabel('1D flame domain')
plt.title('Maximum eigenvalues (from lowest to highest)')
plt.show()

plt.figure()
plt.subplot(2,1,1)
for spec in spec_dictionary:

	# if spec in manual_spec:
	plt.plot(f.grid,global_expl_indices[spec_dictionary[spec],:],label=spec)
plt.legend()
plt.title('Explosive index for selected species along 1D flame domain')

plt.subplot(2,1,2)
for spec in manual_spec:
	y_spec = f.Y[gas.species_index(spec),:]
	plt.semilogy(f.grid,y_spec, label=spec)
plt.legend()	
plt.title('Species mass fraction along 1D flame domain')
plt.show()

