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

gas.set_equivalence_ratio(phi, 'C2H4', 'O2:1.0, N2:3.76')
gas.TP = T,P

### reorder the gas to match pyJac
# n2_ind = gas.species_index('N2')
# # print(n2_ind)
# specs = gas.species()[:]
# # reorder to have the last entry as N2 -> when NOx not considered?
# gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
#         species=specs[:n2_ind] + specs[n2_ind + 1:] + [specs[n2_ind]],
#         reactions=gas.reactions())

# 1D flame simulation
width = 0.03
initial_grid = np.linspace(0,width, 7) # coarse first
f = ct.FreeFlame(gas, grid=initial_grid)
f.set_refine_criteria(ratio=2, slope=0.05, curve=0.1)

f.solve(loglevel=0, auto=True)
# print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))

# Chemical explosive mode analysis for 1D flame
eigenvalues, global_expl_indices, track_species = solve_eig_flame(f,gas)

spec_dictionary = get_species_names(track_species,gas)



manual_spec = ['CH', 'OH', 'H', 'HCO', 'H2O2']
manual_spec = ['CO', 'HCO', 'O', 'O2', 'C']

# TO SEE WHICH SPECIES ARE IN DICTIONARY OF MOST EXPLOSIVE SPECIES
# print spec_dictionary

## EIG FIGURE ## 

plt.figure()
plt.plot(f.grid,np.log10(eigenvalues[0,:]+1),'--')

pdb.set_trace()
# for i in range(len(eigenvalues[:,0])):
# 	plt.plot(f.grid,np.log10(eigenvalues[i,:]+1),'--')

plt.show()

plt.figure()
plt.subplot(2,1,1)
for spec in spec_dictionary:

	if spec in manual_spec:
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

