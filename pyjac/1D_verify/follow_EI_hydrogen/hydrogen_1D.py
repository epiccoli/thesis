import pdb
import sys
sys.path.append('../../hydrogen/')       # path to pyjacob.so
import pyjacob
sys.path.append('../../functions/')     # path to CEMA.py
from CEMA import *

import os
import cantera as ct
import numpy as np
import timeit
import matplotlib.pyplot as plt
import scipy.linalg as LA 

#create gas from original mechanism file gri30.cti
gas = ct.Solution('Li_2003.cti')
fuel_species = 'H2'

#set the gas state
P = 101325 # Pa
phi = 1.0 # Most reactive mixture fraction
T = 1110 # K

# Eigenvalue to follow
eig2track = -9 # 0 for maximum eig

# flame filename to store it and load it again
flame_filename = 'saved_T={:.0f}_phi={:.4f}.xml'.format(T,phi)
# graph filename
graph_filename = 'MODE{}_fuel={}_T={:.0f}_phi={:.4f}.pdf'.format(str(abs(eig2track)),fuel_species,T,phi)

gas.set_equivalence_ratio(phi, fuel_species, 'O2:1.0, N2:3.76')
gas.TP = T,P

# 1D flame simulation
width = 0.03
initial_grid = np.linspace(0,width, 7) # coarse first

# initial_grid = np.linspace(0,width, 300)
f = ct.FreeFlame(gas, grid=initial_grid)

# CREATE FLAME or RESTORE IT (based on if the flame data for the operating point exists)
###############################################################################################
if os.path.isfile(flame_filename):
    f.restore(filename=flame_filename)
else: 
    ### Create Flame object and save it:
    f.set_refine_criteria(ratio=2, slope=0.02, curve=0.05)
    f.solve(loglevel=1, auto=True)
    f.save(filename=flame_filename) #, name='', description='solution of hot stoech. OLI')

print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))
###############################################################################################

##### CEMA for the flame #####

eig_CEM, global_expl_indices, track_species, max_eig_loc = solve_eig_flame(f,gas,eig2track)

# Check temperature profile
# plt.figure()
# plt.plot(f.grid[max_eig_loc:],f.T[max_eig_loc:])
# plt.show()




# Build vector with names of EI vector [T + species] 
EI_tags = get_names(gas)


forward_grid = f.grid[max_eig_loc:]
backward_grid = f.grid[:max_eig_loc]
CEM_fw=eig_CEM[max_eig_loc:]
CEM_bw=eig_CEM[:max_eig_loc]

fig = plt.subplots(figsize=(16,10))
ax1=plt.subplot(2,1,1)

ax1.plot(forward_grid,np.array(CEM_fw)/1e6,linestyle='--',marker='.',label='fw tracking')
ax1.plot(backward_grid,np.array(CEM_bw)/1e6,linestyle='--',marker='.',label='bw tracking')
plt.legend()
plt.title('Timescale of mode {}'.format(str(abs(eig2track))))

ax2 = plt.subplot(2,1,2,sharex=ax1)

for i in range(len(track_species)):

    ax2.plot(f.grid,global_expl_indices[track_species[i],:],linestyle='--',marker='.',label=EI_tags[track_species[i]])

plt.legend()    
plt.title('Most important EI of mode {}'.format(str(abs(eig2track))))

plt.suptitle('Eigenvalue number {} tracked from flame front'.format(str(abs(eig2track))))
plt.savefig(graph_filename)
plt.show()



## SPECIES PLOT
# manual_spec = ['CH', 'OH', 'H', 'HCO', 'H2O2']
# manual_spec = ['CO', 'HCO', 'O', 'O2', 'C']
# plt.figure()
# for spec in manual_spec:
#   y_spec = f.Y[gas.species_index(spec),:]
#   plt.semilogy(f.grid,y_spec, label=spec)



## EIG FIGURE ## 

# PLOT all eigenvalues stored
# for i in range(len(eigenvalues[:,0])):
#     plt.plot(f.grid,np.log10(1+eigenvalues[i,:]),linestyle='--',marker='.',label=str(i))

# # plt.plot(f.grid,eig_CEM[:],linestyle='--',marker='.',label='patched')

# plt.legend()
# plt.xlabel('1D flame domain')
# plt.title('Maximum eigenvalues (from lowest to highest)')
# plt.show()

# if EI_mode == 'debug':
    
#   fig, ax = plt.subplots(2,1) 
#   plt.subplot(2,1,1)
#   for i in range(len(track_species)):
        
#       plt.plot(np.array(tt)*1e6,expl_indices[track_species[i],:],linestyle='--', marker='o',label=EI_keys[track_species[i]])
#   # plt.xlim((7.5e-5,9e-5))
#   plt.xlabel('Residence time')
#   plt.axvline(x=switch_time*1e6, ymin=0., ymax = 1, linewidth=1, color='k')
#   # plt.xticks([7.5e-5, 8e-5, switch_time, 8.5e-5, 9e-5])

#   plt.legend()
#   # plt.show()

#   plt.subplot(2,1,2)
#   plt.plot(np.array(tt)*1e6,eigenvalues[order,:]/1e6,linestyle='--')
#   plt.axvline(x=switch_time*1e6, ymin=0., ymax = 1, linewidth=1, color='k')

#   titlefig = str(order) + '_provaEI' + '.pdf'
#   plt.savefig(titlefig, bbox_inches='tight')

