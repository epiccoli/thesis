import sys
sys.path.append('../../functions/')

from CEMA import *
from parallelism import *
import pyjacob

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time as t
import pdb


## Conditional ##

mode = 'selectedEig' 
mode = 'n_eig'

 
def stats(vector,keys):

	max_val = np.amax(vector)
	max_idx = np.argmax(vector)

	print "Max value: ", max_val, " at position ", max_idx, " -> ", keys[max_idx]


# Gas properties
phi = 1.0
P = 101325
T = 1400
fuel_spec = 'C2H4'

# Simulation conditions
npoints = 1000
timestep = 1.5e-7
CEMA_interval = 1 # only divisors of npoints
data_length = int(npoints/CEMA_interval)
N_eig = 7
N_EI = 1
#### options 
first_eigenmode = True
first_ei = True

# Plots
plot_eigenvalue = True

# Create gas object
gas = ct.Solution('USC.cti')


## REORDER CANTERA SPECIES AS PYJAC WOULD:
specs = gas.species()[:]
N2_ind = gas.species_index('N2')
gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
        species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
        reactions=gas.reactions())
# print gas.species_name(28) # >> should give N2



EI_keys = ['']*gas.n_species

EI_keys[0] = 'T'
for i in range(1,gas.n_species):
	EI_keys[i] = gas.species_name(i-1)
	print gas.species_name(i)
print EI_keys	


# prova manuale : H (entry 1)

## SET EQUIVALENCE RATIO TO phi, temperature and pressure
gas.set_equivalence_ratio(phi,fuel_spec,'O2:1, N2:3.76')
gas.TP = T, P
gas_cp = gas

# # Create constant pressure reactor
r = ct.IdealGasConstPressureReactor(gas)

# # Create simulation PSR object
sim = ct.ReactorNet([r])

# # Initialize time and data vectors
time = 0.0

tim = np.zeros(npoints,'d')
temp = np.zeros(npoints,'d')
press = np.zeros(npoints,'d')
enth = np.zeros(npoints,'d')

# with open('xmgrace.txt','w') as file1:
file1 = open('xmgrace.txt','w')


val=[]
tt=[]
val_ei = []

eigenvalues = np.zeros((N_eig,data_length))
expl_indices = np.zeros((gas.n_species,data_length))

track_species = []

count=0

most_aligned_eig = []
most_aligned_ei = []


start = t.time()
for n in range(npoints):
	time += timestep
	sim.advance(time)
	
	tim[n] = time
	temp[n]= r.T

	if n % CEMA_interval == 0:
		D, L, R = solve_eig_gas(gas)
		max_idx = np.argmax(D)
		val.append(D[max_idx])
		tt.append(time)


		########### START TEST PARALLEL ###############
		# # print n
		# if first_eigenmode == True:
		# 	CEM = eigenmode(D,R)
		# 	first_eigenmode = False

		# idx_CEM = test_parallel(R,CEM)
		# # print D[idx_CEM]
		# most_aligned_eig.append(D[idx_CEM])

		########### END TEST PARALLEL #############



		########### START TEST FOLLOW VALUE ############

		# if count == 0:
		# 	eigenvalues[:,count] = D[np.argsort(D)[-N_eig:]]		# matrix with N_eig lines
		
		# elif count == 1:
		# 	next_eig = match_eig(D,eigenvalues[0,count-1])
		# 	eigenvalues[:,count] = next_eig

		# else:
		# 	last_values = eigenvalues[0,count-2:count]
		# 	slope = np.diff(last_values)[0]
		# 	guess_val = last_values[-1] + slope

		# 	# print guess_val
			
		# 	next_eig = match_eig(D,guess_val)

		# 	eigenvalues[:,count] = next_eig
		# # modified to not take into account eigenvalues = 0.0000000
		# eigenvalues[:,count] = highest_val_excl_0(D,N_eig)
		
		########### END TEST FOLLOW VALUE ############

		########### START TEST EI PARALLEL ###########

		# if first_eigenmode == True:
		# 	ei_vec0 = EI(D,L,R,max_idx)
		# 	first_eigenmode = False

		# idx_CEM = test_ei_parallel(L,R,ei_vec0)


		# most_aligned_ei.append(D[idx_CEM])


		########### END TEST EI PARALLEL ###########

		eigenvalues[:,count] = D[np.argsort(D)[-N_eig:]]
		
		# if n % 100 == 0:
			
		# 	print D
		# 	# print n
		# 	# print eigenvalues[:,count]
		# 	# print D[np.argsort(D)[-N_eig:]]
		# 	pdb.set_trace()

		##### START track EI similar to previous #####
		
		# if first_ei == True:
		# 	ei_prev = EI(D,L,R,max_idx)
		# 	first_ei = False

		# new_ei, idx_followed = follow_ei(D,L,R,ei_prev)
		# ei_prev = new_ei
		# pdb.set_trace()
		# print idx_followed
		# print max_idx
		# print count

		# expl_indices[:,count] = new_ei
		# val_ei.append(D[idx_followed])

		##### END track EI similar to previous #####

		## MANUALLY CHANGE WHICH EIG on which to base EI calc: -6 is eig 1
		#   6 5 4 3 2 1 0
		# -[1 2 3 4 5 6 7]
		if tt[count] < 9.76e-5:
			max_idx = np.argsort(D)[-1]
		else:
			max_idx = np.argsort(D)[-6]

		# max_idx = np.argsort(D)[-5]

		expl_indices[:,count] = EI(D,L,R,max_idx)
		main_EI = np.argsort(expl_indices[:,count])
		# manual identification of important species indices
		# print main_EI
		# pdb.set_trace()



		track_species = np.union1d(main_EI,track_species)

		# stats(expl_indices,EI_keys)

		# track species	

		
		count += 1

		
		
	# enth[n]= r.thermo.enthalpy_mass
	press[n]= r.thermo.P
	file1.write('%10.6f %10.4f %10.6f \n' % (sim.time,r.T,r.thermo.P))



file1.close()

end = t.time()

print end-start, ' seconds'	

dT=np.diff(temp)
dTdt = dT/np.diff(tim)

selected_species = []

# Time series plot of temperature, maximum eigenvalue, temperature gradient
plt.figure(figsize=(10,5))
plt.subplot(3,1,1)
plt.plot(tim,temp)
plt.title('Temp')

plt.subplot(3,1,2)
plt.plot(tt,np.array(val)/1e6)
plt.title('Eig')

plt.subplot(3,1,3)
plt.plot(np.arange(0,len(dT)),dT)
plt.title('Temperature gradient')

plt.show()

# 






track_entries = ['T', 'H', 'O', 'OH', 'HO2', 'O2']
idx_entries = [0, 4, 5, 6, 7, 2]

# idx_entries = [2, 7, 13, 5]
# track_entries = ['H', 'HO2', 'CH3', 'OH']

#### EI plot ####
plt.figure()
for i in range(len(idx_entries)):
	plt.plot(tt,expl_indices[idx_entries[i],:],linestyle='--', marker='o',label=track_entries[i])
plt.axvline(x=9.76e-5, ymin=0., ymax = 1, linewidth=1, color='k')
plt.legend()
plt.show()




#### EIG PLOT ####

if mode == 'n_eig':

	for i in range(N_eig):
		# plt.subplot(N_eig,1,i+1) # comment to overlay instead of subplot
		plt.plot(tt,eigenvalues[i,:]/1e6,'--',label=str(i))
		plt.legend() 
		# pdb.set_trace()


if mode == 'selectedEig':	
	plt.figure()
	# PLOT ONLY FOR COMPARISON
	# i=6
	# plt.plot(tt,eigenvalues[i,:]/1e6,'.',label='CEM high')
	# i=0
	# plt.plot(tt,eigenvalues[i,:]/1e6,'.',label='first neg')
	# i=1
	# plt.plot(tt,eigenvalues[i,:]/1e6,'.',label='CEM')



	i=5
	plt.plot(tt,eigenvalues[i,:]/1e6,'--',label=str(i))
	i=2
	plt.plot(tt,eigenvalues[i,:]/1e6,'--',label=str(i))
	i=4
	plt.plot(tt,eigenvalues[i,:]/1e6,'--',label=str(i))
	i=3
	plt.plot(tt,eigenvalues[i,:]/1e6,'--',label=str(i))



# post_treatment of most_aligned_eig and tt
# discard values of most_aligned_eig/1e6 < 0.5

most_aligned_eig = np.array(most_aligned_eig)
# pdb.set_trace()
tt = np.delete(tt,np.where(most_aligned_eig < -8e5))
most_aligned_eig = np.delete(most_aligned_eig,np.where(most_aligned_eig < -8e5))

# plt.plot(tt,np.array(most_aligned_eig)/1e6,'.',label='most_aligned_eig')
# plt.plot(tt,np.array(most_aligned_ei)/1e6,'.',label='most_aligned_ei')
plt.legend()
plt.show()


print 'maximum heat release rate at time ', tim[np.argmax(np.diff(temp))]
print max(val), 'is max eigenvalue at time ', tt[val.index(max(val))]


