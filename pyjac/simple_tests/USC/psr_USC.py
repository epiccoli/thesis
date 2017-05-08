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
N_eig = 1
N_EI = 1
first_eigenmode = True

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

eigenvalues = np.zeros((N_eig,data_length))
expl_indices = np.zeros((gas.n_species,data_length))

track_species = []

count=0

max_new = []


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
		print n
		if first_eigenmode == True:
			CEM = eigenmode(D,L,R)
			first_eigenmode = False
			print "numbers of passages ici"

		idx_CEM = test_parallel(L,R,CEM)
		# print D[idx_CEM]
		max_new.append(D[idx_CEM])

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
		# modified to not take into account eigenvalues = 0.0000000
		# eigenvalues[:,count] = highest_val_excl_0(D,N_eig)
		
		########### END TEST FOLLOW VALUE ############


		eigenvalues[:,count] = D[np.argsort(D)[-N_eig:]]
		
		# if n % 100 == 0:
			
		# 	print D
		# 	# print n
		# 	# print eigenvalues[:,count]
		# 	# print D[np.argsort(D)[-N_eig:]]
		# 	pdb.set_trace()

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

idx_entries = [2, 7, 13, 5]
track_entries = ['H', 'HO2', 'CH3', 'OH']

#### EI plot ####
# plt.figure()
# for i in range(len(idx_entries)):
# 	plt.plot(tt,expl_indices[idx_entries[i],:],label=track_entries[i])
# plt.legend()
# plt.show()




#### EIG PLOT ####
if plot_eigenvalue == True:
	plt.figure()

	for i in range(N_eig):
		plt.plot(tt,eigenvalues[i,:]/1e6)

	# post_treatment of max_new and tt
	# discard values of max_new/1e6 < 0.5

	max_new = np.array(max_new)
	pdb.set_trace()
	tt = np.delete(tt,np.where(max_new < -8e5))
	max_new = np.delete(max_new,np.where(max_new < -8e5))

	plt.plot(tt,np.array(max_new)/1e6,'.',label='max_new')
	plt.legend()
	plt.show()


print 'maximum heat release rate at time ', tim[np.argmax(np.diff(temp))]
print max(val), 'is max eigenvalue at time ', tt[val.index(max(val))]


