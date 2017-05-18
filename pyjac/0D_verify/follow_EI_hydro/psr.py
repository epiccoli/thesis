import sys
sys.path.append('../../functions/')
sys.path.append('../../hydrogen')

from CEMA import *
import pyjacob

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time as t
import pdb


## Plot toggle
plot_EI = True
plot_EI = False
plot_eigenvalue = False
# mode = 'selectedEig' 
mode = 'n_eig'


## SELECT MODE FOR EI ##
# EI_mode = 'debug'
# EI_mode = 'patch'
EI_mode = 'follow'


# Gas properties
phi = 1.0
P = 101325
T = 1100
fuel_spec = 'H2'

# Simulation conditions
npoints = 3000
timestep = 0.5e-7
CEMA_interval = 1 # only divisors of npoints
data_length = int(npoints/CEMA_interval)
N_eig = 8
N_EI = 3
#### options 
# first_eigenmode = True
# first_ei = True
first_ei_calculation = True


# Create gas object
gas = ct.Solution('Li_2003.cti')
excluded_species = 'AR'

## REORDER CANTERA SPECIES AS PYJAC WOULD:
specs = gas.species()[:]
last_spec_idx = gas.species_index(excluded_species)
gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
        species=specs[:last_spec_idx] + specs[last_spec_idx + 1:] + [specs[last_spec_idx]],
        reactions=gas.reactions())
# print gas.species_name(28) # >> should give N2


### CREATE EI_keys -> vector with strings 
### representing the names of the components of the Phi vector for the jacobian
EI_keys = []
EI_keys.append('T')
for i in range(gas.n_species):
	if gas.species_name(i) != excluded_species:				# assume that the last species by Pyjac was N2 -> check if it is the case for the scheme used to write pyjacob
		EI_keys.append(gas.species_name(i))



## SET EQUIVALENCE RATIO TO phi, temperature and pressure
gas.set_equivalence_ratio(phi,fuel_spec,'O2:1, N2:3.76')
gas.TP = T, P


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
lambda_patched1 = []
lambda_patched2 = []

eigenvalues = np.zeros((N_eig,data_length))
expl_indices = np.zeros((gas.n_species,data_length))
lambda_selected = np.zeros(data_length)

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

		eigenvalues[:,count] = D[np.argsort(D)[-N_eig:]]

		
		## MANUALLY CHANGE WHICH EIG on which to base EI calc: -6 is eig 1
		#    6 5 4 3 2 1 0			i
		# -[ 1 2 3 4 5 6 7 ] 
		switch_time1 = 9.377e-5
		switch_time2 = 9.392e-5

		if EI_mode == 'patch':
			if tt[count] < switch_time1:
				max_idx = np.argsort(D)[-1]
			elif tt[count] > switch_time1 and tt[count] < switch_time2:
				max_idx = np.argsort(D)[-4]
			elif tt[count] > switch_time2:
				max_idx = np.argsort(D)[-5]

			expl_indices[:,count] = EI(D,L,R,max_idx)
			lambda_selected[count] = D[max_idx]

			main_EI = np.argsort(expl_indices[:,count])[-N_EI:]

			# manual identification of important species indices
			# print main_EI
			# pdb.set_trace()

			track_species = np.union1d(main_EI,track_species)
			
			count += 1

		if EI_mode == 'debug':
		# WHEN LOOKING AT WHICH EI IS WHICH (DEBUG EI FOR PLOTS)
			order=-5
			max_idx = np.argsort(D)[order]

			expl_indices[:,count] = EI(D,L,R,max_idx)
			lambda_selected[count] = D[max_idx]

			main_EI = np.argsort(expl_indices[:,count])[-N_EI:]

			# manual identification of important species indices
			# print main_EI
			# pdb.set_trace()

			track_species = np.union1d(main_EI,track_species)
			
			count += 1

		if EI_mode == 'follow':

			# for now simple case of positive eig at the beginning of time (explosive mixture)
			if first_ei_calculation == True:

				max_idx = np.argmax(D)
				ei_previous = EI(D,L,R,max_idx) 	# false "previous"
				first_ei_calculation = False

			# find maximum alignment between other EI and this ei_previous
			alignment = np.zeros(len(D))
			for idx in range(len(D)):
				alignment[idx] = parallelism(EI(D,L,R,idx), ei_previous)
			# pdb.set_trace()

			best_fit_idx = np.argmax(alignment)
			ei_current = EI(D,L,R,best_fit_idx)

			ei_previous = ei_current

			expl_indices[:,count] = ei_current
			lambda_selected[count] = D[best_fit_idx]
		
			main_EI = np.argsort(ei_current)[-N_EI:]

			# # manual identification of important species indices
			# # print main_EI
			# # pdb.set_trace()

			track_species = np.union1d(main_EI,track_species)
		
			count += 1

	press[n]= r.thermo.P
	file1.write('%10.6f %10.4f %10.6f \n' % (sim.time,r.T,r.thermo.P))


file1.close()

end = t.time()

print end-start, ' seconds'	

dT=np.diff(temp)
dTdt = dT/np.diff(tim)

selected_species = []

#### CHECK if timestep and npoints are sufficient to have ignition resolved
# Time series plot of temperature, maximum eigenvalue, temperature gradient
# plt.figure(figsize=(10,5))
# plt.subplot(3,1,1)
# plt.plot(tim,temp)
# plt.title('Temp')

# plt.subplot(3,1,2)
# plt.plot(tt,np.array(val)/1e6)
# plt.title('Eig')

# plt.subplot(3,1,3)
# plt.plot(np.arange(0,len(dT)),dT)
# plt.title('Temperature gradient')

# plt.show()


#### EI plot with tracked_species ####

track_species=map(int,track_species)

if EI_mode == 'debug':

	fig, ax = plt.subplots(figsize=(16, 10))

	plt.subplot(2,1,1)
	
	for i in range(len(track_species)):
		plt.plot(np.array(tt)*1e6,expl_indices[track_species[i],:],linestyle='--', marker='o',label=EI_keys[track_species[i]])
	# plt.xlim((7.5e-5,9e-5))
	plt.xlabel('Residence time')
	plt.axvline(x=switch_time1*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	plt.axvline(x=switch_time2*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	# plt.xticks([7.5e-5, 8e-5, switch_time, 8.5e-5, 9e-5])

	plt.legend()
	# plt.show()

	plt.subplot(2,1,2)
	plt.plot(np.array(tt)*1e6,eigenvalues[order,:]/1e6,linestyle='--', marker='.')
	plt.axvline(x=switch_time1*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	plt.axvline(x=switch_time2*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')

	titlefig = str(order) + '_provaEI' + '.pdf'
	plt.savefig(titlefig, bbox_inches='tight')

if EI_mode == 'patch':

	fig, ax = plt.subplots(figsize=(16,10))

	plt.subplot(2,1,1)
	plt.axvline(x=switch_time1*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	plt.axvline(x=switch_time2*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	for i in range(len(track_species)):
		plt.plot(np.array(tt)*1e6,expl_indices[track_species[i],:],linestyle='--', marker='o',label=EI_keys[track_species[i]])
	plt.legend()
	# plt.xlim((80,100))
	
	plt.subplot(2,1,2)
	plt.axvline(x=switch_time1*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	plt.axvline(x=switch_time2*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
	plt.plot(np.array(tt)*1e6,lambda_selected/1e6,linestyle='--', marker='o',label='patched eig')

	# plt.xlim((80,100))
	plt.xlabel(r'Residence time [$\mu$ s]')
	# plt.xticks([7.5e-5, 8e-5, switch_time, 8.5e-5, 9e-5])


	plt.savefig('patched-4.pdf',bbox_inches='tight')

if EI_mode == 'follow':

	fig= plt.subplots(figsize=(16,10))
	ax1= plt.subplot(2,1,1)
	ax1.plot(np.array(tt)*1e6,lambda_selected/1e6,linestyle='--',marker='.')
	
	ax2 = plt.subplot(2,1,2,sharex=ax1)
	for i in range(len(track_species)):
		ax2.plot(np.array(tt)*1e6,expl_indices[track_species[i],:],linestyle='--', marker='o',label=EI_keys[track_species[i]])
	plt.legend()
	plt.show()



#### EIG PLOT ####

if mode == 'n_eig':
	legend_entry = ['8th','7th','6th','5th','4th','3rd','2nd','1st']
	plt.figure(figsize=(16,10))
	for i in range(N_eig):

		# plt.subplot(N_eig,1,i+1) # comment to overlay instead of subplot
		plt.plot(np.array(tt)*1e6,eigenvalues[i,:]/1e6,linestyle='--',marker='.',label=legend_entry[i])
		plt.axvline(x=switch_time1*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
		plt.axvline(x=switch_time2*1e6, ymin=0., ymax = 1, linewidth=1,linestyle='--', color='k')
		plt.legend() 
	
	plt.show()


print 'maximum heat release rate at time ', tim[np.argmax(np.diff(temp))]
print max(val), 'is max eigenvalue at time ', tt[val.index(max(val))]

plt.show()