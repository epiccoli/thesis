import sys
sys.path.append('../functions/')
sys.path.append('../ethylene/Skeletal29_N/')
from CEMA import *
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

# Simulation conditions
npoints = 1000
timestep = 1.5e-7

# Create gas object
gas = ct.Solution('Skeletal29_N.cti')

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

track_entries = ['H', 'OH', 'CH', 'HO2', 'HCO', 'H2O2', 'CH3', 'O2', 'C']
idx_entries = [1, 4, 19, 8, 20, 6, 10, 5, 18]
# prova manuale : H (entry 1)

## SET EQUIVALENCE RATIO TO phi, temperature and pressure
gas.set_equivalence_ratio(phi,'C2H4','O2:1, N2:3.76')
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

CEMA_interval = 1 # only divisors of npoints
data_length = int(npoints/CEMA_interval)
N_eig = 2
val=[]
tt=[]
ei_sp=[]

eigenvalues = np.zeros((N_eig,data_length))

count=0

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

		eigenvalues[:,count] = D[np.argsort(D)[-N_eig:]]		# matrix with N_eig lines
		print D[np.argsort(D)[-N_eig:]]

		expl_indices = EI(D,L,R,max_idx)
		ei_sp.append(expl_indices[18])
		stats(expl_indices,EI_keys)
		
		count += 1

		
		
	# enth[n]= r.thermo.enthalpy_mass
	press[n]= r.thermo.P
	file1.write('%10.6f %10.4f %10.6f \n' % (sim.time,r.T,r.thermo.P))

# TAI = findAIT(gas,timestep,npoints)
# gas = gas_cp

end = t.time()
# 
print len(val)
print end-start, ' seconds'	



dT=np.diff(temp)
dTdt = dT/np.diff(tim)


file1.close()

plt.figure(figsize=(10,5))
plt.subplot(4,1,1)
plt.plot(tim,temp)
plt.title('Temp')

plt.subplot(4,1,2)
plt.plot(tt,np.array(val)/1e6)
plt.title('Eig')

plt.subplot(4,1,3)
plt.plot(np.arange(0,len(dT)),dT)
plt.title('Temperature gradient')

plt.subplot(4,1,4)
plt.plot(np.arange(0,len(ei_sp)),ei_sp)

plt.show()

plt.figure()

for i in range(N_eig):
	plt.plot(tt,np.log10(eigenvalues[i,:]+1))

plt.show()


print 'maximum heat release rate at time ', tim[np.argmax(np.diff(temp))]
print max(val), 'is max eigenvalue at time ', tt[val.index(max(val))]


