import pyjacob
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time as t
import pdb
from pyjac_functions import *

gas = ct.Solution('Skeletal29_N.cti')

phi = 1.0
P = 101325
Temperatures = [1400]

npoints = 1000

## REORDER CANTERA SPECIES AS PYJAC WOULD:
specs = gas.species()[:]
n2_ind = gas.species_index('N2')
gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
        species=specs[:n2_ind] + specs[n2_ind + 1:] + [specs[n2_ind]],
        reactions=gas.reactions())
# print gas.species_name(28) # >> should give N2


for T in Temperatures:

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

	# # set timestep # adaptive time step with gradient of temperature?
	timestep = 1e-4*1.5e-3

	# with open('xmgrace.txt','w') as file1:
	file1 = open('xmgrace.txt','w')

	val=[]
	tt=[]


	start = t.time()
	for n in range(npoints):
		time += timestep
		sim.advance(time)
		
		tim[n] = time
		temp[n]= r.T

		if n % 1 == 0:
			val.append(find_max_eigenvalue(gas))
			tt.append(time)
			
			
		# enth[n]= r.thermo.enthalpy_mass
		press[n]= r.thermo.P
		file1.write('%10.6f %10.4f %10.6f \n' % (sim.time,r.T,r.thermo.P))
	
	# TAI = findAIT(gas,timestep,npoints)
	# gas = gas_cp
	# eigv, tt = PSR_eig(gas,timestep,npoints,TAI)
	end = t.time()
	# pdb.set_trace()

	print end-start, ' seconds'	



	dT=np.diff(temp)
	dTdt = dT/np.diff(tim)


	file1.close()

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


	print 'maximum heat release rate at time ', tim[np.argmax(np.diff(temp))]
	print max(val), 'is max eigenvalue at time ', tt[val.index(max(val))]




