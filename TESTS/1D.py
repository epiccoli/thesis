import pyjacob
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time as t
import pdb
from pyjac_functions import *

gas = ct.Solution('Skeletal29_N.cti')
gas.transport_model = 'Mix' # 'Multi' or 'Mix'
# phi = 1.0
# T0 = 700
# P0 = 101325

# ADD option to set phi and get Z to input in setSolutionProperties

# gas.set_equivalence_ratio(phi,'C2H4','O2:1, N2:3.76')
# gas.TP = T0, P0

Zv = [0.03, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
Zv = [0.25, 0.35]

Zv = [0.1]
for Z in Zv:

	gas, phi = setSolutionProperties(gas,Z)
	# gas()
	T = gas.T
	print 'phi >>>>> ', phi


	width = 0.03*7

	## REORDER CANTERA SPECIES AS PYJAC WOULD:
	# specs = gas.species()[:]
	# n2_ind = gas.species_index('N2')
	# gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
	#         species=specs[:n2_ind] + specs[n2_ind + 1:] + [specs[n2_ind]],
	#         reactions=gas.reactions())
	# print gas.species_name(28) # >> should give N2


	# Flame object
	initial_grid = np.linspace(0,width , 7) # coarse first
	f = ct.FreeFlame(gas, grid=initial_grid)
	f.set_refine_criteria(ratio=2, slope=0.05, curve=0.1)


	f.solve(loglevel=0, auto=True)
	print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))



	eigenvalues = flame_eig(f,gas)
	# pdb.set_trace()

	# Plot strings
	titleStr = 'Z = ' + str(Z) + ' for Phi_j = 0.8 (Phi_mix = ' + str(phi) + ')' 
	figName = './figs/Z='+'{:.2f}'.format(Z)+'_LLLONG_PhiJet=0.8.png'

	plt.figure()
	plt.plot(f.grid,np.log10(np.array(eigenvalues)+1))

	plt.title(titleStr)
	plt.xlabel('1D flame domain')
	plt.ylabel(r'$log(\lambda_{max} + 1)$')
	# plt.savefig(figName,format='png')
	plt.show()


# HRR = np.sum(f.standard_enthalpies_RT * f.net_production_rates, 0) * ct.gas_constant * f.T *(-1)

# plt.figure()
# plt.plot(f.grid,HRR)
# plt.show()

