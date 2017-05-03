# Finds useful properties for thickened flame model avbp
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import pdb

def main():

	### GAS SETTINGS ###
	phi = 1.2
	T = 300		# Kelvin
	P = 101325 	# Pascal

	### FLAME SETTINGS ###
	width = 0.03	# meter
	initial_grid = np.linspace(0,width , 7) # coarse first

	### Gas object ###
	gas = ct.Solution('Skeletal29_N.cti')
	gas.transport_model = 'Mix' # 'Multi' or 'Mix'
	gas.set_equivalence_ratio(phi,'C2H4','O2:1, N2:3.76')
	gas.TP = T, P

	### Flame object ###
	f = ct.FreeFlame(gas, grid=initial_grid)
	f.set_refine_criteria(ratio=2, slope=0.05, curve=0.1)

	### SOLVE FLAME ###
	f.solve(loglevel=0, auto=True)


	### EXTRACT DATA ###
	temperatures = f.T

	Tu = temperatures[0]
	Tad = temperatures[-1]

	dT = np.diff(temperatures)
	dx = np.diff(f.grid)

	dTdx = dT/dx

	# flame thickness based on max temperature gradient definition
	flame_thick = (Tad-Tu)/np.amax(dTdx)

	# maximum dilatation
	u_vect = f.u
	du = np.diff(u_vect)
	dudx = du/dx

	max_dilatation = np.amax(dudx)


	# maximum heat release rate
	HRR = np.sum(f.standard_enthalpies_RT * f.net_production_rates, 0) * ct.gas_constant * f.T *(-1)
	max_HRR = np.amax(HRR)

	# check that idx are similar (must be around flame)
	dil_idx = np.argmax(dudx)
	tem_idx = np.argmax(dTdx)

	

	print('\nadiabatic flame temperature = \t \t \t {:4f} kelvin'.format(Tad))
	print('mixture-averaged flamespeed = \t \t \t {:7f} m/s'.format(f.u[0]))
	print('flame thickness based on max Temp gradient = \t {:5f} mm'.format(flame_thick*1e3))
	print('maximum heat release rate = \t \t \t {:1.6e} unit'.format(max_HRR))
	




if __name__ == '__main__':
	main()