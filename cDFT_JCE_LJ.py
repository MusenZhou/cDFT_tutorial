# This program is used to illustrate the very basic idea of classic density
# functional theory (cDFT) in the ideal limit , which is without the
# consideration of the term of excess Helmholtz energy
# For more details, please refer to paper (DOI: 10.1021/ed500049m)

# import some module
import math
import numpy as np
import matplotlib.pyplot as plt
import time


# fxn used to calculate external potential in LJ form
def v_ext_cal(epsilon_cal, sigma_cal, r, v_ext):
	for i in range(len(r)):
		v_ext[i] = 4*epsilon_cal*( (sigma_cal/r[i])**12 - (sigma_cal/r[i])**6 )
		if (v_ext[i] >= 100):
			v_ext = 100


# Description: variable rho here represents density.
# dr is the step size of system and r is the length of the system and
# subject to appropriate changes from the user
dr = 0.5e-1
r = np.arange(0.1, 10, dr)
# Picard iteration parameters, which are subject to appropriate changes
# from the user (torr is the tolearance)
mixing_parameter = 0.5
torr = 1e-5



# constant define
# Boltzmann constant, unit: J/Klevin
kb = 1.38064852e-23
# Temperature, unit: Klevin
T = 300
# 1/kBT, unit: 1./J
beta = 1/(kb*T)
# plank constant, unit: J*m
plank_constant = 6.626070040e-34
# Avogadro constant, unit: mol-1
NA = 6.022140857e23
# mass of neon, unit: kg
mass_neon = 20.18e-3/NA
# thermal wavelength, unit: m
const_lambda = math.sqrt((beta*(plank_constant**2)) / (2*math.pi*mass_neon))
# LJ energy parameter, unit: J
epsilon = 1.8e3/NA
# LJ size parameter, unit: Angstrom
sigma = 2.6
# bulk density
rho_bulk = 0.033;
# chemical potential
mu = math.log(rho_bulk*(const_lambda**3)) / beta;

# analytic solution
v_ext_ana = np.empty(r.shape)
v_ext_cal(epsilon, sigma, r, v_ext_ana)
rho_eq_ana = np.empty(r.shape)
for i in range(len(rho_eq_ana)):
	rho_eq_ana[i] = (math.exp(beta*mu)/(const_lambda**3)) * math.exp(-beta*v_ext_ana[i]);


# numerical initial guess
rho_initial = rho_bulk * np.ones(rho_eq_ana.shape)
rho_eq = np.empty(r.shape)
rho_new = np.empty(r.shape)



# plot set up for analytical solution
plt.show()
figure = plt.gca()
figure.set_xlim(-0.5, 10.5)
figure.set_ylim(0, 3)
plt.xlabel('r', fontsize=18)
plt.ylabel(r'$\rho/\rho_{bulk}$', fontsize=18)
line_ana = figure.plot(r, rho_eq_ana/rho_bulk, 'r-', label='analytical solution')
figure.legend()





# iteration till convergence
step=0
while (1):

	# plot for the numerical solution
	if (step==0):
		rho_prev = rho_initial
		line_num, = figure.plot(r, rho_prev/rho_bulk, 'b-', label='numerical solution')
		figure.legend()
	else:
		figure.lines.remove(line_num)
		line_num, = figure.plot(r, rho_new/rho_bulk, 'b-', label='numerical solution')
		figure.legend()



	# calculate the new equilibrium density
	for i in range(len(rho_eq_ana)):
		rho_eq[i] = (math.exp(beta*mu)/(const_lambda**3)) * math.exp(-beta*v_ext_ana[i])
	# picard iteration
	for i in range(len(rho_new)):
		rho_new[i] = (1-mixing_parameter) * rho_prev[i] + mixing_parameter * rho_eq[i]
	# check tolerance
	tot_err = 0
	for i in range(len(rho_new)):
		tot_err = tot_err + (rho_new[i] - rho_prev[i])**2
	tot_err = math.sqrt(tot_err)
	if (tot_err < torr):
		break
	else:
		for i in range(len(rho_new)):
			rho_prev[i] = rho_new[i]

	step = step + 1
	plt.pause(0.5)

plt.show()