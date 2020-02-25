# This program is used to illustrate the very basic idea of classic density
# functional theory (cDFT) in the ideal limit with a python version code.

# import some module
import math
import numpy as np
import matplotlib.pyplot as plt
import time


# fxn to calculate lj potential
#calculate external potential
def v_ext_cal(epsilon_cal, sigma_cal, r, v_ext):
	for i in range(len(r)):
		v_ext[i] = 4*epsilon_cal*( (sigma_cal/r[i])**12 - (sigma_cal/r[i])**6 )
		if (v_ext[i] >= 100):
			v_ext = 100


# Description: rho here is density


# constant define
kb = 1.38064852e-23                                                                 # unit: J/Klevin
T = 300                                                                             # unit: Klevin
beta = 1/(kb*T)                                                                     # unit: 1./J
plank_constant = 6.626070040e-34                                                    # unit: J*m
NA = 6.022140857e23                                                                 # unit: mol-1
mass_neon = 20.18e-3/NA                                                             # unit: kg
const_lambda = math.sqrt((beta*(plank_constant**2)) / (2*math.pi*mass_neon))        # unit: m
epsilon = 1.8e3/NA                                                                 # unit: J
sigma = 2.6                                                                        # unit: Angstrom
dr = 0.5e-1
r = np.arange(0.1, 10, dr)
delta_coeff = 0.1;
rho_bulk = 0.033;

mu = math.log(rho_bulk*(const_lambda**3)) / beta;

#print(r.shape)

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

# Picard iteration parameters
mixing_parameter = 0.5
torr = 1e-5

#plt.scatter(r, rho_eq_ana)
plt.show()
figure = plt.gca()

#plt.ion()
figure.set_xlim(-0.5, 10.5)
figure.set_ylim(0, 3)
line_ana = figure.plot(r, rho_eq_ana/rho_bulk, 'r-')





# iteration till convergence
step=0
while (1):

	# visualization of the code
	if (step==0):
		rho_prev = rho_initial
		line_num, = figure.plot(r, rho_prev/rho_bulk, 'b-')
	else:
		figure.lines.remove(line_num)
		line_num, = figure.plot(r, rho_new/rho_bulk, 'b-')



	# calculate the new equilibrium density
	for i in range(len(rho_eq_ana)):
		rho_eq[i] = (math.exp(beta*mu)/(const_lambda**3)) * math.exp(-beta*v_ext_ana[i])

	# picard iteration
	for i in range(len(rho_new)):
		rho_new[i] = (1-mixing_parameter) * rho_prev[i] + mixing_parameter * rho_eq[i]

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