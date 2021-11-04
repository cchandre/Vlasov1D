###################################################################################################
##              Dictionary of parameters: https://github.com/cchandre/Vlasov1D                   ##
###################################################################################################

import numpy as xp

kappa = 7

Tf = 20
integrator_kinetic = 'position-Verlet'
nsteps = 30
integrator_fluid = 'RK45'
precision = 1e-11

n_moments = 4
n_casimirs = 3

Lx = 2 * xp.pi
Lv = 4
Nx = 2**10
Nv = 2**10

A = 1e-2
k = 0.5
f_init = lambda x, v: (1 - A * xp.cos(k * x)) * v**2 * xp.exp(-v**2 / 2) / xp.sqrt(2 * xp.pi)

output_var = 'E'
output_modes = 'real'

## 'Compute', 'Plot' and/or 'Save'
Kinetic = []
Fluid = ['Compute', 'Plot']

darkmode = True

###################################################################################################
##                               DO NOT EDIT BELOW                                               ##
###################################################################################################
dict = {'kappa': kappa}
dict.update({
        'Tf': Tf,
		'integrator_kinetic': integrator_kinetic,
        'nsteps': nsteps,
		'integrator_fluid': integrator_fluid,
        'precision': precision,
		'n_moments': n_moments,
        'n_casimirs': n_casimirs,
        'Lx': Lx,
		'Lv': Lv,
		'Nx': Nx,
		'Nv': Nv,
		'f_init': f_init,
        'output_var': output_var,
        'output_modes': output_modes,
		'Kinetic': Kinetic,
        'Fluid': Fluid,
        'darkmode': darkmode})
###################################################################################################
