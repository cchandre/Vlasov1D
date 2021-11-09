###################################################################################################
##              Dictionary of parameters: https://github.com/cchandre/Vlasov1D                   ##
###################################################################################################

import numpy as xp
from scipy.special import erf

<<<<<<< Updated upstream
kappa = 1

Tf = 200
=======
kappa = 1.001

Tf = 100
>>>>>>> Stashed changes
integrator_kinetic = 'position-Verlet'
nsteps = 16
integrator_fluid = 'DOP853'
precision = 1e-11

n_moments = 4
n_casimirs = 3

Lx = 2 * xp.pi
Lv = 4
Nx = 2**11
Nv = 2**11

<<<<<<< Updated upstream
A = 1e-4
k = 0.5
alpha = 0.1
# f_init = lambda x, v: (1 - A * xp.cos(k * x)) * v**2 * xp.exp(-v**2 / 2) / xp.sqrt(2 * xp.pi)
f_init = lambda x, v: (1 - A * xp.cos(k * x)) * xp.exp(-v**2 / 2) / xp.sqrt(2 * xp.pi) * (1 + erf(alpha * v / xp.sqrt(2)))
=======
A = 1e-3
k = 0.5
# f_init = lambda x, v: (1 - A * xp.cos(k * x)) * v**2 * xp.exp(-v**2 / 2) / xp.sqrt(2 * xp.pi)
f_init = lambda x, v: (1 - A * xp.cos(k * x)) * xp.exp(-v**2 / 2) / xp.sqrt(2 * xp.pi)
>>>>>>> Stashed changes

output_var = 'E'
output_modes = 'real'

## 'Compute', 'Plot' and/or 'Save'
Kinetic = []
<<<<<<< Updated upstream
Fluid = ['Compute', 'Plot', 'Save']
=======
Fluid = ['Compute', 'Plot']
>>>>>>> Stashed changes

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
