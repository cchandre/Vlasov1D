########################################################################################################################
##                   Dictionary of parameters: https://github.com/cchandre/Vlasov1D                                   ##
########################################################################################################################

import numpy as xp
from scipy.fft import rfft

kappa = 1.35

Tf = 100
integrator_kinetic = 'position-Verlet'
nsteps = 20
integrator_fluid = 'RK45'
precision_fluid = 1e-12

n_moments = 4
n_casimirs = 3

Lx = 2 * xp.pi
Lv = 8
Nx = 2**12
Nv = 2**12

A = 1e-6
k = 0.5
f_init = lambda x, v: (1 - A * xp.cos(k * x[:, None])) * v[None, :]**2 * xp.exp(-v[None, :]**2 / 2) / xp.sqrt(2 * xp.pi)
## ATTENTION: for the fluid approach to work as implemented in this code, f_init should be with S3=0

output_E_modes = 10

ComputeKinetic = False
ComputeFluid = True
SaveKinetic = False
SaveFluid = True
PlotKinetic = False
PlotFluid = True

darkmode = True

########################################################################################################################
##                                                DO NOT EDIT BELOW                                                   ##
########################################################################################################################
dict = {'kappa': kappa}
dict.update({
        'Tf': Tf,
		'integrator_kinetic': integrator_kinetic,
        'nsteps': nsteps,
		'integrator_fluid': integrator_fluid,
        'precision_fluid': precision_fluid,
		'n_moments': n_moments,
        'n_casimirs': n_casimirs,
        'Lx': Lx,
		'Lv': Lv,
		'Nx': Nx,
		'Nv': Nv,
		'f_init': f_init,
        'output_E_modes': output_E_modes,
		'ComputeKinetic': ComputeKinetic,
		'ComputeFluid': ComputeFluid,
		'SaveKinetic': SaveKinetic,
		'SaveFluid': SaveFluid,
		'PlotKinetic': PlotKinetic,
		'PlotFluid': PlotFluid})
########################################################################################################################
