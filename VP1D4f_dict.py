########################################################################################################################
##                   Dictionary of parameters: https://github.com/cchandre/Vlasov1D                                   ##
########################################################################################################################

import numpy as xp

kappa = 1.31

Tf = 48
integrator_kinetic = 'position-Verlet'
nsteps = 32
integrator_fluid = 'RK45'
precision_fluid = 1e-12

n_moments = 4
n_casimirs = 3

Lx = 2 * xp.pi
Lv = 4
Nx = 2**11
Nv = 2**11

A = 1e-6
k = 0.5
f_init = lambda x, v: (1 - A * xp.cos(k * x[:, None])) * v[None, :]**2 * xp.exp(-v[None, :]**2 / 2) / xp.sqrt(2 * xp.pi)

output_E_modes = 8

ComputeKinetic = True
PlotKinetic = True
SaveKinetic = False

ComputeFluid = True
PlotFluid = True
SaveFluid = False

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
        'PlotKinetic': PlotKinetic,
		'SaveKinetic': SaveKinetic,
        'ComputeFluid': ComputeFluid,
        'PlotFluid': PlotFluid,
		'SaveFluid': SaveFluid})
########################################################################################################################
