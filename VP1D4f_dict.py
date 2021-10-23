########################################################################################################################
##                   Dictionary of parameters: https://github.com/cchandre/Vlasov1D                                   ##
########################################################################################################################

import numpy as xp

qe = -1

Tf = 5
TimeStep = 1e-2
integrator_kinetic = 'position-Verlet'
integrator_fluid = 'RK45'
precision_fluid = 1e-11

n_moments = 4
frames = 100

Lx = 2 * xp.pi
Lv = 6
Nx = 2**10
Nv = 2**10

tail_indx = xp.index_exp[Nx//4:3*Nx//4+1]

epsilon = 1e-1
f_init = lambda x, v: (1 - epsilon * xp.cos(xp.pi * x[:, None] / Lx)) * v[None, :]**2 * xp.exp(-v[None, :]**2 / 2)
#f_init = lambda x, v: (1 - epsilon * xp.cos(xp.pi * x[:, None] / Lx)) * xp.exp(-v[None, :]**2 / 2)
## ATTENTION: for the fluid approach to work, f_init should be with S3=0
kappa = 2 * (2*Lx)**(2/3)

ComputeKinetic = False
ComputeFluid = True
SaveKinetic = False
SaveFluid = False
PlotKinetic = False
PlotFluid = True

########################################################################################################################
##                                                DO NOT EDIT BELOW                                                   ##
########################################################################################################################
dict = {'qe': qe}
dict.update({
        'kappa': kappa,
        'Tf': Tf,
		'TimeStep': TimeStep,
		'integrator_kinetic': integrator_kinetic,
		'integrator_fluid': integrator_fluid,
        'precision_fluid': precision_fluid,
		'n_moments': n_moments,
        'frames': frames,
        'Lx': Lx,
		'Lv': Lv,
		'Nx': Nx,
		'Nv': Nv,
        'tail_indx': tail_indx,
		'f_init': f_init,
		'ComputeKinetic': ComputeKinetic,
		'ComputeFluid': ComputeFluid,
		'SaveKinetic': SaveKinetic,
		'SaveFluid': SaveFluid,
		'PlotKinetic': PlotKinetic,
		'PlotFluid': PlotFluid})
########################################################################################################################
