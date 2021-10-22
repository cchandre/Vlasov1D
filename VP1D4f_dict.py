########################################################################################################################
##                   Dictionary of parameters: https://github.com/cchandre/Vlasov1D                                   ##
########################################################################################################################

import numpy as xp

qe = -1

Tf = 45
TimeStep = 1e-3
integrator_kinetic = 'position-Verlet'
integrator_fluid = 'RK45'

n_moments = 4
frames = 100

Lx = 2 * xp.pi
Lp = 6
Nx = 2**10
Np = 2**10

epsilon = 1e-1
f_init = lambda x, p: p[None, :]**2 * (1 - epsilon * xp.cos(x[:, None] / 2)) * xp.exp(-p[None, :]**2 / 2)
## ATTENTION: for the fluid approach to work, f_init should be with S3=0
kappa = 3

ComputeKinetic = True
ComputeFluid = False
SaveKinetic = False
SaveFluid = False
PlotKinetic = True
PlotFluid = False

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
		'n_moments': n_moments,
        'frames': frames,
        'Lx': Lx,
		'Lp': Lp,
		'Nx': Nx,
		'Np': Np,
		'f_init': f_init,
		'ComputeKinetic': ComputeKinetic,
		'ComputeFluid': ComputeFluid,
		'SaveKinetic': SaveKinetic,
		'SaveFluid': SaveFluid,
		'PlotKinetic': PlotKinetic,
		'PlotFluid': PlotFluid})
########################################################################################################################
