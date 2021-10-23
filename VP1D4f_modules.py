#
# BSD 2-Clause License
#
# Copyright (c) 2021, Cristel Chandre
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as xp
from scipy.integrate import solve_ivp
from tqdm import trange
from scipy.io import savemat
import time
from datetime import date
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({
	'text.usetex': True,
	'font.family': 'serif',
	'font.sans-serif': ['Palatino'],
	'font.size': 20,
	'axes.labelsize': 26,
	'figure.figsize': [8, 8],
	'image.cmap': 'Blues'})

def integrate(case):
	f = case.f.copy()
	moments = case.compute_moments(f, max(5, case.n_moments))
	fs = moments[:4, :].reshape(4*case.Nx).copy()
	###########################################
	## ATTENTION: works only if f has S3=0
	if case.ComputeFluid and (xp.max(xp.abs(fs[3*case.Nx:4*case.Nx])) >= case.precision_fluid):
		print('\033[33m        Warning: S3 is not zero \033[00m')
	else:
		fs[3*case.Nx:4*case.Nx] = 0
	moments[3, :] = 0
	fs[2*case.Nx:3*case.Nx] = fs[2*case.Nx:3*case.Nx]**(1/3)
	if case.ComputeFluid and (case.kappa <= xp.min(fs[2*case.Nx:3*case.Nx])):
		print('\033[33m        Warning: the value of kappa may be too small  (kappa < S2^(1/3)) \033[00m')
	###########################################
	Ef = case.E_fluid(fs[:case.Nx])
	Ek = case.E_kinetic(f)
	H0k = case.kinetic_energy(f, Ek)
	H0f = case.fluid_energy(fs)
	plt.ion()
	if case.PlotFluid:
		fig_f = plt.figure(figsize=(8, 13))
		fig_f.canvas.manager.set_window_title('Fluid simulation')
		axs_f = fig_f.subplots(case.n_moments, 1, sharex=True)
		line_Ef, = axs_f[-1].plot(case.x, Ef, 'r', label=r'$E$')
		line_rho_f, = axs_f[0].plot(case.x, moments[0, :], 'k', label=r'$\rho$')
		line_Sf = []
		for m in range(2, min(case.n_moments, 6)):
			line_temp, = axs_f[m-1].plot(case.x, moments[m, :], linewidth=2, label=r'$S_{{{}}}$'.format(m))
			line_Sf.append(line_temp)
		for ax in axs_f:
			ax.set_xlim((-case.Lx, case.Lx))
			ax.legend(loc='upper right')
		axs_f[-1].set_xlabel('$x$')
		plt.draw()
		plt.pause(1e-4)
	if case.PlotKinetic:
		fig = plt.figure(figsize=(7, 5))
		fig.canvas.manager.set_window_title(r'Distribution function f(x,v,t)')
		im = plt.imshow(f.transpose(), interpolation='gaussian', origin='lower', aspect='auto', extent=(-case.Lx, case.Lx, -case.Lv, case.Lv), vmin=xp.min(f), vmax=xp.max(f))
		plt.gca().set_ylabel('$v$')
		plt.gca().set_xlabel('$x$')
		plt.colorbar()
		fig_k = plt.figure(figsize=(8, 13))
		fig_k.canvas.manager.set_window_title('Kinetic simulation')
		axs_k = fig_k.subplots(case.n_moments, 1, sharex=True)
		line_Ek, = axs_k[-1].plot(case.x, Ek, 'r', linewidth=2, label=r'$E$')
		line_rho_k, = axs_k[0].plot(case.x, moments[0, :], 'k', linewidth=2, label=r'$\rho$')
		line_Sk = []
		for m in range(2, case.n_moments):
			line_temp, = axs_k[m-1].plot(case.x, moments[m, :], linewidth=2, label=r'$S_{{{}}}$'.format(m))
			line_Sk.append(line_temp)
		for ax in axs_k:
			ax.set_xlim((-case.Lx, case.Lx))
			ax.legend(loc='upper right')
		axs_k[-1].set_xlabel('$x$')
		plt.draw()
		plt.pause(1e-4)
	nsteps = xp.int32(xp.ceil(case.Tf / case.TimeStep))
	n_int = nsteps // case.frames
	dt = case.TimeStep * (nsteps // case.frames)
	start = time.time()
	for _ in trange(case.frames):
		if case.ComputeFluid:
			sol = solve_ivp(case.eqn_4f, (0, dt), fs, t_eval=(0, dt), method=case.integrator_fluid, max_step=case.TimeStep, atol=case.precision_fluid, rtol=case.precision_fluid)
			fs = sol.y[:, -1]
			if xp.min(fs[2*case.Nx:3*case.Nx]) <= 1e-14:
				print('\033[90m        Error in the fluid simulation (S2<0) \033[00m')
				break
			if case.PlotFluid:
				rho, u, G2, G3 = xp.split(fs, 4)
				line_rho_f.set_ydata(rho)
				line_Ef.set_ydata(case.E_fluid(rho))
				S = case.compute_S(G2, G3)[:min(case.n_moments-2, 4)]
				for m in range(min(case.n_moments-2, 4)):
					line_Sf[m].set_ydata(S[m])
				for ax in axs_f:
					ax.relim()
					ax.autoscale()
					ax.set_xlim((-case.Lx, case.Lx))
				plt.draw()
				plt.pause(1e-4)
		if case.ComputeKinetic:
			for t in range(n_int):
				for coeff, type in zip(case.integr_coeff, case.integr_type):
					if type == 1:
						f, Ek = case.L1(f, Ek, coeff * case.TimeStep)
					elif type == 2:
						f, Ek = case.L2(f, Ek, coeff * case.TimeStep)
				f[f<=1e-14] = 0
				f_ = xp.pad(f, ((0, 1),), mode='wrap')
				f_ /= xp.trapz(xp.trapz(f_, case.v_, axis=1), case.x_)
				f = f_[:-1, :-1]
			if case.PlotKinetic:
				moments = case.compute_moments(f, case.n_moments)
				im.set_data(f.transpose())
				line_Ek.set_ydata(Ek)
				line_rho_k.set_ydata(moments[0, :])
				for m in range(2, case.n_moments):
					line_Sk[m - 2].set_ydata(moments[m, :])
				for ax in axs_k:
					ax.relim()
					ax.autoscale()
					ax.set_xlim((-case.Lx, case.Lx))
				plt.draw()
				plt.pause(1e-4)
	print('\033[90m        Computation finished in {} seconds \033[00m'.format(int(time.time() - start)))
	data_kinetic, data_fluid = [], []
	if case.SaveFluid:
		data_fluid = fs
	if case.SaveKinetic:
		data_kinetic = f
	if case.ComputeKinetic:
		Hk = case.kinetic_energy(f, Ek)
		print('\033[90m        Error in energy (kinetic) = {:.2e}'.format(xp.abs((Hk - H0k) / H0k)))
	if case.ComputeFluid:
		Hf = case.fluid_energy(fs)
		print('\033[90m        Error in energy (fluid) = {:.2e}'.format(xp.abs((Hf - H0f) / H0f)))
	plt.ioff()
	plt.show()

def save_data(data_fluid, data_kinetic, timestr, case, display=True):
	if case.SaveFluid or case.SaveKinetic:
		mdic = case.DictParams.copy()
		mdic.update({'data_k': data_kinetic, 'data_f': data_fluid})
		date_today = date.today().strftime(" %B %d, %Y")
		mdic.update({'date': date_today, 'author': 'cristel.chandre@univ-amu.fr'})
		name_file = type(case).__name__ + '_' + timestr + '.mat'
		savemat(name_file, mdic)
		if display:
			print('\033[90m        Results saved in {} \033[00m'.format(name_file))
