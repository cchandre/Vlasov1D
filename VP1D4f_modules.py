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
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import root
from tqdm import trange
from scipy.io import savemat
import time
from datetime import date
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

def integrate(case):
	timestr = time.strftime("%Y%m%d_%H%M")
	state_k = case.f.copy()
	moments = case.compute_moments(state_k, max(5, case.n_moments))
	moments[xp.abs(moments) <= case.precision] = 0
	state_f = moments[:4, :].reshape(4*case.Nx).copy()
	compute_G = lambda G: case.compute_S0(G) - state_f[2*case.Nx:]
	sol = root(compute_G, case.compute_G0(state_f[2*case.Nx:]), tol=case.precision, method='krylov')
	if ('Compute' in case.Fluid) and (not sol.success):
		print('\033[31m        Error: a correct determination for G2 and G3 was not found \033[00m')
		print('\033[31m               -> fluid computation is skipped \033[00m')
		case.Fluid = []
	else:
		state_f[2*case.Nx:] = sol.x
		state_f[xp.abs(state_f) <= case.precision] = 0
		E_f = case.E(state_f[:case.Nx])
		H0_f = case.energy_fluid(state_f, E_f)
		C0_f = case.casimirs_fluid(state_f)
	E_k = case.E(moments[0, :])
	H0_k = case.energy_kinetic(state_k, E_k)
	C0_k = case.casimirs_kinetic(state_k, case.n_casimirs)
	plt.ion()
	if case.darkmode:
		cs = ['k', 'w']
	else:
		cs = ['w', 'k']
	plt.rc('figure', facecolor=cs[0], titlesize=30)
	plt.rc('text', usetex=True, color=cs[1])
	plt.rc('font', family='sans-serif', size=20)
	plt.rc('axes', facecolor=cs[0], edgecolor=cs[1], labelsize=26, labelcolor=cs[1], titlecolor=cs[1])
	plt.rc('xtick', color=cs[1], labelcolor=cs[1])
	plt.rc('ytick', color=cs[1], labelcolor=cs[1])
	plt.rc('lines', linewidth=3)
	plt.rc('image', cmap='bwr')
	if 'Plot' in case.Fluid:
		fig_f = plt.figure(figsize=(8, 10))
		fig_f.canvas.manager.set_window_title('Fluid simulation')
		axs_f = fig_f.add_gridspec(case.n_moments, hspace=0.2).subplots(sharex=True)
		axs_f[0].set_title('$\omega_p t = 0 $', loc='right', pad=20)
		axs_f[-1].plot(case.x, E_f, 'r--', linewidth=1, label=r'$E(0)$')
		line_Ef, = axs_f[-1].plot(case.x, E_f, 'r', label=r'$E(t)$')
		axs_f[0].plot(case.x, moments[0, :], cs[1], linestyle='--', linewidth=1, label=r'$\rho(0)$')
		line_rho_f, = axs_f[0].plot(case.x, moments[0, :], cs[1], label=r'$\rho(t)$')
		line_Sf = []
		for m in range(2, min(case.n_moments, 6)):
			axs_f[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c--', linewidth=1, label=r'$\rho^{{{}}} S_{{{}}}(0)$'.format(m+1, m))
			line_temp, = axs_f[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c', label=r'$\rho^{{{}}} S_{{{}}}(t)$'.format(m+1, m))
			line_Sf.append(line_temp)
		for ax in axs_f:
			ax.set_xlim((-case.Lx, case.Lx))
			ax.legend(loc='upper right', labelcolor='linecolor')
		axs_f[-1].set_xlabel('$x$')
		plt.draw()
		plt.pause(1e-4)
	if 'Save' in case.Fluid:
		data_f = case.output(0, [E_f, state_f], modes=case.output_modes)
	if 'Plot' in case.Kinetic:
		fig = plt.figure(figsize=(7, 6.5))
		fig.canvas.manager.set_window_title(r'Distribution function f(x,v,t)')
		ax_fxvt = plt.gca()
		ax_fxvt.set_title('$\omega_p t = 0 $', loc='right', pad=-10)
		im = plt.imshow(state_k.transpose(), interpolation='gaussian', origin='lower', aspect='auto', extent=(-case.Lx, case.Lx, -case.Lv, case.Lv), vmin=xp.min(state_k), vmax=xp.max(state_k))
		plt.gca().set_ylabel('$v$')
		plt.gca().set_xlabel('$x$')
		plt.colorbar()
		fig_k = plt.figure(figsize=(8, 10))
		fig_k.canvas.manager.set_window_title('Kinetic simulation')
		axs_k = fig_k.add_gridspec(case.n_moments, hspace=0.2).subplots(sharex=True)
		axs_k[0].set_title('$\omega_p t = 0 $', loc='right', pad=20)
		axs_k[-1].plot(case.x, E_k, 'r--', linewidth=1, label=r'$E(0)$')
		line_Ek, = axs_k[-1].plot(case.x, E_k, 'r', label=r'$E(t)$')
		axs_k[0].plot(case.x, moments[0, :], cs[1], linestyle='--', linewidth=1, label=r'$\rho(0)$')
		line_rho_k, = axs_k[0].plot(case.x, moments[0, :], cs[1], label=r'$\rho(t)$')
		line_Sk = []
		for m in range(2, case.n_moments):
			axs_k[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c--', linewidth=1, label=r'$\rho^{{{}}} S_{{{}}}(0)$'.format(m+1, m))
			line_temp, = axs_k[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c', label=r'$\rho^{{{}}} S_{{{}}}(t)$'.format(m+1, m))
			line_Sk.append(line_temp)
		for ax in axs_k:
			ax.set_xlim((-case.Lx, case.Lx))
			ax.legend(loc='upper right', labelcolor='linecolor')
		axs_k[-1].set_xlabel('$x$')
		plt.draw()
		plt.pause(1e-4)
	if 'Save' in case.Kinetic:
		data_k = case.output(0, [E_k, moments[0:4,:].reshape(4*case.Nx)], modes=case.output_modes)
	TimeStep = 1 / case.nsteps
	t_eval = xp.linspace(1/case.nsteps, 1, case.nsteps)
	start = time.time()
	stop_fluid = False
	for _ in trange(xp.int32(case.Tf)):
		if 'Compute' in case.Fluid and not stop_fluid:
			sol = solve_ivp(case.eqn_4f, (0, 1), state_f, t_eval=t_eval, method=case.integrator_fluid, atol=case.precision, rtol=case.precision)
			if sol.status!=0:
				print('\033[33m        Warning: fluid simulation stopped before the end \033[00m')
				stop_fluid = True
			else:
				state_f = sol.y[:, -1]
				if xp.min(state_f[2*case.Nx:3*case.Nx]) <= case.precision:
					print('\033[31m        Error: fluid simulation with S2<0 \033[00m')
					stop_fluid = True
				if 'Save' in case.Fluid:
					for t in range(case.nsteps):
						E_f = case.E(sol.y[:case.Nx, t])
						data_f = xp.vstack((data_f, case.output(_ +  t_eval[t], [E_f, sol.y[:, t]], modes=case.output_modes)))
			if 'Plot' in case.Fluid:
				axs_f[0].set_title('$\omega_p t = {{{}}}$'.format(_ + 1), loc='right', pad=20)
				rho, u, G2, G3 = xp.split(state_f, 4)
				E_f = case.E(rho)
				line_rho_f.set_ydata(rho)
				line_Ef.set_ydata(E_f)
				S = case.compute_S(G2, G3)[:min(case.n_moments-2, 4)]
				for m in range(min(case.n_moments-2, 4)):
					line_Sf[m].set_ydata(rho**(m+1) * S[m])
				for ax in axs_f:
					ax.relim()
					ax.autoscale()
					ax.set_xlim((-case.Lx, case.Lx))
				plt.draw()
				plt.pause(1e-4)
		if 'Compute' in case.Kinetic:
			for t in range(case.nsteps):
				for coeff, type in zip(case.integr_coeff, case.integr_type):
					if type == 1:
						state_k, E_k = case.L1(state_k, E_k, coeff * TimeStep)
					elif type == 2:
						state_k, E_k = case.L2(state_k, E_k, coeff * TimeStep)
				if 'Save' in case.Kinetic:
					moments = case.compute_moments(state_k, 4)
					data_k = xp.vstack((data_k, case.output(_ + (t+1) * TimeStep, [E_k, moments[0:4,:].reshape(4*case.Nx)], modes=case.output_modes)))
				state_k[state_k<=case.precision] = 0
				f_ = xp.pad(state_k, ((0, 1),), mode='wrap')
				f_ *= case.f0 / simpson(simpson(f_, case.v_, axis=1), case.x_)
				state_k = f_[:-1, :-1]
			if 'Plot' in case.Kinetic:
				ax_fxvt.set_title('$\omega_p t = {{{}}}$'.format(_ + 1), loc='right', pad=-10)
				axs_k[0].set_title('$\omega_p t = {{{}}}$'.format(_ + 1), loc='right', pad=20)
				moments = case.compute_moments(state_k, case.n_moments)
				im.set_data(state_k.transpose())
				line_Ek.set_ydata(E_k)
				line_rho_k.set_ydata(moments[0, :])
				for m in range(2, case.n_moments):
					line_Sk[m - 2].set_ydata(moments[0, :]**(m+1) * moments[m, :])
				for ax in axs_k:
					ax.relim()
					ax.autoscale()
					ax.set_xlim((-case.Lx, case.Lx))
				plt.draw()
				plt.pause(1e-4)
	print('\033[90m        Computation finished in {} seconds \033[00m'.format(int(time.time() - start)))
	if 'Save' in case.Fluid:
		save_data(state_f, data_f, timestr, case, model='Fluid')
	if 'Save' in case.Kinetic:
		save_data(state_k, data_k, timestr, case, model='Kinetic')
	if 'Compute' in case.Kinetic:
		H_k = case.energy_kinetic(state_k, E_k)
		print('\033[90m        Error in energy (kinetic) = {:.2e}'.format(xp.abs(H_k - H0_k)))
		for indx, C in enumerate(case.casimirs_kinetic(state_k, case.n_casimirs)):
			print('\033[90m        Error in Casimir C{:d} (kinetic) = {:.2e}'.format(indx + 1, xp.abs(C - C0_k[indx])))
	if 'Compute' in case.Fluid:
		H = case.energy_fluid(state_f, E_f)
		print('\033[90m        Error in energy (fluid) = {:.2e}'.format(xp.abs(H - H0_f)))
		for indx, C in enumerate(case.casimirs_fluid(state_f)):
			print('\033[90m        Error in Casimir C{:d} (fluid) = {:.2e}'.format(indx + 1, xp.abs(C - C0_f[indx])))
	plt.ioff()
	plt.show()

def save_data(state, data, timestr, case, model=[]):
	mdic = case.DictParams.copy()
	mdic.update({'final': state, 'data': data})
	date_today = date.today().strftime(" %B %d, %Y")
	mdic.update({'date': date_today, 'author': 'cristel.chandre@cnrs.fr'})
	name_file = type(case).__name__ + '_' + model + '_' + timestr + '.mat'
	savemat(name_file, mdic)
	print('\033[90m        {} results saved in {} \033[00m'.format(model, name_file))
