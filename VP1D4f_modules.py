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
	f = case.f.copy()
	moments = case.compute_moments(f, max(5, case.n_moments))
	moments[xp.abs(moments) <= case.precision] = 0
	fs = moments[:4, :].reshape(4*case.Nx).copy()
	compute_G = lambda G: case.compute_S0(G) - fs[2*case.Nx:]
	sol = root(compute_G, case.compute_G0(fs[2*case.Nx:]), tol=case.precision, method='krylov')
	if ('Compute' in case.Fluid) and (not sol.success):
		print('\033[31m        Error: a correct determination for G2 and G3 was not found \033[00m')
		print('\033[31m               -> fluid computation is skipped \033[00m')
		case.Fluid = []
	else:
		fs[2*case.Nx:] = sol.x
		fs[xp.abs(fs) <= case.precision] = 0
		Ef = case.E(fs[:case.Nx])
		H0f = case.energy_fluid(fs, Ef)
		C0f = case.casimirs_fluid(fs)
	Ek = case.E(moments[0, :])
	H0k = case.energy_kinetic(f, Ek)
	C0k = case.casimirs_kinetic(f, case.n_casimirs)
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
		axs_f[-1].plot(case.x, Ef, 'r--', linewidth=1, label=r'$E(0)$')
		line_Ef, = axs_f[-1].plot(case.x, Ef, 'r', label=r'$E(t)$')
		axs_f[0].plot(case.x, moments[0, :], cs[1], linestyle='--', linewidth=1, label=r'$\rho(0)$')
		line_rho_f, = axs_f[0].plot(case.x, moments[0, :], cs[1], label=r'$\rho(t)$')
		line_Sf = []
		for m in range(2, min(case.n_moments, 6)):
			axs_f[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c--', linewidth=1, label=r'$\rho^{{{}}} S_{{{}}}(0)$'.format(m+1, m))
			line_temp, = axs_f[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c', label=r'$\rho^{{{}}} S_{{{}}}(t)$'.format(m+1, m))
			line_Sf.append(line_temp)
		for ax in axs_f:
			ax.set_xlim((-case.Lx, case.Lx))
			ax.legend(loc='upper right')
		axs_f[-1].set_xlabel('$x$')
		plt.draw()
		plt.pause(1e-4)
	if 'Save' in case.Fluid:
		suppl_f = case.output(0, Ef)
	if 'Plot' in case.Kinetic:
		fig = plt.figure(figsize=(7, 6.5))
		fig.canvas.manager.set_window_title(r'Distribution function f(x,v,t)')
		ax_fxvt = plt.gca()
		ax_fxvt.set_title('$\omega_p t = 0 $', loc='right', pad=-10)
		im = plt.imshow(f.transpose(), interpolation='gaussian', origin='lower', aspect='auto', extent=(-case.Lx, case.Lx, -case.Lv, case.Lv), vmin=xp.min(f), vmax=xp.max(f))
		plt.gca().set_ylabel('$v$')
		plt.gca().set_xlabel('$x$')
		plt.colorbar()
		fig_k = plt.figure(figsize=(8, 10))
		fig_k.canvas.manager.set_window_title('Kinetic simulation')
		axs_k = fig_k.add_gridspec(case.n_moments, hspace=0.2).subplots(sharex=True)
		axs_k[0].set_title('$\omega_p t = 0 $', loc='right', pad=20)
		axs_k[-1].plot(case.x, Ek, 'r--', linewidth=1, label=r'$E(0)$')
		line_Ek, = axs_k[-1].plot(case.x, Ek, 'r', label=r'$E(t)$')
		axs_k[0].plot(case.x, moments[0, :], cs[1], linestyle='--', linewidth=1, label=r'$\rho(0)$')
		line_rho_k, = axs_k[0].plot(case.x, moments[0, :], cs[1], label=r'$\rho(t)$')
		line_Sk = []
		for m in range(2, case.n_moments):
			axs_k[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c--', linewidth=1, label=r'$\rho^{{{}}} S_{{{}}}(0)$'.format(m+1, m))
			line_temp, = axs_k[m-1].plot(case.x, moments[0, :]**(m+1) * moments[m, :], 'c', label=r'$\rho^{{{}}} S_{{{}}}(t)$'.format(m+1, m))
			line_Sk.append(line_temp)
		for ax in axs_k:
			ax.set_xlim((-case.Lx, case.Lx))
			ax.legend(loc='upper right')
		axs_k[-1].set_xlabel('$x$')
		plt.draw()
		plt.pause(1e-4)
	if 'Save' in case.Kinetic:
		suppl_k = case.output(0, Ek)
	TimeStep = 1 / case.nsteps
	t_eval = xp.linspace(1/case.nsteps, 1, case.nsteps)
	start = time.time()
	for _ in trange(xp.int32(case.Tf)):
		if 'Compute' in case.Fluid:
			sol = solve_ivp(case.eqn_4f, [0, 1], fs, t_eval=t_eval, method=case.integrator_fluid, max_step=TimeStep, atol=case.precision, rtol=case.precision)
			if sol.status!=0:
				print('\033[33m        Fluid simulation stopped before the end \033[00m')
				break
			if xp.min(fs[2*case.Nx:3*case.Nx]) <= case.precision:
				print('\033[31m        Error: fluid simulation with S2<0 \033[00m')
				break
			if 'Save' in case.Fluid:
				for t in range(case.nsteps):
					Ef = case.E(sol.y[:case.Nx, t])
					suppl_f = xp.vstack((suppl_f, case.output(_ +  t_eval[t], Ef)))
			fs = sol.y[:, -1]
			if 'Plot' in case.Fluid:
				axs_f[0].set_title('$\omega_p t = {{{}}}$'.format(_ + 1), loc='right', pad=20)
				rho, u, G2, G3 = xp.split(fs, 4)
				Ef = case.E(rho)
				line_rho_f.set_ydata(rho)
				line_Ef.set_ydata(Ef)
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
						f, Ek = case.L1(f, Ek, coeff * TimeStep)
					elif type == 2:
						f, Ek = case.L2(f, Ek, coeff * TimeStep)
				if 'Save' in case.Kinetic:
					suppl_k = xp.vstack((suppl_k, case.output(_ + (t+1) * TimeStep, Ek)))
				f[f<=case.precision] = 0
				f_ = xp.pad(f, ((0, 1),), mode='wrap')
				f_ *= case.f0 / simpson(simpson(f_, case.v_, axis=1), case.x_)
				f = f_[:-1, :-1]
			if 'Plot' in case.Kinetic:
				ax_fxvt.set_title('$\omega_p t = {{{}}}$'.format(_ + 1), loc='right', pad=-10)
				axs_k[0].set_title('$\omega_p t = {{{}}}$'.format(_ + 1), loc='right', pad=20)
				moments = case.compute_moments(f, case.n_moments)
				im.set_data(f.transpose())
				line_Ek.set_ydata(Ek)
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
		save_data(fs, suppl_f, timestr, case, model='Fluid')
	if 'Save' in case.Kinetic:
		save_data(f, suppl_k, timestr, case, model='Kinetic')
	if 'Compute' in case.Kinetic:
		Hk = case.energy_kinetic(f, Ek)
		print('\033[90m        Error in energy (kinetic) = {:.2e}'.format(xp.abs(Hk - H0k)))
		for indx, Ck in enumerate(case.casimirs_kinetic(f, case.n_casimirs)):
			print('\033[90m        Error in Casimir C{:d} (kinetic) = {:.2e}'.format(indx + 1, xp.abs(Ck - C0k[indx])))
	if 'Compute' in case.Fluid:
		Hf = case.energy_fluid(fs, Ef)
		print('\033[90m        Error in energy (fluid) = {:.2e}'.format(xp.abs(Hf - H0f)))
		for indx, Cf in enumerate(case.casimirs_fluid(fs)):
			print('\033[90m        Error in Casimir C{:d} (fluid) = {:.2e}'.format(indx + 1, xp.abs(Cf - C0f[indx])))
	plt.ioff()
	plt.show()

def save_data(data, suppl, timestr, case, model=[]):
	mdic = case.DictParams.copy()
	mdic.update({'data': data, 'suppl': suppl})
	date_today = date.today().strftime(" %B %d, %Y")
	mdic.update({'date': date_today, 'author': 'cristel.chandre@cnrs.fr'})
	name_file = type(case).__name__ + '_' + model + '_' + timestr + '.mat'
	savemat(name_file, mdic)
	print('\033[90m        {} results saved in {} \033[00m'.format(model, name_file))
