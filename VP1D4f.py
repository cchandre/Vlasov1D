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
from scipy.fft import rfft, irfft, rfftfreq
from scipy.integrate import simpson
from VP1D4f_modules import integrate
from VP1D4f_dict import dict

def main():
	integrate(VP1D4f(dict))

class VP1D4f:
	def __repr__(self):
		return '{self.__class__.__name__}({self.DictParams})'.format(self=self)

	def __str__(self):
		return 'One-dimensional Vlasov-Poisson equation and its 4-moment Hamiltonian fluid closure ({self.__class__.__name__}) for kappa = {self.kappa}'.format(self=self)

	def __init__(self, dict):
		for key in dict:
			setattr(self, key, dict[key])
		self.DictParams = dict
		self.x = xp.linspace(-self.Lx, self.Lx, self.Nx, endpoint=False, dtype=xp.float64)
		self.v = xp.linspace(-self.Lv, self.Lv, self.Nv, endpoint=False, dtype=xp.float64)
		self.x_ = xp.linspace(-self.Lx, self.Lx, self.Nx+1, dtype=xp.float64)
		self.v_ = xp.linspace(-self.Lv, self.Lv, self.Nv+1, dtype=xp.float64)
		self.kx = xp.pi / self.Lx * rfftfreq(self.Nx, d=1/self.Nx)
		div = xp.divide(1, 1j * self.kx, where=self.kx!=0)
		div[0] = 0
		self.kv = xp.pi / self.Lv * rfftfreq(self.Nv, d=1/self.Nv)
		self.tail_indx = [(xp.s_[3*self.Nx//8:], xp.s_[:],), (xp.s_[:], xp.s_[3*self.Nv//8:],)]
		f_ = self.f_init(self.x_[:, xp.newaxis], self.v_[xp.newaxis, :])
		self.f = f_[:-1, :-1]
		self.f0 = simpson(simpson(f_, self.v_, axis=1), self.x_)
		self.E = lambda rho: irfft(div * self.rfft_(rho))
		if self.integrator_kinetic == 'position-Verlet':
			self.integr_coeff = [0.5, 1, 0.5]
			self.integr_type = [1, 2, 1]
		elif self.integrator_kinetic == 'velocity-Verlet':
			self.integr_coeff = [0.5, 1, 0.5]
			self.integr_type = [2, 1, 2]
		elif self.integrator_kinetic == 'Forest-Ruth':
			theta = 1/(2 - 2**(1/3))
			self.integr_coeff = [theta/2, theta, (1-theta)/2, 1-2*theta, (1-theta)/2, theta, theta/2]
			self.integr_type = [1, 2, 1, 2, 1, 2, 1]
		elif self.integrator_kinetic == 'Yoshida6':
			b = xp.array([1 / (2 - 2**(1/(2*n+1))) for n in range(1, 3)])
			a = 1 - 2 * b
			self.integr_coeff = [b[1]*b[0]/2, b[1]*b[0], b[1]*(b[0]+a[0])/2, b[1]*a[0], b[1]*(b[0]+a[0])/2, b[1]*b[0], (b[1]+a[1])*b[0]/2, a[1]*b[0], a[1]*(b[0]+a[0])/2, a[1]*a[0], a[1]*(b[0]+a[0])/2, a[1]*b[0], (b[1]+a[1])*b[0]/2, b[1]*b[0], b[1]*(b[0]+a[0])/2, b[1]*a[0], b[1]*(b[0]+a[0])/2, b[1]*b[0], b[1]*b[0]/2]
			self.integr_type = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1]
		elif self.integrator_kinetic == 'PEFRL':
			xi = 0.1786178958448091
			lam = -0.2123418310626054
			chi = -0.06626458266981849
			self.integr_coeff = [xi, (1-2*lam)/2, chi, lam, 1-2*(chi+xi), lam, chi, (1-2*lam)/2, xi]
			self.integr_type = [1, 2, 1, 2, 1, 2, 1, 2, 1]
		elif self.integrator_kinetic in ['BM4', 'BM6']:
			if self.integrator_kinetic == 'BM4':
				a = [0.0792036964311957, 0.353172906049774, -0.0420650803577195]
				b = [0.209515106613362, -0.143851773179818, 0.434336666566456]
			elif self.integrator_kinetic == 'BM6':
				a = [0.0502627644003922, 0.413514300428344, 0.0450798897943977, -0.188054853819569, 0.541960678450780]
				b = [0.148816447901042, -0.132385865767784, 0.067307604692185, 0.432666402578175, -0.016404589403618]
			c = xp.empty((len(a) + len(b)), dtype=xp.float64)
			c[0::2] = a
			c[1::2] = b
			self.integr_coeff = xp.hstack((c, 1 - 2 * sum(a), c[::-1]))
			c = xp.empty((len(a) + len(b)), dtype=xp.int64)
			c[0::2] = 1
			c[1::2] = 2
			self.integr_type = xp.hstack((c, 1, c[::-1]))

	def L1(self, f, E, dt):
		ft = irfft(xp.exp(-1j * self.kx[:, None] * self.v[None, :] * dt) * self.rfft_(f, axis=0), axis=0)
		Et = self.E(simpson(xp.pad(ft, ((0, 0), (0, 1)), mode='wrap'), self.v_, axis=1))
		return ft, Et

	def L2(self, f, E, dt):
		ft = irfft(xp.exp(-1j * E[:, None] * self.kv[None, :] * dt) * self.rfft_(f, axis=1), axis=1)
		return ft, E

	def eqn_4f(self, t, f):
		rho, u, G2, G3 = xp.split(f, 4)
		E = self.E(rho)
		S2 = G2**3 + G2 * (self.kappa - G2) * G3**2
		DS2DG3 = 2 * G2 * (self.kappa - G2) * G3
		DS2DG2 = 3 * G2**2 + (self.kappa - 2 * G2) * G3**2
		rho_dot = - irfft(1j * self.kx * self.rfft_(rho * u))
		u_dot = - u * irfft(1j * self.kx * self.rfft_(u)) + E - irfft(1j * self.kx * self.rfft_(rho**3 * S2)) / rho
		G2_dot = - u * irfft(1j * self.kx * self.rfft_(G2)) - irfft(1j * self.kx * self.rfft_(rho**2 * DS2DG3)) / (2 * rho)
		G3_dot = - u * irfft(1j * self.kx * self.rfft_(G3)) - irfft(1j * self.kx * self.rfft_(rho**2 * DS2DG2)) / (2 * rho)
		return xp.hstack((rho_dot, u_dot, G2_dot, G3_dot))

	def rfft_(self, f, axis=0):
		fft_f = rfft(f, axis=axis)
		fft_f[xp.abs(fft_f) <= self.precision] = 0
		fft_f[self.tail_indx[axis][:f.ndim]] = 0
		return fft_f

	def compute_S(self, G2, G3):
		return [G2**3 + G2 * (self.kappa - G2) * G3**2,
			G2 * (self.kappa - G2) * G3 * (3 * G2**2 + (self.kappa - 2 * G2) * G3**2),
				9 * self.kappa * G2**5 / 5 + 6 * G2**3 * (self.kappa - G2)**2 * G3**2 + G2 * (self.kappa - G2) * (self.kappa**2 - 3 * G2 * (self.kappa - G2)) * G3**4,
				9 * self.kappa * G2**5 * (self.kappa - G2) * G3 + 10 * G2**3 * (self.kappa - G2)**3 * G3**3 + G2 * (self.kappa-G2) * (self.kappa - 2 * G2) * (self.kappa**2 - 2 * self.kappa * G2 + 2 * G2**2) * G3**5]

	def compute_S0(self, G):
		G2, G3 = xp.split(G, 2)
		return xp.array([G2**3 + G2 * (self.kappa - G2) * G3**2,
			G2 * (self.kappa - G2) * G3 * (3 * G2**2 + (self.kappa - 2 * G2) * G3**2)]).flatten()

	def compute_G0(self, S):
		S2, S3 = xp.split(S, 2)
		sigma = self.kappa - S2**(1/3)
		return xp.array([S2**(1/3) - S3**2 / (27 * S2**(7/3) * sigma),
				S3 / (3 * S2 * sigma)]).flatten()

	def compute_moments(self, f, n):
		f_ = xp.pad(f, ((0, 1),), mode='wrap')
		rho = simpson(f_, self.v_, axis=1)
		u = simpson(self.v_[None, :] * f_, self.v_, axis=1) / rho
		table_moments = xp.vstack((rho[:-1], u[:-1]))
		for m in range(2, n+1):
			Sm = simpson((self.v_[None, :] - u[:, None])**m * f_, self.v_, axis=1) / rho**(m+1)
			table_moments = xp.vstack((table_moments, Sm[:-1]))
		return table_moments

	def energy_kinetic(self, f, E):
		f_ = xp.pad(f, ((0, 1),), mode='wrap')
		E_ = xp.pad(E, (0, 1), mode='wrap')
		return (simpson(simpson(self.v_[None, :]**2 * f_, self.v_, axis=1), self.x_) + simpson(E_**2, self.x_)) / 2

	def energy_fluid(self, f, E):
		rho_, u_, G2_, G3_ = [xp.pad(_, (0, 1), mode='wrap') for _ in xp.split(f, 4)]
		S2_ = self.compute_S(G2_, G3_)[0]
		E_ = xp.pad(E, (0, 1), mode='wrap')
		return simpson(rho_ * u_**2 + rho_**3 * S2_ + E_**2, self.x_) / 2

	def casimirs_kinetic(self, f, n):
		f_ = xp.pad(f, ((0, 1),), mode='wrap')
		return [simpson(simpson(f_**m, self.v_, axis=1), self.x_) for m in range(1, n+1)]

	def casimirs_fluid(self, f):
		rho_, u_, G2_, G3_ = [xp.pad(_, (0, 1), mode='wrap') for _ in xp.split(f, 4)]
		return [simpson(_, self.x_) for _ in [u_ - rho_ * G2_ * G3_, rho_ * G2_, rho_ * G3_]]

	def output(self, t, data, modes=1):
		if self.output_var=='E':
			var = data[0]
		elif self.output_var=='rho':
			var = data[1][:self.Nx]
		elif self.output_var=='u':
			var = data[1][self.Nx:2*self.Nx]
		elif self.output_var=='P':
			var = data[1][:self.Nx]**3 * data[1][2*self.Nx:3*self.Nx]
		elif self.output_var=='q':
			var = 2 * data[1][:self.Nx]**4 * data[1][3*self.Nx:4*self.Nx]
		else:
			var = data[0]
			print('\033[33m        Warning: the output variable is not valid; used E instead \033[00m')
		if isinstance(modes, int):
			return xp.append(t, self.rfft_(var)[0:modes] / var.size)
		elif modes=='real':
			return xp.append(t, var)
		elif modes=='all':
			return xp.append(t, self.rfft_(var) / var.size)

if __name__ == "__main__":
	main()
