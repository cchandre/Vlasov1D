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
from numpy.fft import fft, ifft, fftfreq
from VP1D4f_modules import integrate
from VP1D4f_dict import dict

def main():
	integrate(VP1D4f(dict))

class VP1D4f:
	def __repr__(self):
		return '{self.__class__.__name__}({self.DictParams})'.format(self=self)

	def __str__(self):
		return 'One-dimensional Vlasov-Poisson equation and its 4-moment fluid Hamiltonian closure ({self.__class__.__name__}) for kappa = {self.kappa} and e = {self.qe}'.format(self=self)

	def __init__(self, dict):
		for key in dict:
			setattr(self, key, dict[key])
		self.DictParams = dict
		self.x = xp.linspace(-self.Lx, self.Lx, self.Nx, endpoint=False, dtype=xp.float64)
		self.v = xp.linspace(-self.Lv, self.Lv, self.Nv, endpoint=False, dtype=xp.float64)
		self.x_ = xp.linspace(-self.Lx, self.Lx, self.Nx+1, dtype=xp.float64)
		self.v_ = xp.linspace(-self.Lv, self.Lv, self.Nv+1, dtype=xp.float64)
		self.kx = xp.pi / self.Lx * fftfreq(self.Nx, d=1/self.Nx)
		div = xp.divide(1, 1j * self.kx, where=self.kx!=0)
		div[0] = 0
		self.kv = xp.pi / self.Lv * fftfreq(self.Nv, d=1/self.Nv)
		f_ = self.f_init(self.x_, self.v_)
		self.f = f_[:-1, :-1] / xp.trapz(xp.trapz(f_, self.v_, axis=1), self.x_)
		self.E_fluid = lambda rho: 4 * xp.pi * self.qe * ifft(div * fft(rho)).real
		self.E_kinetic = lambda f: self.E_fluid(xp.trapz(xp.pad(f, ((0, 1),), mode='wrap'), self.v_, axis=1)[:-1])
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
		elif self.integrator_kinetic == 'PEFRL':
			xi = 0.1786178958448091
			lam = -0.2123418310626054
			chi = -0.06626458266981849
			self.integr_coeff = [xi, (1-2*lam)/2, chi, lam, 1-2*(chi+xi), lam, chi, (1-2*lam)/2, xi]
			self.integr_type = [1, 2, 1, 2, 1, 2, 1, 2, 1]

	def L1(self, f, E, dt):
		ft = ifft(xp.exp(-1j * self.kx[:, None] * self.v[None, :] * dt) * fft(f, axis=0), axis=0).real
		Et = self.E_kinetic(ft)
		return ft, Et

	def L2(self, f, E, dt):
		ft = ifft(xp.exp(-1j * self.qe * E[:, None] * self.kv[None, :] * dt) * fft(f, axis=1), axis=1).real
		return ft, E

	def eqn_4f(self, t, fs):
		rho, u, G2, G3 = xp.split(fs, 4)
		E = self.E_fluid(rho)
		S2 = G2**3 + G2 * (self.kappa - G2) * G3**2
		DS2DG3 = 2 * G2 * (self.kappa - G2) * G3
		DS2DG2 = 3 * G2**2 + (self.kappa - 2 * G2) * G3**2
		rho_dot = - ifft(1j * self.kx * fft(rho * u)).real
		u_dot = - u * ifft(1j * self.kx * fft(u)).real + self.qe * E - ifft(1j * self.kx * fft(rho**3 * S2)).real / rho
		G2_dot = - u * ifft(1j * self.kx * fft(G2)).real - ifft(1j * self.kx * fft(rho**2 * DS2DG3)).real / (2 * rho)
		G3_dot = - u * ifft(1j * self.kx * fft(G3)).real - ifft(1j * self.kx * fft(rho**2 * DS2DG2)).real / (2 * rho)
		return xp.hstack((rho_dot, u_dot, G2_dot, G3_dot))

	def compute_S(self, G2, G3):
		return [G2**3 + G2 * (self.kappa - G2) * G3**2,
        		G2 * G3 * (self.kappa - G2) * (3 * G2**2 + (self.kappa - 2 * G2) * G3**2),
				9 * self.kappa / 5 * G2**5 + 6 * G2**3 * (self.kappa - G2)**2 * G3**2 + G2 * (self.kappa - G2) * (self.kappa**2 - 3 * G2 * (self.kappa - G2)) * G3**4,
        		9 * self.kappa * G2**5 * (self.kappa - G2) * G3 + 10 * G2**3 * (self.kappa - G2)**3 * G3**3 + G2 * (self.kappa-G2) * (self.kappa - 2 * G2) * (self.kappa**2 - 2 * self.kappa * G2 + 2 * G2**2) * G3**5]

	def compute_moments(self, f, n):
		f_ = xp.pad(f, ((0, 1),), mode='wrap')
		rho = xp.trapz(f_, self.v_, axis=1)
		u = xp.trapz(self.v_[None, :] * f_, self.v_, axis=1) / rho
		table_moments = xp.vstack((rho[:-1], u[:-1]))
		for m in range(2, n+1):
			Sm = xp.trapz((self.v_[None, :] - u[:, None])**m * f_, self.v_, axis=1) / rho**(m+1)
			table_moments = xp.vstack((table_moments, Sm[:-1]))
		return table_moments

	def kinetic_energy(self,f, E):
		f_ = xp.pad(f, ((0, 1),), mode='wrap')
		E_ = xp.pad(E, (0, 1), mode='wrap')
		return (xp.trapz(xp.trapz(self.v_[None, :]**2 * f_, self.v_, axis=1), self.x_) + xp.trapz(E_**2, self.x_) / (4 * xp.pi)) / 2

	def fluid_energy(self, fs):
		rho, u, G2, G3 = xp.split(fs, 4)
		rho_ = xp.pad(rho, (0, 1), mode='wrap')
		u_ = xp.pad(u, (0, 1), mode='wrap')
		S2_ = xp.pad(self.compute_S(G2, G3)[0], (0, 1), mode='wrap')
		E_ = xp.pad(self.E_fluid(rho), (0, 1), mode='wrap')
		return xp.trapz(rho_ * u_**2 + rho_**3 * S2_ + E_**2 / (4 * xp.pi), self.x_) / 2

if __name__ == "__main__":
	main()
