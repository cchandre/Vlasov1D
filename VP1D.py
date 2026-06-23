#
# BSD 2-Clause License
#
# Copyright (c) 2026, Cristel Chandre
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

import numpy as np
import sympy as sp
from scipy.fft import rfft, irfft, rfftfreq
# from pyfftw.interfaces.scipy_fft import rfft, irfft, rfftfreq
# import pyfftw
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from dataclasses import dataclass
from pyhamsys import HamSys

# pyfftw.config.PLANNER_EFFORT = 'FFTW_ESTIMATE'

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif', size=20)
plt.rc('axes', labelsize=26)
plt.rc('lines', linewidth=3)
plt.rc('image', cmap='bwr')

@dataclass
class VP1Dparams:
    Lx: float
    Lv: float
    Nx: int
    Nv: int
    precision: float
    f_init: callable
    show_distribution: bool = False

def _rfft(f, precision=None, tail_indx=None, axis=0):
    fft_f = rfft(f, axis=axis)
    if precision is not None:
        fft_f[np.abs(fft_f) <= precision] = 0
    if tail_indx is not None:
        idx = [slice(None)] * f.ndim
        idx[axis] = slice(tail_indx, None)
        fft_f[tuple(idx)] = 0
    return fft_f

class VP1D(HamSys):
	
    def __init__(self, params: VP1Dparams, *, type='kinetic', N_moments=None, fluid_metric=None, fluid_0=None, fluid_2=None):
        super().__init__()
        self.Lx, self.Lv, self.Nx, self.Nv = params.Lx, params.Lv, params.Nx, params.Nv
        self.x = np.linspace(-self.Lx, self.Lx, self.Nx, endpoint=False, dtype=np.float64)
        self.v = np.linspace(-self.Lv, self.Lv, self.Nv, endpoint=False, dtype=np.float64)
        self.kx = np.pi / self.Lx * rfftfreq(self.Nx, d=1/self.Nx)
        self.kv = np.pi / self.Lv * rfftfreq(self.Nv, d=1/self.Nv)
        self.precision = params.precision
        self.tail_indx_x = 3 * self.Nx // 8
        self.tail_indx_v = 3 * self.Nv // 8
        self.div = np.divide(1, 1j * self.kx, where=self.kx!=0, out=None)
        self.div[0] = 0
        self.f0 = params.f_init(self.x[:, None], self.v[None, :])
        if params.show_distribution:
            self.ax, self.im = self.plot_f(self.f0)
        if type == 'fluid':
            if not isinstance(N_moments, int) or N_moments < 1:
                raise ValueError("N_moments must be a positive integer for fluid closure.")
            self.N_moments = N_moments
            self.fluid_metric = fluid_metric
            self.fluid_0 = fluid_0
            self.fluid_2 = fluid_2
            
    def fluid_0(self, t, nu):
        return self.fluid

    def compute_E(self, rho):
        return irfft(self.div * _rfft(rho, precision=self.precision, tail_indx=self.tail_indx_x, axis=0), axis=0)

    def chi(self, h, t, f):
        f_ = irfft(np.exp(-1j * self.kx[:, None] * self.v[None, :] * h) * _rfft(f, precision=self.precision, tail_indx=self.tail_indx_x, axis=0), axis=0)
        rho = trapezoid(f_, self.v, axis=1)
        E = self.compute_E(rho)
        return irfft(np.exp(-1j * E[:, None] * self.kv[None, :] * h) * _rfft(f_, precision=self.precision, tail_indx=self.tail_indx_v,axis=1), axis=1)

    def chi_star(self, h, t, f):
        rho = trapezoid(f, self.v, axis=1)
        E = self.compute_E(rho)
        f_ = irfft(np.exp(-1j * E[:, None] * self.kv[None, :] * h) * _rfft(f, precision=self.precision, tail_indx=self.tail_indx_v, axis=1), axis=1)
        return irfft(np.exp(-1j * self.kx[:, None] * self.v[None, :] * h) * _rfft(f_, precision=self.precision, tail_indx=self.tail_indx_x, axis=0), axis=0)
    
    def compute_moments(self, f, n, centered=True):
        v = self.v
        rho = trapezoid(f, v, axis=1)
        table_moments = [rho]
        if not centered:
            for m in range(1, n + 1):
                Pm = trapezoid(v[None, :]**m * f, v, axis=1)
                table_moments.append(Pm)
        else:
            u = trapezoid(v[None, :] * f, v, axis=1) / rho
            table_moments.append(u)
            _v_ = v[None, :] - u[:, None] 
            for m in range(2, n+1):
                Sm = trapezoid(_v_**m * f, v, axis=1) / rho**(m+1)
                table_moments.append(Sm)
        return np.array(table_moments)
    
    def hamiltonian(self, t, f):
        f_ = f if f.shape == (self.Nx, self.Nv) else f.T
        E = self.compute_E(f_)
        kinetic = trapezoid(trapezoid(self.v[None, :]**2 * f_, self.v, axis=1), self.x)
        potential = trapezoid(E**2, self.x)
        return (kinetic + potential) / 2
    
    def casimirs_kinetic(self, f, n):
        return [trapezoid(trapezoid(f**m, self.v, axis=1), self.x) for m in range(1, n+1)]
    
    def visualize_f(self, t, f, t_eval):
        if hasattr(self, 'ax') and np.any(np.isclose(t, t_eval, atol=1e-5)):
            self.im.set_data(f.transpose())
            self.ax.set_title(f"$t = {t:.1f}$", loc='right', pad=-10)
            self.ax.figure.canvas.draw_idle()
            plt.pause(0.01)
    
    def plot_f(self, f):
        fig, ax = plt.subplots(figsize=(10, 6.5))
        fig.canvas.manager.set_window_title(r'Distribution function f(x,v,t)')
        ax.set_title('$t = 0.0$', loc='right', pad=-10)
        im = ax.imshow(f.transpose(), interpolation='gaussian', origin='lower', aspect='auto', extent=(-self.Lx, self.Lx, -self.Lv, self.Lv))
        ax.set_xlabel('$x$')
        ax.set_ylabel('$v$')
        fig.colorbar(im, ax=ax)
        plt.pause(0.01)
        return ax, im
    
    def eqn_fluid(self, t, nu):
        E = self.compute_E(nu[0])
