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
from pyfftw.interfaces.scipy_fft import rfft, irfft, rfftfreq
from scipy.integrate import trapezoid
from dataclasses import dataclass
from pyhamsys import HamSys

@dataclass
class VP1Dparams:
    Lx: float
    Lv: float
    Nx: int
    Nv: int
    precision: float
    f_init: callable

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
	
    def __init__(self, params: VP1Dparams):
        self.x = np.linspace(-params.Lx, params.Lx, params.Nx, endpoint=False, dtype=np.float64)
        self.v = np.linspace(-params.Lv, params.Lv, params.Nv, endpoint=False, dtype=np.float64)
        self.kx = np.pi / params.Lx * rfftfreq(params.Nx, d=1/params.Nx)
        self.kv = np.pi / params.Lv * rfftfreq(params.Nv, d=1/params.Nv)
        self.precision = params.precision
        self.tail_indx_x = 3 * params.Nx // 8
        self.tail_indx_v = 3 * params.Nv // 8
        div = np.divide(1, 1j * self.kx, where=self.kx!=0)
        div[0] = 0
        self.f_init = params.f_init(self.x[:, None], self.v[None, :])
        rho0 = trapezoid(self.f_init, self.v, axis=1)
        self.E_init = irfft(div * _rfft(rho0, precision=self.precision, tail_indx=self.tail_indx_x, axis=0), axis=0)

    def chi(self, h, t, y):
        f, E = y
        f_ = irfft(np.exp(-1j * self.kx[:, None] * self.v[None, :] * h) * _rfft(f, precision=self.precision, tail_indx=self.tail_indx_x, axis=0), axis=0)
        E_ = self.E(trapezoid(f_, self.v, axis=1))
        f = irfft(np.exp(-1j * E_[:, None] * self.kv[None, :] * h) * _rfft(f_, precision=self.precision, tail_indx=self.tail_indx_v,axis=1), axis=1)
        return f, E_

    def chi_star(self, h, t, y):
        f, E = y
        f_ = irfft(np.exp(-1j * E[:, None] * self.kv[None, :] * h) * _rfft(f, precision=self.precision, tail_indx=self.tail_indx_v, axis=1), axis=1)
        f_ = irfft(np.exp(-1j * self.kx[:, None] * self.v[None, :] * h) * _rfft(f_, precision=self.precision, tail_indx=self.tail_indx_x, axis=0), axis=0)
        E_ = self.E(trapezoid(f_, self.v, axis=1))
        return f_, E_
    
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
    
    def energy_kinetic(self, f, E):
        return (trapezoid(trapezoid(self.v[None, :]**2 * f, self.v, axis=1), self.x) + trapezoid(E**2, self.x)) / 2
    
    def casimirs_kinetic(self, f, n):
        return [trapezoid(trapezoid(f**m, self.v, axis=1), self.x) for m in range(1, n+1)]