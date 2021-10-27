# One-dimensional Vlasov-Poisson equation and its Hamiltonian fluid reductions

- [`Vlasov1D_4field.mlx`](https://github.com/cchandre/Vlasov1D/blob/main/Vlasov1D_4field.mlx): Matlab live script to compute and represent the series expansion of the explicit closure S<sub>4</sub>=S<sub>4</sub>(S<sub>2</sub>,S<sub>3</sub>) and S<sub>5</sub>=S<sub>5</sub>(S<sub>2</sub>,S<sub>3</sub>); also compute the series expansions of the three Casimir invariants C<sub>1</sub>, C<sub>2</sub> and C<sub>3</sub>
- [`ParametricClosure4.nb`](https://github.com/cchandre/Vlasov1D/blob/main/ParametricClosure4.nb): Mathematica notebook to check and represent the parametric closure S<sub>2</sub>=S<sub>2</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>), S<sub>3</sub>=S<sub>3</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>),  S<sub>4</sub>=S<sub>4</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>) and S<sub>5</sub>=S<sub>5</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>)

- `VP1D4f python code`
  - [`VP1D4f_dict.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_dict.py): to be edited to change the parameters of the VP1D4f computation (see below for a dictionary of parameters)

  - [`VP1D4f.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f.py): contains the VP1D4f class and main functions (not to be edited)

  - [`VP1D4f_modules.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_modules.py): contains the methods to run VP1D4f (not to be edited)

  - Once [`VP1D4f_dict.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_dict.py) has been edited with the relevant parameters, run the file as 
  ```sh
  python3 VP1D4f.py
  ```

___
##  Parameter dictionary for VP1D4f

- *kappa*: double; value of &kappa; defining the fluid reduction 
- *Tf*: double; duration of the integration (in units of *&omega;<sub>p</sub><sup>-1</sup>*)
- *integrator_kinetic*: string ('position-Verlet', 'velocity-Verlet', 'Forest-Ruth', 'PEFRL'); choice of solver for the integration of the Vlasov equation
- *nsteps*: integer; number of steps in one period of plasma oscillations (2 *&pi;*/*&omega;<sub>p</sub><sup>*) for the integration of the Vlasov equation
- *integrator_fluid*: string ('RK45', ‘RK23’, ‘DOP853’, ‘Radau’, ‘BDF’, ‘LSODA’); choice of solver for the integration of the fluid equation
- *precision_fluid*: double; numerical precision of the integrator for the fluid equations
- *n_moments*: integer; number of moments to be represented
- *n_casimirs*: integrer; number of Casimir invariants (of the Vlasov model) to be monitored 

- *Lx*: double; the *x*-axis is (-*Lx*, *Lx*)
- *Lv*: double; the *v*-axis is (-*Lv*, *Lv*)
- *Nx*: integer; number of points in *x* to represent the field variables
- *Nv*: integer; number of points in *v* to represent the field variables
- *f_init*: lambda function; initial distribution *f*(*x*,*v*,t=0)
  
- *output_E_modes*: integer; number of Fourier modes of *E(t)* to be saved in the `.mat` file

- *ComputeKinetic*: boolean; if True, integrate the Vlasov-Poisson equation 
- *ComputeFluid*: boolean; if True, integrate the fluid equations
- *SaveKinetic*: boolean; if True, save the results of the integration of the Vlasov-Poisson equation in a `.mat` file
- *SaveFluid*: boolean; if True, save the results of the integration of the fluid equations in a `.mat` file
- *PlotKinetic*: boolean; if True, plot the results of the integration of the Vlasov-Poisson equation
- *PlotFluid*: boolean; if True, plot the results of the integration of the fluid equations
- *darkmode*: boolean; if True, plots are done in dark mode

---
Reference: C. Chandre, B.A. Shadwick, *Four-field Hamiltonian closure of the one-dimensional Vlasov-Poisson equation*, preprint (2021)

For more information: <cristel.chandre@univ-amu.fr>

<p align="center">
  <img src="https://github.com/cchandre/Vlasov1D/blob/main/fig_VP1D4f.png" alt="Example" width="400"/>
</p>
