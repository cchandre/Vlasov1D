# One-dimensional Vlasov-Poisson equation and its Hamiltonian fluid reductions

- [`Vlasov1D_4field.mlx`](https://github.com/cchandre/Vlasov1D/blob/main/Vlasov1D_4field.mlx): MATLAB live script to compute and represent the series expansion of the explicit closure S<sub>4</sub>=S<sub>4</sub>(S<sub>2</sub>,S<sub>3</sub>) and S<sub>5</sub>=S<sub>5</sub>(S<sub>2</sub>,S<sub>3</sub>); also compute the series expansions of the three Casimir invariants C<sub>1</sub>, C<sub>2</sub> and C<sub>3</sub>
- [`ParametricClosure4.nb`](https://github.com/cchandre/Vlasov1D/blob/main/ParametricClosure4.nb): Mathematica notebook to check and represent the parametric closure S<sub>2</sub>=S<sub>2</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>), S<sub>3</sub>=S<sub>3</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>),  S<sub>4</sub>=S<sub>4</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>) and S<sub>5</sub>=S<sub>5</sub>(&Gamma;<sub>2</sub>,&Gamma;<sub>3</sub>)

- `VP1D4f python code`
  - [`VP1D4f_dict.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_dict.py): to be edited to change the parameters of the VP1D4f computation (see below for a dictionary of parameters)

  - [`VP1D4f.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f.py): contains the VP1D4f class and main functions (not to be edited)

  - [`VP1D4f_modules.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_modules.py): contains the methods to run VP1D4f (not to be edited)
  
  - [`VP1D4f_AnalyzeData.m`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_AnalyzeData.m): MATLAB script to analyze the output saved in the `.mat` file  

  - Once [`VP1D4f_dict.py`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_dict.py) has been edited with the relevant parameters, run the file as 
  ```sh
  python3 VP1D4f.py
  ```
  or 
  ```sh
  nohup python3 -u VP1D4f.py &>VP1D4f.out < /dev/null &
  ```
  The list of Python packages and their version are specified in [`modules_version.txt`](https://github.com/cchandre/Vlasov1D/blob/main/modules_version.txt)
___
##  Parameter dictionary for VP1D4f

- *kappa*: double; value of &kappa; defining the fluid reduction 
- *Tf*: double; duration of the integration (in units of *&omega;<sub>p</sub><sup>-1</sup>*)
- *integrator_kinetic*: string ('position-Verlet', 'velocity-Verlet', 'Forest-Ruth', 'PEFRL', 'BM4', 'BM6'); choice of solver for the integration of the Vlasov equation
   - 'Forest-Ruth' from [Forest, Ruth, Physica D 43, 105 (1990)](https://doi.org/10.1016/0167-2789(90)90019-L)
   - 'PEFRL' from [Omelyan, Mryglod, Folk, Comput. Phys. Commun. 146, 188 (2002)](https://doi.org/10.1016/S0010-4655(02)00451-4)
   - 'BM4' and 'BM6' refer respectively to BM<sub>6</sub>4 and BM<sub>10</sub>6 from [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7)
- *nsteps*: integer; number of steps in one period of plasma oscillations (1/*&omega;<sub>p</sub><sup>*) for the integration of the Vlasov equation
- *integrator_fluid*: string ('RK45', ‘RK23’, ‘DOP853’, ‘BDF’, ‘LSODA’); choice of solver for the integration of the fluid equation (see [ivp_solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) for more details)
- *precision*: double; numerical precision of the integrator for the fluid equations; threshold for the Fourier transforms
- *n_casimirs*: integer; number of Casimir invariants to be monitored 

- *Lx*: double; the *x*-axis is (-*Lx*, *Lx*)
- *Lv*: double; the *v*-axis is (-*Lv*, *Lv*)
- *Nx*: integer; number of points in *x* to represent the field variables
- *Nv*: integer; number of points in *v* to represent the field variables
- *f_init*: lambda function; initial distribution *f*(*x*,*v*,*t*=0)
  
- *output_var*: string in ['E', 'rho', 'u', 'P', 'q']; variable to be saved in a `.mat` file
- *output_modes*: integer or string in ['all', 'real']; number of Fourier modes of the variable to be saved in the `.mat` file; if 'all', all Fourier modes are saved; if 'real', the variable is saved in real space; use the MATLAB script [`VP1D4f_AnalyzeData.m`](https://github.com/cchandre/Vlasov1D/blob/main/VP1D4f_AnalyzeData.m) to plot the output 

- *Kinetic*: list of strings in ['Compute', 'Plot', 'Save']; list of instructions for the Vlasov-Poisson simulation; if contains 'Save', the results are saved in a `.mat` file
- *Fluid*: list of strings in ['Compute', 'Plot', 'Save']; list of instructions for the fluid simulation; if contains 'Save', the results are saved in a `.mat` file

- *darkmode*: boolean; if True, plots are done in dark mode

---
Reference: C. Chandre, B.A. Shadwick, *Four-field Hamiltonian fluid closures of the one-dimensional Vlasov-Poisson equation*, [Physics of Plasmas](https://aip.scitation.org/journal/php), to appear (2022); [arXiv:2206.06850](https://arxiv.org/abs/2206.06850)

```bibtex
@article{chandre2022,
         title = {Four-field Hamiltonian fluid closures of the one-dimensional Vlasov-Poisson equation},
         author = {Chandre, C.  and Shadwick, B.A. },
         journal = {Physics of Plasmas},
         volume = {29},
         number = {},
         pages = {},
         year = {2022},
         doi = {10.1063/5.0102418},
         URL = {https://doi.org/10.1063/5.0102418}
}
```

For more information: <cristel.chandre@cnrs.fr>

<p align="center">
  <img src="https://github.com/cchandre/Vlasov1D/blob/main/fig_VP1D4f.png" alt="Example" width="400"/>
</p>
