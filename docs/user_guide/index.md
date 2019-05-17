title: Documentation

![FEST3D](|media|/FEST3D.png)
{: style="text-align: center" }

#Introduction
FEST-3D is a finite-volume solver build to compute compressible fluid flow problems on a structured grid.
## Governing Equations:
\begin{align} \label{GoverningEqn}
\frac{\partial \rho}{\partial t} + \frac{\partial}{\partial x_j}
  (\rho u_j) &= 0\\ \nonumber
\frac{\partial (\rho  u_i)}{\partial t} +
  \frac{\partial}{\partial x_j}(u_j \rho  u_i) &=
  -\frac{\partial p}{\partial x_i} + \frac{\partial}{\partial x_j}\left(2 \mu_{eff} \left( S_{ij} - \frac{1}{3}
  \frac{\partial u_k}{\partial x_k} \delta_{ij} \right) - \frac{2}{3} \rho k \delta_{ij}\right)\\  \nonumber
\frac{\partial (\rho E)}{\partial t} +
  \frac{\partial}{\partial x_j}(u_j \rho H) &=
  \frac{\partial }{\partial x_j}\left(2 \mu_{eff} \left( S_{ij} - \frac{1}{3} \frac{\partial u_k}{\partial x_k} \delta_{ij} \right) - \frac{2}{3} \rho k \delta_{ij}\right)u_i \\ \nonumber
  &+\frac{\partial}{\partial x_j} \left(
  \left(\frac{c_p \hat \mu}{Pr}+\frac{c_p \hat \mu_t}{Pr_t}\right) \frac{\partial T}{\partial x_j}
  + \left(  \mu + \frac{\mu_t}{\sigma_k} \right) \frac{\partial k}{\partial x_j}\right)
\end{align}

The solver provides a tremendous amount of choices to the user in terms of inviscid flux calculation
scheme, higher face-state reconstruction scheme, and time-integration scheme.

##Schemes
### Inviscid Face-flux calculation
* __AUSM__
* __LDFSS__
* __AUSM+__
* __AUSM+-UP__
* __SLAU__

### Face-state reconstruction
* __None__ 1rst order accurate in space
* __MUSCL__ 3rd order accurate in space
* __PPM__ 4th order accurate in space
* __WENO__ 5th order accurate in space
* __WENO-NM__ 5th order accurate in space (specifically for non-uniform grid)


### Time-integration
#### Explicit
* __Euler Explicit__ First order accurate in time
* __RK2__ 2nd order accurate in time, Runge-Kutta method
* __RK4__ 4th order accurate in time, Runge-Kutta method
* __TVDRK2__ Total variation diminishing RK2 method for Weno scheme
* __TVDRK3__ Total variation diminishing RK3 method for Weno scheme
#### Implicit
* __implicit__ Matrix free LU-SGS method, first order accurate in time.
* __PLUSGS__ Preconditioned Matrix free LU-SGS method for very low speed flow; first order accuate in time

### Turbulence
* __SA__
* __SST__
* __SST2003__
* __k-kL__

### Transition
* __Gamma LCTM2015__
* __BC__


## FEST-3D Team
Over the last five years, many individuals have contributed to the development of FEST-3D. 

* ```Jatinder Pal Singh Sandhu```  
  Ph.D. Student (Current)  
  _Added turbulence and transition models; implicit time-integration method and weno scheme_   

* ```Raskesh Ramakrishnan```  
   Dual Degree Student (2016)   
  _Developed FEST-3D into a three-dimensional laminar flow solver_   

* ```Anant Girdhar```  
    B.Tech Student (2015)   
  _Initiated the FEST-3D code as two-dimensional inviscid flow solver_  

* ```R. D. Teja```  
   B.Tech Student (2016)  
  _Introudued the MPI capability to the solver_  

All the above individuals were guided by ```Dr. Santanu Ghosh```.
