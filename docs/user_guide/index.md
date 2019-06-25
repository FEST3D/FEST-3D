title: Documentation

![FEST-3D](|media|/FEST3D.png)
{: style="text-align: center" }

#Introduction
FEST-3D is a finite-volume solver build to compute incompressible/compressible, and inviscid/laminar/transitonal/turbulent fluid flow problems on structured grids.

Highlights: The solver provides multiple choices to the user in terms of inviscid flux calculation
scheme, higher-order face-state reconstruction scheme, time-integration scheme, and turbulence and transition models.

##Schemes
### Inviscid flux calculation
* __AUSM__   
_Liou, M.-S. and Steffen, C., “A New Flux Splitting Scheme,” J. Comput. Phys., vol. 107, no. 1, pp. 23-39, 1993._
* __LDFSS__  
_Edwards, J.R., A low-diffusion flux-splitting scheme for Navier-Stokes calculations. Computers & Fluids, vol. 26, no. 6, pp.635-659, 1997._
* __AUSM+__  
_Liou, M. S., “A sequel to AUSM: AUSM+,” Journal of Computational Physics, vol. 129, no. 2, pp. 364–382, 1996._
* __AUSM+-UP__  
_Liou, M. S., “A sequel to AUSM, Part II: AUSM+-up for all speeds,” Journal of Computational Physics, vol. 214, no. 1, pp. 137–170, 2006._
* __SLAU__  
_Shima, E., and Kitamura, K., “Parameter-Free Simple Low-Dissipation AUSM-Family Scheme for All Speeds,” AIAA Journal, vol. 49, no. 8, pp. 1693–1709, 2011._
### Higher-order spatial reconstruction
* __None__ 1rst order accurate in space
* __MUSCL__ 3rd order accurate in space  
_van Leer, B., Towards the Ultimate Conservative Difference Scheme, V. A Second Order Sequel to Godunov's Method, J. Com. Phys., vol. 32, no. 1, pp. 101–136, 1979_
* __PPM__ 4th order accurate in space  
_Colella, P., and Woodward, P. R., “The Piecewise Parabolic Method (PPM) for gas-dynamical simulations,” Journal of Computational Physics, vol. 54, no. 1, pp. 174–201, 1984._
* __WENO__ 5th order accurate in space  
_Shu, C.-W., “High-order Finite Difference and Finite Volume WENO Schemes and Discontinuous Galerkin Methods for CFD,” International Journal of Computational Fluid Dynamics, vol. 17, no. 2, pp. 107–118, 2003._
* __WENO-NM__ 5th order accurate in space (specifically for non-uniform grid)  
_Huang, W.-F., Ren, Y.-X., and Jiang, X., “A simple algorithm to improve the performance of the WENO scheme on non-uniform grids,” Acta Mechanica Sinica, vol. 34, no. 1, pp. 37–47, 2018._

### Temporal integration
#### Explicit
* __Euler Explicit__ First order accurate in time
* __RK2__ 2nd order accurate in time, Runge-Kutta method
* __RK4__ 4th order accurate in time, Runge-Kutta method
* __TVDRK2__ Total variation diminishing RK2 method for Weno scheme
* __TVDRK3__ Total variation diminishing RK3 method for Weno scheme  
_Hoffmann, Klaus A., and Steve T. Chiang. "Computational fluid dynamics volume I." Engineering Education System, 2000._

#### Implicit
* __implicit__ Matrix free LU-SGS method, first order accurate in time.  
_Chen, R. F., and Wang, Z. J., “Fast , Block Lower-Upper Symmetric Gauss – Seidel Scheme Introduction,” AIAA Journal, vol. 38, no. 12, pp. 2238–2245, 2000._  
* __PLUSGS__ Preconditioned Matrix free LU-SGS method for very low speed flow; first order accuate in time  
_Kitamura, K., Shima, E., Fujimoto, K., and Wang, Z. J., “Performance of Low-Dissipation Euler Fluxes and Preconditioned LU-SGS at Low Speeds,” Communications in Computational Physics, vol. 10, no. 1, pp. 90–119, 2011._

### Turbulence model
* __SA__   
_Allmaras, S. R., Johnson, F. T., and Spalart, P. R., “Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model,” Seventh International Conference on Computational Fluid Dynamics (ICCFD7), 2012._<br>
_Spalart, P. R., and Allmaras, S., “A one-equation turbulence model for aerodynamic flows,” 30th Aerospace Sciences Meeting and Exhibit, 1992._
* __SST__   
_Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications," AIAA Journal, vol. 32, no. 8, pp. 1598-1605, 1994._
* __SST2003__   
_Menter, F. R., Kuntz, M., and Langtry, R., "Ten Years of Industrial Experience with the SST Turbulence Model," Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano, and M. Tummers, Begell House, Inc., pp. 625 - 632, 2003._
* __k-kL__   
_Menter, F. R., and Egorov, Y., “The scale-adaptive simulation method for unsteady turbulent flow predictions. part 1: Theory and model description,” Flow, Turbulence and Combustion, vol. 85, no. 1, pp. 113–138, 2010._

### Transition model
* __Gamma LCTM2015__   
_Menter, F. R., Smirnov, P. E., Liu, T., and Avancha, R., “A One-Equation Local Correlation-Based Transition Model,” Flow, Turbulence and Combustion, vol. 95, no. 4, pp. 583–619, 2015._
* __SA-BC__   
_Cakmakcioglu, S. C., Bas, O., and Kaynak, U., “A correlation-based algebraic transition model,” Proceedings of the Institution of Mechanical Engineers, Part C: Journal of Mechanical Engineering Science, vol. 232, no. 21, pp. 3915–3929, 2018._


## FEST-3D Team
Over the last five years, many individuals have contributed to the development of FEST-3D. The team includes:

* ```Jatinder Pal Singh Sandhu```  
  Ph.D. Student (Current)  
  _Added turbulence and transition models:SST, SA, k-kL; implicit time-integration method: LU-SGS and PLU-SGS; approximate Reimann solver: SLAU, AUSM+-UP, AUSM+; and 5th order weno scheme_   

* ```R. D. Teja```  
   B.Tech Student (2016)  
  _Parallelized FEST-3D using MPI routines_  

* ```Raskesh Ramakrishnan```  
   Dual Degree Student (2016)   
  Modified FEST-3D into a three-dimensional laminar flow solver_   

* ```Anant Girdhar```  
    B.Tech Student (2015)   
  _Developed the FEST-3D code as a modular two-dimensional inviscid flow solver for strutured grids_  

All the above individuals were guided by ```Dr. Santanu Ghosh```.
