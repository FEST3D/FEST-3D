title: Documentation

![FEST-3D](|media|/FEST3D.png)
{: style="text-align: center" }

#Introduction
FEST-3D is a finite-volume solver build to compute compressible fluid flow problems on a structured grid.

The solver provides a tremendous amount of choices to the user in terms of inviscid flux calculation
scheme, higher face-state reconstruction scheme, and time-integration scheme.

## Algorithms
Details about the algorithms used in the FEST-3D code can be found in the publication by Jatinder Pal Singh Sandhu et al. (_Singh Sandhu, J. P., Girdhar, A., Ramakrishnan, R., Teja, R. D., & Ghosh, S., **A convergence study of solutions using two two-equation RANS turbulence models on a finite volume solver for structured grids**, AIAA 2018-3859_).

##Schemes
### Inviscid flux calculation
* __AUSM__   
_Liou, M.-S. and Steffen, C., “A New Flux Splitting Scheme,” J. Comput. Phys., Vol. 107, 23-39, 1993._
* __LDFSS__  
_Edwards, J.R., 1997. A low-diffusion flux-splitting scheme for Navier-Stokes calculations. Computers & Fluids, 26(6), pp.635-659._
* __AUSM+__  
_Liou, M. S., “A sequel to AUSM: AUSM+,” Journal of Computational Physics, vol. 129, 1996, pp. 364–382._
* __AUSM+-UP__  
_Liou, M. S., “A sequel to AUSM, Part II: AUSM+-up for all speeds,” Journal of Computational Physics, vol. 214, 2006, pp. 137–170._
* __SLAU__  
_Shima, E., and Kitamura, K., “Parameter-Free Simple Low-Dissipation AUSM-Family Scheme for All Speeds,” AIAA Journal, vol. 49, Aug. 2011, pp. 1693–1709._
### Higher-order reconstruction
* __None__ 1rst order accurate in space
* __MUSCL__ 3rd order accurate in space  
_van Leer, B. (1979), Towards the Ultimate Conservative Difference Scheme, V. A Second Order Sequel to Godunov's Method, J. Com. Phys.., 32, 101–136_
* __PPM__ 4th order accurate in space  
_Colella, P., and Woodward, P. R., “The Piecewise Parabolic Method (PPM) for gas-dynamical simulations,” Journal of Computational Physics, vol. 54, Apr. 1984, pp. 174–201._
* __WENO__ 5th order accurate in space  
_Britain, G., and Press, P., “Brief of finite volume WENO method,” vol. v, 1995, pp. 1–4._
* __WENO-NM__ 5th order accurate in space (specifically for non-uniform grid)  
_Huang, W. F., Ren, Y. X., and Jiang, X., “A simple algorithm to improve the performance of the WENO scheme on non-uniform grids,” Acta Mechanica Sinica/Lixue Xuebao, 2017, pp. 1–11._


### Time-integration method
#### Explicit
* __Euler Explicit__ First order accurate in time
* __RK2__ 2nd order accurate in time, Runge-Kutta method
* __RK4__ 4th order accurate in time, Runge-Kutta method
* __TVDRK2__ Total variation diminishing RK2 method for Weno scheme
* __TVDRK3__ Total variation diminishing RK3 method for Weno scheme  
_Computational Fluid Dynamics, Book by Klaus A. Hoffmann and Steve T. Chiang_  
#### Implicit
* __implicit__ Matrix free LU-SGS method, first order accurate in time.  
_Chen, R. F., and Wang, Z. J., “Fast , Block Lower-Upper Symmetric Gauss – Seidel Scheme Introduction,” AIAA Journal, vol. 38, 2000, pp. 2238–2245._  
* __PLUSGS__ Preconditioned Matrix free LU-SGS method for very low speed flow; first order accuate in time  
_Kitamura, K., Shima, E., Fujimoto, K., and Wang, Z. J., “Performance of low-dissipation euler fluxes and preconditioned LU-SGS at low speeds,” Communications in Computational Physics, vol. 10, 2011, pp. 90–119._

### Turbulence model
* __SA__   
_Allmaras, S. R., Johnson, F. T., and Spalart, P. R., “Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model,” Seventh International Conference on Computational Fluid Dynamics (ICCFD7), 2012, pp. 1–11._
* __SST__   
_Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications," AIAA Journal, Vol. 32, No. 8, August 1994, pp. 1598-1605_
* __SST2003__   
_Menter, F. R., Kuntz, M., and Langtry, R., "Ten Years of Industrial Experience with the SST Turbulence Model," Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano, and M. Tummers, Begell House, Inc., 2003, pp. 625 - 632._
* __k-kL__   
_Menter, F. R., and Egorov, Y., “The scale-adaptive simulation method for unsteady turbulent flow predictions. part 1: Theory and model description,” Flow, Turbulence and Combustion, vol. 85, 2010, pp. 113–138._

### Transition model
* __Gamma LCTM2015__   
_Menter, F. R., Smirnov, P. E., Liu, T., and Avancha, R., “A One-Equation Local Correlation-Based Transition Model,” Flow, Turbulence and Combustion, vol. 95, 2015, pp. 583–619._
* __BC__   
_Cakmakcioglu, S. C., Bas, O., and Kaynak, U., “A correlation-based algebraic transition model,” Proceedings of the Institution of Mechanical Engineers, Part C: Journal of Mechanical Engineering Science, vol. 232, Nov. 2018, pp. 3915–3929._


## FEST-3D Team
Over the last five years, many individuals have contributed to the development of FEST-3D. 

* ```Jatinder Pal Singh Sandhu```  
  Ph.D. Student (Current)  
  _Added turbulence and transition models; implicit time-integration method and weno scheme_   

* ```R. D. Teja```  
   B.Tech Student (2016)  
  _Introudued the MPI capability to the solver_  
 
* ```Raskesh Ramakrishnan```  
   Dual Degree Student (2016)   
  _Developed FEST-3D into a three-dimensional laminar flow solver_   

* ```Anant Girdhar```  
    B.Tech Student (2015)   
  _Initiated the FEST-3D code as two-dimensional inviscid flow solver_  

All the above individuals were guided by ```Dr. Santanu Ghosh```.
