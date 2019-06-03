---
title: 'FEST-3D: Finite-volume Explicit STructured 3-Dimensional solver'
authors:
- affiliation: '1'
  name: Jatinder Pal Singh Sandhu
  orcid: 0000-0002-1546-3855
- affiliation: 1
  name: Anant Girdhar
- affiliation: 1
  name: Rakesh Ramakrishnan
- affiliation: 1
  name: R. D. Teja
- affiliation: 1
  name: Dr. Santanu Ghosh
date: "28 May 2019"
output: pdf_document
bibliography: paper.bib
tags:
- CFD
- Turbulence-Transition modeling
- Structured grid
- RANS modeling
- Higher order method
- FORTRAN 90
affiliations:
- index: 1
  name: Department of Aerospace engineering, IIT Madras, Chennai
---

# Summary

The perpetual tryst of automobile, aerospace, and turbo-machinery industries to 
increase the efficiency of their product and decrease its cost of operation 
has lead to the development of large number software tools. 
The well-know open-Source software like OpenFoam and SU2 are such tool which 
helps in improving the design of industrial products and also aide in academic research.

``FEST-3D`` is another CFD software in same league, written in FORTRAN 90 language, to solve 
[compressible Favre-averaged Navier-Stokes equation](https://turbmodels.larc.nasa.gov/implementrans.html) 
using finite-volume method. 
It provides large number of options for schemes and models, 
along with latest turbulence and transition model not available 
in any other Open-source CFD software. 
FEST-3D provides latest one-equation $\gamma$ transition model [@2015Menter] 
and zero-equation BC transition model [@2017Caka]. 
Along with standard turbulence models: SA [@SAmodel] and SST [@menter1994], 
it also provides k-kL [@kkl2015] turbulence model. 
Since ``FEST-3D`` uses structured grid to solve fluid flow problem, a higher order
method of 3rd [@muscl], 4th [@ppm], and 5th [@weno] order accuracy in space 
can be employed; which is difficult obtain for unstructured grid and data-structure

A python script is provided to ease the user interface with the main FEST-3D code.
A large number inputs can easily be set in the first few line of run.py script:

| Variable |  Expected Input | Description|
|:---------|:----------------|:----------------------|
|RunDir    |String           | Give any name to the Run directory |
|GridDir   |String (path)    | Directory name which contains only grid files|
|NumberOfBlocks| Integer | Total number of blocks|
|AbsBinaryPath| String (path)| provide the absolute path to the binary|
|CFL | Real Number greater than zero | Courant–Friedrichs–Lewy number.|
|LoadLevel | Integer | Restart folder number in the time_directories/ directory |
|MaxIterations | Integer greater than SaveIterations | Maximum number of iteration |
|SaveIterations| Integer lesser than MaxIterations | Save solution state after every these many iteration |
|OutputFileFormat | 'vtk' or 'tecplot' | Format of the solution output file |
|OutputDataFormat | 'ASCII'  | Type of the data in the output folder. Only ASCII for now.|
|InputFileFormat | 'vtk' or 'tecplot' | Format of the solution file from which solution will be restarted |
|InputDataFormat | 'ASCII'| Type of the data in the restart file. Only ASCII is supported for now.|
|Precision | Integer, lesser than 14 and greater than 1 | Data precision for residual output, not used for solution output.|
|Purge | Integer |  Number of recent solution folder to keep and delete others. __0__ input will keep all the folders|
|ResidualWriteInterval | Integer greater than zero |  Number of iteration after which to save the residual output in the file|
|Tolerance | Real number and ["Mass_abs", "Viscous_abs", "Mass_abs", "Resnorm_abs", "Viscous_abs", "Turbulent_abs", "Continuity_abs", "X-mom_abs", "Y-mom_abs", "Z-mom_abs", "Energy_abs", "Mass_rel", "Resnorm_rel", "Viscous_rel", "Turublent_rel", "Continuity_abs", "X-mom_rel", "Y-mom_rel", "Z-mom_rel", "Energy_rel", "TKE_abs", "Tv_abs", "Dissipation_abs", "Omega_abs", "Kl_abs", "TKE_rel", "Tv_rel", "Dissipation_rel", "Omega_rel", "Kl_rel"] | Tolerance value and variable. The solver will stop once this value is achived. List of tolerace variable that can be used is given below|
|DebugLevel | 1, 2, 3, 4, or 5 | 5-Only important information is logged, 1-All the information is logged which helps in debuging. Will be remove in later release|
|InviscidFlux | 'ausm', 'slau', 'ausmUP', or 'ldfss0' | Scheme to calcualte inviscid fluxes through cell faces|
|FaceState | 'none', 'muscl', 'ppm', or 'weno' | Scheme for higher-order face-state reconstuction|
|Limiter | '1 1 1  0 0 0' or '0 0 0 0 0 0'|Switch for limiters and pressure based switching when using higher order face-state reconstuction. Three value for I,j, and k direction 1->on  and 0-> off|
|TurbulenceLimiter | '1 1 1' or '0 0 0' | Switch for limiters when used for higher-order face-state reconstuction of turbulent variables. 1->on  and 0-> off|
|TurbulenceModel | 'none', 'sa', 'sst', or 'sst2003' |  Turbulence model|
|TransitionModel | 'none', 'bc', 'lctm2015' |  Transition model|
|TimeStep | 'l' or 'g <optional time step>' | Time-step for time-integration. 'l' for local and 'g' for global. In case for global method you can provide exact value to time-step here.|
|TimeIntegration | 'none', 'RK2', 'RK4', 'TVDRK2', 'TVDRK3', 'implicit', or 'plusgs' | Method for time-integration|
|HigherOrderBC | 0 or 1 |  Higher order boundary conditions.  1->on  and 0-> off.|
|NumberOfVariables | 5 | Total number of variables to solve |
|DensityInf | Real Number | Free-stream density |
|UInf | Real Number | Free-stream x-component of velocity|
|VInf | Real Number | Free-stream y-component of velocity|
|WInf | Real Number | Free-stream z-component of velocity|
|PressureInf | Real Number | Free-stream pressure|
|TurbulenceIntensity | Real Number | Free-stream trubulent intensity in (percentage|
|ViscosityRatio | Real Number | Free-stream viscosity ratio|
|Intermitency | Real Number | Free-stream turbulent intermittency|
|ReferenceViscosity | Real Number | Reference viscosity|
|ViscosityLaw | 'sutherland_law' or 'constant' |  Law used for viscosity variation|
|ReferenceTemp | Real Number | Reference temperature for viscosity variation|
|SutherlandTemp | Real Number | Sutherland temperature|
|PrandtlNumbers | Two real numbers | Prandtl number and turbulent prandtl number|
|SpecificHeatRatio | Real number | Specific heat ratio|
|GasConstant | Real | Gas Constant|
|OutputControl['Out'] | [ "Velocity" , "Density" , "Pressure" , "Mu" , "Mu_t" , "TKE" , "Omega" , "kL" , "tv" , "Wall_distance" , "resnorm" , "TKE_residue" , "Mass_residue" , "X_mom_residue" , "Y_mom_residue" , "Z_mom_residue" , "energy_residue" , "DuDx",   "Dudy",   "DuDz" , "DvDx",   "DvDy",   "DvDz" , "DwDx",   "DWDy",   "DwDz" , "DTDx",   "DTDy",   "DTDz" , "DtkDx",  "DtkDy",  "DtkDz" , "DtwDx",  "DtwDy",  "DtwDz" , "DtvDx",  "DtvDy",  "DtvDz" , "DtkLDx", "DtkLDy", "DtkLDz"]| Variables to write in the output file|
|OutputControl['In'] | ["Velocity" ,"Density" ,"Pressure" ,"viscosity" ,"TKE" ,"Omega" ,"kL" ,"tv"]  | Variables to read in case of restrart|
|ResidualControl['Out'] | ["Mass_abs", "Viscous_abs", "Mass_abs", "Resnorm_abs", "Viscous_abs", "Turbulent_abs", "Continuity_abs", "X-mom_abs", "Y-mom_abs", "Z-mom_abs", "Energy_abs", "Mass_rel", "Resnorm_rel", "Viscous_rel", "Turublent_rel", "Continuity_abs", "X-mom_rel", "Y-mom_rel", "Z-mom_rel", "Energy_rel", "TKE_abs", "Tv_abs", "Dissipation_abs", "Omega_abs", "Kl_abs", "TKE_rel", "Tv_rel", "Dissipation_rel", "Omega_rel", "Kl_rel"] | Residual to write in the resnorm file|
|BoundaryConditions | [-3, -4, -5, -8, -6, -6]  where <-1:'SUPERSONIC INLET', -2:'SUPERSONIC OUTFLOW', -3:'SUBSONIC INFLOW', -4:'SUBSONIC OUTFLOW', -5:'WALL', -6:'SYMMETRY', -7:'Pole', -8:'Far-field', -11:'Total inlet'> | Boundary conditions to used for the six face of the domain|

# Higher order methods
Most of the modern CFD software are based on unstructured-grid data-structure 
and limit themself to maximum 3rd order of accuracy in space, 
as it is computationally expensive and difficult to implement such higher
order method using unstructured data structure.
But, ``FEST-3D`` used structured-grid data-structure 
and provides higher order methods like MUSCL (3rd order
accurate in space), PPM (4th order accurate in space) 
and WENO (5th order accurate in space).
Such higher order methods may proves essential in academic research field.


# Past and current application
``FEST-3D`` is suitable for academic research as well as industrial research. 
It was part of the study performed on the effect of limiter on the solution 
when using higher order methods[@2018jatinder].
Currently, ``FEST-3D`` is being used for the development of a new 
transition model and a 3D immersed boundary method for compressible flows.
The FEST-3D is also being used for teaching purpose in the department of Aerospace Engineering, IIT Madras.


# Acknowledgement
We acknowledge the open source projects which are part of the FEST-3D code. 
[Seth Mortron's](https://github.com/SethMMorton) template project, available at GitHub 
as [cmake_fortran_tempalte](https://github.com/SethMMorton/cmake_fortran_template), 
for Fortran using CMake as the build system. 
To automatically generate documentation, 
[Fortran F/OSS Programmers Group's](https://github.com/Fortran-FOSS-Programmers) 
document generator, called [FORD](https://github.com/Fortran-FOSS-Programmers/ford), is used. 


# References

