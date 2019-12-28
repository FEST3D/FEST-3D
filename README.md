FEST-3D: Finite-volume Explicit STructured 3-Dimensional
========
![FEST-3D](docs/media/FEST3D.png)  

[![Build Status](https://travis-ci.com/FEST3D/FEST-3D.svg?branch=master)](https://travis-ci.com/FEST3D/FEST-3D)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Join the chat at https://gitter.im/FEST-3D/community](https://badges.gitter.im/FEST-3D/community.svg)](https://gitter.im/FEST-3D/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


## AIM
FEST3D (Finite-volume Explicit STructured 3-Dimensional) is computational fluid dynamic code written in Fortran 90 for solving Navier-Stokes equations on a structured grid using state of the art finite-volume numerical methods. It is a modular, multiblock finite-volume code developed to solve compressible flow problems encountered in the field of aerodynamics.

## Documentation
[Full Documentation](https://fest3d.github.io/index.html) for the FEST-3D solver is generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford) documentation generator. Documentation is available at [website](https://fest3d.github.io/index.html) .

## Installation
For installation instruction check out [Documentation](https://fest3d.github.io/page/01_install.html) guide.

### Dependencies
FEST-3D has the following dependencies:

 * Fortran compiler
 * MPI library
 * CMAKE and MAKE

So, all the dependencies can be installed using following commands:
```
$sudo apt-get update
$sudo apt-get install build-essential
$sudo apt-get install cmake
$sudo apt-get install gfortran
$sudo apt-get install mpich
```
The CMAKE implementation in FEST-3D is adapted from the template developed by [Seth Morton](https://github.com/SethMMorton) available at [GitHub](https://github.com/SethMMorton/cmake_fortran_template).


### Building
To build the FEST-3D code, we first generate a new "build" directory/folder in the root directory of the FEST-3D code. Inside the build directory, we execute Cmake with specific FORTRAN compiler. <br>
Following list of comnands will build the executalbe
```
$mkdir build && cd build
$FC=mpif90 cmake ..
$make -j 4
```
These commands will create a binary file named FEST3D in the ```bin``` folder; located in the root directory of the FEST-3D code.


## Tutorials
For tutorials check out description of [test cases](https://fest3d.github.io/page/05_tutorials/index.html).

## Expected Inputs
A python script is provided in the run folder to ease the user interface with FEST-3D. Following inputs are
required in order to setup any case in FEST-3D.   

| Variable |  Expected Input | Description|
|:---------|:----------------|:----------------------|
|RunDir    |String           | Give any name to the Run directory |
|GridDir   |String (path)    | Directory name which contains only grid files. If the domain is decomposed into 4 blocks, then 4 separate files are required. Although, you can give any unique name to the grid files, for simplicity the python code expect grid file name in __grid_nn.txt__, where the nn is the _block-number-1_. So, for 4 blocks, we will use the following: grid_00.txt, grid_01.txt, grid_02.txt, and grid_03.txt. |
|NumberOfBlocks| Integer | Total number of blocks|
|AbsBinaryPath| String (path)| provide the absolute path to the binary. In order to use the executable build in the binary folder of the FEST-3D code, a soft link is created between the FEST3D executable ```bin``` folder and FEST3D in ```bin``` folder of Run folder|
|Control['CFL'] | Real Number greater than zero | Courant–Friedrichs–Lewy number. Low value (less than 1) for explicit scheme and high value for implicit scheme|
|Control['LoadLevel'] | Integer | Restart folder number in the time_directories/ directory |
|Control['MaxIterations'] | Integer greater than SaveIterations | Maximum number of iteration |
|Control['SaveIterations']| Integer lesser than MaxIterations | Save solution state after every these many iteration |
|Control['OutputFileFormat'] | 'vtk' or 'tecplot' | Format of the solution output file |
|Control['OutputDataFormat'] | 'ASCII'  | Type of the data in the output folder. Only ASCII for now. _BINARY_ will be added in later release |
|Control['InputFileFormat'] | 'vtk' or 'tecplot' | Format of the solution file from which solution will be restarted |
|Control['InputDataFormat'] | 'ASCII'| Type of the data in the restart file <br>. Similar to output data type, only ASCII is supported for now.|
|Control['Precision'] | Integer, lesser than 14 and greater than 1 | Data precision for residual output, not used for solution output.|
|Control['Purge'] | Integer |  Number of recent solution folder to keep and delete others. __0__ input will keep all the folders|
|Control['ResidualWriteInterval'] | Integer greater than zero |  Number of iteration after which to save the residual output in the file|
|Control['Tolerance'] | Real number and ["Mass_abs", "Viscous_abs", "Mass_abs", "Resnorm_abs", "Viscous_abs", "Turbulent_abs", "Continuity_abs", "X-mom_abs", "Y-mom_abs", "Z-mom_abs", "Energy_abs", "Mass_rel", "Resnorm_rel", "Viscous_rel", "Turublent_rel", "Continuity_abs", "X-mom_rel", "Y-mom_rel", "Z-mom_rel", "Energy_rel", "TKE_abs", "Tv_abs", "Dissipation_abs", "Omega_abs", "Kl_abs", "TKE_rel", "Tv_rel", "Dissipation_rel", "Omega_rel", "Kl_rel"] | Tolerance value and variable. The solver will stop once this value is achived. List of tolerace variable that can be used is given below|
|Control['DebugLevel'] | 1, 2, 3, 4, or 5 | 5-Only important information is logged, 1-All the information is logged which helps in debuging. Will be remove in later release|
|Scheme['InviscidFlux'] | 'ausm', 'slau', 'ausmUP', or 'ldfss0' | Scheme to calcualte inviscid fluxes through cell faces|
|Scheme['FaceState'] | 'none', 'muscl', 'ppm', or 'weno' | Scheme for higher-order face-state reconstuction|
|Scheme['Limiter'] | '1 1 1  0 0 0' or '0 0 0 0 0 0'|Switch for limiters and pressure based switching when using higher order face-state reconstuction. Three value for I,j, and k direction 1->on  and 0-> off|
|Scheme['TurbulenceLimiter'] | '1 1 1' or '0 0 0' | Switch for limiters when used for higher-order face-state reconstuctio of turbulent variables. 1->on  and 0-> off|
|Scheme['TurbulenceModel'] | 'none', 'sa', 'sst', or 'sst2003' |  Turbulence model|
|Scheme['TransitionModel'] | 'none', 'bc', 'lctm2015' |  Transition model|
|Scheme['TimeStep'] | 'l' or 'g <optional time step>' | Time-step for time-integration. 'l' for local and 'g' for global. In case for global method you can provide exact value to time-step here.|
|Scheme['TimeIntegration'] | 'none', 'RK2', 'RK4', 'TVDRK2', 'TVDRK3', 'implicit', or 'plusgs' | Method for time-integration|
|Scheme['HigherOrderBC'] | 0 or 1 |  Higher order boundary conditions.  1->on  and 0-> off.|
|Flow["NumberOfVariables"] | 5 | Total number of variables to solver. Reduntant and will be depricated in later release. |
|Flow["DensityInf"] | Real Number | Free-stream density |
|Flow["UInf"] | Real Number | Free-stream x-component of velocity|
|Flow["VInf"] | Real Number | Free-stream y-component of velocity|
|Flow["WInf"] | Real Number | Free-stream z-component of velocity|
|Flow["PressureInf"] | Real Number | Free-stream pressure|
|Flow["TurbulenceIntensity"] | Real Number | Free-stream trubulent intensity in (percentage|
|Flow["ViscosityRatio"] | Real Number | Free-stream ratio of turbulence viscosity to molecular viscosity|
|Flow["Intermitency"] | Real Number | Free-stream turbulent intermittency|
|Flow["ReferenceViscosity"] | Real Number | Reference viscosity|
|Flow["ViscosityLaw"] | 'sutherland_law' or 'constant' |  Law used for viscosity variation|
|Flow["ReferenceTemp"] | Real Number | Reference temperature for viscosity variation|
|Flow["SutherlandTemp"] | Real Number | Sutherland temperature|
|Flow["PrandtlNumbers"] | Two real numbers | Prandtl number and turbulent prandtl number|
|Flow["SpecificHeatRatio"] | Real number | Specific heat ratio|
|Flow["GasConstant"] | Real | Gas Constant|
|OutputControl['Out'] | [ "Velocity" , "Density" , "Pressure" , "Mu" , "Mu_t" , "TKE" , "Omega" , "kL" , "tv" , "Wall_distance" , "resnorm" , "TKE_residue" , "Mass_residue" , "X_mom_residue" , "Y_mom_residue" , "Z_mom_residue" , "energy_residue" , "DuDx",   "Dudy",   "DuDz" , "DvDx",   "DvDy",   "DvDz" , "DwDx",   "DWDy",   "DwDz" , "DTDx",   "DTDy",   "DTDz" , "DtkDx",  "DtkDy",  "DtkDz" , "DtwDx",  "DtwDy",  "DtwDz" , "DtvDx",  "DtvDy",  "DtvDz" , "DtkLDx", "DtkLDy", "DtkLDz"]| Variables to write in the output file|
|OutputControl['In'] | ["Velocity" ,"Density" ,"Pressure" ,"viscosity" ,"TKE" ,"Omega" ,"kL" ,"tv"]  | Variables to read in case of restrart|
|ResidualControl['Out'] | ["Mass_abs", "Viscous_abs", "Mass_abs", "Resnorm_abs", "Viscous_abs", "Turbulent_abs", "Continuity_abs", "X-mom_abs", "Y-mom_abs", "Z-mom_abs", "Energy_abs", "Mass_rel", "Resnorm_rel", "Viscous_rel", "Turublent_rel", "Continuity_abs", "X-mom_rel", "Y-mom_rel", "Z-mom_rel", "Energy_rel", "TKE_abs", "Tv_abs", "Dissipation_abs", "Omega_abs", "Kl_abs", "TKE_rel", "Tv_rel", "Dissipation_rel", "Omega_rel", "Kl_rel"] | Residual to write in the resnorm file|
|BoundaryConditions | [-3, -4, -5, -8, -6, -6]  where <-1:'SUPERSONIC INLET', -2:'SUPERSONIC OUTFLOW', -3:'SUBSONIC INFLOW', -4:'SUBSONIC OUTFLOW', -5:'WALL', -6:'SYMMETRY', -7:'Pole', -8:'Far-field', -11:'Total inlet'> | Boundary conditions to used for the six face of the domain|

## Reference and Citation
Details about the algorithms used in the FEST-3D code can be found in the publication by Jatinder Pal Singh Sandhu et al. (_Singh Sandhu, J. P., Girdhar, A., Ramakrishnan, R., Teja, R. D., & Ghosh, S., **A convergence study of solutions using two two-equation RANS turbulence models on a finite volume solver for structured grids**, AIAA 2018-3859_). If you want to cite the FEST-3D code in you work, use the same reference.


## FEST-3D Team

### Advisors:
- [Dr. Santanu Ghosh](https://sites.google.com/view/santanu-ghosh-ae-iitm/home?authuser=2)  
  Assistant Professor  
  Department of Aerospace Engineering  
  Indian Institute of Technology Madras

### Code contributors:
- [Jatinder Pal Singh Sandhu](https://github.com/jayten)   
  Ph.D. Student (Current)  
  Department of Aerospace Engineering  
  Indian Institute of Technology Madras

- Anant Girdhar  
  B.Tech Student (2015)  
  Department of Aerospace Engineering  
  Indian Institute of Technology Madras

- Rakesh Ramakrishnan   
  Dual Degree Student (2016)  
  Department of Aerospace Engineering   
  Indian Institute of Technology Madras

- R D Teja  
  B.Tech Student (2016)    
  Department of Aerospace Engineering  
  Indian Institute of Technology Madras
