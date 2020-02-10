FEST-3D: Finite-volume Explicit STructured 3-Dimensional
========
![FEST-3D](docs/media/FEST3D.png)  

[![DOI](https://joss.theoj.org/papers/10.21105/joss.01555/status.svg)](https://doi.org/10.21105/joss.01555)
[![Build Status](https://travis-ci.com/FEST3D/FEST-3D.svg?branch=master)](https://travis-ci.com/FEST3D/FEST-3D)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Join the chat at https://gitter.im/FEST-3D/community](https://badges.gitter.im/FEST-3D/community.svg)](https://gitter.im/FEST-3D/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/187151232.svg)](https://zenodo.org/badge/latestdoi/187151232)

## AIM
FEST3D (Finite-volume Explicit STructured 3-Dimensional) is computational fluid dynamic code written in Fortran 90 for solving Navier-Stokes equations on a structured grid using state of the art finite-volume numerical methods. It is a modular, multiblock finite-volume code developed to solve compressible flow problems encountered in the field of aerodynamics.

## Documentation
[Full Documentation](https://fest3d.github.io/index.html) for the FEST-3D solver is generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford) documentation generator. Documentation is available at [website](https://fest3d.github.io/index.html) .

## Installation
An install.sh bash script has been provided for ubuntu, CentOS and MacOS platforms. This script will take care of dependencies, build executable, run unit and integrated test and set `FEST3D` enviornment variable. You can run install.sh using the following command:
```
$ . install.sh
```
You might be prompted for password to install dependencies. This scirpt will add a following line at the end of either `.bashrc` or `.bash_profile` file
```
export FEST3D="<edit-path-to-FEST-3D-folder>/bin/FEST"
```
For manual installation instruction check out [Documentation](https://fest3d.github.io/page/01_install.html) guide.

### Dependencies
FEST-3D has the following dependencies:

 * Fortran and C++11 compiler
 * MPI library
 * CMAKE and MAKE
 * Python3-numpy for integrated test

So, all the dependencies can be installed using following commands on ubuntu:
```
$sudo apt-get update
$sudo apt-get install build-essential
$sudo apt-get install cmake
$sudo apt-get install gfortran
$sudo apt-get install mpich
$sudo apt-get install python3-numpy
```
The CMAKE implementation in FEST-3D is adapted from the template developed by [Seth Morton](https://github.com/SethMMorton) available at [GitHub](https://github.com/SethMMorton/cmake_fortran_template).


### Building
To build the FEST-3D code, we first generate a new "build" directory/folder in the root directory of the FEST-3D code. Inside the build directory, we execute Cmake with specific FORTRAN compiler. <br>
Following list of comnands will build the executalbe
```
$mkdir build && cd build
$cmake ..
$make -j 4
```
These commands will create a binary file named FEST3D in the ```bin``` folder; located in the root directory of the FEST-3D code.

### Check
Make sure the `FEST3D` enviornment variable is set properly by using following command:
```
echo $FEST3D
```
If output is string of zero length then execute following and additionaly add it to `.bashrc` or `.bash_profile` for permanent effect.
```
export FEST3D="<edit-path-to-FEST-3D-folder>/bin/FEST"
```


## Tutorials
For tutorials check out description of [test cases](https://fest3d.github.io/page/05_tutorials/index.html). A step by step instruction for Lid-Driven cavity test case is given [here](https://fest3d.github.io/page/04_Steps_to_run_FEST3D.html).


## Citation
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01555/status.svg)](https://doi.org/10.21105/joss.01555)

User are requested to cite the following JOSS article for any research publications made using the FEST-3D solver.

Jatinder Pal Singh Sandhu, Anant Girdhar, Rakesh Ramakrishnan, R. Teja, Santanu Ghosh, **FEST-3D: Finite-volume Explicit STructured 3-Dimensional solver**, *Journal of Open Source Software*, 5(46), 1555, [https://doi.org/10.21105/joss.01555](https://joss.theoj.org/papers/10.21105/joss.01555).


## Reference
Details about the algorithms used in the FEST-3D code can be found in the publication by Jatinder Pal Singh Sandhu et al. (_Singh Sandhu, J. P., Girdhar, A., Ramakrishnan, R., Teja, R. D., & Ghosh, S., **A convergence study of solutions using two two-equation RANS turbulence models on a finite volume solver for structured grids**, AIAA 2018-3859_).


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
