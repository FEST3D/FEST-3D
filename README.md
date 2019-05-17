FEST-3D
========

## AIM
FEST3D (Finite-volume Explicit STructured 3-Dimensional) is computational fluid dynamic code written in Fortran 90 for solving Navier-Stokes equations on a structured grid using state of the art finite-volume numerical methods. It is a modular, multiblock finite-volume code developed to solve compressible flow problems encountered in the field of aerodynamics.

## Documentation
[Full Documentation](https://fest3d.github.io/index.html) for the FEST-3D solver is available at [website](https://fest3d.github.io/index.html) .

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

