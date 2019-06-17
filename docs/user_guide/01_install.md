title: Installation

#Installing FEST-3D
The installation instructions are listed here for Linux distribution. Specifically, examples are given for UBUNTU distribution.

## Dependencies
FEST-3D has the following dependencies:

 * Fortran compiler
 * MPI library
 * CMAKE and MAKE

### Compiler
To build the FEST-3D code written in FORTRAN 90 code, a FORTRAN compiler is required. Either open-source gfortran or commercial ifort FORTRAN compilers can be used.
To install gfortran on Ubuntu, use following command:<br>
```$sudo apt-get install gfortran```
### Distributive MPI Library
To have a faster computation, the FEST-3D code uses the MPICH library to distribute computation load to different processors within one machine or across multiple machines. You can install MPICH on Ubuntu machine using:<br>
```$sudo apt-get install mpich```
### CMAKE and Make
CMake is required to generate a Makefile, which in turn generate executable. The following command is useful to install cmake:<br>
```$sudo apt-get install build-essential``` <br>
```$sudo apt-get install cmake```

So, in summary, all the dependencies can be installed using following commands:
```
$sudo apt-get update
$sudo apt-get install build-essential
$sudo apt-get install cmake
$sudo apt-get install gfortran
$sudo apt-get install mpich
```

## Building
To build the FEST-3D code, we first generate a new "build" directory/folder in the root directory of the FEST-3D code. Inside the build directory, we execute Cmake with specific FORTRAN compiler. <br>
```$FC=mpif90 cmake ..```<br>
```$make -j 4 ```<br>
Here -j is a parallel building option which speeds up the building process. You can replace the number after "-j" with a number of processors you want to use for the building. 
@note if you have installed OPENMPI also along with MPICH, then the mpif90 might be pointing to the `mpif90.openmpi` instead of `mpif90.mpich`. In that case use following command<br>
```$FC=mpif90.mpich cmake ..```<br>



so, in summary, following list of comnands will build the executalbe
```
$mkdir build && cd build
$FC=mpif90 cmake ..
$make -j 4
```

These commands will create a binary file named FEST3D in the ```bin``` folder; located in the root directory of the FEST-3D code.

