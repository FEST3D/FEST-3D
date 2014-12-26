Solver Structure
================

This document describes the different components of the solver.

This solver is a basic solver that incorporates several CFD schemes to solve
various problems. There are several planned uses for this basic solver:
- [ ] To find metrics that indicate the level of artifical / numerical  
  dissipation present in a scheme. 
- [ ] To solve turbulence related problems. 

To facilitate this, the solver has been partitioned into the following
components that can then be reused as, when and how they are required:

#### State Module

This module contains the state of the system as well as certain procedures that
operate on them. 

The state of the system is defined in terms of the primitive variables:
pressure, velocity and density. 

#### Grid Module

This module contains the specification of the grid as well as procedures to
read in the grid from a file. 

#### Geometry Module

The geometry module consists of various geometrical parameters corresponding to
the grid like normals, areas and volumes of the grid cells. The related
procedures to compute these are also bundled in.

#### The Scheme Files

Various files will contain the various schemes that will be used in the solver.
There will be one file per scheme. 

#### Utils

Various utility procedures are also provided to carry out various routine
operations and to improve the code readability.

#### Main Program (Solver)

The main program, the solver, will reference all these modules to carry out the
required task. Several add-ons can be added to effect the final problems that
are to be tackled. 


There is separate documentation describing the mathematics and physics of these
components. 


Programming Language
--------------------

Most of the solver is developed in FORTRAN 90. Several utilities have been
developed in python and bash. 

A styleguide for FORTRAN has been compiled. For python, refer to PEP 8. 


Development Status
------------------

- [ ] State module
  - [ ] Initialization
  - [ ] State save
  - [ ] Output results

- [ ] Grid module
  - [ ] Grid read

- [ ] Geometry module
  - [ ] Normal calculation
  - [ ] Area calculation
  - [ ] Volume calculation

- [ ] Scheme files
  - [ ] Van Leer scheme
  - [ ] LDFSS(0) scheme
  - [ ] HLLE / CUSP

- [ ] Utils
  - [ ] Memory deallocation

- [ ] Basic solver

