title: Flat plate (Transition)

# Subsonic transition flow of Mach 0.015, Reynolds number 0.36 million over a flat plate

<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/T3ADomain.png" alt="Contours" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.1 Domain of computation for transition flow over flat plate test case.</figcaption>
  </div>
</figure>


## Problem Statement
In this test case, a subsonic transition flow of Mach 0.015 over a flat plate will be simulated. The domain used, and boundary
condition applied to the domain are illustrated in Fig. 1. The case definition and grid used were obtained from
[ Turbulence Modeling Resource](https://turbmodels.larc.nasa.gov/flatplate.html). 


## Mesh
A structured grid of size 285 × 161 x 2 will be used as shown in Fig. 2. The grid is available in the 
tutorial folder __|rootFolder|/run/Tutorial/TurbulentFlatPlate/CreateBlocks/__. 
@note
If `run` folder is empty, please download the content from [Github](https://github.com/FEST3D/run) 
direcotry or download the zip file [here](https://github.com/FEST3D/run/archive/master.zip).
@endnote
In the ```blocking_point.f90```
edit the number of blocks in the I-direction and J-direction. In this test case, four blocks will
be used. <br>
```integer, parameter :: xblocks = 4```<br>
```integer, parameter :: yblocks = 1```<br>
To compile the `blocking_point.f90` code with gfortran and execute it, commands are written 
in the makefile. Just use the following command to generate grid files in the __grid/__ folder.<br>
```
$make
```
Run the previous command in the __|rootFolder|/run/Tutorial/TurbulentFlatPlate/CreateBlocks/__ directory.
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/T3AGrid.png" alt="mesh" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.2 The Mesh for turbulent flow over flat plate test case.</figcaption>
  </div>
</figure>

## Setup
To set up the case directory, a python automation script is provided: ```automaton.py```. First, set up the most important
parameters, the paths to the grid files and main executable binary to FEST-3D<br>
```RunDir = 'T3A'```   Name of the run directory to create for the current case. <br>
```GridDir= 'CreateBlocks/grid/'```  Path to the folder which only contains the grid files. <br>
```NumberOfBlocks = 4``` Total number of blocks. It should match with the number of grid files available in the __GridDir__ folder <br>
```
AbsBinaryPath="/home/usr/FEST3D/bin/FEST3D"
``` 
Absolute path to the FEST-3D binary. Should be in __|rootfolder|/bin/FEST3D__<br>
Now provide the common control parameters.<br>
```Control['CFL'] = 30.0``` High CFL since implicit time-integration method will be used.<br>
```Control['LoadLevel'] = 0 ```Since simulation will be started from scratch, __0__ is specified.<br>
```Control['MaxIterations'] = 400000``` Maximum number of iteration to perform. <br>
```Control['SaveIterations'] = 10000``` Solution folder will be written every 1000 iteration. <br>
```Control['OutputFileFormat'] = 'tecplot'``` Type of solution file to write. If you have vtk file viewer you can user __vtk__ instead of __tecplot__<br>
```Control['Purge'] = 1``` Only one latest solution folder will be kept in the __time_directories__ and rest will be deleted. <br>
```Control['ResidualWriteInterval'] = 20 ```Write residual in __time_directories/aux/resnorm__ after every 5 iteration. <br>
```Control['Tolerance'] = "1e-13 Continuity_abs"``` Stop the iteration if the absolute residual value of _continuity_ equation is less than 1e-13. <br>
Few scheme parameters for inviscid flow:<br>
```Scheme['InviscidFlux'] = 'ausmUP'``` Inviscid flux-reconstruction shceme. You can use: __slau__, __ldfss0__, __ausmP__ , and __ausm__ instead of __ausmUP__.<br>
```Scheme['FaceState'] = 'muscl'``` Higher-order  face-state reconstruction method. you can use: __none__, __ppm__, and __weno__. <br>
```Scheme['Limiter'] = '0 0 0  0 0 0'``` Switch off the limiter for I,J,and K direction and switch of pressure based switching for all direction.<br>
```Scheme['TurbulenceLimiter'] = '1 1 1'``` Using limiter only for turbulent variables face-state reconstruction improve convergence.<br>
```Scheme['TurbulenceModel']='sst2003'``` using SST2003 turbulence model. Other models __kkl__, __sa__, and __sst__ can be user also. <br>
```Scheme['TransitionModel']='lctm2015'```  __$\gamma$__ one-equation transition model is metioned as transition model. <br>
```Scheme['TimeStep']='g'``` Global time stepping method. You can use local time-stepping method also __l__.<br>
```Scheme['TimeIntegration']='plusgs'``` Preconditioned LU-SGS matrix-free time integration method.<br>
Now, lets define the flow feature of test case:<br>
```Flow["DensityInf"] = 1.2``` Free-stream density.<br>
```Flow["UInf"] = 5.18``` Free-stream x-component of velcotiy vector. (Mach = 0.015 at T_ref=300K)<br> 
```Flow["VInf"] = 0.0``` Free-stream y-component of velcotiy vector. <br>
```Flow["WInf"] = 0.0``` Free-stream z-component of velcotiy vector. <br>
```Flow["PressureInf"] = 103320.0``` Free-stream pressure. <br>
```Flow["ReferenceViscosity"] = 1.8e-5``` Reference viscosity is set such that the Reynolds number (L_ref=1m) is 0.36 million<br>

```OutputControl['Out'] = ["Velocity", "Density", "Pressure", "TKE", "Omega", "Mu", "Mu_t"]``` Variables to write in the output file.<br>
```ResidualControl['Out'] = ["Mass_abs", "Viscous_abs", "Continuity_abs", "TKE_abs", "Omega_abs"] ``` Residual to write in the resnorm file. <br>
```BoundaryConditions = [-3, -4, -5, -8, -6, -6]``` Broad boundary condition.[Subsonic inlet, subsonic outlet, no-slip wall, far-field, and slip wall for rest ]. Since the bottom faces
of the full domain have two different boundary conditions, as shown in Fig. 1, and the blocking provided here is such that first out the four blocks
have inviscid region as a boundary condition at the bottom face and other three have the no-slip adiabatic wall as a boundary condition at
the bottom faces. The boundary condition at the bottom face (__Jmin__) will be changed manually in the layout.md file after the automatic 
case setup. <br>

Rest of the variables should be left to their default value. To execute this script, use the following command:
```
$python automaton.py
```
Now you will see a new folder created with ```RunDir``` name (_FlatPlate_Turb_). Switch to that directory to run the test case. To make sure setup is correct, check
the ```layout.md``` file located in __system/mesh/layout/layout.md__. The file should look like the following:
```
## BLOCK LAYOUT FILE
## ==========================
## NUMBER OF PROCESSES
4
## NUMBER OF ENTRIES PER PROCESS
9
## PROCESS_NO GRID BC_FILE IMIN IMAX JMIN JMAX KMIN KMAX
## ===================================
## PROCESS 0
00  grid_00.txt  bc_00.md  -003  0001  -005  -006  -006  -006
## PROCESS 1
01  grid_01.txt  bc_01.md  0000  0002  -005  -006  -006  -006
## PROCESS 2
02  grid_02.txt  bc_02.md  0001  0003  -005  -006  -006  -006
## PROCESS 3
03  grid_03.txt  bc_03.md  0002  -004  -005  -006  -006  -006
```
Change the Jmin boundary condition for the first block from no-slip wall (-5) to slip wall (-6).
```
## BLOCK LAYOUT FILE
## ==========================
## NUMBER OF PROCESSES
4
## NUMBER OF ENTRIES PER PROCESS
9
## PROCESS_NO GRID BC_FILE IMIN IMAX JMIN JMAX KMIN KMAX
## ===================================
## PROCESS 0
00  grid_00.txt  bc_00.md  -003  0001  -006  -006  -006  -006
## PROCESS 1
01  grid_01.txt  bc_01.md  0000  0002  -005  -006  -006  -006
## PROCESS 2
02  grid_02.txt  bc_02.md  0001  0003  -005  -006  -006  -006
## PROCESS 3
03  grid_03.txt  bc_03.md  0002  -004  -005  -006  -006  -006
```

Finally, to run the simulation use the following command:
```
$nohup bash run.sh &
```
`nohup` helps in avoiding any output on the screen, __&__ execute the last command in the background, allowing you to keep using the terminal.

## Results

<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/T3ACf.png" alt="Contour of v speed" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.3 Coefficient of friction along the surface of the plate for SST turbulence model. </figcaption>
  </div>
</figure>
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/T3ATu.png" alt="Contour of v speed" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.4 Coefficient of friction along the surface of the plate for k-kL turbulence model. </figcaption>
  </div>
</figure>
