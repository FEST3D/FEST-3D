title: 2D Smooth Bump

# Subsonic inviscid flow of Mach 0.5 past a smooth 2D-bump

<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/2DBumpPressureContour.png" alt="Contour of u speed" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.1 Pressure contour over smoooth bump.</figcaption>
  </div>
</figure>

## Problem Statement
In this test case, a subsonic inviscid flow of Mach 0.5 past a smooth 2D-bump will be simulated. The domain used, and boundary
condition applied to the domain are illustrated in Fig. 2. The case definition and grid used were obtained from
[4th International Workshop on High-Order CFD Methods ](https://how4.cenaero.be/content/bi2-inviscid-flow-over-bump). 
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/2DSmoothBumpDomain.png" alt="Domain" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.2 Domain and boundary conditions for subsonic flow over 2D smooth bump.</figcaption>
  </div>
</figure>


## Mesh
A structured grid of size 97 x 49 × 2 will be used, as shown in Fig. 3. 
The grid is available in the tutorial folder __|rootFolder|/run/Tutorial/2DBump/CreateBlocks/__. 
@note
If `run` folder is empty, please download the content from [Github](https://github.com/FEST3D/run) 
direcotry or download the zip file [here](https://github.com/FEST3D/run/archive/master.zip).
@endnote
In the ```bump.py```
edit the number of blocks in the I-direction . In this test case, only two blocks will
be used. <br>
```imax = 97``` Maximum number of grid points in the I-direction<br>
```jmax = 49``` Maximum number of grid points in the J-direction<br>
```kmax = 2```  Maximum number of grid points in the K-direction<br>
```blocks = 2``` Number of blocks in the I-direction<br>
In order to execute the ```bump.py``` script with python, all commands are written 
in the makefile. Just use the following command to generate grid files in the __grid/__ folder.<br>
```
$make
```
Run the previous command in the __|rootFolder|/run/Tutorial/2DBump/CreateBlocks/__ directory.
<figure>
  <div style="display: flex; justify-content: center;">
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/2DSmoothBumpMesh.png" alt="mesh" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.3 The Mesh for subsonic flow over 2D smooth bump.</figcaption>
  </div>
</figure>

## Setup
To setup the case directory, a python automation script is provided: ```automaton.py```. First setup the most important
parameters, the paths to the grid files and main executable binary to FEST-3D<br>
```RunDir = 'GiveAnyName'```   Name of the run directory to create for current case. <br>
```GridDir= 'CreateBlocks/grid/'```  Path to the folder which only contains grid file. <br>
```NumberOfBlocks = 2``` Total number of blocks. It should match with number of gridfiles avaiable in the __GridDir__ folder <br>
```
AbsBinaryPath="/home/usr/FEST3D/bin/FEST3D"
``` 
Absolute path to the FEST-3D binary. Should be in __|rootfolder|/bin/FEST3D__<br>
Now provide the common control parameters.<br>
```Control['CFL'] = 10.0``` High CFL since implicit time-integration method will be used.<br>
```Control['LoadLevel'] = 0 ```Since simulation will be started from scratch, __0__ is specified.<br>
```Control['MaxIterations'] = 10000``` Maximum number of iteration to perform. <br>
```Control['SaveIterations'] = 1000``` Solution folder will be written every 1000 iteration. <br>
```Control['OutputFileFormat'] = 'tecplot'``` Type of solution file to write. If you have vtk file viewer you can user __vtk__ instead of __tecplot__<br>
```Control['Purge'] = 1``` Only one latest solution folder will be kept in the __time_directories__ and rest will be deleted. <br>
```Control['ResidualWriteInterval'] = 5 ```Write residual in __time_directories/aux/resnorm__ after every 5 iteration. <br>
```Control['Tolerance'] = "1e-13 Continuity_abs"``` Stop the iteration if the absolute residual value of _continuity_ equation is less than 1e-13. <br>
Few scheme parameters for inviscid flow:<br>
```Scheme['InviscidFlux'] = 'slau'``` Inviscid flux-reconstruction shceme. You can use: __ausm__, __ldfss0__, __ausmP__ , and __ausmUP instead of __slau__.<br>
```Scheme['FaceState'] = 'muscl'``` Higher-order  face-state reconstruction method. you can use: __none__, __ppm__, and __weno__. <br>
```Scheme['Limiter'] = '0 0 0  0 0 0'``` Switch on the limiter for I,J,and K direction and switch of pressure based switching for all direction.<br>
```Scheme['TimeStep']='l'``` Local time-stepping method. You can use global time-stepping method also __g__.<br>
```Scheme['TimeIntegration']='implicit'``` LU-SGS matrix-free time integration method.<br>
Now, lets define the flow feature of test case:<br>
```Flow["DensityInf"] = 1.225``` Free-stream density.<br>
```Flow["UInf"] = 170.14``` Free-stream x-component of velcotiy vector. <br>
```Flow["VInf"] = 0.0``` Free-stream y-component of velcotiy vector. <br>
```Flow["WInf"] = 0.0``` Free-stream z-component of velcotiy vector. <br>
```Flow["PressureInf"] = 101325.0``` Free-stream pressure. <br>
```Flow["ReferenceViscosity"] = 0.0``` Set reference viscosity to zero for inviscid flow.<br>

```OutputControl['Out'] = ["Velocity", "Density", "Pressure"]``` Variables to write in the output file.<br>
```ResidualControl['Out'] = ["Mass_abs", "Viscous_abs", "Continuity_abs"]``` Residual to write in the resnorm file. <br>
```BoundaryConditions = [-8, -4, -6, -6, -6, -6]``` Broad boundary condition.[Riemann-inlet, subsonic outlet, rest are Slip-walls]. Free-stream pressure is used to fix the pressure outlet value.<br>

Rest of the variables should be left to their default value. In order to execute this script use following command:
```
$python automaton.py
```
Now you will see a new folder created with ```RunDir``` name. Switch to that directory to run the test case. To make sure setup is correct, check
the ```layout.md``` file located in __system/mesh/layout/layout.md__. The file should look like the following:
```
## BLOCK LAYOUT FILE
## ==========================
## NUMBER OF PROCESSES
2
## NUMBER OF ENTRIES PER PROCESS
9
## PROCESS_NO GRID BC_FILE IMIN IMAX JMIN JMAX KMIN KMAX
## ===================================
## PROCESS 0
00  grid_00.txt  bc_00.md  -008  0001  -006  -006  -006  -006
## PROCESS 1
01  grid_01.txt  bc_01.md  0000  -004  -006  -006  -006  -006
```

The last two uncommented line, in sequence, indicates following:

```
Block_Number    GridFile   Boundary_condition_file   Imin_boundary_condition_number    Imax_boundary_condition_number Jmin_boundary_condition_number Jmax_boundary_condition_number Kmin_boundary_condition_number  Kmax_boundary_condition_number
```

Since there are two blocks, there are two rows of entries, one for each block. All the lines with `#` as the first character are skipped while reading by the FEST-3D solver.
All the positive numbers define the interface boundary condition, and negative numbers define the physical boundary conditions.
Finally, to run the simulation use following command:
```
$nohup bash run.sh &
```
`nohup` helps in avoiding any output on the screen, __&__ execute the last command in the background, allowing you to keep using the terminal.

## Results
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/tutorials/2DSmoothBumpContour.png" alt="Contour of v speed" style="width:350px">
  </div>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <figcaption> Fig.4 Pressure contour for different higher order methods. </figcaption>
  </div>
</figure>
