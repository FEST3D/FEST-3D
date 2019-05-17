title: Run

# How to run the code to solve a real problem
To set up the run folder, for solving a problem, few 
input files are needed by the FEST-3D solver. To facilitate the process of setting up
all these input files, a Python code is provided at [**Github**](https://github.com/FEST3D/run.git).


## Dependencies
 * **C++11 compiler**
 * **Python 2.7**
 * **bash**

## Inputs
First, we need to create a directory to perform simulation and save all the input/output data.<br>
```RunDir = 'Test'``` Give any name to the Run directory <br>
One of the important input is the grid/mesh files. FEST-3D code requires a separate file for each process.
If the domain is decomposed into 4 blocks, then 4 separate files are required. Although, you can give
any unique name to the grid files, for simplicity the python code expect grid file name in __grid_nn.txt__, where
the nn is the _block-number-1_. So, for 4 blocks, we will use the following: grid_00.txt, grid_01.txt, grid_02.txt, and
grid_03.txt. You should keep all the grid files in a separate folder and mention name of that folder 
at<br>
```GridDir='Mesh'```, here Mesh is the folder in which all the grid files are kept.<br>
```NumberOfBlocks = 1``` Total number of blocks<br>
In order to use the executable build in the binary folder of the FEST-3D code, a soft link is created
between the FEST3D executable ```bin``` folder and FEST3D in ```bin``` folder of Run folder.<br>
```AbsBinaryPath="/home/jatinder/solver/FEST3D/bin/FEST3D"``` provide the absolute path to the binary.<br>
Now, you are required to fix different parameter of the solver based on the problem you are simulating.
Meaning of the most input is self-explanatory from the name.
```Control['CFL'] = 1.0``` Courant–Friedrichs–Lewy number<br>
```Control['LoadLevel'] = 0``` Restart folder number <br>
```Control['MaxIterations'] = 100000``` Maximum number of iteration <br>
```Control['SaveIterations'] = 1000```  Save solution state after every these many iteration <br>
```Control['OutputFileFormat'] = 'vtk'``` Format of the solution output file <br>
```Control['OutputDataFormat'] = 'ASCII'``` Type of the data in the output folder. Only __ASCII__ for now. _BINARY_ will be added in later release <br>
```Control['InputFileFormat'] = 'vtk'``` Format of the solution file from which solution will be restarted <br>
```Control['InputDataFormat'] = 'ASCII'``` Type of the data in the restart file <br>. Similar to output data type, only ASCII is supported for now.<br>
```Control['Precision'] = 6``` Data precision for residual output, not used for solution output.<br>
```Control['Purge'] = 2``` Number of recent solution folder to keep and delete others. __0__ input will keep all the folders<br>
```Control['ResidualWriteInterval'] = 100``` Number of iteration after which to save the residual output in the file<br>
```Control['Tolerance'] = "1e-13 Continuity_abs"``` Tolerance value and variable. The solver will stop once this value is achived. List of tolerace variable that can be used is given below. <br>
```Control['DebugLevel'] = 5```5-Only important information is logged, 1-All the information is logged which helps in debuging. Will be remove in later release.<br>

```Scheme['InviscidFlux'] = 'ausm'``` Scheme to calcualte inviscid fluxes through cell faces.<br>
```Scheme['FaceState'] = 'muscl'``` Scheme for higher-order face-state reconstuction <br>
```Scheme['Limiter'] = '1 1 1  0 0 0'```Switch for limiters and pressure based switching when using higher order face-state reconstuction. Three value for I,j, and k direction.
   1->on  and 0-> off<br>
```Scheme['TurbulenceLimiter'] = '1 1 1'``` Switch for limiters when used for higher-order face-state reconstuctio of turbulent variables.  1->on  and 0-> off<br>
```Scheme['TurbulenceModel']='none'``` Turbulence model<br>
```Scheme['TimeStep']=' g 1e-4'``` Time-step for time-integration. 'l' for local and 'g' for global. In case for global method you can provide exact value to time-step here.<br>
```Scheme['TimeIntegration']='RK4'``` Method for time-integration<br>
```Scheme['HigherOrderBC']='0'``` Higher order boundary conditions.  1->on  and 0-> off.<br>

```Flow["NumberOfVariables"] = 5``` Total number of variables to solver. Reduntant and will be depricated in later release.<br>
```Flow["DensityInf"] = 1.176``` Free-stream density.<br>
```Flow["UInf"] = 173.59``` Free-stream x-component of velocity.<br>
```Flow["VInf"] = 0.0``` Free-stream y-component of velocity.<br>
```Flow["WInf"] = 0.0``` Free-stream z-component of velocity.<br>
```Flow["PressureInf"] = 101325.0``` Free-stream pressure.<br>
```Flow["TurbulenceIntensity"] = 0.0``` Free-stream trubulent intensity in (percentage)<br>
```Flow["ViscosityRatio"] = 0.0``` Free-stream ratio of turbulence viscosity to molecular viscosity.<br>
```Flow["Intermitency"] = 0.0``` Free-stream turbulent intermittency.<br>
```Flow["ReferenceViscosity"] = 0.2049e-3```. Reference viscosity.<br>
```Flow["ViscosityLaw"] = "sutherland_law"``` Law used for viscosity variation<br>
```Flow["ReferenceTemp"] = 300``` Reference temperature for viscosity variation<br>
```Flow["SutherlandTemp"] = 110.5``` Sutherland temperature<br>
```Flow["PrandtlNumbers"] = "0.72 0.9"``` Prandtl number and turbulent prandtl number<br>
```Flow["SpecificHeatRatio"]=1.4``` Specific heat ratio.<br>
```Flow["GasConstant"]=287.0``` Gas Constant.<br>

```OutputControl['Out'] = ["Velocity", "Density", "Pressure", "Mu"]``` Variables to write in the output file<br>
```OutputControl['In'] = ["Velocity", "Density", "Pressure", "Mu"]``` Variables to read in case of restrart<br>
```ResidualControl['Out'] = ["Mass_abs", "Viscous_abs", "Continuity_abs"]``` Residual to write in the resnorm file<br>
```BoundaryConditions = [-3, -4, -5, -8, -6, -6]``` Boundary conditions to used for the six face of the domain<br>

## Expected Inputs
### CFL
 A floating number greater than zero.
### LoadLevel
 Number of any folder present in the time_directories/
### MaxIterations
 Integer, greater than SaveIterations
### SaveIterations
 Integer, lesser than MaxIteratinos
### InputFileFormat 
 *  vtk
 * tecplot
### InputDataFormat
 * ASCII
### OutputFileFormat
 * vtk 
 * tecplot
### OutputDataFormat
 * ASCII
### Precision
 Integer, lesser than 14 and greater than 1
### Purge
 Integer
### ResidualWriteInterval
 Integer greater than 0
### Tolerance
 * Mass_abs
 * Resnorm_abs
 * Viscous_abs
 * Turbulent_abs
 * Continuity_abs
 * X-mom_abs
 * Y-mom_abs
 * Z-mom_abs
 * Energy_abs
 * Mass_rel
 * Resnorm_rel
 * Viscous_rel
 * Turublent_rel
 * Continuity_abs
 * X-mom_rel
 * Y-mom_rel
 * Z-mom_rel
 * Energy_rel
 * TKE_abs
 * Tv_abs
 * Dissipation_abs
 * Omega_abs
 * Kl_abs
 * TKE_rel
 * Tv_rel
 * Dissipation_rel
 * Omega_rel
 * Kl_rel
### DebugLevel 
 1, 2, 3, 4, or 5
### InviscidFlux
 * ausm
 * ldfss0
 * slau
 * van_leer
 * ausmup
 * ausmp
### FaceState 
 * muscl
 * none
 * ppm
 * weno
### Limiter
 * 1 1 1  0 0 0
### TurbulenceLimiter
 1 1 1 
###  TurbulenceModel
 * none
 * sst
 * kkl
 * sa
###  TransitionModel
 * lctm2015
 * bc
### TimeStep
 * g 1e-5
 * l
### TimeIntegration
 * RK4
 * RK2
 * TVDRK2
 * TVDRK3
 * implicit
 * plusgs
 * none
### HigherOrderBC
 0 or 1
### ViscosityLaw
 * sutherland_law
 * constant

### OutputControl[Out]
 * Velocity
 * Density
 * Pressure
 * Mu
 * Mu_t
 * TKE
 * Omega
 * kL
 * tv
 * Wall_distance
 * resnorm
 * TKE_residue
 * Mass_residue
 * X_mom_residue
 * Y_mom_residue
 * Z_mom_residue,
 * energy_residue
 * DuDx,   Dudy,   DuDz
 * DvDx,   DvDy,   DvDz
 * DwDx,   DWDy,   DwDz
 * DTDx,   DTDy,   DTDz
 * DtkDx,  DtkDy,  DtkDz
 * DtwDx,  DtwDy,  DtwDz
 * DtvDx,  DtvDy,  DtvDz
 * DtkLDx, DtkLDy, DtkLDz
    
### OutputControl[In]
 * Velocity
 * Density
 * Pressure
 * viscosity
 * TKE
 * Omega
 * kL
 * tv
    
### ResidualControl['Out'] 
 * Mass_abs
 * Viscous_abs
 * Mass_abs
 * Resnorm_abs
 * Viscous_abs
 * Turbulent_abs
 * Continuity_abs
 * X-mom_abs
 * Y-mom_abs
 * Z-mom_abs
 * Energy_abs
 * Mass_rel
 * Resnorm_rel
 * Viscous_rel
 * Turublent_rel
 * Continuity_abs
 * X-mom_rel
 * Y-mom_rel
 * Z-mom_rel
 * Energy_rel
 * TKE_abs
 * Tv_abs
 * Dissipation_abs
 * Omega_abs
 * Kl_abs
 * TKE_rel
 * Tv_rel
 * Dissipation_rel
 * Omega_rel
 * Kl_rel
]

### Boundary conditions
 * -1:'SUPERSONIC INLET'
 * -2:'SUPERSONIC OUTFLOW'
 * -3:'SUBSONIC INFLOW',
 * -4:'SUBSONIC OUTFLOW', 
 * -5:'WALL', 
 * -6:'SYMMETRY', 
 * -7:'Pole', 
 * -8:'Far-field', 
 * -9:'Total inlet'

```bash 
$python automaton.py
```
@note
Make sure to provide the absolute path of the FEST3D binary in the automaton.py script before executing. And also
the number of files in the GridDir folder should be equal to the number of blocks as input. 

## Directory structure
Executing automation.py will create a directory with usual directory structure:
 * __system__:all the input files including the mesh files and boundary condition file is located in this directory
 * __time\_directory__: all the output files will be stored in this folder.
 * __bin__: a soft link between original FEST-3D binary is stored here.
 * __run.sh__: bash script to run the solver. This helps to remove log clutter on screen and save it in a log file named `out` in the __time\_directory/aux/__ directory along with `resnorm` file which store the residual values. 


## Check Layout file
Although, automaton.py python script tries to handle the boundary condition by its own; it is still not full-proof. So, always the check the __layout.md__ file in the _system/mesh/layout_ directory. Make sure all the boundary condition number are as you expect. In the case of pole boundary condition, some random number will be mentioned and required to change manually to -007. The layout file in explained in the later section.

## Execute
```bash
$mpiexec.hydra -np 16 bin/FEST3D 
```
On linux os you can use following command to run FEST-3D in the background:<br>
```nohup bash run.sh &```<br>
