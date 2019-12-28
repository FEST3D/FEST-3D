title: Run

# How to run the code to solve a real problem
To set up the run folder, for solving a problem, few 
input files are needed by the FEST-3D solver. To facilitate the process of setting up
all these input files, a Python code is provided at [**Github**](https://github.com/FEST3D/run.git).


## Dependencies
 * **C++11 compiler**
 * **Python 3**
 * **bash**

## Inputs
@note
A python scipt is provided to ease the user interface with FEST-3D solver.
In order to run the script provide all the input variables described below
and run the script  using:
```
python edit-automaton.py
```
Above command will create a new directory. This new directory host all the input files.
Before runing the solver, change directory to newly created directory with name provide in variable *$RunDir*
eg:<br>
`$python edit-automaton.py` <br>
`$cd <New directory>` <br>
If required, tweak the input files, and after that, run the solver using<br>
`$nohup bash run.sh &` <br>
To check the current run status<br>
`$tail -f time_directories/aux/out` <br>
To plot the residual using gnuplot script `gnplt`<br>
`$gnuplot gnplt`<br>
@endnote
First, we need to create a directory to perform simulation and save all the input/output data.<br>
```RunDir = 'Test'``` Give any name to the Run directory <br>
<hr>
One of the important input is the grid/mesh files. FEST-3D code requires a separate file for each process.
If the domain is decomposed into 4 blocks, then 4 separate files are required. Although, you can give
any unique name to the grid files, for simplicity the python code expect grid file name in __grid_nn.txt__, where
the nn is the _block-number-1_. So, for 4 blocks, we will use the following: grid_00.txt, grid_01.txt, grid_02.txt, and
grid_03.txt. You should keep all the grid files in a separate folder and mention name of that folder 
at<br>
```GridDir='Mesh'```, here Mesh is the folder in which all the grid files are kept.<br>
For more details about the grid/mesh read the subsection [Mesh](./02_mesh.html)
<hr>
```NumberOfBlocks = 1``` Total number of blocks<br>
<hr>
Now, you are required to fix different input parameter of the solver based on the problem you are simulating.
Meaning of the most input is self-explanatory from the name.

| Variable |  Expected Input | Description|
|:---------|:----------------|:----------------------|
|<hr>|<hr>|<hr>|
|Control['CFL'] | Real Number greater than zero | Courant–Friedrichs–Lewy number. Low value (less than 1) for explicit scheme and high value for implicit scheme|
|<hr>|<hr>|<hr>|
|Control['LoadLevel'] | Integer | Restart folder number in the time_directories/ directory |
|<hr>|<hr>|<hr>|
|Control['MaxIterations'] | Integer greater than SaveIterations | Maximum number of iteration |
|<hr>|<hr>|<hr>|
|Control['SaveIterations']| Integer lesser than MaxIterations | Save solution state after every these many iteration |
|<hr>|<hr>|<hr>|
|Control['OutputFileFormat'] | 'vtk' or 'tecplot' | Format of the solution output file |
|<hr>|<hr>|<hr>|
|Control['OutputDataFormat'] | 'ASCII'  | Type of the data in the output folder. Only ASCII for now. _BINARY_ will be added in later release |
|<hr>|<hr>|<hr>|
|Control['InputFileFormat'] | 'vtk' or 'tecplot' | Format of the solution file from which solution will be restarted |
|<hr>|<hr>|<hr>|
|Control['InputDataFormat'] | 'ASCII'| Type of the data in the restart file <br>. Similar to output data type, only ASCII is supported for now.|
|<hr>|<hr>|<hr>|
|Control['Precision'] | Integer, lesser than 14 and greater than 1 | Data precision for residual output, not used for solution output.|
|<hr>|<hr>|<hr>|
|Control['Purge'] | Integer |  Number of recent solution folder to keep and delete others. __0__ input will keep all the folders|
|<hr>|<hr>|<hr>|
|Control['ResidualWriteInterval'] | Integer greater than zero |  Number of iteration after which to save the residual output in the file|
|<hr>|<hr>|<hr>|
|Control['Tolerance'] | Real number and ["Mass_abs", "Viscous_abs", "Mass_abs", "Resnorm_abs", "Viscous_abs", "Turbulent_abs", "Continuity_abs", "X-mom_abs", "Y-mom_abs", "Z-mom_abs", "Energy_abs", "Mass_rel", "Resnorm_rel", "Viscous_rel", "Turublent_rel", "Continuity_abs", "X-mom_rel", "Y-mom_rel", "Z-mom_rel", "Energy_rel", "TKE_abs", "Tv_abs", "Dissipation_abs", "Omega_abs", "Kl_abs", "TKE_rel", "Tv_rel", "Dissipation_rel", "Omega_rel", "Kl_rel"] | Tolerance value and variable. The solver will stop once this value is achived. List of tolerace variables is given in expected input column. You can use only one input at a time.|
|<hr>|<hr>|<hr>|
|Control['DebugLevel'] | 1, 2, 3, 4, or 5 | Print the information about the function/subroutines called in the solver. This helps in debugging the code. 5-Only important information is logged, 1-All the information is logged which helps in debuging. Will be remove in later release|
|<hr>|<hr>|<hr>|
|Scheme['InviscidFlux'] | 'ausm', 'slau', 'ausmUP', or 'ldfss0' | Scheme to calcualte inviscid fluxes through cell faces|
|<hr>|<hr>|<hr>|
|Scheme['FaceState'] | 'none', 'muscl', 'ppm', or 'weno' | Scheme for higher-order face-state reconstuction|
|<hr>|<hr>|<hr>|
|Scheme['Limiter'] | '1 1 1  0 0 0' or '0 0 0 0 0 0'|Switch for limiters and pressure based switching when using higher order face-state reconstuction. Three value for i,j, and k direction 1->on  and 0-> off. Koren limiter is being used with MUSCL scheme.|
|<hr>|<hr>|<hr>|
|Scheme['TurbulenceLimiter'] | '1 1 1' or '0 0 0' | Switch for limiters when used for higher-order face-state reconstuctio of turbulent variables.  1->on  and 0-> off|
|<hr>|<hr>|<hr>|
|Scheme['TurbulenceModel'] | 'none', 'sa', 'sst', or 'sst2003' |  Turbulence model|
|<hr>|<hr>|<hr>|
|Scheme['TransitionModel'] | 'none', 'bc', 'lctm2015' |  Transition model|
|<hr>|<hr>|<hr>|
|Scheme['TimeStep'] | 'l' or 'g [optional time step]' | Time-step for time-integration. 'l' for local and 'g' for global. In case for global method you can provide exact value to time-step here.|
|<hr>|<hr>|<hr>|
|Scheme['TimeIntegration'] | 'none', 'RK2', 'RK4', 'TVDRK2', 'TVDRK3', 'implicit', or 'plusgs' | Method for time-integration|
|<hr>|<hr>|<hr>|
|Scheme['HigherOrderBC'] | 0 or 1 |  Higher order boundary conditions.  1->on  and 0-> off.|
|<hr>|<hr>|<hr>|
|Flow["NumberOfVariables"] | 5 | Total number of variables to solver. Reduntant and will be depricated in later release. |
|<hr>|<hr>|<hr>|
|Flow["DensityInf"] | Real Number | Free-stream density |
|<hr>|<hr>|<hr>|
|Flow["UInf"] | Real Number | Free-stream x-component of velocity|
|<hr>|<hr>|<hr>|
|Flow["VInf"] | Real Number | Free-stream y-component of velocity|
|<hr>|<hr>|<hr>|
|Flow["WInf"] | Real Number | Free-stream z-component of velocity|
|<hr>|<hr>|<hr>|
|Flow["PressureInf"] | Real Number | Free-stream pressure|
|<hr>|<hr>|<hr>|
|Flow["TurbulenceIntensity"] | Real Number | Free-stream trubulent intensity in (percentage)|
|<hr>|<hr>|<hr>|
|Flow["ViscosityRatio"] | Real Number | Free-stream ratio of turbulence viscosity to molecular viscosity|
|<hr>|<hr>|<hr>|
|Flow["Intermitency"] | Real Number | Free-stream turbulent intermittency|
|<hr>|<hr>|<hr>|
|Flow["ReferenceViscosity"] | Real Number | Reference viscosity|
|<hr>|<hr>|<hr>|
|Flow["ViscosityLaw"] | 'sutherland_law' or 'constant' |  Law used for viscosity variation|
|<hr>|<hr>|<hr>|
|Flow["ReferenceTemp"] | Real Number | Reference temperature for viscosity variation|
|<hr>|<hr>|<hr>|
|Flow["SutherlandTemp"] | Real Number | Sutherland temperature|
|<hr>|<hr>|<hr>|
|Flow["PrandtlNumbers"] | Two real numbers | Prandtl number and turbulent prandtl number|
|<hr>|<hr>|<hr>|
|Flow["SpecificHeatRatio"] | Real number | Specific heat ratio|
|<hr>|<hr>|<hr>|
|Flow["GasConstant"] | Real | Gas Constant|
|<hr>|<hr>|<hr>|
|OutputControl['Out'] | [ "Velocity" , "Density" , "Pressure" , "Mu" , "Mu_t" , "TKE" , "Omega" , "kL" , "tv" , "Wall_distance" , "DuDx",   "Dudy",   "DuDz" , "DvDx",   "DvDy",   "DvDz" , "DwDx",   "DWDy",   "DwDz" , "DTDx",   "DTDy",   "DTDz" , "DtkDx",  "DtkDy",  "DtkDz" , "DtwDx",  "DtwDy",  "DtwDz" , "DtvDx",  "DtvDy",  "DtvDz" , "DtkLDx", "DtkLDy", "DtkLDz"]| Variables to write in the output file|
|<hr>|<hr>|<hr>|
|OutputControl['In'] | ["Velocity" ,"Density" ,"Pressure" ,"viscosity" ,"TKE" ,"Omega" ,"kL" ,"tv"]  | Variables to read in case of restrart|
|<hr>|<hr>|<hr>|
|ResidualControl['Out'] | ["Mass_abs", "Viscous_abs", "Mass_abs", "Resnorm_abs", "Viscous_abs", "Turbulent_abs", "Continuity_abs", "X-mom_abs", "Y-mom_abs", "Z-mom_abs", "Energy_abs", "Mass_rel", "Resnorm_rel", "Viscous_rel", "Turublent_rel", "Continuity_abs", "X-mom_rel", "Y-mom_rel", "Z-mom_rel", "Energy_rel", "TKE_abs", "Tv_abs", "Dissipation_abs", "Omega_abs", "Kl_abs", "TKE_rel", "Tv_rel", "Dissipation_rel", "Omega_rel", "Kl_rel"] | Residual to write in the resnorm file|
|<hr>|<hr>|<hr>|
|BoundaryConditions | [-3, -4, -5, -8, -6, -6]  where <-1:'SUPERSONIC INLET', -2:'SUPERSONIC OUTFLOW', -3:'SUBSONIC INFLOW', -4:'SUBSONIC OUTFLOW', -5:'WALL', -6:'SYMMETRY', -7:'Pole', -8:'Far-field', -11:'Total inlet'> | Boundary conditions to used for the six face of the domain|
|<hr>|<hr>|<hr>|


```bash 
$python edit-automaton.py
```
@note
Make sure to provide the absolute path of the FEST3D binary in the edit-automaton.py script before executing. And also
the number of files in the GridDir folder should be equal to the number of blocks as input. 

## Directory structure
Executing edit-automaton.py will create a directory with usual directory structure:
###
 *  **system**:all the input files including the mesh files and boundary condition file is located in this directory
 *  **time\_directory**: all the output files will be stored in this folder.
 *  **bin**: a soft link between original FEST-3D binary is stored here.
 *  **run.sh**: bash script to run the solver. This helps to remove log clutter on screen and save it in a log file named `out` in the __time\_directory/aux/__ directory along with `resnorm` file which store the residual values. 


## Check Layout file
Although, edit-automaton.py python script tries to handle the boundary condition by its own; it is still not full-proof. So, always the check the __layout.md__ file in the _system/mesh/layout_ directory. Make sure all the boundary condition number are as you expect. In the case of pole boundary condition, some random number will be mentioned and required to change manually to -007. The layout file in explained in the later section.

## Execute
```bash
$mpiexec.hydra -np 16 bin/FEST3D 
```
On linux os you can use following command to run FEST-3D in the background:<br>
```nohup bash run.sh &```<br>

## Check Status
The screen output is directed to the text file: ```[RunDir]/time_directories/aux/out```<br>
The residual are stored in the text file: ```[RunDir]/time_directories/aux/resnorm```<br>
You can output the current status or Current iteration number of the run on screen using
```
tail -f time_directories/aux/out
```
or you can check the residual by using Gnuplot or other similar software:
```
gnuplot gnplt
```
```gnplt``` is the script provided for the Gnuplot software. You can download Gnuplot using:
```
sudo apt-get install gnuplot
```
Use ctl-c to stop the Gnuplot.


## Post Processing
Since FEST-3D output solution file in Tecplot and VTK format, you can use either commerical TECPLOT software or open-source softwares like: <a href="https://wci.llnl.gov/simulation/computer-codes/visit/downloads" target="_blank">Visit</a>, <a href="https://www.paraview.org/download/" target="_blank">Paraview</a>, <a href="https://docs.enthought.com/mayavi/mayavi/installation.html" target="_blank">Mayavi</a>.
