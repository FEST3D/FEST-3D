CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
van_leer

## Higher order extension: none or ppm or muscl
none

## CFL
0.5

## Time-stepping method: global (g) or local (l)
l

## Higher Order Time Integration: none or RK4
none

## Tolerance for residue norm comparisons
1e-6

## Grid file
bumpgrid-3D-1row.txt

## State Load File ('~' for no load file)
##state.fvtk
~

## Max Iterations
60000

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
500

## Debug level: Most detail (1) to least detail (5)
5

# FLOW SPECIFIC

## gamma (ratio of specific heats)
1.4

## R\_gas (specific gas constant)
287.

## FREE STREAM PROPERTIES

### Number of variables
5

### Free Stream Density
0.53584

### Free Stream X Speed
403.796

### Free Stream Y Speed
0.

### Free Stream Z Speed
0.

### Free Stream Pressure
31840.4567

### Viscous effects
### Give mu reference = 0 for inviscid
### Using Sutherlands law for coefficient of viscosity
### mu reference or mu0 (in kg/ms)
0.0

### T reference or T0 (in K)
273.15

### Sutherland Temparature (K)
110.4

### Prandtl Number
0.7

### Post Processing

### Plot type
Pseudocolor

### Plot variable
Mach_No
