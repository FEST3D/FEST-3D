CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
ausm

## Higher order extension: none or ppm or muscl
ppm

## CFL
0.5

## Time-stepping method: global (g) or local (l)
l

## Higher Order Time extension: none or RK4
RK4

## Tolerance for residue norm comparisons
1e-6

## Grid file
Rchannel.txt

## State Load File ('~' for no load file)
###load.fvtk
~

## Max Iterations
20000

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
1000

## Debug level: Most detail (1) to least detail (5)
5

# FLOW SPECIFIC

## gamma (ratio of specific heats)
1.4

## R\_gas (specific gas constant)
287.

## FREE STREAM PROPERTIES

### Free Stream Density
0.1613

### Free Stream X Speed (at M = 2.4 and T = 300)
567.157

### Free Stream Y Speed
0.

### Free Stream Pressure
5930.316

### Viscous effects
### Give mu reference = 0 for inviscid
### Using Sutherlands law for coefficient of viscosity
### mu reference or mu0 (in kg/ms)
0.0

### T reference or T0 (in K)
273.15

### Sutherland Temperature (K)
110.4

### Prandlt Number
0.7

### Post Processing

### Plot type
Pseudocolor

### Plot variable
Mach_No
