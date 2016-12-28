CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
ausm

## Higher order extension: none or ppm or muscl
muscl

## Higher order extension  switch:
### 0 -> off
### 1 -> on
### limiter /&/  pressure based switch
     1              0

## CFL
3.0

## Time-stepping method: global (g) or local (l)
l

## Higher Order Time Integration: none or RK4
RK4

## Tolerance for residue norm comparisons
1e-6

## Grid file
grid_bump_97x25x4.dat

## State Load File ('0' for no load file)
0

## Max Iterations
1000

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
100

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
1.225

### Free Stream X Speed
170.146995272

### Free Stream Y Speed
0.

### Free Stream Z Speed
0.

### Free Stream Pressure
101325

### Viscous effects
### Give mu reference = 0 for inviscid
### Using Sutherlands law for coefficient of viscosity
### mu reference or mu0 (in kg/ms)
0.0

### T reference or T0 (in K)
273.15

### Sutherland Temparature (K)
120

### Prandtl Number
0.7

### Tubulence: none / sst / 
none

### Post Processing

### Plot type
Pseudocolor

### Plot variable
Mach_No
