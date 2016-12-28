
SCHEME CONTROL FILE
===================

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

### Tubulence: none / sst / 
none

## Time-stepping method: global (g) or local (l)
l

## Higher Order Time Integration: none or RK4
RK4
