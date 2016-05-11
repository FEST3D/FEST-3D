import numpy as np

mu = 0.002

V = 100

d = 0.002

rho = 1.13

T = 300

P = rho * 287 * T

Re = rho * V * d / mu

print 'Re = ', Re
print 'T = ', T
print 'P = ', P

x = 0.04

Re_x = rho * V * x / mu
blt = 4.91 * x / np.sqrt(Re_x)

print 'BLT = ', blt
