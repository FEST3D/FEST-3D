# This script calculates the free stream pressure, density, etc., using
# the stagnation conditions

import numpy as np

P0 = 101325     # In Pascals
Rho0 = 1.225    # In Kg/m^3

M_inf = 2.5

gamma = 1.4
R_gas = 288.68

T0 = P0/(Rho0 * R_gas)



quant_in_bracket = 1 +(((gamma-1)/2)*(M_inf**2))

T_inf = T0/quant_in_bracket
Rho_inf = Rho0/(quant_in_bracket**(1/(gamma-1)))
P_inf = P0/(quant_in_bracket**(gamma/(gamma-1)))

speed_sound = np.sqrt(gamma*R_gas*T_inf)
X_speed_inf = M_inf * speed_sound

#P_inf = 1013.25
#T_inf = 300
#Rho_inf = P_inf / (R_gas * T_inf)

print 'Rho_inf = ', Rho_inf
print 'X_speed_inf = ', X_speed_inf
print 'P_inf = ', P_inf
print 'T_inf = ', T_inf
