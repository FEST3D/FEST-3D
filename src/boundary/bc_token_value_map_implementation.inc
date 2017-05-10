if (token == 'FIX_DENSITY') then
    val = BC_FIX_DENSITY
else if (token == 'FIX_X_SPEED') then
    val = BC_FIX_X_SPEED
else if (token == 'FIX_Y_SPEED') then
    val = BC_FIX_Y_SPEED
else if (token == 'FIX_Z_SPEED') then
    val = BC_FIX_Z_SPEED
else if (token == 'FIX_PRESSURE') then
    val = BC_FIX_PRESSURE
else if (token == 'COPY_DENSITY') then
    val = BC_COPY_DENSITY
else if (token == 'COPY_X_SPEED') then
    val = BC_COPY_X_SPEED
else if (token == 'COPY_Y_SPEED') then
    val = BC_COPY_Y_SPEED
else if (token == 'COPY_Z_SPEED') then
    val = BC_COPY_Z_SPEED
else if (token == 'COPY_PRESSURE') then
    val = BC_COPY_PRESSURE
else if (token == 'FLOW_TANGENCY') then
    val = BC_FLOW_TANGENCY
else if (token == 'NO_SLIP') then
    val = BC_NO_SLIP
else if (token == 'INTERFACE') then
    val = BC_INTERFACE
else if (token == 'PERIODIC') then
    val = BC_PERIODIC
    include "turbulence_models/sst/bc/bc_token_value_map_implementation.inc"
else
    val = 0
end if

