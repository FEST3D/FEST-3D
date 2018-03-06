!should be included in state file, in subroutine name same as this file name 
         
    call dmsg(1, 'turbulence_models/sa/state', 'init_infinity_values')
        
    tv_inf = 0.005*3*mu_ref/density_inf
