!should be included in state file, in subroutine name same as this file name 
         
        call dmsg(1, 'turbulence_models/kkl/state', 'init_infinity_values')
        tk_inf = 9*(1e-9)*(sound_speed_inf()**2)
        tkl_inf = 1.5589*(1e-6)*(mu_ref*sound_speed_inf())/density_inf
