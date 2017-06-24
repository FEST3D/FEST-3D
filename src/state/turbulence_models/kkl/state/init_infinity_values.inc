!should be included in state file, in subroutine name same as this file name 
         
              call dmsg(1, 'turbulence_models/kkl/state', 'init_infinity_values')
!              tk_inf = 1.5*((0.05**2)*&
!                            (x_speed_inf**2+&
!                             y_speed_inf**2+&
!                             z_speed_inf**2))
!              tw_inf = tk_inf/(mu_ref*0.001)
!
        tk_inf = 9*(1e-9)*sound_speed_inf()**2
        tkl_inf = 1.5589*(1e-6)*(mu_ref*sound_speed_inf()**2)/density_inf
        
!        tk_inf = 1.125*(vel_mag**2)/(vel_mag*density_inf/mu_ref)
!        tw_inf = 125*vel_mag
