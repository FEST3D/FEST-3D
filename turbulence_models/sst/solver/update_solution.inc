!include in solver.f90
            
!            call dmsg(1, 'turbulence_models/sst/solver', 'update_solution')


               tk_temp = tk(i, j, k) - &
                            (( (-1 * tk(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( TKE_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               tw_temp = tw(i, j, k) - &
                            (( (-1 * tw(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( omega_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               tk(i, j, k) = tk_temp
               tw(i, j, k) = tw_temp
