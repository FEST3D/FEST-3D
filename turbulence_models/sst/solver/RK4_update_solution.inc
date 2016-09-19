!include in solver.f90
                tk_temp = qp_n(i, j, k, 6) - &
                               (((dEdx_1(i, j, k, 6) / 6.0) + &
                                 (dEdx_2(i, j, k, 6) / 3.0) + &
                                 (dEdx_3(i, j, k, 6) / 3.0) + &
                                 (( (-1 * tk(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( TKE_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                tw_temp = qp_n(i, j, k, 7) - &
                               (((dEdx_1(i, j, k, 7) / 6.0) + &
                                 (dEdx_2(i, j, k, 7) / 3.0) + &
                                 (dEdx_3(i, j, k, 7) / 3.0) + &
                                 (( (-1 * tw(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( omega_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
            
                tk(i, j, k) = tk_temp
                tw(i, j, k) = tw_temp
