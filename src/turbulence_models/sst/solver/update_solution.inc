!include in solver.f90
            
!            call dmsg(1, 'turbulence_models/sst/solver', 'update_solution')


!               tk_temp = tk(i, j, k) - &
!                            (( (-1 * tk(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( TKE_residue(i, j, k) / density(i, j, k)) ) * &
!                            delta_t(i, j, k) / volume(i, j, k))
!
!               tw_temp = tw(i, j, k) - &
!                            (( (-1 * tw(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( omega_residue(i, j, k) / density(i, j, k)) ) * &
!                            delta_t(i, j, k) / volume(i, j, k))

      beta = beta1*sst_F1(i,j,k) + (1. - sst_F1(i,j,k))*beta2
              tk_temp = tk(i,j,k) - &
                        (1.*delta_t(i,j,k)/volume(i,j,k))&
                       *(-1.*tk(i,j,k)*mass_residue(i,j,k)/density_temp &
                       + 1./(1.+bstar*tw(i,j,k)*delta_t(i,j,k))&
                       * (TKE_residue(i,j,k))/density_temp)

              tw_temp = tw(i,j,k) - &
                        (1.*delta_t(i,j,k)/volume(i,j,k))&
                       *(-1.*tw(i,j,k)*mass_residue(i,j,k)/density_temp&
                       + 1./(1. + 2.*beta*tw(i,j,k)*delta_t(i,j,k))&
                       *(omega_residue(i,j,k))/density_temp)
               tk(i, j, k) = tk_temp
               tw(i, j, k) = tw_temp

               if (tk(i,j,k)<0.) tk(i,j,k) =0.
