!include in solver.f90
            
!            call dmsg(1, 'turbulence_models/sst/solver', 'update_solution')


              eta = density(i,j,k)*dist(i,j,k)*(sqrt(0.3*tk(i,j,k))/(20*mu(i,j,k)))
              fphi = (1+cd1*eta)/(1+eta**4)

              pressure_temp = pressure_temp - &
                              (-(gm-1)*TKE_residue(i,j,k))*&
                              (delta_t(i,j,k)/volume(i,j,k))
              tk_temp = tk(i,j,k) - &
                        (1.*delta_t(i,j,k)/volume(i,j,k))&
                       *(-1.*tk(i,j,k)*mass_residue(i,j,k)/density(i,j,k) &
                       + 1./(1.+(2.5*((cmu**0.75)*density(i,j,k)*(tk(i,j,k)**1.5)/tkl(i,j,k))&
                     -(2*mu(i,j,k)/dist(i,j,k))*delta_t(i,j,k)))&
                       * (TKE_residue(i,j,k))/density(i,j,k))

              tkl_temp = tw(i,j,k) - &
                        (1.*delta_t(i,j,k)/volume(i,j,k))&
                       *(-1.*tw(i,j,k)*mass_residue(i,j,k)/density(i,j,k)&
                       + 1./(1. -(6*mu(i,j,k)*fphi)*delta_t(i,j,k))&
                       *(kl_residue(i,j,k))/density(i,j,k))

               if (tk_temp>0.) tk(i, j, k) = tk_temp
               if (tkl_temp>0.) tw(i, j, k) = tkl_temp
