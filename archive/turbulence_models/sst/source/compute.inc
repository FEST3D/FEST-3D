! add in source module-  within add_source term residue subroutine
    
!        call dmsg(1, "turbulence_models/sst/source/", "add_source _term_residue")    

            ! ! __VORTICITY__ ! !
            vort = sqrt(2*(&
                            + (gradu_y(i, j, k) - gradv_x(i, j, k))**2 &
                            + (gradv_z(i, j, k) - gradw_y(i, j, k))**2 &
                            + (gradw_x(i, j, k) - gradu_z(i, j, k))**2 &
                          )&
                       )

            ! ! __FACE value of rho turb_kinectic_energy, omega__ ! !
            rho = density(i, j, k)
            KE  = tk(i, j, k)
            omega = tw(i, j, k)

            ! ! Viscosity   ! !
            T = pressure(i, j, k)/(  density(i, j, k) * R_gas )
            mu = mu_ref * (T / T_ref)**1.5 * ((T_ref + & 
                Sutherland_temp) / (T + Sutherland_temp)) 
            nu = mu/rho !TODO change it to accurate form

            CD = (2*rho*sigma_omega2*(1./omega)*(&
                    gradtk_x(i, j,  k)  * gradtw_x(i,  j,  k) &
                  + gradtk_y(i, j,  k)  * gradtw_y(i,  j,  k)))

            arg1 = min(max(sqrt(KE)/(bs*omega*dist(i,j,k)), 500*nu/(omega*dist(i,j,k)**2)),&
                      4*rho*sigma_omega2*KE/(max(CD,1e-20)*dist(i,j,k)**2))

            F1(i,j,k) = tanh(arg1**4)

            ! ! __ FINDING sst_mu__  ! !
            arg2 = max(2.0*sqrt(KE)/(bs*omega*dist(i,j,k)), 500.0*nu/(dist(i,j,k)*omega))
            F2 = tanh(arg2**2)
            sst_mu(i,j,k) = rho*a1*KE/(max(a1*omega,vort*F2))

            ! ! ____BLEND CONSTANT____  ! ! 
            g1_w = (beta1/bs) - (sigma_omega1*kappa**2)/sqrt(bs)
            g2_w = (beta2/bs) - (sigma_omega2*kappa**2)/sqrt(bs)
            sigma_k =         F1(i,j,k)*sigma_k1 + (1. - F1(i,j,k))*sigma_k2
            sigma_omega = F1(i,j,k)*sigma_omega1 + (1. - F1(i,j,k))*sigma_omega2
            g_w  =               F1(i,j,k)*g1_w  + (1. - F1(i,j,k))*g2_w
            beta =               F1(i,j,k)*beta1 + (1. - F1(i,j,k))*beta2


            !! !! ____SOURCE TERM ___ !! !!
                ! ! ____DISSIPATION TERM____ ! !
                D_k = bs*rho*omega*KE
                D_w = beta*rho*omega**2

                ! ! ____ PRODUCTION term____ ! !
                P_k = sst_mu(i,j,k)*vort**2
                P_k = min(P_k,20*bs*rho*omega*KE)
                P_w = min(g_w*rho*P_k*sst_mu(i,j,k), 10.0*D_w)

            lamda = (1. - F1(i,j,k))*CD      !last blend term in omega equation

            S_k = P_k - D_k           !Source term TKE
            S_w = P_w - D_w  +lamda   !source term omega

            S_k = S_k * vol(i, j, k)
            S_w = S_w * vol(i, j, k)

            TKE_residue(i, j, k)   = TKE_residue(i, j, k) + S_k
            omega_residue(i, j, k) = omega_residue(i, j, k) + S_w

