!include in viscous (transport) module


           
            ! Gradients at face as average of gradients at cell centres
            dtkdx = 0.5 * (gradtk_x(i, j-1, k) + gradtk_x(i, j, k))
            dtkdy = 0.5 * (gradtk_y(i, j-1, k) + gradtk_y(i, j, k))
            dtkdz = 0.5 * (gradtk_z(i, j-1, k) + gradtk_z(i, j, k))
            dtwdx = 0.5 * (gradtw_x(i, j-1, k) + gradtw_x(i, j, k))
            dtwdy = 0.5 * (gradtw_y(i, j-1, k) + gradtw_y(i, j, k))
            dtwdz = 0.5 * (gradtw_z(i, j-1, k) + gradtw_z(i, j, k))


            ! normal_comp is the component along r_ij
            ! W_j - W_i = W_right_cell - W_left_cell
            !           = W(i, j, k) - W(i, j, k-1)
            ! For this, the values of state variables at the cell
            ! centres are needed
            normal_comp = ( (tk(i, j, k) - &
                             tk(i, j, k-1)) - &
                          ((dtkdx * (xc_R - xc_L)) + &
                           (dtkdy * (yc_R - yc_L)) + &
                           (dtkdz * (zc_R - zc_L))) &
                          ) / d_LR
            dtkdx = dtkdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dtkdy = dtkdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dtkdz = dtkdz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (tw(i, j, k) - &
                             tw(i, j, k-1)) - &
                          ((dtwdx * (xc_R - xc_L)) + &
                           (dtwdy * (yc_R - yc_L)) + &
                           (dtwdz * (zc_R - zc_L))) &
                          ) / d_LR
            dtwdx = dtwdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dtwdy = dtwdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dtwdz = dtwdz + (normal_comp * (zc_R - zc_L) / d_LR)

            sst_f_mu = 0.5*(sst_mu(i, j, k-1) + sst_mu(i, j, k))
            F1_f = 0.5*(F1(i, j, k-1) + F1(i, j, k))

            ! ! ____BLEND CONSTANT____  ! ! 
            g1_w = (beta1/bs) - (sigma_omega1*kappa**2)/sqrt(bs)
            g2_w = (beta2/bs) - (sigma_omega2*kappa**2)/sqrt(bs)
            sigma_k =         F1_f*sigma_k1 + (1. - F1_f)*sigma_k2
            sigma_omega = F1_f*sigma_omega1 + (1. - F1_f)*sigma_omega2
            g_w  =               F1_f*g1_w  + (1. - F1_f)*g2_w
            beta =               F1_f*beta1 + (1. - F1_f)*beta2


            ! Using lambda = -2 * mu / 3
            ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
            ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
            Tau_xx = 2. * sst_f_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_yy = 2. * sst_f_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_zz = 2. * sst_f_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_xy = sst_f_mu * (dvdx + dudy)
            Tau_xz = sst_f_mu * (dwdx + dudz)
            Tau_yz = sst_f_mu * (dwdy + dvdz)

            ! ! ____ TURBULENT THERMAL CONDUCTIVITY ____ ! !
            !TODO different prandtl number for turbulence 
            ! Pr: Prandtl Number
            ! Qx, Qy, Qz: Conduction fluxes
            K_heat = sst_f_mu * gm * R_gas / ((gm - 1) * Pr)
            Qx = K_heat*dTdx
            Qy = K_heat*dTdy
            Qz = K_heat*dTdz


            ! Note that the xi-direction faces only need the following quantities:
            ! Tau_xx, Tau_xy, Tau_xz -> dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz
            ! Qx -> dTdx
            ! The mass flux has no viscous component
            ! momentum for xi-face:
            H(i, j, k, 2) = H(i, j, k, 2) - ( ((Tau_xx * znx(i, j, k)) + &
                            (Tau_xy * zny(i, j, k)) + (Tau_xz * znz(i, j, k))) * &
                            zA(i, j, k))
            H(i, j, k, 3) = H(i, j, k, 3) - ( ((Tau_xy * znx(i, j, k)) + &
                            (Tau_yy * zny(i, j, k)) + (Tau_yz * znz(i, j, k))) * &
                            zA(i, j, k))
            H(i, j, k, 4) = H(i, j, k, 4) - ( ((Tau_xz * znx(i, j, k)) + &
                            (Tau_yz * zny(i, j, k)) + (Tau_zz * znz(i, j, k))) * &
                            zA(i, j, k))
           
            ! Energy flux
            ! (TijVj + Qi)ni
            H(i, j, k, 5) = H(i, j, k, 5) - (zA(i, j, k) * ( &
                            ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                              Qx) * znx(i, j, k)) + &
                            ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                              Qy) * zny(i, j, k)) + &
                            ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                              Qz) * znz(i, j, k)) ) )
            H(i, j, k, 6) = H(i, j, k, 6) - (zA(i, j, k)*( &
                            (mu + sigma_k*sst_f_mu)*(dtkdx * znx(i, j, k)&
                              +dtkdy * zny(i, j, k) + dtkdz * znz(i, j, k))))
            
            H(i, j, k, 7) = H(i, j, k, 7) - (zA(i, j, k)*( &
                            (mu + sigma_omega*sst_f_mu)*(dtwdx * znx(i, j, k)&
                              +dtwdy * zny(i, j, k) + dtwdz * znz(i, j, k))))
