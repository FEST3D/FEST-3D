module sst_turbulent_flux

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx
    use global_vars, only : grid_x
    use global_vars, only : grid_y
    use global_vars, only : grid_z

    use global_vars, only : xnx, xny, xnz !face unit normal x
    use global_vars, only : ynx, yny, ynz !face unit normal y
    use global_vars, only : znx, zny, znz !face unit normal z
    use global_vars, only : xA, yA, zA    !face area
    use global_vars, only : vol => volume
    use global_vars, only :   left_ghost_centroid
    use global_vars, only :  right_ghost_centroid
    use global_vars, only :  front_ghost_centroid
    use global_vars, only :   back_ghost_centroid
    use global_vars, only :    top_ghost_centroid
    use global_vars, only : bottom_ghost_centroid
    
    use global, only: FILE_NAME_LENGTH
    use global_vars, only : gm
    use global_vars, only : n_var
    use global_vars, only : R_gas
    use global_vars, only : mu_ref
    use global_vars, only : T_ref
    use global_vars, only : Pr
    use global_vars, only : Sutherland_temp
    use global_vars, only : density
    use global_vars, only : x_speed
    use global_vars, only : y_speed
    use global_vars, only : z_speed
    use global_vars, only : pressure
    use global_vars, only : tk
    use global_vars, only : tw
    use global_vars  ,only :   gradu_x
    use global_vars  ,only :   gradu_y
    use global_vars  ,only :   gradu_z
    use global_vars  ,only :   gradv_x
    use global_vars  ,only :   gradv_y
    use global_vars  ,only :   gradv_z
    use global_vars  ,only :   gradw_x
    use global_vars  ,only :   gradw_y
    use global_vars  ,only :   gradw_z
    use global_vars  ,only :   gradT_x
    use global_vars  ,only :   gradT_y
    use global_vars  ,only :   gradT_z
    use global_vars  ,only :   gradtk_x
    use global_vars  ,only :   gradtk_y
    use global_vars  ,only :   gradtk_z
    use global_vars  ,only :   gradtw_x
    use global_vars  ,only :   gradtw_y
    use global_vars  ,only :   gradtw_z
    use global_vars  ,only :   mu
    use global_vars  ,only :   sst_mu
    use global_vars  ,only :   dist
    use global_vars, only : turbulence
    use global_vars, only : process_id
    use global_sst
    use utils, only: alloc, dealloc, dmsg
    use string
    use face_interpolant, only: x_tk_left, x_tk_right, &
                                y_tk_left, y_tk_right, &
                                z_tk_left, z_tk_right, &
                                x_tw_left, x_tw_right, &
                                y_tw_left, y_tw_right, &
                                z_tw_left, z_tw_right

    implicit none
    private

    public :: compute_sst_fluxes

    contains

    subroutine compute_xi_turbulent_fluxes(F)
        implicit none
        
        real :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
        real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
        real :: mu_f
        real :: sst_f_mu
        real :: omega
        real :: KE
        real :: rho

        real ::  F1
        real ::  sigma_k
        real ::  sigma_w


        integer :: i, j, k
        real, dimension(:, :, :, :), pointer :: F


        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx

           
            ! Gradients at face as average of gradients at cell centres
            dtkdx = 0.5 * (gradtk_x(i-1, j, k) + gradtk_x(i, j, k))
            dtkdy = 0.5 * (gradtk_y(i-1, j, k) + gradtk_y(i, j, k))
            dtkdz = 0.5 * (gradtk_z(i-1, j, k) + gradtk_z(i, j, k))
            dtwdx = 0.5 * (gradtw_x(i-1, j, k) + gradtw_x(i, j, k))
            dtwdy = 0.5 * (gradtw_y(i-1, j, k) + gradtw_y(i, j, k))
            dtwdz = 0.5 * (gradtw_z(i-1, j, k) + gradtw_z(i, j, k))

            if (i .eq. 1) then
                xc_L = left_ghost_centroid(j, k, 1)
                yc_L = left_ghost_centroid(j, k, 2)
                zc_L = left_ghost_centroid(j, k, 3)
            else
                ! Coordinate of left cell centre: element (i-1, j, k)
                xc_L = (grid_x(i-1, j, k) + grid_x(i, j, k) + &
                        grid_x(i, j+1, k) + grid_x(i-1, j+1, k) + &
                        grid_x(i-1, j, k+1) + grid_x(i, j, k+1) + &
                        grid_x(i, j+1, k+1) + grid_x(i-1, j+1, k+1) &
                        ) * 0.125
                yc_L = (grid_y(i-1, j, k) + grid_y(i, j, k) + &
                        grid_y(i, j+1, k) + grid_y(i-1, j+1, k) + &
                        grid_y(i-1, j, k+1) + grid_y(i, j, k+1) + &
                        grid_y(i, j+1, k+1) + grid_y(i-1, j+1, k+1) &
                        ) * 0.125
                zc_L = (grid_z(i-1, j, k) + grid_z(i, j, k) + &
                        grid_z(i, j+1, k) + grid_z(i-1, j+1, k) + &
                        grid_z(i-1, j, k+1) + grid_z(i, j, k+1) + &
                        grid_z(i, j+1, k+1) + grid_z(i-1, j+1, k+1) &
                        ) * 0.125
            end if

            if (i .eq. imx) then
                xc_R = right_ghost_centroid(j, k, 1)
                yc_R = right_ghost_centroid(j, k, 2)
                zc_R = right_ghost_centroid(j, k, 3)
            else           
                ! Correcting the odd-even problem
                ! Coordinate of right cell centre: element (i, j, k)
                xc_R = (grid_x(i, j, k) + grid_x(i+1, j, k) + &
                        grid_x(i+1, j+1, k) + grid_x(i, j+1, k) + &
                        grid_x(i, j, k+1) + grid_x(i+1, j, k+1) + &
                        grid_x(i+1, j+1, k+1) + grid_x(i, j+1, k+1) &
                        ) * 0.125
                yc_R = (grid_y(i, j, k) + grid_y(i+1, j, k) + &
                        grid_y(i+1, j+1, k) + grid_y(i, j+1, k) + &
                        grid_y(i, j, k+1) + grid_y(i+1, j, k+1) + &
                        grid_y(i+1, j+1, k+1) + grid_y(i, j+1, k+1) &
                        ) * 0.125
                zc_R = (grid_z(i, j, k) + grid_z(i+1, j, k) + &
                        grid_z(i+1, j+1, k) + grid_z(i, j+1, k) + &
                        grid_z(i, j, k+1) + grid_z(i+1, j, k+1) + &
                        grid_z(i+1, j+1, k+1) + grid_z(i, j+1, k+1) &
                        ) * 0.125
            end if

            d_LR = sqrt((xc_R - xc_L)**2 + (yc_R - yc_L)**2 + &
                        (zc_R - zc_L)**2)

            ! normal_comp is the component along r_ij
            ! W_j - W_i = W_right_cell - W_left_cell
            !           = W(i, j, k) - W(i-1, j, k)
            ! For this, the values of state variables at the cell
            ! centres are needed

            normal_comp = ( (tk(i, j, k) - &
                             tk(i-1, j, k)) - &
                          ((dtkdx * (xc_R - xc_L)) + &
                           (dtkdy * (yc_R - yc_L)) + &
                           (dtkdz * (zc_R - zc_L))) &
                          ) / d_LR
            dtkdx = dtkdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dtkdy = dtkdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dtkdz = dtkdz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (tw(i, j, k) - &
                             tw(i-1, j, k)) - &
                          ((dtwdx * (xc_R - xc_L)) + &
                           (dtwdy * (yc_R - yc_L)) + &
                           (dtwdz * (zc_R - zc_L))) &
                          ) / d_LR
            dtwdx = dtwdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dtwdy = dtwdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dtwdz = dtwdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! __ blending function calculation __
            mu_f   = 0.5*(     mu(i-1,j,k) +      mu(i,j,k))
            omega = 0.5*(     tw(i-1,j,k) +      tw(i,j,k))
            KE    = 0.5*(     tk(i-1,j,k) +      tk(i,j,k))
            rho   = 0.5*(density(i-1,j,k) + density(i,j,k)) 
            F1 = 0.5*(sst_F1(i-1,j,k) + sst_F1(i,j,k))

            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w2*(1. - F1)
            sst_f_mu = 0.5*(sst_mu(i-1, j, k) + sst_mu(i, j, k))

            F(i, j, k, 6) = F(i, j, k, 6) - (xA(i, j, k)*( &
                            (mu_f + sigma_k*sst_f_mu)*(dtkdx * xnx(i, j, k)&
                              +dtkdy * xny(i, j, k) + dtkdz * xnz(i, j, k))))
            
            F(i, j, k, 7) = F(i, j, k, 7) - (xA(i, j, k)*( &
                            (mu_f + sigma_w*sst_f_mu)*(dtwdx * xnx(i, j, k)&
                              +dtwdy * xny(i, j, k) + dtwdz * xnz(i, j, k))))
           
          end do
         end do
        end do

    end subroutine compute_xi_turbulent_fluxes

    subroutine compute_eta_turbulent_fluxes(G)
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the eta-face
    ! using a different scheme instead of taking the average of the 
    ! cell centre gradients. The latter leads to odd-even decoupling
    ! problem, which the different scheme corrects.
    !-----------------------------------------------------------------

        implicit none
        
        real :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
        real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
        real :: mu_f
        real :: sst_f_mu
        real :: omega
        real :: KE
        real :: rho

        real ::  F1
        real ::  sigma_k
        real ::  sigma_w

        integer :: i, j, k
        real, dimension(:, :, :, :), pointer :: G

        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do k = 1, kmx - 1
         do j = 1, jmx
          do i = 1, imx - 1

            ! Gradients at face as average of gradients at cell centres
            dtkdx = 0.5 * (gradtk_x(i, j-1, k) + gradtk_x(i, j, k))
            dtkdy = 0.5 * (gradtk_y(i, j-1, k) + gradtk_y(i, j, k))
            dtkdz = 0.5 * (gradtk_z(i, j-1, k) + gradtk_z(i, j, k))
            dtwdx = 0.5 * (gradtw_x(i, j-1, k) + gradtw_x(i, j, k))
            dtwdy = 0.5 * (gradtw_y(i, j-1, k) + gradtw_y(i, j, k))
            dtwdz = 0.5 * (gradtw_z(i, j-1, k) + gradtw_z(i, j, k))

            if (j .eq. 1) then
                xc_L = front_ghost_centroid(i, k, 1)
                yc_L = front_ghost_centroid(i, k, 2)
                zc_L = front_ghost_centroid(i, k, 3)
            else
                ! Coordinate of back cell centre: element (i, j-1, k)
                xc_L = (grid_x(i, j-1, k) + grid_x(i+1, j-1, k) + &
                        grid_x(i+1, j, k) + grid_x(i, j, k) + &
                        grid_x(i, j-1, k+1) + grid_x(i+1, j-1, k+1) + &
                        grid_x(i+1, j, k+1) + grid_x(i, j, k+1) &
                        ) * 0.125
                yc_L = (grid_y(i, j-1, k) + grid_y(i+1, j-1, k) + &
                        grid_y(i+1, j, k) + grid_y(i, j, k) + &
                        grid_y(i, j-1, k+1) + grid_y(i+1, j-1, k+1) + &
                        grid_y(i+1, j, k+1) + grid_y(i, j, k+1) &
                        ) * 0.125
                zc_L = (grid_z(i, j-1, k) + grid_z(i+1, j-1, k) + &
                        grid_z(i+1, j, k) + grid_z(i, j, k) + &
                        grid_z(i, j-1, k+1) + grid_z(i+1, j-1, k+1) + &
                        grid_z(i+1, j, k+1) + grid_z(i, j, k+1) &
                        ) * 0.125
            end if
            
            if (j .eq. jmx) then
                xc_R = back_ghost_centroid(i, k, 1)
                yc_R = back_ghost_centroid(i, k, 2)
                zc_R = back_ghost_centroid(i, k, 3)
            else
                ! Correcting the odd-even problem
                ! Coordinate of forward cell centre: element (i, j, k)
                xc_R = (grid_x(i, j, k) + grid_x(i+1, j, k) + &
                        grid_x(i+1, j+1, k) + grid_x(i, j+1, k) + &
                        grid_x(i, j, k+1) + grid_x(i+1, j, k+1) + &
                        grid_x(i+1, j+1, k+1) + grid_x(i, j+1, k+1) &
                        ) * 0.125
                yc_R = (grid_y(i, j, k) + grid_y(i+1, j, k) + &
                        grid_y(i+1, j+1, k) + grid_y(i, j+1, k) + &
                        grid_y(i, j, k+1) + grid_y(i+1, j, k+1) + &
                        grid_y(i+1, j+1, k+1) + grid_y(i, j+1, k+1) &
                        ) * 0.125
                zc_R = (grid_z(i, j, k) + grid_z(i+1, j, k) + &
                        grid_z(i+1, j+1, k) + grid_z(i, j+1, k) + &
                        grid_z(i, j, k+1) + grid_z(i+1, j, k+1) + &
                        grid_z(i+1, j+1, k+1) + grid_z(i, j+1, k+1) &
                        ) * 0.125
            end if

            d_LR = sqrt((xc_R - xc_L)**2 + (yc_R - yc_L)**2 + &
                        (zc_R - zc_L)**2)

            ! normal_comp is the component along r_ij
            ! W_j - W_i = W_front_cell - W_back_cell
            !           = W(i, j, k) - W(i, j-1, k)
            ! For this, the values of state variables at the cell
            ! centres are needed

            normal_comp = ( (tk(i, j, k) - &
                             tk(i, j-1, k)) - &
                          ((dtkdx * (xc_R - xc_L)) + &
                           (dtkdy * (yc_R - yc_L)) + &
                           (dtkdz * (zc_R - zc_L))) &
                          ) / d_LR
            dtkdx = dtkdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dtkdy = dtkdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dtkdz = dtkdz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (tw(i, j, k) - &
                             tw(i, j-1, k)) - &
                          ((dtwdx * (xc_R - xc_L)) + &
                           (dtwdy * (yc_R - yc_L)) + &
                           (dtwdz * (zc_R - zc_L))) &
                          ) / d_LR
            dtwdx = dtwdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dtwdy = dtwdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dtwdz = dtwdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! __ blending function calculation __
        
            mu_f    = 0.5*(     mu(i,j-1,k) +      mu(i,j,k))
            omega = 0.5*(     tw(i,j-1,k) +      tw(i,j,k))
            KE    = 0.5*(     tk(i,j-1,k) +      tk(i,j,k))
            rho   = 0.5*(density(i,j-1,k) + density(i,j,k)) 
            F1 = 0.5*(sst_F1(i,j-1,k)+sst_F1(i,j,k))

            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w2*(1. - F1)
            sst_f_mu = 0.5*(sst_mu(i, j-1, k) + sst_mu(i, j, k))


            G(i, j, k, 6) = G(i, j, k, 6) - (yA(i, j, k)*( &
                            (mu_f + sigma_k*sst_f_mu)*(dtkdx * ynx(i, j, k)&
                              +dtkdy * yny(i, j, k) + dtkdz * ynz(i, j, k))))
            
            G(i, j, k, 7) = G(i, j, k, 7) - (yA(i, j, k)*( &
                            (mu_f + sigma_w*sst_f_mu)*(dtwdx * ynx(i, j, k)&
                              +dtwdy * yny(i, j, k) + dtwdz * ynz(i, j, k))))
          end do
         end do
        end do

    end subroutine compute_eta_turbulent_fluxes

    subroutine compute_zeta_turbulent_fluxes(H)
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the zeta-face
    ! using a different scheme instead of taking the average of the 
    ! cell centre gradients. The latter leads to odd-even decoupling
    ! problem, which the different scheme corrects.
    !-----------------------------------------------------------------

        implicit none
        
        real :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
        real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
        real :: mu_f
        integer :: i, j, k
        real :: sst_f_mu
        real :: omega
        real :: KE
        real :: rho

        real ::  F1
        real ::  sigma_k
        real ::  sigma_w

        real, dimension(:, :, :, :), pointer :: H

        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do k = 1, kmx
         do j = 1, jmx - 1
          do i = 1, imx - 1

            ! Gradients at face as average of gradients at cell centres
            dtkdx = 0.5 * (gradtk_x(i, j, k-1) + gradtk_x(i, j, k))
            dtkdy = 0.5 * (gradtk_y(i, j, k-1) + gradtk_y(i, j, k))
            dtkdz = 0.5 * (gradtk_z(i, j, k-1) + gradtk_z(i, j, k))
            dtwdx = 0.5 * (gradtw_x(i, j, k-1) + gradtw_x(i, j, k))
            dtwdy = 0.5 * (gradtw_y(i, j, k-1) + gradtw_y(i, j, k))
            dtwdz = 0.5 * (gradtw_z(i, j, k-1) + gradtw_z(i, j, k))

            if (k .eq. 1) then
                xc_L = bottom_ghost_centroid(i, j, 1)
                yc_L = bottom_ghost_centroid(i, j, 2)
                zc_L = bottom_ghost_centroid(i, j, 3)
            else
                ! Coordinate of bottom cell centre: element (i, j, k-1)
                xc_L = (grid_x(i, j, k-1) + grid_x(i+1, j, k-1) + &
                        grid_x(i+1, j+1, k-1) + grid_x(i, j+1, k-1) + &
                        grid_x(i, j, k) + grid_x(i+1, j, k) + &
                        grid_x(i+1, j+1, k) + grid_x(i, j+1, k) &
                        ) * 0.125
                yc_L = (grid_y(i, j, k-1) + grid_y(i+1, j, k-1) + &
                        grid_y(i+1, j+1, k-1) + grid_y(i, j+1, k-1) + &
                        grid_y(i, j, k) + grid_y(i+1, j, k) + &
                        grid_y(i+1, j+1, k) + grid_y(i, j+1, k) &
                        ) * 0.125
                zc_L = (grid_z(i, j, k-1) + grid_z(i+1, j, k-1) + &
                        grid_z(i+1, j+1, k-1) + grid_z(i, j+1, k-1) + &
                        grid_z(i, j, k) + grid_z(i+1, j, k) + &
                        grid_z(i+1, j+1, k) + grid_z(i, j+1, k) &
                        ) * 0.125
            end if
            
            if (k .eq. kmx) then
                xc_R = top_ghost_centroid(i, j, 1)
                yc_R = top_ghost_centroid(i, j, 2)
                zc_R = top_ghost_centroid(i, j, 3)
            else
                ! Correcting the odd-even problem
                ! Coordinate of top cell centre: element (i, j, k)
                xc_R = (grid_x(i, j, k) + grid_x(i+1, j, k) + &
                        grid_x(i+1, j+1, k) + grid_x(i, j+1, k) + &
                        grid_x(i, j, k+1) + grid_x(i+1, j, k+1) + &
                        grid_x(i+1, j+1, k+1) + grid_x(i, j+1, k+1) &
                        ) * 0.125
                yc_R = (grid_y(i, j, k) + grid_y(i+1, j, k) + &
                        grid_y(i+1, j+1, k) + grid_y(i, j+1, k) + &
                        grid_y(i, j, k+1) + grid_y(i+1, j, k+1) + &
                        grid_y(i+1, j+1, k+1) + grid_y(i, j+1, k+1) &
                        ) * 0.125
                zc_R = (grid_z(i, j, k) + grid_z(i+1, j, k) + &
                        grid_z(i+1, j+1, k) + grid_z(i, j+1, k) + &
                        grid_z(i, j, k+1) + grid_z(i+1, j, k+1) + &
                        grid_z(i+1, j+1, k+1) + grid_z(i, j+1, k+1) &
                        ) * 0.125
            end if

            d_LR = sqrt((xc_R - xc_L)**2 + (yc_R - yc_L)**2 + &
                        (zc_R - zc_L)**2)

            ! normal_comp is the component along r_ij
            ! W_j - W_i = W_top_cell - W_bottom_cell
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


            ! __ blending function calculation __
            mu_f  = 0.5*(     mu(i,j,k-1) +      mu(i,j,k))
            omega = 0.5*(     tw(i,j,k-1) +      tw(i,j,k))
            KE    = 0.5*(     tk(i,j,k-1) +      tk(i,j,k))
            rho   = 0.5*(density(i,j,k-1) + density(i,j,k)) 
            F1 = 0.5*(sst_F1(i,j,k-1)+sst_F1(i,j,k))

            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w2*(1. - F1)
            sst_f_mu = 0.5*(sst_mu(i, j, k-1) + sst_mu(i, j, k))


            H(i, j, k, 6) = H(i, j, k, 6) - (zA(i, j, k)*( &
                            (mu_f + sigma_k*sst_f_mu)*(dtkdx * znx(i, j, k)&
                              +dtkdy * zny(i, j, k) + dtkdz * znz(i, j, k))))
            
            H(i, j, k, 7) = H(i, j, k, 7) - (zA(i, j, k)*( &
                            (mu_f + sigma_w*sst_f_mu)*(dtwdx * znx(i, j, k)&
                              +dtwdy * zny(i, j, k) + dtwdz * znz(i, j, k))))


          end do
         end do
        end do

    end subroutine compute_zeta_turbulent_fluxes

    subroutine compute_sst_fluxes(F, G, H)

        implicit none

        real, dimension(:, :, :, :), pointer :: F, G, H
       
        call compute_xi_turbulent_fluxes(F)
        call compute_eta_turbulent_fluxes(G)
        call compute_zeta_turbulent_fluxes(H)

    end subroutine compute_sst_fluxes

end module sst_turbulent_flux
 
