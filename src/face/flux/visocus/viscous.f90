module viscous
  !-----------------------------------------------------------------
  ! The viscous module contains the viscous flux calculations and 
  ! the boundary conditions to be imposed
  !-----------------------------------------------------------------
  !TODO: Viscous: Change to single subroutine for all directions  

  use global     , only: FILE_NAME_LENGTH
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
  use global_vars, only : gradu_x
  use global_vars, only : gradu_y
  use global_vars, only : gradu_z
  use global_vars, only : gradv_x
  use global_vars, only : gradv_y
  use global_vars, only : gradv_z
  use global_vars, only : gradw_x
  use global_vars, only : gradw_y
  use global_vars, only : gradw_z
  use global_vars, only : gradT_x
  use global_vars, only : gradT_y
  use global_vars, only : gradT_z
  use global_vars, only : gradtk_x
  use global_vars, only : gradtk_y
  use global_vars, only : gradtk_z
  use global_vars, only : gradtw_x
  use global_vars, only : gradtw_y
  use global_vars, only : gradtw_z
  use global_vars, only : mu_v=>mu
  use global_vars, only : mu_t
  use global_vars, only : turbulence
  use utils      , only : alloc, dealloc, dmsg
  use string
  use face_interpolant, only: x_density_left 
  use face_interpolant, only: y_density_left
  use face_interpolant, only: z_density_left
  use face_interpolant, only: x_density_right
  use face_interpolant, only: y_density_right
  use face_interpolant, only: z_density_right
  use face_interpolant, only: x_x_speed_left
  use face_interpolant, only: x_y_speed_left
  use face_interpolant, only: x_z_speed_left
  use face_interpolant, only: x_x_speed_right
  use face_interpolant, only: x_y_speed_right
  use face_interpolant, only: x_z_speed_right
  use face_interpolant, only: y_x_speed_left
  use face_interpolant, only: y_y_speed_left
  use face_interpolant, only: y_z_speed_left
  use face_interpolant, only: y_x_speed_right
  use face_interpolant, only: y_y_speed_right
  use face_interpolant, only: y_z_speed_right
  use face_interpolant, only: z_x_speed_left
  use face_interpolant, only: z_y_speed_left
  use face_interpolant, only: z_z_speed_left
  use face_interpolant, only: z_x_speed_right
  use face_interpolant, only: z_y_speed_right
  use face_interpolant, only: z_z_speed_right
  use face_interpolant, only: x_pressure_left
  use face_interpolant, only: y_pressure_left
  use face_interpolant, only: z_pressure_left
  use face_interpolant, only: x_pressure_right
  use face_interpolant, only: y_pressure_right
  use face_interpolant, only: z_pressure_right
  implicit none
  private

  public :: compute_viscous_fluxes

  contains


    subroutine compute_xi_viscous_fluxes(F)
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the xi-face
    ! using a different scheme instead of taking the average of the 
    ! cell centre gradients. The latter leads to odd-even decoupling
    ! problem, which the different scheme corrects.
    !-----------------------------------------------------------------

        implicit none
        
        real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
        real :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                dTdx, dTdy, dTdz
        real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
        real :: T_LE, T_RE
!        real :: T_L, T_R, T_face, K_heat, mu, Qx, Qy, Qz
        real :: K_heat, mu, Qx, Qy, Qz
        real :: uface, vface, wface
        integer :: i, j, k
        real, dimension(:, :, :, :), pointer :: F

!        include "turbulence_models/include/transport/Ftransport_init.inc" 

        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx

           
            ! Gradients at face as average of gradients at cell centres
            dudx = 0.5 * (gradu_x(i-1, j, k) + gradu_x(i, j, k))
            dudy = 0.5 * (gradu_y(i-1, j, k) + gradu_y(i, j, k))
            dudz = 0.5 * (gradu_z(i-1, j, k) + gradu_z(i, j, k))
            dvdx = 0.5 * (gradv_x(i-1, j, k) + gradv_x(i, j, k))
            dvdy = 0.5 * (gradv_y(i-1, j, k) + gradv_y(i, j, k))
            dvdz = 0.5 * (gradv_z(i-1, j, k) + gradv_z(i, j, k))
            dwdx = 0.5 * (gradw_x(i-1, j, k) + gradw_x(i, j, k))
            dwdy = 0.5 * (gradw_y(i-1, j, k) + gradw_y(i, j, k))
            dwdz = 0.5 * (gradw_z(i-1, j, k) + gradw_z(i, j, k))
            dTdx = 0.5 * (gradT_x(i-1, j, k) + gradT_x(i, j, k))
            dTdy = 0.5 * (gradT_y(i-1, j, k) + gradT_y(i, j, k))
            dTdz = 0.5 * (gradT_z(i-1, j, k) + gradT_z(i, j, k))

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
            normal_comp = ( (x_speed(i, j, k) - &
                             x_speed(i-1, j, k)) - &
                          ((dudx * (xc_R - xc_L)) + &
                           (dudy * (yc_R - yc_L)) + &
                           (dudz * (zc_R - zc_L))) &
                          ) / d_LR
            dudx = dudx + (normal_comp * (xc_R - xc_L) / d_LR)
            dudy = dudy + (normal_comp * (yc_R - yc_L) / d_LR)
            dudz = dudz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (y_speed(i, j, k) - &
                             y_speed(i-1, j, k)) - &
                          ((dvdx * (xc_R - xc_L)) + &
                           (dvdy * (yc_R - yc_L)) + &
                           (dvdz * (zc_R - zc_L))) &
                          ) / d_LR
            dvdx = dvdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dvdy = dvdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dvdz = dvdz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (z_speed(i, j, k) - &
                             z_speed(i-1, j, k)) - &
                          ((dwdx * (xc_R - xc_L)) + &
                           (dwdy * (yc_R - yc_L)) + &
                           (dwdz * (zc_R - zc_L))) &
                          ) / d_LR
            dwdx = dwdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dwdy = dwdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dwdz = dwdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! Finding the temperature of left and right element to the
            ! face i, j, k
            T_LE = pressure(i-1, j, k) / (density(i-1, j, k) * R_gas)
            T_RE = pressure(i, j, k) / (density(i, j, k) * R_gas)
            normal_comp = ( (T_RE - T_LE) - &
                          ((dTdx * (xc_R - xc_L)) + &
                           (dTdy * (yc_R - yc_L)) + &
                           (dTdz * (zc_R - zc_L))) &
                          ) / d_LR
            dTdx = dTdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dTdy = dTdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dTdz = dTdz + (normal_comp * (zc_R - zc_L) / d_LR)
            
            ! mu requires T at the face. Hence:
            ! T_L and T_R are the left and right states of the face i,j,k
            ! The values at face used instead of element values
!            T_L = x_pressure_left(i, j, k) / (x_density_left(i, j, k) * R_gas)
!            T_R = x_pressure_right(i, j, k) / (x_density_right(i, j, k) * R_gas)
!            T_face = 0.5 * (T_L + T_R)
!            mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &
!                          Sutherland_temp) / (T_face + Sutherland_temp))
        !   mu = 0.5 * ( (mu_ref * (T_L / T_ref)**1.5 * (T_ref + &
        !                 Sutherland_temp) / (T_L + Sutherland_temp)) + & 
        !                (mu_ref * (T_R / T_ref)**1.5 * (T_ref + &
        !                 Sutherland_temp) / (T_R + Sutherland_temp)) ) 
        !   mu = 2.0

            ! Using lambda = -2 * mu / 3
            ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
            ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
            if (turbulence/='none') then
              mu = 0.5*(mu_v(i-1,j,k) + mu_v(i,j,k)) &
                 + 0.5*(mu_t(i-1,j,k) + mu_t(i,j,k)) 
            else
              mu = 0.5*(mu_v(i-1,j,k) + mu_v(i,j,k))
            end if
            Tau_xx = 2. * mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_yy = 2. * mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_zz = 2. * mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_xy = mu * (dvdx + dudy)
            Tau_xz = mu * (dwdx + dudz)
            Tau_yz = mu * (dwdy + dvdz)

            ! Pr: Prandtl Number
            ! Qx, Qy, Qz: Conduction fluxes
            K_heat = mu * gm * R_gas / ((gm - 1) * Pr)
            Qx = K_heat*dTdx
            Qy = K_heat*dTdy
            Qz = K_heat*dTdz

            ! Note that the xi-direction faces only need the following quantities:
            ! Tau_xx, Tau_xy, Tau_xz -> dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz
            ! Qx -> dTdx
            ! The mass flux has no viscous component
            ! momentum for xi-face:
            F(i, j, k, 2) = F(i, j, k, 2) - ( ((Tau_xx * xnx(i, j, k)) + &
                            (Tau_xy * xny(i, j, k)) + (Tau_xz * xnz(i, j, k))) * &
                            xA(i, j, k))
            F(i, j, k, 3) = F(i, j, k, 3) - ( ((Tau_xy * xnx(i, j, k)) + &
                            (Tau_yy * xny(i, j, k)) + (Tau_yz * xnz(i, j, k))) * &
                            xA(i, j, k))
            F(i, j, k, 4) = F(i, j, k, 4) - ( ((Tau_xz * xnx(i, j, k)) + &
                            (Tau_yz * xny(i, j, k)) + (Tau_zz * xnz(i, j, k))) * &
                            xA(i, j, k))
           
            ! Energy flux
            uface = 0.5 * (x_x_speed_left(i, j, k) + x_x_speed_right(i, j, k))
            vface = 0.5 * (x_y_speed_left(i, j, k) + x_y_speed_right(i, j, k))
            wface = 0.5 * (x_z_speed_left(i, j, k) + x_z_speed_right(i, j, k))
            ! (TijVj + Qi)ni
            F(i, j, k, 5) = F(i, j, k, 5) - (xA(i, j, k) * ( &
                            ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                              Qx) * xnx(i, j, k)) + &
                            ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                              Qy) * xny(i, j, k)) + &
                            ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                              Qz) * xnz(i, j, k)) ) )
            
           ! gradu_y_f(i, j, k) = dTdy 
!        include "turbulence_models/include/transport/Ftransport_find.inc" 
           
          end do
         end do
        end do

    end subroutine compute_xi_viscous_fluxes

    subroutine compute_eta_viscous_fluxes(G)
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the eta-face
    ! using a different scheme instead of taking the average of the 
    ! cell centre gradients. The latter leads to odd-even decoupling
    ! problem, which the different scheme corrects.
    !-----------------------------------------------------------------

        implicit none
        
        real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
        real :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                dTdx, dTdy, dTdz
        real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
        real :: T_LE, T_RE
!        real :: T_L, T_R, T_face, K_heat, mu, Qx, Qy, Qz
        real :: K_heat, mu, Qx, Qy, Qz
        real :: uface, vface, wface
        integer :: i, j, k
        real, dimension(:, :, :, :), pointer :: G

        !include "turbulence_models/include/transport/Gtransport_init.inc" 
        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do k = 1, kmx - 1
         do j = 1, jmx
          do i = 1, imx - 1

            ! Gradients at face as average of gradients at cell centres
            dudx = 0.5 * (gradu_x(i, j-1, k) + gradu_x(i, j, k))
            dudy = 0.5 * (gradu_y(i, j-1, k) + gradu_y(i, j, k))
            dudz = 0.5 * (gradu_z(i, j-1, k) + gradu_z(i, j, k))
            dvdx = 0.5 * (gradv_x(i, j-1, k) + gradv_x(i, j, k))
            dvdy = 0.5 * (gradv_y(i, j-1, k) + gradv_y(i, j, k))
            dvdz = 0.5 * (gradv_z(i, j-1, k) + gradv_z(i, j, k))
            dwdx = 0.5 * (gradw_x(i, j-1, k) + gradw_x(i, j, k))
            dwdy = 0.5 * (gradw_y(i, j-1, k) + gradw_y(i, j, k))
            dwdz = 0.5 * (gradw_z(i, j-1, k) + gradw_z(i, j, k))
            dTdx = 0.5 * (gradT_x(i, j-1, k) + gradT_x(i, j, k))
            dTdy = 0.5 * (gradT_y(i, j-1, k) + gradT_y(i, j, k))
            dTdz = 0.5 * (gradT_z(i, j-1, k) + gradT_z(i, j, k))

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
            normal_comp = ( (x_speed(i, j, k) - &
                             x_speed(i, j-1, k)) - &
                          ((dudx * (xc_R - xc_L)) + &
                           (dudy * (yc_R - yc_L)) + &
                           (dudz * (zc_R - zc_L)) ) & 
                          ) / d_LR
            dudx = dudx + (normal_comp * (xc_R - xc_L) / d_LR)
            dudy = dudy + (normal_comp * (yc_R - yc_L) / d_LR)
            dudz = dudz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (y_speed(i, j, k) - &
                             y_speed(i, j-1, k)) - &
                          ((dvdx * (xc_R - xc_L)) + &
                           (dvdy * (yc_R - yc_L)) + &
                           (dvdz * (zc_R - zc_L))) &
                          ) / d_LR
            dvdx = dvdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dvdy = dvdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dvdz = dvdz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (z_speed(i, j, k) - &
                             z_speed(i, j-1, k)) - &
                          ((dwdx * (xc_R - xc_L)) + &
                           (dwdy * (yc_R - yc_L)) + &
                           (dwdz * (zc_R - zc_L))) &
                          ) / d_LR
            dwdx = dwdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dwdy = dwdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dwdz = dwdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! Finding the temperature of forward and backward element to the
            ! face i, j, k
            T_LE = pressure(i, j-1, k) / (density(i, j-1, k) * R_gas)
            T_RE = pressure(i, j, k) / (density(i, j, k) * R_gas)
            normal_comp = ( (T_RE - T_LE) - &
                          ((dTdx * (xc_R - xc_L)) + &
                           (dTdy * (yc_R - yc_L)) + &
                           (dTdz * (zc_R - zc_L))) &
                          ) / d_LR
            dTdx = dTdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dTdy = dTdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dTdz = dTdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! mu requires T at the face. Hence:
            ! T_L and T_R are the left and right states of the face i,j,k
            ! The values at face used instead of element values
!            T_L = y_pressure_left(i, j, k) / (y_density_left(i, j, k) * R_gas)
!            T_R = y_pressure_right(i, j, k) / (y_density_right(i, j, k) * R_gas)
!            T_face = 0.5 * (T_L + T_R)
!            mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &
!                          Sutherland_temp) / (T_face + Sutherland_temp))
        !   mu = 0.5 * ( (mu_ref * (T_L / T_ref)**1.5 * (T_ref + &
        !                 Sutherland_temp) / (T_L + Sutherland_temp)) + & 
        !                (mu_ref * (T_R / T_ref)**1.5 * (T_ref + &
        !                 Sutherland_temp) / (T_R + Sutherland_temp)) ) 
        !   mu = 2.00

            ! Using lambda = -2 * mu / 3
            ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
            ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
            if (turbulence/='none') then
              mu = 0.5*(mu_v(i,j-1,k) + mu_v(i,j,k)) &
                 + 0.5*(mu_t(i,j-1,k) + mu_t(i,j,k)) 
            else 
              mu = 0.5*(mu_v(i,j-1,k) + mu_v(i,j,k))
            end if

            Tau_xx = 2. * mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_yy = 2. * mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_zz = 2. * mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_xy = mu * (dvdx + dudy)
            Tau_xz = mu * (dwdx + dudz)
            Tau_yz = mu * (dwdy + dvdz)

            
            ! Pr: Prandtl Number
            ! Qx, Qy, Qz: Conduction fluxes
            K_heat = mu * gm * R_gas / ((gm - 1) * Pr)
            Qx = K_heat*dTdx
            Qy = K_heat*dTdy
            Qz = K_heat*dTdz

            ! Note that the xi-direction faces only need the following quantities:
            ! Tau_xx, Tau_xy, Tau_xz -> dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz
            ! Qx -> dTdx
            ! The mass flux has no viscous component
            ! momentum for xi-face:
            G(i, j, k, 2) = G(i, j, k, 2) - ( ((Tau_xx * ynx(i, j, k)) + &
                            (Tau_xy * yny(i, j, k)) + (Tau_xz * ynz(i, j, k))) * &
                            yA(i, j, k))
            G(i, j, k, 3) = G(i, j, k, 3) - ( ((Tau_xy * ynx(i, j, k)) + &
                            (Tau_yy * yny(i, j, k)) + (Tau_yz * ynz(i, j, k))) * &
                            yA(i, j, k))
            G(i, j, k, 4) = G(i, j, k, 4) - ( ((Tau_xz * ynx(i, j, k)) + &
                            (Tau_yz * yny(i, j, k)) + (Tau_zz * ynz(i, j, k))) * &
                            yA(i, j, k))
          
            ! Wall boundary condition
            !TODO MAKE BOUNDARY CONSTION GENERAL NOT ONY AT Y AND Z FACE
!            if ((j .eq. 1) .or. (j .eq. jmx)) then
!            if ((j .eq. 1)) then
!                uface = 0.0
!                vface = 0.0
!                wface = 0.0
!            else
                uface = 0.5 * (y_x_speed_left(i, j, k) + y_x_speed_right(i, j, k))
                vface = 0.5 * (y_y_speed_left(i, j, k) + y_y_speed_right(i, j, k))
                wface = 0.5 * (y_z_speed_left(i, j, k) + y_z_speed_right(i, j, k))
!            end if

            ! Energy flux
            ! (TijVj - Qi)ni
            G(i, j, k, 5) = G(i, j, k, 5) - (yA(i, j, k) * ( &
                            ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                              Qx) * ynx(i, j, k)) + &
                            ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                              Qy) * yny(i, j, k)) + &
                            ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                              Qz) * ynz(i, j, k)) ) )
          !include "turbulence_models/include/transport/Gtransport_find.inc" 
          end do
         end do
        end do


    end subroutine compute_eta_viscous_fluxes

    subroutine compute_zeta_viscous_fluxes(H)
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the zeta-face
    ! using a different scheme instead of taking the average of the 
    ! cell centre gradients. The latter leads to odd-even decoupling
    ! problem, which the different scheme corrects.
    !-----------------------------------------------------------------

        implicit none
        
        real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
        real :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                dTdx, dTdy, dTdz
        real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
        real :: T_LE, T_RE
!        real :: T_L, T_R, T_face, K_heat, mu, Qx, Qy, Qz
        real :: K_heat, mu, Qx, Qy, Qz
        real :: uface, vface, wface
        integer :: i, j, k
        real, dimension(:, :, :, :), pointer :: H

!        include "turbulence_models/include/transport/Htransport_init.inc" 
        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do k = 1, kmx
         do j = 1, jmx - 1
          do i = 1, imx - 1

            ! Gradients at face as average of gradients at cell centres
            dudx = 0.5 * (gradu_x(i, j, k-1) + gradu_x(i, j, k))
            dudy = 0.5 * (gradu_y(i, j, k-1) + gradu_y(i, j, k))
            dudz = 0.5 * (gradu_z(i, j, k-1) + gradu_z(i, j, k))
            dvdx = 0.5 * (gradv_x(i, j, k-1) + gradv_x(i, j, k))
            dvdy = 0.5 * (gradv_y(i, j, k-1) + gradv_y(i, j, k))
            dvdz = 0.5 * (gradv_z(i, j, k-1) + gradv_z(i, j, k))
            dwdx = 0.5 * (gradw_x(i, j, k-1) + gradw_x(i, j, k))
            dwdy = 0.5 * (gradw_y(i, j, k-1) + gradw_y(i, j, k))
            dwdz = 0.5 * (gradw_z(i, j, k-1) + gradw_z(i, j, k))
            dTdx = 0.5 * (gradT_x(i, j, k-1) + gradT_x(i, j, k))
            dTdy = 0.5 * (gradT_y(i, j, k-1) + gradT_y(i, j, k))
            dTdz = 0.5 * (gradT_z(i, j, k-1) + gradT_z(i, j, k))

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
            normal_comp = ( (x_speed(i, j, k) - &
                             x_speed(i, j, k-1)) - &
                          ((dudx * (xc_R - xc_L)) + &
                           (dudy * (yc_R - yc_L)) + &
                           (dudz * (zc_R - zc_L))) &
                          ) / d_LR
            dudx = dudx + (normal_comp * (xc_R - xc_L) / d_LR)
            dudy = dudy + (normal_comp * (yc_R - yc_L) / d_LR)
            dudz = dudz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (y_speed(i, j, k) - &
                             y_speed(i, j, k-1)) - &
                          ((dvdx * (xc_R - xc_L)) + &
                           (dvdy * (yc_R - yc_L)) + &
                           (dvdz * (zc_R - zc_L))) &
                          ) / d_LR
            dvdx = dvdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dvdy = dvdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dvdz = dvdz + (normal_comp * (zc_R - zc_L) / d_LR)

            normal_comp = ( (z_speed(i, j, k) - &
                             z_speed(i, j, k-1)) - &
                          ((dwdx * (xc_R - xc_L)) + &
                           (dwdy * (yc_R - yc_L)) + &
                           (dwdz * (zc_R - zc_L))) &
                          ) / d_LR
            dwdx = dwdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dwdy = dwdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dwdz = dwdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! Finding the temperature of top and bottom element to the
            ! face i, j, k
            T_LE = pressure(i, j, k-1) / (density(i, j, k-1) * R_gas)
            T_RE = pressure(i, j, k) / (density(i, j, k) * R_gas)
            normal_comp = ( (T_RE - T_LE) - &
                          ((dTdx * (xc_R - xc_L)) + &
                           (dTdy * (yc_R - yc_L)) + &
                           (dTdz * (zc_R - zc_L))) &
                          ) / d_LR
            dTdx = dTdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dTdy = dTdy + (normal_comp * (yc_R - yc_L) / d_LR)
            dTdz = dTdz + (normal_comp * (zc_R - zc_L) / d_LR)

            ! mu is required at the face. Hence:
            ! T_L and T_R are the left and right states of the face i,j,k
            ! The values at face used instead of element values
!            T_L = z_pressure_left(i, j, k) / (z_density_left(i, j, k) * R_gas)
!            T_R = z_pressure_right(i, j, k) / (z_density_right(i, j, k) * R_gas)
!            T_face = 0.5 * (T_L + T_R)
!            mu = mu_ref * (T_face / T_ref)**1.5 * ((T_ref + &
!                          Sutherland_temp) / (T_face + Sutherland_temp))
        !   mu = 0.5 * ( (mu_ref * (T_L / T_ref)**1.5 * (T_ref + &
        !                 Sutherland_temp) / (T_L + Sutherland_temp)) + & 
        !                (mu_ref * (T_R / T_ref)**1.5 * (T_ref + &
        !                 Sutherland_temp) / (T_R + Sutherland_temp)) ) 
        !   mu = 2.000

            ! Using lambda = -2 * mu / 3
            ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
            ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
            if (turbulence/='none') then
              mu = 0.5*(mu_v(i,j,k-1) + mu_v(i,j,k)) &
                 + 0.5*(mu_t(i,j,k-1) + mu_t(i,j,k)) 
            else
              mu = 0.5*(mu_v(i,j,k-1) + mu_v(i,j,k))
            end if

            Tau_xx = 2. * mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_yy = 2. * mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_zz = 2. * mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
            Tau_xy = mu * (dvdx + dudy)
            Tau_xz = mu * (dwdx + dudz)
            Tau_yz = mu * (dwdy + dvdz)

            ! Pr: Prandtl Number
            ! Qx, Qy, Qz: Conduction fluxes
            K_heat = mu * gm * R_gas / ((gm - 1) * Pr)
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
            
            ! Wall boundary condition
!            if ((k .eq. 1) .or. (k .eq. kmx)) then
!                uface = 0.0
!                vface = 0.0
!                wface = 0.0
!            else
                uface = 0.5 * (z_x_speed_left(i, j, k) + z_x_speed_right(i, j, k))
                vface = 0.5 * (z_y_speed_left(i, j, k) + z_y_speed_right(i, j, k))
                wface = 0.5 * (z_z_speed_left(i, j, k) + z_z_speed_right(i, j, k))
!            end if

            ! Energy flux
            ! (TijVj - Qi)ni
            H(i, j, k, 5) = H(i, j, k, 5) - (zA(i, j, k) * ( &
                            ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                              Qx) * znx(i, j, k)) + &
                            ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                              Qy) * zny(i, j, k)) + &
                            ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                              Qz) * znz(i, j, k)) ) )
!            include "turbulence_models/include/transport/Htransport_find.inc" 
          end do
         end do
        end do

    end subroutine compute_zeta_viscous_fluxes

    subroutine compute_viscous_fluxes(F, G, H)

        implicit none

        real, dimension(:, :, :, :), pointer :: F, G, H
        
        call compute_xi_viscous_fluxes(F)
        call compute_eta_viscous_fluxes(G)
        call compute_zeta_viscous_fluxes(H)

    end subroutine compute_viscous_fluxes

end module viscous
