module viscous
  !-----------------------------------------------------------------
  ! The viscous module contains the viscous flux calculations and 
  ! the boundary conditions to be imposed
  !-----------------------------------------------------------------
  !TODO: Viscous: Change to single subroutine for all directions  
#include "error.inc"
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
  
  use global_vars, only : process_id
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
  use global_vars, only : tkl
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
  use global_vars, only : gradtkl_x
  use global_vars, only : gradtkl_y
  use global_vars, only : gradtkl_z
  use global_vars, only : mu
  use global_vars, only : mu_t
  use global_vars, only : kkl_mu
  use global_vars, only : sst_mu
  use global_vars, only : turbulence
  use global_sst , only : sst_F1
  use global_sst , only : sigma_k1
  use global_sst , only : sigma_k2
  use global_sst , only : sigma_w1
  use global_sst , only : sigma_w2
  use global_kkl , only : sigma_k
  use global_kkl , only : sigma_phi
  use utils      , only : alloc, dealloc, dmsg
  use string
  implicit none
  private

  public :: compute_viscous_fluxes

  contains

    subroutine compute_viscous_fluxes(F, G, H)

        implicit none

        real, dimension(:, :, :, :), pointer :: F, G, H
        
        select case(trim(turbulence))
          case('none')
            call compute_xi_viscous_fluxes_laminar(F)
            call compute_eta_viscous_fluxes_laminar(G)
            if(kmx==2)then
              continue
            else
              call compute_zeta_viscous_fluxes_laminar(H)
            end if
          case('sst', 'sst2003')
            call compute_xi_viscous_fluxes_sst(F)
            call compute_eta_viscous_fluxes_sst(G)
            if(kmx==2)then
              continue
            else
              call compute_zeta_viscous_fluxes_sst(H)
            end if
          case('kkl')
            call compute_xi_viscous_fluxes_kkl(F)
            call compute_eta_viscous_fluxes_kkl(G)
            if(kmx==2)then
              continue
            else
              call compute_zeta_viscous_fluxes_kkl(H)
            end if
          case DEFAULT
            Fatal_error
        end select


            if (any(isnan(G))) then
              Fatal_error
            end if
            if (any(isnan(F))) then
              Fatal_error
            end if
            if (any(isnan(H))) then
              Fatal_error
            end if

    end subroutine compute_viscous_fluxes

    subroutine compute_xi_viscous_fluxes_laminar(F)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: F

      !---------------------------------------------------------------------
      ! Calculating the fluxes at the faces
      ! Calculating for the interior xi-faces
      !--------------------------------------------------------------------
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx

          !--- FACE Gradients ---!
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
          !--- For ODD-EVEN coupling error ---!
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

          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz

              mu_f = 0.5*(mu(i-1,j,k) + mu(i,j,k))
              mut_f = 0. 

          !--- end of ODD-EVEN coupling correction ---!
          total_mu = mu_f + mut_f
          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9)* gm * R_gas / (gm - 1)
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
          uface = 0.5 * (x_speed(i-1, j, k) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i-1, j, k) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i-1, j, k) + z_speed(i, j, k))
          ! (TijVj + Qi)ni
          F(i, j, k, 5) = F(i, j, k, 5) - (xA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * xnx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * xny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * xnz(i, j, k)) ) )
          
         
        end do
       end do
      end do
    end subroutine compute_xi_viscous_fluxes_laminar


    subroutine compute_xi_viscous_fluxes_sst(F)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: F

      !--- sst variable requirement ---!
      real :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
      real ::  F1
      real ::  sigma_kf
      real ::  sigma_wf

      !---------------------------------------------------------------------
      ! Calculating the fluxes at the faces
      ! Calculating for the interior xi-faces
      !--------------------------------------------------------------------
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx

          !--- FACE Gradients ---!
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
          dtkdx = 0.5 * (gradtk_x(i-1, j, k) + gradtk_x(i, j, k))
          dtkdy = 0.5 * (gradtk_y(i-1, j, k) + gradtk_y(i, j, k))
          dtkdz = 0.5 * (gradtk_z(i-1, j, k) + gradtk_z(i, j, k))
          dtwdx = 0.5 * (gradtw_x(i-1, j, k) + gradtw_x(i, j, k))
          dtwdy = 0.5 * (gradtw_y(i-1, j, k) + gradtw_y(i, j, k))
          dtwdz = 0.5 * (gradtw_z(i-1, j, k) + gradtw_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
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

          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz

          mu_f     = 0.5*(     mu(i-1,j,k) +      mu(i,j,k))
          F1       = 0.5*(sst_F1(i-1,j,k) + sst_F1(i,j,k))
          sigma_kf  =    sigma_k1*F1  +    sigma_k2*(1. - F1)
          sigma_wf  =    sigma_w1*F1  +    sigma_w2*(1. - F1)
          mut_f    = 0.5*(sst_mu(i-1, j, k) + sst_mu(i, j, k))

          !--- end of ODD-EVEN coupling correction ---!
          total_mu = mu_f + mut_f
          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! k in reynolds stress
          Tau_xx = Tau_xx-((density(i-1,j,k)+density(i,j,k))*(tk(i-1,j,k)+tk(i,j,k)))/6.
          Tau_yy = Tau_yy-((density(i-1,j,k)+density(i,j,k))*(tk(i-1,j,k)+tk(i,j,k)))/6.
          Tau_zz = Tau_zz-((density(i-1,j,k)+density(i,j,k))*(tk(i-1,j,k)+tk(i,j,k)))/6.
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9)* gm * R_gas / (gm - 1)
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
          uface = 0.5 * (x_speed(i-1, j, k) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i-1, j, k) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i-1, j, k) + z_speed(i, j, k))
          ! (TijVj + Qi)ni
          F(i, j, k, 5) = F(i, j, k, 5) - (xA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * xnx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * xny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * xnz(i, j, k)) ) )
          F(i, j, k, 6) = F(i, j, k, 6) - (xA(i, j, k)*( &
                          (mu_f + sigma_kf*mut_f)*(dtkdx * xnx(i, j, k)&
                            +dtkdy * xny(i, j, k) + dtkdz * xnz(i, j, k))))
          
          F(i, j, k, 7) = F(i, j, k, 7) - (xA(i, j, k)*( &
                          (mu_f + sigma_wf*mut_f)*(dtwdx * xnx(i, j, k)&
                            +dtwdy * xny(i, j, k) + dtwdz * xnz(i, j, k))))
         
        end do
       end do
      end do
    end subroutine compute_xi_viscous_fluxes_sst


    subroutine compute_xi_viscous_fluxes_kkl(F)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: F
      !--- kkl variable requirement  ---!
      real :: dtkdx, dtkdy, dtkdz
      real :: dtkldx, dtkldy, dtkldz


      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx

          !--- FACE Gradients ---!
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
          dtkdx  = 0.5 * (gradtk_x(i-1, j, k)  + gradtk_x(i, j, k))
          dtkdy  = 0.5 * (gradtk_y(i-1, j, k)  + gradtk_y(i, j, k))
          dtkdz  = 0.5 * (gradtk_z(i-1, j, k)  + gradtk_z(i, j, k))
          dtkldx = 0.5 * (gradtkl_x(i-1, j, k) + gradtkl_x(i, j, k))
          dtkldy = 0.5 * (gradtkl_y(i-1, j, k) + gradtkl_y(i, j, k))
          dtkldz = 0.5 * (gradtkl_z(i-1, j, k) + gradtkl_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
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
          normal_comp = ( (tk(i, j, k) - &
                           tk(i-1, j, k)) - &
                        ((dtkdx * (xc_R - xc_L)) + &
                         (dtkdy * (yc_R - yc_L)) + &
                         (dtkdz * (zc_R - zc_L))) &
                        ) / d_LR
          dtkdx = dtkdx + (normal_comp * (xc_R - xc_L) / d_LR)
          dtkdy = dtkdy + (normal_comp * (yc_R - yc_L) / d_LR)
          dtkdz = dtkdz + (normal_comp * (zc_R - zc_L) / d_LR)

          normal_comp = ( (tkl(i, j, k) - &
                           tkl(i-1, j, k)) - &
                        ((dtkldx * (xc_R - xc_L)) + &
                         (dtkldy * (yc_R - yc_L)) + &
                         (dtkldz * (zc_R - zc_L))) &
                        ) / d_LR
          dtkldx = dtkldx + (normal_comp * (xc_R - xc_L) / d_LR)
          dtkldy = dtkldy + (normal_comp * (yc_R - yc_L) / d_LR)
          dtkldz = dtkldz + (normal_comp * (zc_R - zc_L) / d_LR)

          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz

          mu_f   = 0.5*(     mu(i-1,j,k) +      mu(i,j,k))
          mut_f  = 0.5*(kkl_mu(i-1, j, k) + kkl_mu(i, j, k))

          !--- end of ODD-EVEN coupling correction ---!
          total_mu = mu_f + mut_f
          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! k in reynolds stress
          Tau_xx = Tau_xx-((density(i-1,j,k)+density(i,j,k))*(tk(i-1,j,k)+tk(i,j,k)))/6.
          Tau_yy = Tau_yy-((density(i-1,j,k)+density(i,j,k))*(tk(i-1,j,k)+tk(i,j,k)))/6.
          Tau_zz = Tau_zz-((density(i-1,j,k)+density(i,j,k))*(tk(i-1,j,k)+tk(i,j,k)))/6.
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9)* gm * R_gas / (gm - 1)
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
          uface = 0.5 * (x_speed(i-1, j, k) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i-1, j, k) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i-1, j, k) + z_speed(i, j, k))
          ! (TijVj + Qi)ni
          F(i, j, k, 5) = F(i, j, k, 5) - (xA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * xnx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * xny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * xnz(i, j, k)) ) )
          F(i, j, k, 6) = F(i, j, k, 6) - (xA(i, j, k)*( &
                          (mu_f + sigma_k*mut_f)*(dtkdx * xnx(i, j, k)&
                            +dtkdy * xny(i, j, k) + dtkdz * xnz(i, j, k))))
          
          F(i, j, k, 7) = F(i, j, k, 7) - (xA(i, j, k)*( &
                          (mu_f + sigma_phi*mut_f)*(dtkldx * xnx(i, j, k)&
                            +dtkldy * xny(i, j, k) + dtkldz * xnz(i, j, k))))
         
        end do
       end do
      end do
    end subroutine compute_xi_viscous_fluxes_kkl

  !---------- xi face calculation end ----------------!


  !---------- eta face calculation begin ----------------!

    subroutine compute_eta_viscous_fluxes_laminar(G)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz 
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: G


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
          !--- For ODD-EVEN coupling error ---!
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
          !--- end of ODD-EVEN coupling correction ---!


          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz

          mu_f = 0.5*(mu(i,j-1,k) + mu(i,j,k))
          mut_f = 0. 
          total_mu = mu_f + mut_f

          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9) * gm * R_gas / (gm - 1)
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
        
          uface = 0.5 * (x_speed(i, j-1, k) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i, j-1, k) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i, j-1, k) + z_speed(i, j, k))

          ! Energy flux
          ! (TijVj - Qi)ni
          G(i, j, k, 5) = G(i, j, k, 5) - (yA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * ynx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * yny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * ynz(i, j, k)) ) )
        end do
       end do
      end do
    end subroutine compute_eta_viscous_fluxes_laminar


    subroutine compute_eta_viscous_fluxes_sst(G)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz 
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: G

      !--- sst variable requirement ---!
      real :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
      real ::  F1
      real ::  sigma_kf
      real ::  sigma_wf

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
          dtkdx = 0.5 * (gradtk_x(i, j-1, k) + gradtk_x(i, j, k))
          dtkdy = 0.5 * (gradtk_y(i, j-1, k) + gradtk_y(i, j, k))
          dtkdz = 0.5 * (gradtk_z(i, j-1, k) + gradtk_z(i, j, k))
          dtwdx = 0.5 * (gradtw_x(i, j-1, k) + gradtw_x(i, j, k))
          dtwdy = 0.5 * (gradtw_y(i, j-1, k) + gradtw_y(i, j, k))
          dtwdz = 0.5 * (gradtw_z(i, j-1, k) + gradtw_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
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

          !--- end of ODD-EVEN coupling correction ---!


          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz

          mu_f     = 0.5*(     mu(i,j-1,k) +      mu(i,j,k))
          F1       = 0.5*(sst_F1(i,j-1,k) + sst_F1(i,j,k))
          sigma_kf  =    sigma_k1*F1  +    sigma_k2*(1. - F1)
          sigma_wf  =    sigma_w1*F1  +    sigma_w2*(1. - F1)
          mut_f    = 0.5*(sst_mu(i, j-1, k) + sst_mu(i, j, k))

          total_mu = mu_f + mut_f

          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! k in reynolds stress
          Tau_xx = Tau_xx-((density(i,j-1,k)+density(i,j,k))*(tk(i,j-1,k)+tk(i,j,k)))/6.
          Tau_yy = Tau_yy-((density(i,j-1,k)+density(i,j,k))*(tk(i,j-1,k)+tk(i,j,k)))/6.
          Tau_zz = Tau_zz-((density(i,j-1,k)+density(i,j,k))*(tk(i,j-1,k)+tk(i,j,k)))/6.
          
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9) * gm * R_gas / (gm - 1)
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
        
          uface = 0.5 * (x_speed(i, j-1, k) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i, j-1, k) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i, j-1, k) + z_speed(i, j, k))

          ! Energy flux
          ! (TijVj - Qi)ni
          G(i, j, k, 5) = G(i, j, k, 5) - (yA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * ynx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * yny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * ynz(i, j, k)) ) )
          G(i, j, k, 6) = G(i, j, k, 6) - (yA(i, j, k)*( &
                          (mu_f + sigma_kf*mut_f)*(dtkdx * ynx(i, j, k)&
                            +dtkdy * yny(i, j, k) + dtkdz * ynz(i, j, k))))
          
          G(i, j, k, 7) = G(i, j, k, 7) - (yA(i, j, k)*( &
                          (mu_f + sigma_wf*mut_f)*(dtwdx * ynx(i, j, k)&
                            +dtwdy * yny(i, j, k) + dtwdz * ynz(i, j, k))))
        end do
       end do
      end do
    end subroutine compute_eta_viscous_fluxes_sst


    subroutine compute_eta_viscous_fluxes_kkl(G)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz 
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: G

      !--- kkl variable requirement  ---!
      real :: dtkdx, dtkdy, dtkdz
      real :: dtkldx, dtkldy, dtkldz

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
          dtkdx = 0.5 * (gradtk_x(i, j-1, k) + gradtk_x(i, j, k))
          dtkdy = 0.5 * (gradtk_y(i, j-1, k) + gradtk_y(i, j, k))
          dtkdz = 0.5 * (gradtk_z(i, j-1, k) + gradtk_z(i, j, k))
          dtkldx = 0.5 * (gradtkl_x(i, j-1, k) + gradtkl_x(i, j, k))
          dtkldy = 0.5 * (gradtkl_y(i, j-1, k) + gradtkl_y(i, j, k))
          dtkldz = 0.5 * (gradtkl_z(i, j-1, k) + gradtkl_z(i, j, k))
          !--- For ODD-EVEN coupling error ---!
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
          normal_comp = ( (tk(i, j, k) - &
                           tk(i, j-1, k)) - &
                        ((dtkdx * (xc_R - xc_L)) + &
                         (dtkdy * (yc_R - yc_L)) + &
                         (dtkdz * (zc_R - zc_L))) &
                        ) / d_LR
          dtkdx = dtkdx + (normal_comp * (xc_R - xc_L) / d_LR)
          dtkdy = dtkdy + (normal_comp * (yc_R - yc_L) / d_LR)
          dtkdz = dtkdz + (normal_comp * (zc_R - zc_L) / d_LR)

          normal_comp = ( (tkl(i, j, k) - &
                           tkl(i, j-1, k)) - &
                        ((dtkldx * (xc_R - xc_L)) + &
                         (dtkldy * (yc_R - yc_L)) + &
                         (dtkldz * (zc_R - zc_L))) &
                        ) / d_LR
          dtkldx = dtkldx + (normal_comp * (xc_R - xc_L) / d_LR)
          dtkldy = dtkldy + (normal_comp * (yc_R - yc_L) / d_LR)
          dtkldz = dtkldz + (normal_comp * (zc_R - zc_L) / d_LR)

          !--- end of ODD-EVEN coupling correction ---!


          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz

          mu_f    = 0.5*(     mu(i,j-1,k) +      mu(i,j,k))
          mut_f = 0.5*(kkl_mu(i, j-1, k) + kkl_mu(i, j, k))
          total_mu = mu_f + mut_f

          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! k in reynolds stress
          Tau_xx = Tau_xx-((density(i,j-1,k)+density(i,j,k))*(tk(i,j-1,k)+tk(i,j,k)))/6.
          Tau_yy = Tau_yy-((density(i,j-1,k)+density(i,j,k))*(tk(i,j-1,k)+tk(i,j,k)))/6.
          Tau_zz = Tau_zz-((density(i,j-1,k)+density(i,j,k))*(tk(i,j-1,k)+tk(i,j,k)))/6.
          
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9) * gm * R_gas / (gm - 1)
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
        
          uface = 0.5 * (x_speed(i, j-1, k) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i, j-1, k) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i, j-1, k) + z_speed(i, j, k))

          ! Energy flux
          ! (TijVj - Qi)ni
          G(i, j, k, 5) = G(i, j, k, 5) - (yA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * ynx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * yny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * ynz(i, j, k)) ) )
          G(i, j, k, 6) = G(i, j, k, 6) - (yA(i, j, k)*( &
                          (mu_f + sigma_k*mut_f)*(dtkdx * ynx(i, j, k)&
                            +dtkdy * yny(i, j, k) + dtkdz * ynz(i, j, k))))
          
          G(i, j, k, 7) = G(i, j, k, 7) - (yA(i, j, k)*( &
                          (mu_f + sigma_phi*mut_f)*(dtkldx * ynx(i, j, k)&
                            +dtkldy * yny(i, j, k) + dtkldz * ynz(i, j, k))))
        end do
       end do
      end do
    end subroutine compute_eta_viscous_fluxes_kkl

  !---------- eta face calculation end ----------------!


  !---------- zeta face calculation begin ----------------!

    subroutine compute_zeta_viscous_fluxes_laminar(H)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz 
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: H

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

          !--- For ODD-EVEN coupling error ---!
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

          !--- end of ODD-EVEN coupling correction ---!

          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
          mu_f = 0.5*(mu(i,j,k-1) + mu(i,j,k))
          mut_f = 0. 
          total_mu= mu_f + mut_f

          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9) * gm * R_gas / (gm - 1)
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
          
          uface = 0.5 * (x_speed(i, j, k-1) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i, j, k-1) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i, j, k-1) + z_speed(i, j, k))

          ! Energy flux
          ! (TijVj - Qi)ni
          H(i, j, k, 5) = H(i, j, k, 5) - (zA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * znx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * zny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * znz(i, j, k)) ) )
        end do
       end do
      end do
    end subroutine compute_zeta_viscous_fluxes_laminar

    subroutine compute_zeta_viscous_fluxes_sst(H)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz 
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: H

      !--- sst variable requirement ---!
      real :: dtkdx, dtkdy, dtkdz, dtwdx, dtwdy, dtwdz
      real ::  F1
      real ::  sigma_kf
      real ::  sigma_wf

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
          dtkdx = 0.5 * (gradtk_x(i, j, k-1) + gradtk_x(i, j, k))
          dtkdy = 0.5 * (gradtk_y(i, j, k-1) + gradtk_y(i, j, k))
          dtkdz = 0.5 * (gradtk_z(i, j, k-1) + gradtk_z(i, j, k))
          dtwdx = 0.5 * (gradtw_x(i, j, k-1) + gradtw_x(i, j, k))
          dtwdy = 0.5 * (gradtw_y(i, j, k-1) + gradtw_y(i, j, k))
          dtwdz = 0.5 * (gradtw_z(i, j, k-1) + gradtw_z(i, j, k))

          !--- For ODD-EVEN coupling error ---!
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
          !turbulence variable sst
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

          !--- end of ODD-EVEN coupling correction ---!

          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
          mu_f     = 0.5*(     mu(i,j,k-1) +      mu(i,j,k))
          F1       = 0.5*(sst_F1(i,j,k-1) + sst_F1(i,j,k))
          sigma_kf  =    sigma_k1*F1  +    sigma_k2*(1. - F1)
          sigma_wf  =    sigma_w1*F1  +    sigma_w2*(1. - F1)
          mut_f    = 0.5*(sst_mu(i, j, k-1) + sst_mu(i, j, k))

          total_mu= mu_f + mut_f

          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! k in reynolds stress
          Tau_xx = Tau_xx-((density(i,j,k-1)+density(i,j,k))*(tk(i,j,k-1)+tk(i,j,k)))/6.
          Tau_yy = Tau_yy-((density(i,j,k-1)+density(i,j,k))*(tk(i,j,k-1)+tk(i,j,k)))/6.
          Tau_zz = Tau_zz-((density(i,j,k-1)+density(i,j,k))*(tk(i,j,k-1)+tk(i,j,k)))/6.
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9) * gm * R_gas / (gm - 1)
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
          
          uface = 0.5 * (x_speed(i, j, k-1) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i, j, k-1) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i, j, k-1) + z_speed(i, j, k))

          ! Energy flux
          ! (TijVj - Qi)ni
          H(i, j, k, 5) = H(i, j, k, 5) - (zA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * znx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * zny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * znz(i, j, k)) ) )
          H(i, j, k, 6) = H(i, j, k, 6) - (zA(i, j, k)*( &
                          (mu_f + sigma_kf*mut_f)*(dtkdx * znx(i, j, k)&
                            +dtkdy * zny(i, j, k) + dtkdz * znz(i, j, k))))
          
          H(i, j, k, 7) = H(i, j, k, 7) - (zA(i, j, k)*( &
                          (mu_f + sigma_wf*mut_f)*(dtwdx * znx(i, j, k)&
                            +dtwdy * zny(i, j, k) + dtwdz * znz(i, j, k))))
        end do
       end do
      end do
    end subroutine compute_zeta_viscous_fluxes_sst


    subroutine compute_zeta_viscous_fluxes_kkl(H)
      implicit none
      
      real :: Tau_xx, Tau_xy, Tau_xz, Tau_yy, Tau_yz, Tau_zz
      real :: dudx, dudy, dudz
      real :: dvdx, dvdy, dvdz
      real :: dwdx, dwdy, dwdz 
      real :: dTdx, dTdy, dTdz
      real :: xc_L, yc_L, zc_L, xc_R, yc_R, zc_R, d_LR, normal_comp
      real :: T_LE, T_RE
      real :: K_heat, Qx, Qy, Qz
      real :: mu_f
      real :: mut_f
      real :: total_mu
      real :: uface, vface, wface
      integer :: i, j, k
      real, dimension(:, :, :, :), pointer :: H

      !--- kkl variable requirement  ---!
      real :: dtkdx, dtkdy, dtkdz
      real :: dtkldx, dtkldy, dtkldz

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
          dtkdx = 0.5 * (gradtk_x(i, j, k-1) + gradtk_x(i, j, k))
          dtkdy = 0.5 * (gradtk_y(i, j, k-1) + gradtk_y(i, j, k))
          dtkdz = 0.5 * (gradtk_z(i, j, k-1) + gradtk_z(i, j, k))
          dtkldx = 0.5 * (gradtkl_x(i, j, k-1) + gradtkl_x(i, j, k))
          dtkldy = 0.5 * (gradtkl_y(i, j, k-1) + gradtkl_y(i, j, k))
          dtkldz = 0.5 * (gradtkl_z(i, j, k-1) + gradtkl_z(i, j, k))

          !--- For ODD-EVEN coupling error ---!
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
          !turbulence variables kkl
          normal_comp = ( (tk(i, j, k) - &
                           tk(i, j, k-1)) - &
                        ((dtkdx * (xc_R - xc_L)) + &
                         (dtkdy * (yc_R - yc_L)) + &
                         (dtkdz * (zc_R - zc_L))) &
                        ) / d_LR
          dtkdx = dtkdx + (normal_comp * (xc_R - xc_L) / d_LR)
          dtkdy = dtkdy + (normal_comp * (yc_R - yc_L) / d_LR)
          dtkdz = dtkdz + (normal_comp * (zc_R - zc_L) / d_LR)

          normal_comp = ( (tkl(i, j, k) - &
                           tkl(i, j, k-1)) - &
                        ((dtkldx * (xc_R - xc_L)) + &
                         (dtkldy * (yc_R - yc_L)) + &
                         (dtkldz * (zc_R - zc_L))) &
                        ) / d_LR
          dtkldx = dtkldx + (normal_comp * (xc_R - xc_L) / d_LR)
          dtkldy = dtkldy + (normal_comp * (yc_R - yc_L) / d_LR)
          dtkldz = dtkldz + (normal_comp * (zc_R - zc_L) / d_LR)

          !--- end of ODD-EVEN coupling correction ---!

          ! Using lambda = -2 * mu / 3
          ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
          ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
          mu_f   = 0.5*(     mu(i,j,k-1) +      mu(i,j,k))
          mut_f  = 0.5*(kkl_mu(i, j, k-1) + kkl_mu(i, j, k))
          total_mu= mu_f + mut_f

          Tau_xx = 2. * total_mu * (dudx - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_yy = 2. * total_mu * (dvdy - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_zz = 2. * total_mu * (dwdz - ((dudx + dvdy + dwdz) / 3.)) 
          Tau_xy = total_mu * (dvdx + dudy)
          Tau_xz = total_mu * (dwdx + dudz)
          Tau_yz = total_mu * (dwdy + dvdz)

          ! k in reynolds stress
          Tau_xx = Tau_xx-((density(i,j,k-1)+density(i,j,k))*(tk(i,j,k-1)+tk(i,j,k)))/6.
          Tau_yy = Tau_yy-((density(i,j,k-1)+density(i,j,k))*(tk(i,j,k-1)+tk(i,j,k)))/6.
          Tau_zz = Tau_zz-((density(i,j,k-1)+density(i,j,k))*(tk(i,j,k-1)+tk(i,j,k)))/6.
          ! Pr: Prandtl Number
          ! Qx, Qy, Qz: Conduction fluxes
          K_heat = (mu_f/Pr + mut_f/0.9) * gm * R_gas / (gm - 1)
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
          
          uface = 0.5 * (x_speed(i, j, k-1) + x_speed(i, j, k))
          vface = 0.5 * (y_speed(i, j, k-1) + y_speed(i, j, k))
          wface = 0.5 * (z_speed(i, j, k-1) + z_speed(i, j, k))

          ! Energy flux
          ! (TijVj - Qi)ni
          H(i, j, k, 5) = H(i, j, k, 5) - (zA(i, j, k) * ( &
                          ((Tau_xx*uface + Tau_xy*vface + Tau_xz*wface + &
                            Qx) * znx(i, j, k)) + &
                          ((Tau_xy*uface + Tau_yy*vface + Tau_yz*wface + &
                            Qy) * zny(i, j, k)) + &
                          ((Tau_xz*uface + Tau_yz*vface + Tau_zz*wface + &
                            Qz) * znz(i, j, k)) ) )
          H(i, j, k, 6) = H(i, j, k, 6) - (zA(i, j, k)*( &
                          (mu_f + sigma_k*mut_f)*(dtkdx * znx(i, j, k)&
                            +dtkdy * zny(i, j, k) + dtkdz * znz(i, j, k))))
          
          H(i, j, k, 7) = H(i, j, k, 7) - (zA(i, j, k)*( &
                          (mu_f + sigma_phi*mut_f)*(dtkldx * znx(i, j, k)&
                            +dtkldy * zny(i, j, k) + dtkldz * znz(i, j, k))))
        end do
       end do
      end do
    end subroutine compute_zeta_viscous_fluxes_kkl

  !---------- eta face calculation end ----------------!
end module viscous
