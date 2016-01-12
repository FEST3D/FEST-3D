module viscous
    !-----------------------------------------------------------------
    ! The viscous module contains the viscous flux calculations and 
    ! the boundary conditions to be imposed
    !-----------------------------------------------------------------

    use global, only: FILE_NAME_LENGTH
    use utils, only: alloc, dealloc, dmsg
    use string
    use grid, only: imx, jmx, grid_x, grid_y
    use geometry, only: xnx, xny, ynx, yny, &
                        xA, yA, left_ghost_centroid, &
        right_ghost_centroid, top_ghost_centroid, bottom_ghost_centroid
    use geometry, only: vol => volume
    use state, only: gm, R_gas, mu_ref, T_ref, Pr, Sutherland_temp, &
                     density, x_speed, y_speed, pressure
    use face_interpolant, only: x_density_left, x_density_right, &
        y_density_left, y_density_right, &
        x_x_speed_left, x_x_speed_right, x_y_speed_left, x_y_speed_right, &
        y_x_speed_left, y_x_speed_right, y_y_speed_left, y_y_speed_right, &
        x_pressure_left, x_pressure_right, y_pressure_left, y_pressure_right

    implicit none
    private

!   real :: mu
    
    real, public, dimension(:, :), allocatable, target :: gradu_x, gradu_y, &
                                        gradv_x, gradv_y, gradT_x, gradT_y
    real, dimension(:, :), allocatable, target :: x_T_face, y_T_face

    public :: setup_viscous
    public :: compute_viscous_fluxes
    public :: destroy_viscous

    contains

    subroutine setup_viscous()

        implicit none

        call dmsg(1, 'viscous', 'setup_viscous')

        call alloc(gradu_x, 0, imx, 0, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradu_x - viscous')
        call alloc(gradu_y, 0, imx, 0, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradu_y - viscous')

        call alloc(gradv_x, 0, imx, 0, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradv_x - viscous')
        call alloc(gradv_y, 0, imx, 0, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'Gradv_y - viscous')

        call alloc(gradT_x, 0, imx, 0, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'GradT_x - viscous')
        call alloc(gradT_y, 0, imx, 0, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'GradT_y - viscous')

        call alloc(x_T_face, 1, imx, 1, jmx-1, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'x_T_face - viscous')
        call alloc(y_T_face, 1, imx-1, 1, jmx, &
                errmsg='Error: Unable to allocate memory for ' // &
                    'x_T_face - viscous')
    
    end subroutine setup_viscous

    subroutine destroy_viscous()

        implicit none

        call dmsg(1, 'viscous', 'destroy_viscous')

        call dealloc(gradu_x)
        call dealloc(gradu_y)

        call dealloc(gradv_x)
        call dealloc(gradv_y)
        
        call dealloc(gradT_x)
        call dealloc(gradT_y)
            
    end subroutine destroy_viscous
    
    subroutine compute_gradients_cell_centre()
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the cell
    ! centre
    !-----------------------------------------------------------------

        implicit none
        call dmsg(1, 'viscous', 'compute_gradients_cell_centre')

        x_T_face = 0.5 * (x_pressure_left/(x_density_left * R_gas)) + &
                   (x_pressure_right/(x_density_right * R_gas)) 
        y_T_face = 0.5 * (y_pressure_left/(y_density_left * R_gas)) + &
                   (y_pressure_right/(y_density_right * R_gas)) 

        gradu_x = 0.0
        gradu_y = 0.0
        gradv_x = 0.0
        gradv_y = 0.0
        gradT_x = 0.0
        gradT_y = 0.0

        gradu_x(1:imx-1, 1:jmx-1) = ( &
          - ((x_x_speed_left(1:imx-1, 1:jmx-1) + x_x_speed_right(1:imx-1, 1:jmx-1)) * &
               xnx(1:imx-1,1:jmx-1) * xA(1:imx-1,1:jmx-1)) &
          + ((x_x_speed_left(2:imx, 1:jmx-1) + x_x_speed_right(2:imx, 1:jmx-1)) * &
               xnx(2:imx,1:jmx-1) * xA(2:imx,1:jmx-1)) &
          - ((y_x_speed_left(1:imx-1, 1:jmx-1) + y_x_speed_right(1:imx-1, 1:jmx-1)) * &
               ynx(1:imx-1,1:jmx-1) * yA(1:imx-1,1:jmx-1)) &
          + ((y_x_speed_left(1:imx-1, 2:jmx) + y_x_speed_right(1:imx-1, 2:jmx)) * &
               ynx(1:imx-1,2:jmx) * yA(1:imx-1,2:jmx)) &
          ) * 0.5 / vol

        gradu_y(1:imx-1, 1:jmx-1) = ( &
          - ((x_x_speed_left(1:imx-1, 1:jmx-1) + x_x_speed_right(1:imx-1, 1:jmx-1)) * &
               xny(1:imx-1,1:jmx-1) * xA(1:imx-1,1:jmx-1)) &
          + ((x_x_speed_left(2:imx, 1:jmx-1) + x_x_speed_right(2:imx, 1:jmx-1)) * &
               xny(2:imx,1:jmx-1) * xA(2:imx,1:jmx-1)) &
          - ((y_x_speed_left(1:imx-1, 1:jmx-1) + y_x_speed_right(1:imx-1, 1:jmx-1)) * &
               yny(1:imx-1,1:jmx-1) * yA(1:imx-1,1:jmx-1)) &
          + ((y_x_speed_left(1:imx-1, 2:jmx) + y_x_speed_right(1:imx-1, 2:jmx)) * &
               yny(1:imx-1,2:jmx) * yA(1:imx-1,2:jmx)) &
          ) * 0.5 / vol

        gradv_x(1:imx-1, 1:jmx-1) = ( &
          - ((x_y_speed_left(1:imx-1, 1:jmx-1) + x_y_speed_right(1:imx-1, 1:jmx-1)) * &
               xnx(1:imx-1,1:jmx-1) * xA(1:imx-1,1:jmx-1)) &
          + ((x_y_speed_left(2:imx, 1:jmx-1) + x_y_speed_right(2:imx, 1:jmx-1)) * &
               xnx(2:imx,1:jmx-1) * xA(2:imx,1:jmx-1)) &
          - ((y_y_speed_left(1:imx-1, 1:jmx-1) + y_y_speed_right(1:imx-1, 1:jmx-1)) * &
               ynx(1:imx-1,1:jmx-1) * yA(1:imx-1,1:jmx-1)) &
          + ((y_y_speed_left(1:imx-1, 2:jmx) + y_y_speed_right(1:imx-1, 2:jmx)) * &
               ynx(1:imx-1,2:jmx) * yA(1:imx-1,2:jmx)) &
          ) * 0.5 / vol

        gradv_y(1:imx-1, 1:jmx-1) = ( &
          - ((x_y_speed_left(1:imx-1, 1:jmx-1) + x_y_speed_right(1:imx-1, 1:jmx-1)) * &
               xny(1:imx-1,1:jmx-1) * xA(1:imx-1,1:jmx-1)) &
          + ((x_y_speed_left(2:imx, 1:jmx-1) + x_y_speed_right(2:imx, 1:jmx-1)) * &
               xny(2:imx,1:jmx-1) * xA(2:imx,1:jmx-1)) &
          - ((y_y_speed_left(1:imx-1, 1:jmx-1) + y_y_speed_right(1:imx-1, 1:jmx-1)) * &
               yny(1:imx-1,1:jmx-1) * yA(1:imx-1,1:jmx-1)) &
          + ((y_y_speed_left(1:imx-1, 2:jmx) + y_y_speed_right(1:imx-1, 2:jmx)) * &
               yny(1:imx-1,2:jmx) * yA(1:imx-1,2:jmx)) &
          ) * 0.5 / vol

        gradT_x(1:imx-1, 1:jmx-1) = ( &
          - (x_T_face(1:imx-1, 1:jmx-1) * xnx(1:imx-1,1:jmx-1) * xA(1:imx-1,1:jmx-1)) &
          + (x_T_face(2:imx, 1:jmx-1) * xnx(2:imx,1:jmx-1) * xA(2:imx,1:jmx-1)) &
          - (y_T_face(1:imx-1, 1:jmx-1) * ynx(1:imx-1,1:jmx-1) * yA(1:imx-1,1:jmx-1)) &
          + (y_T_face(1:imx-1, 2:jmx) * ynx(1:imx-1,2:jmx) * yA(1:imx-1,2:jmx)) &
          ) / vol

        gradT_y(1:imx-1, 1:jmx-1) = ( &
          - (x_T_face(1:imx-1, 1:jmx-1) * xny(1:imx-1,1:jmx-1) * xA(1:imx-1,1:jmx-1)) &
          + (x_T_face(2:imx, 1:jmx-1) * xny(2:imx,1:jmx-1) * xA(2:imx,1:jmx-1)) &
          - (y_T_face(1:imx-1, 1:jmx-1) * yny(1:imx-1,1:jmx-1) * yA(1:imx-1,1:jmx-1)) &
          + (y_T_face(1:imx-1, 2:jmx) * yny(1:imx-1,2:jmx) * yA(1:imx-1,2:jmx)) &
          ) / vol

        ! Ghost cells same as interior
        gradu_x(0, 1:jmx-1) = gradu_x(1, 1:jmx-1)
        gradu_y(0, 1:jmx-1) = gradu_y(1, 1:jmx-1)
        gradv_x(0, 1:jmx-1) = gradv_x(1, 1:jmx-1)
        gradv_y(0, 1:jmx-1) = gradv_y(1, 1:jmx-1)
        gradT_x(0, 1:jmx-1) = gradT_x(1, 1:jmx-1)
        gradT_y(0, 1:jmx-1) = gradT_y(1, 1:jmx-1)

        gradu_x(imx, 1:jmx-1) = gradu_x(imx - 1, 1:jmx-1)
        gradu_y(imx, 1:jmx-1) = gradu_y(imx - 1, 1:jmx-1)
        gradv_x(imx, 1:jmx-1) = gradv_x(imx - 1, 1:jmx-1)
        gradv_y(imx, 1:jmx-1) = gradv_y(imx - 1, 1:jmx-1)
        gradT_x(imx, 1:jmx-1) = gradT_x(imx - 1, 1:jmx-1)
        gradT_y(imx, 1:jmx-1) = gradT_y(imx - 1, 1:jmx-1)

        gradu_x(1:imx-1, 0) = gradu_x(1:imx-1, 1)
        gradu_y(1:imx-1, 0) = gradu_y(1:imx-1, 1)
        gradv_x(1:imx-1, 0) = gradv_x(1:imx-1, 1)
        gradv_y(1:imx-1, 0) = gradv_y(1:imx-1, 1)
        gradT_x(1:imx-1, 0) = gradT_x(1:imx-1, 1)
        gradT_y(1:imx-1, 0) = gradT_y(1:imx-1, 1)

        gradu_x(1:imx-1, jmx) = gradu_x(1:imx-1, jmx - 1)
        gradu_y(1:imx-1, jmx) = gradu_y(1:imx-1, jmx - 1)
        gradv_x(1:imx-1, jmx) = gradv_x(1:imx-1, jmx - 1)
        gradv_y(1:imx-1, jmx) = gradv_y(1:imx-1, jmx - 1)
        gradT_x(1:imx-1, jmx) = gradT_x(1:imx-1, jmx - 1)
        gradT_y(1:imx-1, jmx) = gradT_y(1:imx-1, jmx - 1)

    end subroutine compute_gradients_cell_centre

    subroutine compute_viscous_face_fluxes(Flux_p, f_dir)
    !-----------------------------------------------------------------
    ! Computes the gradients of velocity and temperature at the xi-face
    ! using a different scheme instead of taking the average of the 
    ! cell centre gradients. The latter leads to odd-even decoupling
    ! problem, which the different scheme corrects.
    !-----------------------------------------------------------------

        implicit none
        
        real :: Tau_xx, Tau_xy, Tau_yy
        real :: dudx, dudy, dvdx, dvdy, dTdx, dTdy 
        real :: xc_L, yc_L, xc_R, yc_R, d_LR, normal_comp
        real :: T_LE, T_RE
        real :: K_heat, mu, Qx, Qy
        real :: uface, vface
        integer :: i, j, i_f, j_f
        real, dimension(:, :, :), pointer :: Flux_p
        real, dimension(:, :), pointer :: fnx, fny, fA, f_T_face
        real, dimension(1:imx, 1:jmx) :: f_x_speed, f_y_speed
        
        character, intent(in) :: f_dir

        call dmsg(1, 'viscous', 'compute_viscous_face_fluxes')

        select case (f_dir)
            case ('x')
                i_f = 1
                j_f = 0
                f_T_face => x_T_face
                fnx => xnx
                fny => xny
                fA => xA
                f_x_speed(1:imx, 1:jmx-1) = 0.5 * (x_x_speed_left + &
                                                   x_x_speed_right)
                f_y_speed(1:imx, 1:jmx-1) = 0.5 * (x_y_speed_left + &
                                                   x_y_speed_right)
            case ('y')
                i_f = 0
                j_f = 1
                f_T_face => y_T_face
                fnx => ynx
                fny => yny
                fA => yA
                f_x_speed(1:imx-1, 1:jmx) = 0.5 * (y_x_speed_left + &
                                                   y_x_speed_right)
                f_y_speed(1:imx-1, 1:jmx) = 0.5 * (y_y_speed_left + &
                                                   y_y_speed_right)
            case default
                call dmsg(5, 'ppm', 'pressure_based_switching', &
                        'Direction not recognised')
                stop
        end select


        ! Calculating the fluxes at the faces
        ! A different calculation is to be done for interior faces as compared
        ! to the bounday
        ! Calculating for the interior xi-faces
        do j = 1, jmx - 1 + j_f
         do i = 1, imx - 1 + i_f
            ! Gradients at face as average of gradients at cell centres
            dudx = 0.5 * (gradu_x(i- i_f, j - j_f) + gradu_x(i, j))
            dudy = 0.5 * (gradu_y(i- i_f, j - j_f) + gradu_y(i, j))
            dvdx = 0.5 * (gradv_x(i- i_f, j - j_f) + gradv_x(i, j))
            dvdy = 0.5 * (gradv_y(i- i_f, j - j_f) + gradv_y(i, j))
            dTdx = 0.5 * (gradT_x(i- i_f, j - j_f) + gradT_x(i, j))
            dTdy = 0.5 * (gradT_y(i- i_f, j - j_f) + gradT_y(i, j))

            ! Correcting the odd-even problem
            if (i_f .eq. 1) then
              if (i .eq. 1) then
                xc_L = left_ghost_centroid(j, 1)
                yc_L = left_ghost_centroid(j, 2)
              else
                ! Coordinate of left cell centre: element (i-1, j)
                xc_L = (grid_x(i-1, j) + grid_x(i, j) + &
                        grid_x(i, j+1) + grid_x(i-1, j+1)) * 0.25
                yc_L = (grid_y(i-1, j) + grid_y(i, j) + &
                        grid_y(i, j+1) + grid_y(i-1, j+1)) * 0.25
              end if

              if (i .eq. imx) then
                xc_R = right_ghost_centroid(j, 1)
                yc_R = right_ghost_centroid(j, 2)
              else           
                ! Coordinate of right cell centre: element (i, j)
                xc_R = (grid_x(i, j) + grid_x(i+1, j) + &
                        grid_x(i+1, j+1) + grid_x(i, j+1)) * 0.25
                yc_R = (grid_y(i, j) + grid_y(i+1, j) + &
                        grid_y(i+1, j+1) + grid_y(i, j+1)) * 0.25
              end if
            else 
              if (j .eq. 1) then
                xc_L = bottom_ghost_centroid(i, 1)
                yc_L = bottom_ghost_centroid(i, 2)
              else
                ! Coordinate of bottom cell centre: element (i, j-1)
                xc_L = (grid_x(i, j-1) + grid_x(i+1, j-1) + &
                        grid_x(i+1, j) + grid_x(i, j)) * 0.25
                yc_L = (grid_y(i, j-1) + grid_y(i+1, j-1) + &
                        grid_y(i+1, j) + grid_y(i, j)) * 0.25
              end if

              if (j .eq. jmx) then
                xc_R = top_ghost_centroid(i, 1)
                yc_R = top_ghost_centroid(i, 2)
              else           
                ! Coordinate of top cell centre: element (i, j)
                xc_R = (grid_x(i, j) + grid_x(i+1, j) + &
                        grid_x(i+1, j+1) + grid_x(i, j+1)) * 0.25
                yc_R = (grid_y(i, j) + grid_y(i+1, j) + &
                        grid_y(i+1, j+1) + grid_y(i, j+1)) * 0.25
              end if
            end if

            d_LR = sqrt((xc_R - xc_L)**2 + (yc_R - yc_L)**2)

            ! normal_comp is the component along r_ij
            ! W_j - W_i = W_right_cell - W_left_cell
            !           = W(i, j, k) - W(i-1, j, k)
            ! For this, the values of state variables at the cell
            ! centres are needed
            normal_comp = ( (x_speed(i, j) - &
                             x_speed(i-i_f, j-j_f)) - &
                          ((dudx * (xc_R - xc_L)) + &
                           (dudy * (yc_R - yc_L)) &
                          )) / d_LR
            dudx = dudx + (normal_comp * (xc_R - xc_L) / d_LR)
            dudy = dudy + (normal_comp * (yc_R - yc_L) / d_LR)

            normal_comp = ( (y_speed(i, j) - &
                             y_speed(i-i_f, j-j_f)) - &
                          ((dvdx * (xc_R - xc_L)) + &
                           (dvdy * (yc_R - yc_L)) &
                          )) / d_LR
            dvdx = dvdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dvdy = dvdy + (normal_comp * (yc_R - yc_L) / d_LR)

            ! Finding the temperature of left and right element to the
            ! face i, j, k
            T_LE = pressure(i-i_f, j-j_f) / (density(i-i_f, j-j_f) * R_gas)
            T_RE = pressure(i, j) / (density(i, j) * R_gas)
            normal_comp = ( (T_RE - T_LE) - &
                          ((dTdx * (xc_R - xc_L)) + &
                           (dTdy * (yc_R - yc_L)) &
                          )) / d_LR
            dTdx = dTdx + (normal_comp * (xc_R - xc_L) / d_LR)
            dTdy = dTdy + (normal_comp * (yc_R - yc_L) / d_LR)
            
            ! mu requires T at the face. Hence:
            ! T_L and T_R are the left and right states of the face i,j,k
            ! The values at face used instead of element values
            mu = mu_ref * (f_T_face(i,j) / T_ref)**1.5 * ((T_ref + &
                          Sutherland_temp) / (f_T_face(i,j) + Sutherland_temp))

            ! Using lambda = -2 * mu / 3
            ! For the xi direction fluxes, only Tau_xx, Tau_yx, 
            ! Tau_zx is needed. Tau_yx = Tau_xy and T_zx = Tau_xz
            Tau_xx = 2. * mu * (dudx - ((dudx + dvdy) / 3.)) 
            Tau_yy = 2. * mu * (dvdy - ((dudx + dvdy) / 3.)) 
            Tau_xy = mu * (dvdx + dudy)

            ! Pr: Prandtl Number
            ! Qx, Qy, Qz: Conduction fluxes
            K_heat = mu * gm * R_gas / ((gm - 1) * Pr)
            Qx = K_heat*dTdx
            Qy = K_heat*dTdy

            ! Note that the xi-direction faces only need the following quantities:
            ! Tau_xx, Tau_xy, Tau_xz -> dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz
            ! Qx -> dTdx
            ! The mass flux has no viscous component
            ! momentum for xi-face:
            Flux_p(i, j, 2) = Flux_p(i, j, 2) - ( ((Tau_xx * fnx(i, j)) + &
                         (Tau_xy * fny(i, j))) * fA(i, j))
            Flux_p(i, j, 3) = Flux_p(i, j, 3) - ( ((Tau_xy * fnx(i, j)) + &
                         (Tau_yy * fny(i, j))) * fA(i, j))
           
            ! Energy flux
            uface = f_x_speed(i, j)
            vface = f_y_speed(i, j)
            ! (TijVj + Qi)ni
            Flux_p(i, j, 4) = Flux_p(i, j, 4) - (fA(i, j) * ( &
                            ((Tau_xx*uface + Tau_xy*vface + Qx) * &
                             fnx(i, j)) + &
                            ((Tau_xy*uface + Tau_yy*vface + Qy) * &
                             fny(i, j)) ) )
         end do
        end do

    end subroutine compute_viscous_face_fluxes

    subroutine compute_viscous_fluxes(F, G)

        implicit none

        real, dimension(:, :, :), pointer :: F, G
        call dmsg(1, 'viscous', 'compute_viscous_fluxes')
        
        call compute_gradients_cell_centre()
        call compute_viscous_face_fluxes(F, 'x')
        call compute_viscous_face_fluxes(G, 'y')

    end subroutine compute_viscous_fluxes

end module viscous
