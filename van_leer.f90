module van_leer
    !-------------------------------------------------------------------
    ! The Van-Leer scheme is a type of flux-splitting scheme
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx
    use geometry, only: xnx, xny, ynx, yny, xA, yA
    use state, only: gm
    use face_interpolant, only: x_qp_left, x_qp_right, y_qp_left, y_qp_right, &
            x_density_left, x_x_speed_left, x_y_speed_left, x_pressure_left, &
            x_density_right, x_x_speed_right, x_y_speed_right, &
                x_pressure_right, &
            y_density_left, y_x_speed_left, y_y_speed_left, y_pressure_left, &
            y_density_right, y_x_speed_right, y_y_speed_right, &
                y_pressure_right, &
            x_sound_speed_left, x_sound_speed_right, &
            y_sound_speed_left, y_sound_speed_right

    implicit none
    private

    real, dimension(:, :), allocatable :: x_sound_speed_avg
    real, dimension(:, :), allocatable :: x_M_perp_left, x_M_perp_right
    real, dimension(:, :), allocatable :: y_sound_speed_avg
    real, dimension(:, :), allocatable :: y_M_perp_left, y_M_perp_right
    real, dimension(:, :), allocatable :: x_alpha_plus, x_alpha_minus
    real, dimension(:, :), allocatable :: y_alpha_plus, y_alpha_minus
    real, dimension(:, :), allocatable :: x_beta_left, x_beta_right
    real, dimension(:, :), allocatable :: y_beta_left, y_beta_right
    real, dimension(:, :), allocatable :: x_c_plus, x_c_minus
    real, dimension(:, :), allocatable :: y_c_plus, y_c_minus
    real, dimension(:, :), allocatable :: x_scrD_plus, x_scrD_minus
    real, dimension(:, :), allocatable :: y_scrD_plus, y_scrD_minus

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: get_residue

    ! Public members which can be used in schemes derived from Van Leer
    public :: compute_xi_face_quantities
    public :: compute_eta_face_quantities
    public :: compute_residue
    public :: x_M_perp_left, x_M_perp_right
    public :: y_M_perp_left, y_M_perp_right
    public :: x_beta_left, x_beta_right
    public :: y_beta_left, y_beta_right
    public :: x_c_plus, x_c_minus
    public :: y_c_plus, y_c_minus

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'van_leer', 'setup_scheme')

            call alloc(x_M_perp_left, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_M_perp_left.')
            call alloc(x_M_perp_right, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_M_perp_right.')
            call alloc(y_M_perp_left, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_M_perp_left.')
            call alloc(y_M_perp_right, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_M_perp_right.')

            call alloc(x_alpha_plus, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_alpha_plus.')
            call alloc(x_alpha_minus, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_alpha_minus.')
            call alloc(y_alpha_plus, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_alpha_plus.')
            call alloc(y_alpha_minus, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_alpha_minus.')

            call alloc(x_beta_left, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for x_beta_left.')
            call alloc(x_beta_right, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_beta_right.')
            call alloc(y_beta_left, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for y_beta_left.')
            call alloc(y_beta_right, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_beta_right.')

            call alloc(x_c_plus, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for x_c_plus.')
            call alloc(x_c_minus, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for x_c_minus.')
            call alloc(y_c_plus, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for y_c_plus.')
            call alloc(y_c_minus, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for y_c_minus.')

            call alloc(x_scrD_plus, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for x_scrD_plus.')
            call alloc(x_scrD_minus, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_scrD_minus.')
            call alloc(y_scrD_plus, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for y_scrD_plus.')
            call alloc(y_scrD_minus, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_scrD_minus.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'van_leer', 'destroy_scheme')

            call dealloc(x_M_perp_left)
            call dealloc(x_M_perp_right)
            call dealloc(y_M_perp_left)
            call dealloc(y_M_perp_right)

            call dealloc(x_alpha_plus)
            call dealloc(x_alpha_minus)
            call dealloc(y_alpha_plus)
            call dealloc(y_alpha_minus)

            call dealloc(x_beta_left)
            call dealloc(x_beta_right)
            call dealloc(y_beta_left)
            call dealloc(y_beta_right)

            call dealloc(x_c_plus)
            call dealloc(x_c_minus)
            call dealloc(y_c_plus)
            call dealloc(y_c_minus)

            call dealloc(x_scrD_plus)
            call dealloc(x_scrD_minus)
            call dealloc(y_scrD_plus)
            call dealloc(y_scrD_minus)

        end subroutine destroy_scheme

        subroutine compute_xi_face_quantities()
            !-----------------------------------------------------------
            ! Compute xi direction quantities used in F flux computation
            !
            ! x_c_plus and x_scrD_plus are used to calculate F_plus. 
            ! Similarly for F_minus. This subroutine computes these.
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx, jmx-1) :: x_face_normal_speeds
            real, dimension(imx, jmx-1) :: M_plus, M_minus
            real, dimension(imx, jmx-1) :: D_plus, D_minus

            call dmsg(1, 'van_leer', 'compute_xi_face_quantities')

            x_sound_speed_avg = 0.5 * &
                    (x_sound_speed_left() + x_sound_speed_right())

            ! Compute the '+' direction quantities
            x_face_normal_speeds = x_x_speed_left * xnx + x_y_speed_left * xny
            x_M_perp_left = x_face_normal_speeds / x_sound_speed_avg
            x_alpha_plus = 0.5 * (1.0 + sign(1.0, x_M_perp_left))
            x_beta_left = -max(0, 1 - floor(abs(x_M_perp_left)))
            M_plus = 0.25 * ((1. + x_M_perp_left) ** 2.)
            D_plus = 0.25 * ((1. + x_M_perp_left) ** 2.) * (2. - x_M_perp_left)
            x_c_plus = (x_alpha_plus * (1.0 + x_beta_left) * x_M_perp_left) - &
                    x_beta_left * M_plus
            x_scrD_plus = (x_alpha_plus * (1. + x_beta_left)) - &
                    (x_beta_left * D_plus)

            ! Compute the '-' direction quantities
            x_face_normal_speeds = x_x_speed_right * xnx + &
                    x_y_speed_right * xny
            x_M_perp_right = x_face_normal_speeds / x_sound_speed_avg
            x_alpha_minus = 0.5 * (1.0 - sign(1.0, x_M_perp_right))
            x_beta_right = -max(0, 1 - floor(abs(x_M_perp_right)))
            M_minus = - 0.25 * ((1 - x_M_perp_right) ** 2.)
            D_minus = 0.25 * ((1 - x_M_perp_right) ** 2.) * &
                    (2. + x_M_perp_right)
            x_c_minus = (x_alpha_minus * (1.0 + x_beta_right) * &
                    x_M_perp_right) - (x_beta_right * M_minus)
            x_scrD_minus = (x_alpha_minus * (1.0 + x_beta_right)) - &
                    (x_beta_right * D_minus)

            call dmsg(1, 'van_leer', 'compute_xi_face_quantities', 'Ended')

        end subroutine compute_xi_face_quantities

        function compute_F_plus() result(F_plus)

            implicit none
            real, dimension(imx, jmx-1, 4) :: F_plus

            ! First construct the F mass flux
            F_plus(:, :, 1) = x_density_left * x_sound_speed_avg * x_c_plus
            ! Construct other fluxes in terms of the F mass flux
            F_plus(:, :, 2) = (F_plus(:, :, 1) * x_x_speed_left) + &
                    (x_scrD_plus * x_pressure_left * xnx)
            F_plus(:, :, 3) = (F_plus(:, :, 1) * x_y_speed_left) + &
                    (x_scrD_plus * x_pressure_left * xny)
            F_plus(:, :, 4) = F_plus(:, :, 1) * &
                    ((0.5 * (x_x_speed_left ** 2. + x_y_speed_left ** 2.)) + &
                    ((gm / (gm - 1.)) * x_pressure_left / x_density_left))

            ! Multiply in the face areas
            F_plus(:, :, 1) = F_plus(:, :, 1) * xA
            F_plus(:, :, 2) = F_plus(:, :, 2) * xA
            F_plus(:, :, 3) = F_plus(:, :, 3) * xA
            F_plus(:, :, 4) = F_plus(:, :, 4) * xA

        end function compute_F_plus

        function compute_F_minus() result(F_minus)
            
            implicit none
            real, dimension(imx, jmx-1, 4) :: F_minus

            ! First construct the F mass flux
            F_minus(:, :, 1) = x_density_right * x_sound_speed_avg * x_c_minus
            ! Construct other fluxes in terms of the F mass flux
            F_minus(:, :, 2) = (F_minus(:, :, 1) * x_x_speed_right) + &
                    (x_scrD_minus * x_pressure_right * xnx)
            F_minus(:, :, 3) = (F_minus(:, :, 1) * x_y_speed_right) + &
                    (x_scrD_minus * x_pressure_right * xny)
            F_minus(:, :, 4) = F_minus(:, :, 1) * &
                    ((0.5 * (x_x_speed_right ** 2. + x_y_speed_right ** 2.)) + &
                    ((gm / (gm - 1.)) * x_pressure_right / x_density_right))

            ! Multiply in the face areas
            F_minus(:, :, 1) = F_minus(:, :, 1) * xA
            F_minus(:, :, 2) = F_minus(:, :, 2) * xA
            F_minus(:, :, 3) = F_minus(:, :, 3) * xA
            F_minus(:, :, 4) = F_minus(:, :, 4) * xA

        end function compute_F_minus

        subroutine compute_eta_face_quantities()
            !-----------------------------------------------------------
            ! Compute eta direction quantities used in G flux computation
            !
            ! y_c_plus and y_scrD_plus are used to calculate G_plus. 
            ! Similarly for G_minus. This subroutine computes these.
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx-1, jmx) :: y_face_normal_speeds
            real, dimension(imx-1, jmx) :: M_plus, M_minus
            real, dimension(imx-1, jmx) :: D_plus, D_minus

            call dmsg(1, 'van_leer', 'compute_eta_face_quantities')

            y_sound_speed_avg = 0.5 * &
                    (y_sound_speed_left() + y_sound_speed_right())

            ! Compute the '+' direction quantities
            y_face_normal_speeds = y_x_speed_left * ynx + y_y_speed_left * yny
            y_M_perp_left = y_face_normal_speeds / y_sound_speed_avg
            y_alpha_plus = 0.5 * (1.0 + sign(1.0, y_M_perp_left))
            y_beta_left = -max(0, 1 - floor(abs(y_M_perp_left)))
            M_plus = 0.25 * ((1. + y_M_perp_left) ** 2.)
            D_plus = 0.25 * ((1. + y_M_perp_left) ** 2.) * (2. - y_M_perp_left)
            y_c_plus = (y_alpha_plus * (1.0 + y_beta_left) * y_M_perp_left) - &
                    y_beta_left * M_plus
            y_scrD_plus = (y_alpha_plus * (1. + y_beta_left)) - &
                    (y_beta_left * D_plus)

            ! Compute the '-' direction quantities
            y_face_normal_speeds = y_x_speed_right * ynx + &
                    y_y_speed_right(:, :) * yny
            y_M_perp_right = y_face_normal_speeds / y_sound_speed_avg
            y_alpha_minus = 0.5 * (1.0 - sign(1.0, y_M_perp_right))
            y_beta_right = -max(0, 1 - floor(abs(y_M_perp_right)))
            M_minus = - 0.25 * ((1 - y_M_perp_right) ** 2.)
            D_minus = 0.25 * ((1 - y_M_perp_right) ** 2.) * &
                    (2. + y_M_perp_right)
            y_c_minus = (y_alpha_minus * (1.0 + y_beta_right) * &
                    y_M_perp_right) - (y_beta_right * M_minus)
            y_scrD_minus = (y_alpha_minus * (1.0 + y_beta_right)) - &
                    (y_beta_right * D_minus)
            call dmsg(1, 'van_leer', 'compute_eta_face_quantities', 'Ended')

        end subroutine compute_eta_face_quantities

        function compute_G_plus() result(G_plus)

            implicit none
            real, dimension(imx-1, jmx, 4) :: G_plus

            ! First construct the G mass flux
            G_plus(:, :, 1) = y_density_left * y_sound_speed_avg * y_c_plus
            ! Make the convective G mass flux zero at the wall
            G_plus(:, 1, 1) = 0
            G_plus(:, jmx, 1) = 0
            ! Construct other fluxes in terms of the G mass flux
            G_plus(:, :, 2) = (G_plus(:, :, 1) * y_x_speed_left) + &
                    (y_scrD_plus * y_pressure_left * ynx)
            G_plus(:, :, 3) = (G_plus(:, :, 1) * y_y_speed_left) + &
                    (y_scrD_plus * y_pressure_left * yny)
            G_plus(:, :, 4) = G_plus(:, :, 1) * &
                    ((0.5 * (y_x_speed_left ** 2. + y_y_speed_left ** 2.)) + &
                    ((gm / (gm - 1.)) * y_pressure_left / y_density_left))

            ! Multiply in the face areas
            G_plus(:, :, 1) = G_plus(:, :, 1) * yA
            G_plus(:, :, 2) = G_plus(:, :, 2) * yA
            G_plus(:, :, 3) = G_plus(:, :, 3) * yA
            G_plus(:, :, 4) = G_plus(:, :, 4) * yA

        end function compute_G_plus

        function compute_G_minus() result(G_minus)
            
            implicit none
            real, dimension(imx-1, jmx, 4) :: G_minus

            ! First construct the G mass flux
            G_minus(:, :, 1) = y_density_right * y_sound_speed_avg * y_c_minus
            ! Make the convective G mass flux zero at the wall
            G_minus(:, 1, 1) = 0
            G_minus(:, jmx, 1) = 0
            ! Construct other fluxes in terms of the G mass flux
            G_minus(:, :, 2) = (G_minus(:, :, 1) * y_x_speed_right) + &
                    (y_scrD_minus * y_pressure_right * ynx)
            G_minus(:, :, 3) = (G_minus(:, :, 1) * y_y_speed_right) + &
                    (y_scrD_minus * y_pressure_right * yny)
            G_minus(:, :, 4) = G_minus(:, :, 1) * &
                    ((0.5 * (y_x_speed_right ** 2. + &
                        y_y_speed_right ** 2.)) + &
                    ((gm / (gm - 1.)) * y_pressure_right / y_density_right))

            ! Multiply in the face areas
            G_minus(:, :, 1) = G_minus(:, :, 1) * yA
            G_minus(:, :, 2) = G_minus(:, :, 2) * yA
            G_minus(:, :, 3) = G_minus(:, :, 3) * yA
            G_minus(:, :, 4) = G_minus(:, :, 4) * yA

        end function compute_G_minus

        function compute_residue() result(residue)
            !-----------------------------------------------------------
            ! Compute the residue using the Van-Leer scheme
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx, jmx-1, 4) :: F_plus, F_minus
            real, dimension(imx-1, jmx, 4) :: G_plus, G_minus
            real, dimension(imx-1, jmx-1, 4) :: residue

            call dmsg(1, 'van_leer', 'compute_residue')

            F_plus = compute_F_plus()
            F_minus = compute_F_minus()
            !print *, 'F_plus: ', F_plus
            if (any(isnan(F_plus))) stop
            !print *, 'F_minus: ', F_minus
            if (any(isnan(F_minus))) stop

            G_plus = compute_G_plus()
            G_minus = compute_G_minus()
            !print *, 'G_plus: ', G_plus
            if (any(isnan(G_plus))) stop
            !print *, 'G_minus: ', G_minus
            if (any(isnan(G_minus))) stop

            residue = F_plus(2:imx, 1:jmx-1, :) &
                    + F_minus(2:imx, 1:jmx-1, :) &
                    - F_plus(1:imx-1, 1:jmx-1, :) &
                    - F_minus(1:imx-1, 1:jmx-1, :) &
                    + G_plus(1:imx-1, 2:jmx, :) &
                    + G_minus(1:imx-1, 2:jmx, :) &
                    - G_plus(1:imx-1, 1:jmx-1, :) &
                    - G_minus(1:imx-1, 1:jmx-1, :)

        end function compute_residue

        function get_residue() result(residue)
            !-----------------------------------------------------------
            ! Return the VL residue
            !-----------------------------------------------------------
            
            implicit none
            real, dimension(imx-1, jmx-1, 4) :: residue

            call compute_xi_face_quantities()
            call compute_eta_face_quantities()

            residue = compute_residue()

        end function get_residue

end module van_leer
