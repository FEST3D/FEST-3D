module ldfss0
    !-------------------------------------------------------------------
    ! LDFSS is a class of flux-splitting schemes
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx
    use geometry, only: xnx, xny, ynx, yny, xA, yA
    use state, only: qp, pressure, density, x_speed, y_speed, gm, R_gas, &
            x_a, y_a, xi_face_normal_speeds, eta_face_normal_speeds

    implicit none
    private

    real, dimension(:, :), allocatable :: x_c_plus, x_c_minus
    real, dimension(:, :), allocatable :: y_c_plus, y_c_minus
    real, dimension(:, :), allocatable :: x_scrD_plus, x_scrD_minus
    real, dimension(:, :), allocatable :: y_scrD_plus, y_scrD_minus

    ! Public methods
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_residue

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ldfss', 'setup_scheme')

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

            call dmsg(1, 'ldfss', 'destroy_scheme')

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
            real, dimension(imx, jmx-1) :: x_M_perp_left, x_M_perp_right
            real, dimension(imx, jmx-1) :: x_alpha_plus, x_alpha_minus
            real, dimension(imx, jmx-1) :: x_beta_left, x_beta_right
            real, dimension(imx, jmx-1) :: M_plus, M_minus
            real, dimension(imx, jmx-1) :: D_plus, D_minus
            real, dimension(imx, jmx-1) :: M_ldfss

            call dmsg(1, 'ldfss', 'compute_xi_face_quantities')

            ! Begin by applying the Van-Leer scheme

            ! Compute the '+' direction quantities
            x_M_perp_left = xi_face_normal_speeds('+') / x_a
            x_alpha_plus = 0.5 * (1.0 + sign(1.0, x_M_perp_left))
            x_beta_left = -max(0, 1 - floor(abs(x_M_perp_left)))
            M_plus = 0.25 * ((1. + x_M_perp_left) ** 2.)
            D_plus = 0.25 * ((1. + x_M_perp_left) ** 2.) * (2. - x_M_perp_left)
            x_c_plus = (x_alpha_plus * (1.0 + x_beta_left) * x_M_perp_left) - &
                    x_beta_left * M_plus
            x_scrD_plus = (x_alpha_plus * (1. + x_beta_left)) - &
                    (x_beta_left * D_plus)

            !print *, 'x_a: ', x_a
            if (any(isnan(x_a))) stop
            !print *, 'x_M_perp_left: ', x_M_perp_left
            if (any(isnan(x_M_perp_left))) stop
            !print *, 'x_alpha_plus: ', x_alpha_plus
            if (any(isnan(x_alpha_plus))) stop
            !print *, 'x_beta: ', x_beta
            if (any(isnan(x_beta_left))) stop
            !print *, 'M_plus: ', M_plus
            if (any(isnan(M_plus))) stop
            !print *, 'D_plus: ', D_plus
            if (any(isnan(D_plus))) stop
            !print *, 'x_c_plus: ', x_c_plus
            if (any(isnan(x_c_plus))) stop
            !print *, 'x_scrD_plus: ', x_scrD_plus
            if (any(isnan(x_scrD_plus))) stop

            ! Compute the '-' direction quantities
            x_M_perp_right = xi_face_normal_speeds('-') / x_a
            x_alpha_minus = 0.5 * (1.0 - sign(1.0, x_M_perp_right))
            x_beta_right = -max(0, 1 - floor(abs(x_M_perp_right)))
            M_minus = - 0.25 * ((1 - x_M_perp_right) ** 2.)
            D_minus = 0.25 * ((1 - x_M_perp_right) ** 2.) * (2. + x_M_perp_right)
            x_c_minus = (x_alpha_minus * (1.0 + x_beta_right) * x_M_perp_right) - &
                    (x_beta_right * M_minus)
            x_scrD_minus = (x_alpha_minus * (1.0 + x_beta_right)) - &
                    (x_beta_right * D_minus)

            !print *, 'x_M_perp_right: ', x_M_perp_right
            if (any(isnan(x_M_perp_right))) stop
            !print *, 'x_alpha_minus: ', x_alpha_minus
            if (any(isnan(x_alpha_minus))) stop
            !print *, 'x_beta_right: ', x_beta_right
            if (any(isnan(x_beta_right))) stop
            !print *, 'M_minus: ', M_minus
            if (any(isnan(M_minus))) stop
            !print *, 'D_minus: ', D_minus
            if (any(isnan(D_minus))) stop
            !print *, 'x_c_minus: ', x_c_minus
            if (any(isnan(x_c_minus))) stop
            !print *, 'x_scrD_minus: ', x_scrD_minus
            if (any(isnan(x_scrD_minus))) stop

            ! Update the Van-Leer computed speeds (x_c_[plus|minus])
            M_ldfss = 0.25 * x_beta_left * x_beta_right * &
                    (sqrt((x_M_perp_left ** 2. + x_M_perp_right ** 2.) * 0.5) &
                    - 1) ** 2.
            x_c_plus = x_c_plus - M_ldfss
            x_c_minus = x_c_minus + M_ldfss

            call dmsg(1, 'ldfss', 'compute_xi_face_quantities', 'Ended')

        end subroutine compute_xi_face_quantities

        function compute_F_plus() result(F_plus)

            implicit none
            real, dimension(imx, jmx-1, 4) :: F_plus

            ! First construct the F mass flux
            F_plus(:, :, 1) = density(0:imx-1, 1:jmx-1) * x_a * x_c_plus
            ! Construct other fluxes in terms of the F mass flux
            F_plus(:, :, 2) = (F_plus(:, :, 1) * x_speed(0:imx-1, 1:jmx-1)) + &
                    (x_scrD_plus * pressure(0:imx-1, 1:jmx-1) * xnx)
            F_plus(:, :, 3) = (F_plus(:, :, 1) * y_speed(0:imx-1, 1:jmx-1)) + &
                    (x_scrD_plus * pressure(0:imx-1, 1:jmx-1) * xny)
            F_plus(:, :, 4) = F_plus(:, :, 1) * &
                    ((0.5 * (x_speed(0:imx-1, 1:jmx-1) ** 2. + &
                        y_speed(0:imx-1, 1:jmx-1) ** 2.)) + &
                    ((gm / (gm - 1.)) * pressure(0:imx-1, 1:jmx-1) / &
                        density(0:imx-1, 1:jmx-1)))

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
            F_minus(:, :, 1) = density(1:imx, 1:jmx-1) * x_a * x_c_minus
            ! Construct other fluxes in terms of the F mass flux
            F_minus(:, :, 2) = (F_minus(:, :, 1) * x_speed(1:imx, 1:jmx-1)) + &
                    (x_scrD_minus * pressure(1:imx, 1:jmx-1) * xnx)
            F_minus(:, :, 3) = (F_minus(:, :, 1) * y_speed(1:imx, 1:jmx-1)) + &
                    (x_scrD_minus * pressure(1:imx, 1:jmx-1) * xny)
            F_minus(:, :, 4) = F_minus(:, :, 1) * &
                    ((0.5 * (x_speed(1:imx, 1:jmx-1) ** 2. + &
                        y_speed(1:imx, 1:jmx-1) ** 2.)) + &
                    ((gm / (gm - 1.)) * pressure(1:imx, 1:jmx-1) / &
                        density(1:imx, 1:jmx-1)))

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
            real, dimension(imx-1, jmx) :: y_M_perp_left, y_M_perp_right
            real, dimension(imx-1, jmx) :: y_alpha_plus, y_alpha_minus
            real, dimension(imx-1, jmx) :: y_beta_left, y_beta_right
            real, dimension(imx-1, jmx) :: M_plus, M_minus
            real, dimension(imx-1, jmx) :: D_plus, D_minus
            real, dimension(imx-1, jmx) :: M_ldfss

            call dmsg(1, 'ldfss', 'compute_eta_face_quantities')

            ! Begin by applying the Van-Leer scheme

            ! Compute the '+' direction quantities
            y_M_perp_left = eta_face_normal_speeds('+') / y_a
            y_alpha_plus = 0.5 * (1.0 + sign(1.0, y_M_perp_left))
            y_beta_left = -max(0, 1 - floor(abs(y_M_perp_left)))
            M_plus = 0.25 * ((1. + y_M_perp_left) ** 2.)
            D_plus = 0.25 * ((1. + y_M_perp_left) ** 2.) * (2. - y_M_perp_left)
            y_c_plus = (y_alpha_plus * (1.0 + y_beta_left) * y_M_perp_left) - &
                    y_beta_left * M_plus
            y_scrD_plus = (y_alpha_plus * (1. + y_beta_left)) - (y_beta_left * D_plus)

            ! Compute the '-' direction quantities
            y_M_perp_right = eta_face_normal_speeds('-') / y_a
            y_alpha_minus = 0.5 * (1.0 - sign(1.0, y_M_perp_right))
            y_beta_right = -max(0, 1 - floor(abs(y_M_perp_right)))
            M_minus = - 0.25 * ((1 - y_M_perp_right) ** 2.)
            D_minus = 0.25 * ((1 - y_M_perp_right) ** 2.) * (2. + y_M_perp_right)
            y_c_minus = (y_alpha_minus * (1.0 + y_beta_right) * y_M_perp_right) - &
                    (y_beta_right * M_minus)
            y_scrD_minus = (y_alpha_minus * (1.0 + y_beta_right)) - &
                    (y_beta_right * D_minus)

            ! Update the Van-Leer computed speeds (x_c_[plus|minus])
            M_ldfss = 0.25 * y_beta_left * y_beta_right * &
                    (sqrt((y_M_perp_left ** 2. + y_M_perp_right ** 2.) * 0.5) &
                    - 1) ** 2.
            y_c_plus = y_c_plus - M_ldfss
            y_c_minus = y_c_minus + M_ldfss

            call dmsg(1, 'ldfss', 'compute_eta_face_quantities', 'Ended')

        end subroutine compute_eta_face_quantities

        function compute_G_plus() result(G_plus)

            implicit none
            real, dimension(imx-1, jmx, 4) :: G_plus

            ! First construct the G mass flux
            G_plus(:, :, 1) = density(1:imx-1, 0:jmx-1) * y_a * y_c_plus
            ! Make the convective G mass flux zero at the wall
            G_plus(:, 1, 1) = 0
            G_plus(:, jmx, 1) = 0
            ! Construct other fluxes in terms of the G mass flux
            G_plus(:, :, 2) = (G_plus(:, :, 1) * x_speed(1:imx-1, 0:jmx-1)) + &
                    (y_scrD_plus * pressure(1:imx-1, 0:jmx-1) * ynx)
            G_plus(:, :, 3) = (G_plus(:, :, 1) * y_speed(1:imx-1, 0:jmx-1)) + &
                    (y_scrD_plus * pressure(1:imx-1, 0:jmx-1) * yny)
            G_plus(:, :, 4) = G_plus(:, :, 1) * &
                    ((0.5 * (x_speed(1:imx-1, 0:jmx-1) ** 2. + &
                        y_speed(1:imx-1, 0:jmx-1) ** 2.)) + &
                    ((gm / (gm - 1.)) * pressure(1:imx-1, 0:jmx-1) / &
                        density(1:imx-1, 0:jmx-1)))

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
            G_minus(:, :, 1) = density(1:imx-1, 1:jmx) * y_a * y_c_minus
            ! Make the convective G mass flux zero at the wall
            G_minus(:, 1, 1) = 0
            G_minus(:, jmx, 1) = 0
            ! Construct other fluxes in terms of the G mass flux
            G_minus(:, :, 2) = (G_minus(:, :, 1) * x_speed(1:imx-1, 1:jmx)) + &
                    (y_scrD_minus * pressure(1:imx-1, 1:jmx) * ynx)
            G_minus(:, :, 3) = (G_minus(:, :, 1) * y_speed(1:imx-1, 1:jmx)) + &
                    (y_scrD_minus * pressure(1:imx-1, 1:jmx) * yny)
            G_minus(:, :, 4) = G_minus(:, :, 1) * &
                    ((0.5 * (x_speed(1:imx-1, 1:jmx) ** 2. + &
                        y_speed(1:imx-1, 1:jmx) ** 2.)) + &
                    ((gm / (gm - 1.)) * pressure(1:imx-1, 1:jmx) / &
                        density(1:imx-1, 1:jmx)))

            ! Multiply in the face areas
            G_minus(:, :, 1) = G_minus(:, :, 1) * yA
            G_minus(:, :, 2) = G_minus(:, :, 2) * yA
            G_minus(:, :, 3) = G_minus(:, :, 3) * yA
            G_minus(:, :, 4) = G_minus(:, :, 4) * yA

        end function compute_G_minus

        function compute_residue() result(residue)
            !-----------------------------------------------------------
            ! Compute the residue using the LDFSS(0) scheme
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx, jmx-1, 4) :: F_plus, F_minus
            real, dimension(imx-1, jmx, 4) :: G_plus, G_minus
            real, dimension(imx-1, jmx-1, 4) :: residue

            call dmsg(1, 'ldfss', 'compute_residue')

            call compute_xi_face_quantities()
            F_plus = compute_F_plus()
            F_minus = compute_F_minus()
            !print *, 'F_plus: ', F_plus
            if (any(isnan(F_plus))) stop
            !print *, 'F_minus: ', F_minus
            if (any(isnan(F_minus))) stop

            call compute_eta_face_quantities()
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

end module ldfss0
