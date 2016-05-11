module ldfss0
    !-------------------------------------------------------------------
    ! LDFSS is a class of flux-splitting schemes
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx
    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz, xA, yA, zA
    use van_leer, only: setup_scheme_VL => setup_scheme, &
            destroy_scheme_VL => destroy_scheme, &
            compute_residue_VL => compute_residue, &
            compute_xi_face_quantities_VL => compute_xi_face_quantities, &
            compute_eta_face_quantities_VL => compute_eta_face_quantities, &
            compute_tau_face_quantities_VL => compute_tau_face_quantities, &
            x_M_perp_left, x_M_perp_right, &
            y_M_perp_left, y_M_perp_right, &
            z_M_perp_left, z_M_perp_right, &
            x_beta_left, x_beta_right, &
            y_beta_left, y_beta_right, &
            z_beta_left, z_beta_right, &
            x_c_plus, x_c_minus, &
            y_c_plus, y_c_minus
            z_c_plus, z_c_minus

    implicit none
    private

    real, dimension(:, :, :), allocatable :: x_M_ldfss
    real, dimension(:, :, :), allocatable :: y_M_ldfss
    real, dimension(:, :, :), allocatable :: z_M_ldfss

    ! Public methods
    public :: setup_scheme
    public :: destroy_scheme
    public :: get_residue

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ldfss', 'setup_scheme')

            call setup_scheme_VL()

            call alloc(x_M_ldfss, 1, imx, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for x_M_ldfss.')
            call alloc(y_M_ldfss, 1, imx-1, 1, jmx, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for y_M_ldfss.')
            call alloc(z_M_ldfss, 1, imx-1, 1, jmx-1, 1, kmx, &
                    errmsg='Error: Unable to allocate memory for z_M_ldfss.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ldfss', 'destroy_scheme')

            call dealloc(x_M_ldfss)
            call dealloc(y_M_ldfss)
            call dealloc(z_M_ldfss)
            call destroy_scheme_VL()

        end subroutine destroy_scheme

        subroutine ldfss_modify_xi_face_quantities()
            !-----------------------------------------------------------
            ! Update the Van-Leer computed speeds: x_c_plus & x_c_minus
            !-----------------------------------------------------------

            implicit none

            x_M_ldfss = 0.25 * x_beta_left * x_beta_right * &
                    (sqrt((x_M_perp_left ** 2. + x_M_perp_right ** 2.) * 0.5) &
                    - 1) ** 2.
            x_c_plus = x_c_plus - x_M_ldfss
            x_c_minus = x_c_minus + x_M_ldfss

        end subroutine ldfss_modify_xi_face_quantities

        subroutine ldfss_modify_eta_face_quantities()
            !-----------------------------------------------------------
            ! Update the Van-Leer computed speeds: y_c_plus & y_c_minus
            !-----------------------------------------------------------

            implicit none

            y_M_ldfss = 0.25 * y_beta_left * y_beta_right * &
                    (sqrt((y_M_perp_left ** 2. + y_M_perp_right ** 2.) * 0.5) &
                    - 1) ** 2.
            y_c_plus = y_c_plus - y_M_ldfss
            y_c_minus = y_c_minus + y_M_ldfss

        end subroutine ldfss_modify_eta_face_quantities

        function get_residue() result(residue)
            !-----------------------------------------------------------
            ! Return the LDFSS(0) residue
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx-1, jmx-1, 1, kmx-1, n_var) :: residue

            call compute_xi_face_quantities_VL()
            call compute_eta_face_quantities_VL()
            call compute_tau_face_quantities_VL()
            
            ! Update the face variables according to the LDFSS(0) specs
            call ldfss_modify_xi_face_quantities()
            call ldfss_modify_eta_face_quantities()
            call ldfss_modify_tau_face_quantities()

            residue = compute_residue_VL()

        end function get_residue

end module ldfss0
