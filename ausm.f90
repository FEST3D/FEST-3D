module ausm
    !-------------------------------------------------------------------
    ! AUSM (Advection Upstream Splitting Method) is a class of 
    ! flux-splitting schemes derived from the Van Leer scheme. It was 
    ! invented to get rid of the erroneous mass flux that Van Leer 
    ! generates.
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use state, only: n_var
    use grid, only: imx, jmx, kmx
    use geometry, only: xnx, xny, xnz, ynx, yny, ynz, znx, zny, znz, xA, yA, zA
    use van_leer, only: setup_scheme_VL => setup_scheme, &
            destroy_scheme_VL => destroy_scheme, &
            compute_residue_VL => compute_residue, &
            compute_xi_face_quantities_VL => compute_xi_face_quantities, &
            compute_eta_face_quantities_VL => compute_eta_face_quantities, &
            compute_tau_face_quantities_VL => compute_tau_face_quantities, &
            x_c_plus, x_c_minus, &
            y_c_plus, y_c_minus, &
            z_c_plus, z_c_minus

    implicit none
    private

    ! Public methods
    public :: setup_scheme
    public :: destroy_scheme
    public :: get_residue

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ausm', 'setup_scheme')

            call setup_scheme_VL()

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ausm', 'destroy_scheme')

            call destroy_scheme_VL()

        end subroutine destroy_scheme

        subroutine ausm_modify_xi_face_quantities()
            !-----------------------------------------------------------
            ! Update the Van-Leer computed speeds: x_c_plus & x_c_minus
            !-----------------------------------------------------------

            implicit none
            real, dimension(size(x_c_plus, 1), size(x_c_plus, 2)) :: temp

            temp = x_c_plus + x_c_minus
            x_c_plus = max(0., temp)
            x_c_minus = min(0., temp)

        end subroutine ausm_modify_xi_face_quantities

        subroutine ausm_modify_eta_face_quantities()
            !-----------------------------------------------------------
            ! Update the Van-Leer computed speeds: y_c_plus & y_c_minus
            !-----------------------------------------------------------

            implicit none
            real, dimension(size(y_c_plus, 1), size(y_c_plus, 2)) :: temp

            temp = y_c_plus + y_c_minus
            y_c_plus = max(0., temp)
            y_c_minus = min(0., temp)

        end subroutine ausm_modify_eta_face_quantities

        function get_residue() result(residue)
            !-----------------------------------------------------------
            ! Return the AUSM residue
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx-1, jmx-1, kmx-1, n_var) :: residue

            call compute_xi_face_quantities_VL()
            call compute_eta_face_quantities_VL()
            
            ! Update the face variables according to the AUSM specs
            call ausm_modify_xi_face_quantities()
            call ausm_modify_eta_face_quantities()

            residue = compute_residue_VL()

        end function get_residue

end module ausm
