module ausm
    !-------------------------------------------------------------------
    ! AUSM (Advection Upstream Splitting Method) is a class of 
    ! flux-splitting schemes derived from the Van Leer scheme. It was 
    ! invented to get rid of the erroneous mass flux that Van Leer 
    ! generates.
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx
    use van_leer, only: setup_scheme_VL => setup_scheme, &
            destroy_scheme_VL => destroy_scheme, &
            get_residue_VL => get_residue, &
            compute_face_quantities_VL => compute_face_quantities, &
            compute_fluxes_VL => compute_fluxes, &
            x_c_plus, x_c_minus, &
            y_c_plus, y_c_minus, &
            F_van_leer => F, G_van_leer => G

    implicit none
    private

    real, public, dimension(:, :, :), pointer :: F, G

    ! Public methods
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_face_quantities
    public :: compute_fluxes
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
            call dmsg(1, 'ausm', 'modify_xi_face_quantities')

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
            call dmsg(1, 'ausm', 'modify_eta_face_quantities')

            temp = y_c_plus + y_c_minus
            y_c_plus = max(0., temp)
            y_c_minus = min(0., temp)

        end subroutine ausm_modify_eta_face_quantities

        subroutine compute_face_quantities()

            implicit none
            call dmsg(1, 'ausm', 'compute_face_quantities')

            call compute_face_quantities_VL()
            
            ! Update the face variables according to the AUSM specs
            call ausm_modify_xi_face_quantities()
            call ausm_modify_eta_face_quantities()

        end subroutine compute_face_quantities

        subroutine compute_fluxes()

            implicit none

            call dmsg(1, 'ausm', 'compute_fluxes')
            call compute_fluxes_VL()
            F => F_van_leer
            G => G_van_leer

        end subroutine compute_fluxes

        function get_residue() result(residue)
            !-----------------------------------------------------------
            ! Return the AUSM residue
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx-1, jmx-1, 4) :: residue
            call dmsg(1, 'ausm', 'compute_residue')

            residue = get_residue_VL()

        end function get_residue

end module ausm
