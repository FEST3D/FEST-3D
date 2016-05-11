module weno_reconstruction
    !-------------------------------------------------------------------
    ! WENO (Weighted Essentially Non-Oscillatory) is an approximation 
    ! technique. To solve the hyperbolic conservation equations, this 
    ! code implements the WENO reconstruction procedure.
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx, kmx
    use state, only: qp, n_var

    implicit none
    private

    real, dimension(:, :, :), allocatable :: x_qp_face_estimate
    real, dimension(:, :, :), allocatable :: y_qp_face_estimate
    real, dimension(:, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :), allocatable, target :: y_qp_left, y_qp_right

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_weno_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right

    contains

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'weno_reconstruction', 'setup_scheme')

            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'weno_reconstruction', 'destroy_scheme')

            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)

        end subroutine destroy_scheme

        subroutine compute_xi_approximations()

            implicit none

        end subroutine compute_xi_approximations

        subroutine compute_eta_approximations()

            implicit none

        end subroutine compute_eta_approximations

        subroutine compute_weno_states()

            implicit none

            call compute_xi_approximations()
            call compute_eta_approximations()

        end subroutine compute_weno_states

end module weno_reconstruction
