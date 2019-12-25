module face_interpolant

#include "../../debug.h"
#include "../../error.h"
    use vartypes
    use utils  , only: alloc
    use muscl  , only: compute_muscl_states
    use ppm    , only: compute_ppm_states
    use weno   , only: compute_weno_states
    use weno_NM, only: compute_weno_NM_states
    implicit none
    private


    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :, :), allocatable, target :: y_qp_left, y_qp_right
    real, dimension(:, :, :, :), allocatable, target :: z_qp_left, z_qp_right
    real, dimension(:, :, :), allocatable :: pdif
   !< Used for pressure based witch


    ! Public members
    public :: setup_interpolant_scheme
    public :: compute_face_interpolant
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right
    public :: pdif

    contains

!        subroutine allocate_memory()
!            implicit none
!            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'x_qp_left.')
!            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'x_qp_right.')
!            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'y_qp_left.')
!            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'y_qp_right.')
!            call alloc(z_qp_left, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'z_qp_left.')
!            call alloc(z_qp_right, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for ' // &
!                        'z_qp_right.')
!        call alloc(pdif, 0, imx, 0, jmx, 0, kmx, &
!                errmsg='Error: Unable to allocate memory for' // &
!                    'pdif')
!        end subroutine allocate_memory


        subroutine setup_interpolant_scheme(dims)
            implicit none
            !type(controltype), intent(in) :: control
            !type(schemetype), intent(in) :: scheme
            type(extent), intent(in) :: dims

            call alloc(x_qp_left, 0, dims%imx+1, 1, dims%jmx-1, 1, dims%kmx-1, 1, dims%n_var, &
                    errmsg='Error: Unable to allocate memory for x_qp_left.')
            call alloc(x_qp_right, 0, dims%imx+1, 1, dims%jmx-1, 1, dims%kmx-1, 1, dims%n_var, &
                    errmsg='Error: Unable to allocate memory for x_qp_right.')
            call alloc(y_qp_left, 1, dims%imx-1, 0, dims%jmx+1, 1, dims%kmx-1, 1, dims%n_var, &
                    errmsg='Error: Unable to allocate memory for y_qp_left.')
            call alloc(y_qp_right, 1, dims%imx-1, 0, dims%jmx+1, 1, dims%kmx-1, 1, dims%n_var, &
                    errmsg='Error: Unable to allocate memory for y_qp_right.')
            call alloc(z_qp_left, 1, dims%imx-1, 1, dims%jmx-1, 0, dims%kmx+1, 1, dims%n_var, &
                    errmsg='Error: Unable to allocate memory for z_qp_left.')
            call alloc(z_qp_right, 1, dims%imx-1, 1, dims%jmx-1, 0, dims%kmx+1, 1, dims%n_var, &
                    errmsg='Error: Unable to allocate memory for z_qp_right.')
            call alloc(pdif, 0, dims%imx, 0, dims%jmx, 0, dims%kmx, &
                    errmsg='Error: Unable to allocate memory for pdif')

        end subroutine setup_interpolant_scheme


        subroutine extrapolate_cell_averages_to_faces(qp, dims)
            implicit none
            type(extent), intent(in) :: dims
            real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp

            DebugCall('extrapolate_cell_averages_to_faces')

            x_qp_left(:, :, :, :) = qp(-1:dims%imx, 1:dims%jmx-1, 1:dims%kmx-1, 1:dims%n_var)
            x_qp_right(:, :, :, :) = qp(0:dims%imx+1, 1:dims%jmx-1, 1:dims%kmx-1, 1:dims%n_var)
            y_qp_left(:, :, :, :) = qp(1:dims%imx-1, -1:dims%jmx, 1:dims%kmx-1, 1:dims%n_var)
            y_qp_right(:, :, :, :) = qp(1:dims%imx-1, 0:dims%jmx+1, 1:dims%kmx-1, 1:dims%n_var)
            z_qp_left(:, :, :, :) = qp(1:dims%imx-1, 1:dims%jmx-1, -1:dims%kmx, 1:dims%n_var)
            z_qp_right(:, :, :, :) = qp(1:dims%imx-1, 1:dims%jmx-1, 0:dims%kmx+1, 1:dims%n_var)
        end subroutine extrapolate_cell_averages_to_faces


        subroutine compute_face_interpolant(qp, cells, scheme, flow, dims)
            implicit none
            type(extent), intent(in) :: dims
            real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in):: qp
            type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
            type(schemetype), intent(in) :: scheme
            type(flowtype), intent(in) :: flow
            select case (scheme%interpolant)
                case ("none")
                    call extrapolate_cell_averages_to_faces(qp, dims)
                case ("ppm")
                    call compute_ppm_states(qp,  x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, pdif, scheme, flow, dims)
                case ("muscl")
                    call compute_muscl_states(qp, x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, pdif, scheme, flow, dims)
                case ("weno")
                    call compute_weno_states(qp, x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, dims)
                case ("weno_NM")
                    call compute_weno_NM_states(qp, x_qp_left, x_qp_right, y_qp_left, y_qp_right, z_qp_left, z_qp_right, cells, dims)
                case default
                    Fatal_error
            end select
        end subroutine compute_face_interpolant

end module face_interpolant
