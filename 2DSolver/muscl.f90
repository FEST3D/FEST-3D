module muscl
    !-----------------------------------------------------------------
    ! MUSCL (Monotone Upwing Schemes for Scalar Conservation Laws is
    ! a scheme which replaces the piecewise constant approximation by
    ! reconstructing the states at the left and right side of each face.
    ! This is a one parameter upwind scheme which results in at most 3rd
    ! order accuracy.
    !
    ! The MUSCL scheme alone creates non-physical oscillations near 
    ! discontinuities like shocks. Hence, MUSCL is combined with
    ! some TVD (Total Variation Diminishing) to reduce such oscillations.
    ! TVD schemes also ensure that no new extrema of the state variables
    ! is created at the faces.
    !-----------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx
    use state, only: qp, pressure, pressure_inf

    implicit none
    private

    ! Private variables
    real, dimension(:, :, :), allocatable, target :: x_qp_left, &
        x_qp_right, y_qp_left, y_qp_right
    real, dimension(:, :, :), pointer :: f_qp_left, f_qp_right
    real, dimension(:, :), allocatable :: pdif
    real :: phi, kappa

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_muscl_states
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right

    contains
        
        subroutine setup_scheme()

        implicit none

        call dmsg(1, 'muscl', 'setup_muscl')

        phi = 1.0

        kappa = 1.0/3.0

        call alloc(x_qp_left, 1, imx, 1, jmx-1, 1, 4, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_left.')
        call alloc(x_qp_right, 1, imx, 1, jmx-1, 1, 4, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_right.')

        call alloc(y_qp_left, 1, imx-1, 1, jmx, 1, 4, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_left.')
        call alloc(y_qp_right, 1, imx-1, 1, jmx, 1, 4, &
            errmsg='Error: Unable to allocate memory for ' // &
                'x_qp_right.')

        call alloc(pdif, 0, imx, 0, jmx, &
            errmsg='Error: Unable to allocate memory for ' // &
                'pdif.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'muscl', 'destroy_muscl')

            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)

        end subroutine destroy_scheme

        subroutine compute_xi_face_states()

            implicit none

            integer :: i, j, k
            real :: psi1, psi2, fd, bd, r

            phi = 1.0
            kappa = -1.0

            do k = 1, 4
             do j = 1, jmx - 1
              do i = 1, imx - 1
                ! All faces interior only (even at boundaries)
                ! Hence (i=1, left) and (i=imx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                !TODO: Generalise TVD scheme functions
                fd = qp(i+1, j, k) - qp(i, j, k)
                bd = qp(i, j, k) - qp(i-1, j, k)
                r = fd / bd
             !  psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
                r = bd / fd
             !  psi2 = min(1., (3 - kappa) * r / (1 - kappa))
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))

                x_qp_left(i+1, j, k) = qp(i, j, k) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                x_qp_right(i, j, k) = qp(i, j, k) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
              end do
             end do
            end do

            do j = 1, jmx - 1
                ! Exterior boundaries
                x_qp_left(1, j, :) = 0.5 * (qp(0, j, :) + qp(1, j, :))
                x_qp_right(imx, j, :) = 0.5 * (qp(imx-1, j, :) + qp(imx, j, :))
            end do

        end subroutine compute_xi_face_states

        subroutine compute_eta_face_states()

            implicit none

            integer :: i, j, k
            real :: psi1, psi2, fd, bd, r

            do k = 1, 4
             do j = 1, jmx - 1
              do i = 1, imx - 1
                ! All faces interior only (even at boundaries)
                ! Hence (j=1, left) and (j=jmx, right) will be dealt separately
                ! Koren limiter for now
                ! From paper: delta: forward difference 'fd'
                !             nabla: backward difference 'bd'
                !TODO: Generalise TVD scheme functions
                fd = qp(i, j+1, k) - qp(i, j, k)
                bd = qp(i, j, k) - qp(i, j-1, k)
                r = fd / bd
                psi1 = max(0., min(2*r, (2 + r)/3., 2.))
             !  psi1 = min(1., (3 - kappa) * r / (1 - kappa))
                r = bd / fd
                psi2 = max(0., min(2*r, (2 + r)/3., 2.))
             !  psi2 = min(1., (3 - kappa) * r / (1 - kappa))

                y_qp_left(i, j+1, k) = qp(i, j, k) + 0.25*phi* &
                    (((1-kappa) * psi1 * bd) + ((1+kappa) * psi2 * fd))
                y_qp_right(i, j, k) = qp(i, j, k) - 0.25*phi* &
                    (((1+kappa) * psi1 * bd) + ((1-kappa) * psi2 * fd))
              end do
             end do
            end do

            do i = 1, imx - 1
                ! Exterior boundaries
                y_qp_left(i, 1, :) = 0.5 * (qp(i, 0, :) + qp(i, 1, :))
                y_qp_right(i, jmx, :) = 0.5 * (qp(i, jmx-1, :) + qp(i, jmx, :))
            end do

        end subroutine compute_eta_face_states

        subroutine pressure_based_switching(f_dir)

            implicit none
            ! Character can be x or y or z
            character, intent(in) :: f_dir
            integer :: i, j, i_end, j_end
            integer :: i_f, j_f  ! Flags to determine face direction
            real :: pd2

            call dmsg(1, 'muscl', 'pressure_based_switching')

            select case (f_dir)
                case ('x')
                    f_qp_left => x_qp_left
                    f_qp_right => x_qp_right
                    i_f = 1
                    j_f = 0
                    i_end = imx
                    j_end = jmx - 1
                case ('y')
                    f_qp_left => y_qp_left
                    f_qp_right => y_qp_right
                    i_f = 0
                    j_f = 1
                    i_end = imx - 1
                    j_end = jmx 
                case default
                    call dmsg(5, 'ppm', 'pressure_based_switching', &
                            'Direction not recognised')
                    stop
            end select

            ! i_end and j_end denote number of faces
            ! Total number of cells including ghost_cells is
            ! (i_end+1) * j_end for xi faces and i_end*(j_end+1) for
            ! eta faces. 

            ! Loop over cells (physical)
            do j = 1, jmx - 1
             do i = 1, imx - 1
                pd2 = abs(pressure(i + i_f*1, j + j_f*1) - &
                          pressure(i - i_f*1, j - j_f*1))
                pdif(i, j) = 1 - (pd2/(pd2 + pressure_inf))
             end do
            end do

            ! Update at ghost cells
            pdif((1-i_f):(1-i_f)*(imx-1), (1-j_f):(1-j_f)*(jmx-1)) = &
                pdif(1:imx-1 - i_f*(imx-2), 1:jmx-1 - j_f*(jmx-2))

            ! Loop over faces
            do j = 1, jmx - (1 - j_f)
             do i = 1, imx - (1 - i_f)
                f_qp_left(i, j, :) = qp(i - i_f*1, j - j_f*1, :) + (&
                    pdif(i - i_f*1, j - j_f*1) * ( &
                    f_qp_left(i, j, :) - qp(i - i_f*1, j - j_f*1, :)))

                f_qp_right(i, j, :) = qp(i, j, :) - (&
                    pdif(i, j) * ( &
                    qp(i, j, :) - f_qp_right(i, j, :)))
             end do
            end do

        end subroutine pressure_based_switching
        
        subroutine compute_muscl_states()
            !---------------------------------------------------------
            ! Implement MUSCL scheme to get left and right states at
            ! each face
            !---------------------------------------------------------
            
            call compute_xi_face_states()
            call pressure_based_switching('x')

            call compute_eta_face_states()
            call pressure_based_switching('y')

        end subroutine compute_muscl_states

end module muscl
