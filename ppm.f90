module ppm
    !-------------------------------------------------------------------
    ! PPM (Piecewise Parabolic Method) is an interpolation technique
    ! which can be used to create higher order extensions of schemes
    ! like the Van-Leer and LDFSS schemes.
    !-------------------------------------------------------------------

    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx
    use state, only: qp

    implicit none
    private

    real, dimension(:, :, :), allocatable :: x_qp_face_estimate
    real, dimension(:, :, :), allocatable :: y_qp_face_estimate
    real, dimension(:, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :), allocatable, target :: y_qp_left, y_qp_right

    ! Public members
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_ppm_states
    public :: output_data
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right

    contains

        function xmd(a, b)

            implicit none
            real, dimension(:, :), intent(in) :: a, b
            real, dimension(size(a, 1), size(a, 2)) :: xmd

            xmd = 0.5 * (sign(1., a) + sign(1., b)) * min(abs(a), abs(b))

        end function xmd

        subroutine setup_scheme()

            implicit none

            call dmsg(1, 'ppm', 'setup_ppm')

            call alloc(x_qp_face_estimate, 0, imx+1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_face_estimate.')
            call alloc(y_qp_face_estimate, 1, imx-1, 0, jmx+1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_face_estimate.')
            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')

        end subroutine setup_scheme

        subroutine destroy_scheme()

            implicit none

            call dmsg(1, 'ppm', 'destroy_ppm')

            call dealloc(x_qp_face_estimate)
            call dealloc(y_qp_face_estimate)
            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)

        end subroutine destroy_scheme

        subroutine compute_xi_face_estimate()
            !-----------------------------------------------------------
            ! Estimate the value of the state at the xi cell interfaces
            !-----------------------------------------------------------

            implicit none
            integer :: i, j
            integer, dimension(imx - 2) :: ifaces
            integer, dimension(jmx - 1) :: jrows

            ifaces = (/ (i, i = 2, imx-1) /)
            jrows = (/ (j, j = 1, jmx-1) /)

            ! Estimate face values for interior faces
            x_qp_face_estimate(ifaces, jrows, :) = &
                    ((7. * (qp(ifaces, jrows, :) + qp(ifaces-1, jrows, :))) - &
                    (qp(ifaces+1, jrows, :) + qp(ifaces-2, jrows, :))) / 12.

            ! Estimate face values for other faces
            x_qp_face_estimate(1, jrows, :) = &
                    0.5 * (qp(0, jrows, :) + qp(1, jrows, :))
            x_qp_face_estimate(0, jrows, :) = &
                    (4. * qp(0, jrows, :) - qp(1, jrows, :)) / 3.
            x_qp_face_estimate(imx, jrows, :) = &
                    0.5 * (qp(imx, jrows, :) + qp(imx-1, jrows, :))
            x_qp_face_estimate(imx+1, jrows, :) = &
                    (4. * (qp(imx, jrows, :) - qp(imx-1, jrows, :))) / 3.

        end subroutine compute_xi_face_estimate

        subroutine init_left_and_right_xi_estimates()

            implicit none

            ! x_qp_left and x_qp_right are stored at faces.
            x_qp_left = x_qp_face_estimate
            x_qp_right = x_qp_face_estimate

        end subroutine init_left_and_right_xi_estimates

        subroutine remove_xi_extrema()

            implicit none
            integer :: i, j, k
            ! dqrl is the difference (d) in state (qp) between right (r)
            ! and left (l) estimates for a cell
            ! The left estimate for cell 'i' is the right estimate for
            ! face 'i'. The right estimate for cell 'i' is the left
            ! estimate for face 'i+1'.
            real :: dqrl
            ! dq6 is the difference (d) in the state (qp) between the
            ! center and the average of the edges times 6
            real :: dq6
            ! Both these parameters are as defined in the original paper
            ! on PPM by Collela and Woodward

            !TODO: Can finding the number of flow variables (4) in the
            ! parameter vector be made more generic?
            do k = 1, 4
                do j = 1, jmx - 1
                    do i = 0, imx
                        if ((x_qp_right(i, j, k) - qp(i, j, k)) * &
                                (qp(i, j, k) - x_qp_left(i+1, j, k)) <= 0) then
                            x_qp_right(i, j, k) = qp(i, j, k)
                            x_qp_left(i+1, j, k) = x_qp_right(i, j, k)
                        else
                            dqrl = x_qp_left(i+1, j, k) - x_qp_right(i, j, k)
                            dq6 = 6. * (qp(i, j, k) - &
                                    (0.5 * (x_qp_right(i, j, k) + &
                                        x_qp_left(i+1, j, k))))
                            if (dqrl * dq6 > dqrl ** 2.) then
                                x_qp_right(i, j, k) = 3. * qp(i, j, k) - &
                                        2. * x_qp_left(i+1, j, k)
                            else if (-(dqrl ** 2.)  > dqrl * dq6) then
                                x_qp_left(i+1, j, k) = 3. * qp(i, j, k) - &
                                        2. * x_qp_right(i, j, k)
                            end if
                        end if
                    end do
                end do
            end do

            if (any(isnan(x_qp_left))) then
                call dmsg(5, 'ppm', 'remove_xi_extrema', &
                        msg='ERROR: NaN detected in x_qp_left.')
                stop
            else if (any(isnan(x_qp_right))) then
                call dmsg(5, 'ppm', 'remove_xi_extrema', &
                        msg='ERROR: NaN detected in x_qp_right.')
                stop
            end if

        end subroutine remove_xi_extrema

        subroutine reset_xi_edge_values()

            implicit none

            integer :: j
            integer, dimension(jmx - 1) :: jrows
            jrows = (/ (j, j = 1, jmx-1) /)

            ! Reset outer boundaries of domain
            x_qp_left(1, jrows, :) = 0.5 * &
                    (qp(0, jrows, :) + qp(1, jrows, :))
            x_qp_right(imx, jrows, :) = 0.5 * &
                    (qp(imx, jrows, :) + qp(imx-1, jrows, :))

            ! Reset inner boundaries of domain
            x_qp_right(1, jrows, :) = qp(1, jrows, :) - &
                    0.5 * xmd((qp(2, jrows, :) - qp(1, jrows, :)), &
                    (qp(1, jrows, :) - qp(0, jrows, :)))
            x_qp_left(imx, jrows, :) = qp(imx-1, jrows, :) + &
                    0.5 * xmd((qp(imx-1, jrows, :) - qp(imx-2, jrows, :)), &
                    (qp(imx, jrows, :) - qp(imx-1, jrows, :)))

        end subroutine reset_xi_edge_values

        subroutine compute_eta_face_estimate()
            !-----------------------------------------------------------
            ! Estimate the value of the state at the eta cell interfaces
            !-----------------------------------------------------------

            implicit none
            integer :: i, j
            integer, dimension(imx - 1) :: irows
            integer, dimension(jmx - 2) :: jfaces

            irows = (/ (i, i = 1, imx-1) /)
            jfaces = (/ (j, j = 2, jmx-1) /)

            ! Estimate face values for interior faces
            y_qp_face_estimate(irows, jfaces, :) = &
                    ((7. * (qp(irows, jfaces, :) + qp(irows, jfaces-1, :))) - &
                    (qp(irows, jfaces+1, :) + qp(irows, jfaces-2, :))) / 12.

            ! Estimate face values for other faces
            y_qp_face_estimate(irows, 1, :) = &
                    0.5 * (qp(irows, 0, :) + qp(irows, 1, :))
            y_qp_face_estimate(irows, 0, :) = &
                    (4. * qp(irows, 0, :) - qp(irows, 1, :)) / 3.
            y_qp_face_estimate(irows, jmx, :) = &
                    0.5 * (qp(irows, jmx, :) + qp(irows, jmx-1, :))
            y_qp_face_estimate(irows, jmx+1, :) = &
                    (4. * (qp(irows, jmx, :) - qp(irows, jmx-1, :))) / 3.

        end subroutine compute_eta_face_estimate

        subroutine init_left_and_right_eta_estimates()

            implicit none

            ! y_qp_left and y_qp_right are stored at faces.
            y_qp_left = y_qp_face_estimate
            y_qp_right = y_qp_face_estimate

        end subroutine init_left_and_right_eta_estimates

        subroutine remove_eta_extrema()

            implicit none
            integer :: i, j, k
            ! dqrl is the difference (d) in state (qp) between right and
            ! left (rl) estimates
            real :: dqrl
            ! dq6 is the difference (d) in the state (qp) between the
            ! center and the average of the edges times 6
            real :: dq6
            ! Both these parameters are as defined in the original paper
            ! on PPM by Collela and Woodward

            !TODO: Can finding the number of flow variables (4) in the
            ! parameter vector be made more generic?
            do k = 1, 4
                do j = 0, jmx
                    do i = 1, imx - 1
                        if ((y_qp_right(i, j, k) - qp(i, j, k)) * &
                                (qp(i, j, k) - y_qp_left(i, j+1, k)) <= 0) then
                            y_qp_right(i, j, k) = qp(i, j, k)
                            y_qp_left(i, j+1, k) = y_qp_right(i, j, k)
                        else
                            dqrl = y_qp_left(i, j+1, k) - y_qp_right(i, j, k)
                            dq6 = 6. * (qp(i, j, k) - &
                                    (0.5 * (y_qp_right(i, j, k) + &
                                        y_qp_left(i, j+1, k))))
                            if (dqrl * dq6 > dqrl ** 2.) then
                                y_qp_right(i, j, k) = 3. * qp(i, j, k) - &
                                        2. * y_qp_left(i, j+1, k)
                            else if (-(dqrl ** 2.)  > dqrl * dq6) then
                                y_qp_left(i, j+1, k) = 3. * qp(i, j, k) - &
                                        2. * y_qp_right(i, j, k)
                            end if
                        end if
                    end do
                end do
            end do

            if (any(isnan(y_qp_left))) then
                call dmsg(5, 'ppm', 'remove_eta_extrema', &
                        msg='ERROR: NaN detected in y_qp_left.')
                stop
            else if (any(isnan(y_qp_right))) then
                call dmsg(5, 'ppm', 'remove_eta_extrema', &
                        msg='ERROR: NaN detected in y_qp_right.')
                stop
            end if

        end subroutine remove_eta_extrema

        subroutine reset_eta_edge_values()

            implicit none

            integer :: i
            integer, dimension(imx - 1) :: irows
            irows = (/ (i, i = 1, imx-1) /)

            ! Reset outer boundaries of domain
            y_qp_left(irows, 1, :) = 0.5 * &
                    (qp(irows, 0, :) + qp(irows, 1, :))
            y_qp_right(irows, jmx, :) = 0.5 * &
                    (qp(irows, jmx, :) + qp(irows, jmx-1, :))

            ! Reset inner boundaries of domain
            y_qp_right(irows, 1, :) = qp(irows, 1, :) - &
                    0.5 * xmd((qp(irows, 2, :) - qp(irows, 1, :)), &
                    (qp(irows, 1, :) - qp(irows, 0, :)))
            y_qp_left(irows, jmx, :) = qp(irows, jmx-1, :) + &
                    0.5 * xmd((qp(irows, jmx-1, :) - qp(irows, jmx-2, :)), &
                    (qp(irows, jmx, :) - qp(irows, jmx-1, :)))

        end subroutine reset_eta_edge_values

        subroutine compute_ppm_states()

            implicit none

            call compute_xi_face_estimate()
            call init_left_and_right_xi_estimates()
            call remove_xi_extrema

            call compute_eta_face_estimate()
            call init_left_and_right_eta_estimates()
            call remove_eta_extrema

            call reset_xi_edge_values()
            call reset_eta_edge_values()

        end subroutine compute_ppm_states

        subroutine output_data(filename)
            implicit none
            character(len=*), intent(in) :: filename
            integer :: i, j
            open(71, file=filename)
            write(71, *) 'density'
            write(71, *) 'xi_left'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_left(i, j, 1)
                end do
            end do
            write(71, *) 'xi_right'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_right(i, j, 1)
                end do
            end do
            write(71, *) 'eta_left'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_left(i, j, 1)
                end do
            end do
            write(71, *) 'eta_right'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_right(i, j, 1)
                end do
            end do
            write(71, *) 'x_speed'
            write(71, *) 'xi_left'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_left(i, j, 2)
                end do
            end do
            write(71, *) 'xi_right'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_right(i, j, 2)
                end do
            end do
            write(71, *) 'eta_left'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_left(i, j, 2)
                end do
            end do
            write(71, *) 'eta_right'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_right(i, j, 2)
                end do
            end do
            write(71, *) 'y_speed'
            write(71, *) 'xi_left'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_left(i, j, 3)
                end do
            end do
            write(71, *) 'xi_right'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_right(i, j, 3)
                end do
            end do
            write(71, *) 'eta_left'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_left(i, j, 3)
                end do
            end do
            write(71, *) 'eta_right'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_right(i, j, 3)
                end do
            end do
            write(71, *) 'pressure'
            write(71, *) 'xi_left'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_left(i, j, 4)
                end do
            end do
            write(71, *) 'xi_right'
            do j = 1, jmx-1
                do i = 0, imx+1
                    write(71, *) x_qp_right(i, j, 4)
                end do
            end do
            write(71, *) 'eta_left'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_left(i, j, 4)
                end do
            end do
            write(71, *) 'eta_right'
            do j = 0, jmx+1
                do i = 1, imx-1
                    write(71, *) y_qp_right(i, j, 4)
                end do
            end do
            close(71)
        end subroutine output_data

end module ppm
