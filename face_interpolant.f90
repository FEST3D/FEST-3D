module face_interpolant

    use global, only: INTERPOLANT_NAME_LENGTH
    use utils, only: alloc, dealloc, dmsg
    use grid, only: imx, jmx
    use state, only: qp, gm
    use ppm, only: setup_scheme_ppm => setup_scheme, &
            destroy_scheme_ppm => destroy_scheme, &
            compute_ppm_states, &
            x_qp_left_ppm => x_qp_left, &
            x_qp_right_ppm => x_qp_right, &
            y_qp_left_ppm => y_qp_left, &
            y_qp_right_ppm => y_qp_right

    implicit none
    private

    character(len=INTERPOLANT_NAME_LENGTH) :: interpolant

    real, dimension(:, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :), allocatable, target :: y_qp_left, y_qp_right
    real, dimension(:, :), pointer :: x_density_left, x_density_right
    real, dimension(:, :), pointer :: x_x_speed_left, x_x_speed_right
    real, dimension(:, :), pointer :: x_y_speed_left, x_y_speed_right
    real, dimension(:, :), pointer :: x_pressure_left, x_pressure_right
    real, dimension(:, :), pointer :: y_density_left, y_density_right
    real, dimension(:, :), pointer :: y_x_speed_left, y_x_speed_right
    real, dimension(:, :), pointer :: y_y_speed_left, y_y_speed_right
    real, dimension(:, :), pointer :: y_pressure_left, y_pressure_right

    ! Public members
    public :: interpolant
    public :: setup_interpolant_scheme
    public :: destroy_interpolant_scheme
    public :: compute_face_interpolant
    public :: x_sound_speed_left
    public :: x_sound_speed_right
    public :: y_sound_speed_left
    public :: y_sound_speed_right
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: x_density_left, x_density_right
    public :: x_x_speed_left, x_x_speed_right
    public :: x_y_speed_left, x_y_speed_right
    public :: x_pressure_left, x_pressure_right
    public :: y_density_left, y_density_right
    public :: y_x_speed_left, y_x_speed_right
    public :: y_y_speed_left, y_y_speed_right
    public :: y_pressure_left, y_pressure_right

    contains

        subroutine allocate_memory()
            implicit none
            call alloc(x_qp_left, 1, imx, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 1, imx, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 1, jmx, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 1, jmx, 1, 4, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')
        end subroutine allocate_memory

        subroutine link_aliases()
            implicit none

            ! Link xi faces left pointers
            x_density_left(1:imx, 1:jmx-1) => x_qp_left(:, :, 1)
            x_x_speed_left(1:imx, 1:jmx-1) => x_qp_left(:, :, 2)
            x_y_speed_left(1:imx, 1:jmx-1) => x_qp_left(:, :, 3)
            x_pressure_left(1:imx, 1:jmx-1) => x_qp_left(:, :, 4)

            ! Link xi faces right pointers
            x_density_right(1:imx, 1:jmx-1) => x_qp_right(:, :, 1)
            x_x_speed_right(1:imx, 1:jmx-1) => x_qp_right(:, :, 2)
            x_y_speed_right(1:imx, 1:jmx-1) => x_qp_right(:, :, 3)
            x_pressure_right(1:imx, 1:jmx-1) => x_qp_right(:, :, 4)

            ! Link eta faces left pointers
            y_density_left(1:imx-1, 1:jmx) => y_qp_left(:, :, 1)
            y_x_speed_left(1:imx-1, 1:jmx) => y_qp_left(:, :, 2)
            y_y_speed_left(1:imx-1, 1:jmx) => y_qp_left(:, :, 3)
            y_pressure_left(1:imx-1, 1:jmx) => y_qp_left(:, :, 4)

            ! Link eta faces right pointers
            y_density_right(1:imx-1, 1:jmx) => y_qp_right(:, :, 1)
            y_x_speed_right(1:imx-1, 1:jmx) => y_qp_right(:, :, 2)
            y_y_speed_right(1:imx-1, 1:jmx) => y_qp_right(:, :, 3)
            y_pressure_right(1:imx-1, 1:jmx) => y_qp_right(:, :, 4)

        end subroutine link_aliases

        subroutine setup_interpolant_scheme()
            implicit none
            select case (interpolant)
                case ("none")
                    ! Do nothing
                    continue
                case ("ppm")
                    call setup_scheme_ppm()
                case default
                    call dmsg(5, 'state_interpolant', &
                            'setup_interpolant_scheme', &
                            'Interpolant not recognized.')
                    stop
            end select
            call allocate_memory()
            call link_aliases()
        end subroutine setup_interpolant_scheme

        subroutine deallocate_memory()
            implicit none
            call dealloc(x_qp_left)
            call dealloc(x_qp_right)
            call dealloc(y_qp_left)
            call dealloc(y_qp_right)
        end subroutine deallocate_memory

        subroutine unlink_aliases()
            implicit none

            ! Unlink xi faces left pointers
            nullify(x_density_left)
            nullify(x_x_speed_left)
            nullify(x_y_speed_left)
            nullify(x_pressure_left)

            ! Unlink xi faces right pointers
            nullify(x_density_right)
            nullify(x_x_speed_right)
            nullify(x_y_speed_right)
            nullify(x_pressure_right)

            ! Unlink eta faces left pointers
            nullify(y_density_left)
            nullify(y_x_speed_left)
            nullify(y_y_speed_left)
            nullify(y_pressure_left)

            ! Unlink eta faces right pointers
            nullify(y_density_right)
            nullify(y_x_speed_right)
            nullify(y_y_speed_right)
            nullify(y_pressure_right)

        end subroutine unlink_aliases

        subroutine destroy_interpolant_scheme()
            implicit none
            call unlink_aliases()
            call deallocate_memory()
            select case (interpolant)
                case ("none")
                    ! Do nothing
                    continue
                case ("ppm")
                    call destroy_scheme_ppm()
                case default
                    call dmsg(5, 'state_interpolant', &
                            'destroy_interpolant_scheme', &
                            'Interpolant not recognized.')
                    stop
            end select
        end subroutine destroy_interpolant_scheme

        subroutine extrapolate_cell_averages_to_faces()
            implicit none
            x_qp_left(:, :, :) = qp(0:imx-1, 1:jmx-1, 1:4)
            x_qp_right(:, :, :) = qp(1:imx, 1:jmx-1, 1:4)
            y_qp_left(:, :, :) = qp(1:imx-1, 0:jmx-1, 1:4)
            y_qp_right(:, :, :) = qp(1:imx-1, 1:jmx, 1:4)
        end subroutine extrapolate_cell_averages_to_faces

        subroutine compute_face_interpolant()
            implicit none
            select case (interpolant)
                case ("none")
                    call extrapolate_cell_averages_to_faces()
                case ("ppm")
                    call compute_ppm_states()
                    x_qp_left(:, :, :) = x_qp_left_ppm(1:imx, :, :)
                    x_qp_right(:, :, :) = x_qp_right_ppm(1:imx, :, :)
                    y_qp_left(:, :, :) = y_qp_left_ppm(:, 1:jmx, :)
                    y_qp_right(:, :, :) = y_qp_right_ppm(:, 1:jmx, :)
                case default
                    call dmsg(5, 'state_interpolant', &
                            'compute_face_interpolant', &
                            'Interpolant not recognized.')
                    stop
            end select
        end subroutine compute_face_interpolant

        function x_sound_speed_left()
            implicit none
            real, dimension(1:imx, 1:jmx-1) :: x_sound_speed_left
            x_sound_speed_left = sqrt(gm * x_pressure_left / x_density_left)
        end function x_sound_speed_left

        function x_sound_speed_right()
            implicit none
            real, dimension(1:imx, 1:jmx-1) :: x_sound_speed_right
            x_sound_speed_right = sqrt(gm * x_pressure_right / x_density_right)
        end function x_sound_speed_right

        function y_sound_speed_left()
            implicit none
            real, dimension(1:imx-1, 1:jmx) :: y_sound_speed_left
            y_sound_speed_left = sqrt(gm * y_pressure_left / y_density_left)
        end function y_sound_speed_left

        function y_sound_speed_right()
            implicit none
            real, dimension(1:imx-1, 1:jmx) :: y_sound_speed_right
            y_sound_speed_right = sqrt(gm * y_pressure_right / y_density_right)
        end function y_sound_speed_right

end module face_interpolant
