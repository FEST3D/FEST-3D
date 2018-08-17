module face_interpolant

    use global, only: INTERPOLANT_NAME_LENGTH

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx

    use global_vars, only : qp
    use global_vars, only : gm
    use global_vars, only : n_var
    use global_vars, only : turbulence

    use global_vars, only : interpolant

    use utils, only: alloc, dealloc, dmsg
    use muscl, only: setup_scheme_muscl => setup_scheme, &
            destroy_scheme_muscl => destroy_scheme, &
            compute_muscl_states, &
            x_qp_left_muscl => x_qp_left, &
            x_qp_right_muscl => x_qp_right, &
            y_qp_left_muscl => y_qp_left, &
            y_qp_right_muscl => y_qp_right, &
            z_qp_left_muscl => z_qp_left, &
            z_qp_right_muscl => z_qp_right
    use ppm, only: setup_scheme_ppm => setup_scheme, &
            destroy_scheme_ppm => destroy_scheme, &
            compute_ppm_states, &
            x_qp_left_ppm => x_qp_left, &
            x_qp_right_ppm => x_qp_right, &
            y_qp_left_ppm => y_qp_left, &
            y_qp_right_ppm => y_qp_right, &
            z_qp_left_ppm => z_qp_left, &
            z_qp_right_ppm => z_qp_right
    use weno, only: setup_scheme_weno => setup_scheme, &
            destroy_scheme_weno => destroy_scheme, &
            compute_weno_states, &
             x_qp_left_weno => x_qp_left, &
            x_qp_right_weno => x_qp_right, &
             y_qp_left_weno => y_qp_left, &
            y_qp_right_weno => y_qp_right, &
             z_qp_left_weno => z_qp_left, &
            z_qp_right_weno => z_qp_right
    use weno_NM, only: setup_scheme_weno_NM => setup_scheme, &
            destroy_scheme_weno_NM => destroy_scheme, &
            compute_weno_NM_states, &
             x_qp_left_weno_NM => x_qp_left, &
            x_qp_right_weno_NM => x_qp_right, &
             y_qp_left_weno_NM => y_qp_left, &
            y_qp_right_weno_NM => y_qp_right, &
             z_qp_left_weno_NM => z_qp_left, &
            z_qp_right_weno_NM => z_qp_right
    include "turbulence_models/include/face_interpolant/import_module.inc"

    implicit none
    private


    real, dimension(:, :, :, :), allocatable, target :: x_qp_left, x_qp_right
    real, dimension(:, :, :, :), allocatable, target :: y_qp_left, y_qp_right
    real, dimension(:, :, :, :), allocatable, target :: z_qp_left, z_qp_right
    real, dimension(:, :, :), pointer :: x_density_left, x_density_right
    real, dimension(:, :, :), pointer :: x_x_speed_left, x_x_speed_right
    real, dimension(:, :, :), pointer :: x_y_speed_left, x_y_speed_right
    real, dimension(:, :, :), pointer :: x_z_speed_left, x_z_speed_right
    real, dimension(:, :, :), pointer :: x_pressure_left, x_pressure_right
    real, dimension(:, :, :), pointer :: y_density_left, y_density_right
    real, dimension(:, :, :), pointer :: y_x_speed_left, y_x_speed_right
    real, dimension(:, :, :), pointer :: y_y_speed_left, y_y_speed_right
    real, dimension(:, :, :), pointer :: y_z_speed_left, y_z_speed_right
    real, dimension(:, :, :), pointer :: y_pressure_left, y_pressure_right
    real, dimension(:, :, :), pointer :: z_density_left, z_density_right
    real, dimension(:, :, :), pointer :: z_x_speed_left, z_x_speed_right
    real, dimension(:, :, :), pointer :: z_y_speed_left, z_y_speed_right
    real, dimension(:, :, :), pointer :: z_z_speed_left, z_z_speed_right
    real, dimension(:, :, :), pointer :: z_pressure_left, z_pressure_right

    !turbulent variable left and right with public deceleration
    include "turbulence_models/include/face_interpolant/variables_deceleration.inc"

    ! Public members
    public :: interpolant
    public :: setup_interpolant_scheme
    public :: extrapolate_cell_averages_to_faces
    public :: destroy_interpolant_scheme
    public :: compute_face_interpolant
    public :: x_sound_speed_left
    public :: x_sound_speed_right
    public :: y_sound_speed_left
    public :: y_sound_speed_right
    public :: z_sound_speed_left
    public :: z_sound_speed_right
    public :: x_qp_left, x_qp_right
    public :: y_qp_left, y_qp_right
    public :: z_qp_left, z_qp_right
    public :: x_density_left, x_density_right
    public :: x_x_speed_left, x_x_speed_right
    public :: x_y_speed_left, x_y_speed_right
    public :: x_z_speed_left, x_z_speed_right
    public :: x_pressure_left, x_pressure_right
    public :: y_density_left, y_density_right
    public :: y_x_speed_left, y_x_speed_right
    public :: y_y_speed_left, y_y_speed_right
    public :: y_z_speed_left, y_z_speed_right
    public :: y_pressure_left, y_pressure_right
    public :: z_density_left, z_density_right
    public :: z_x_speed_left, z_x_speed_right
    public :: z_y_speed_left, z_y_speed_right
    public :: z_z_speed_left, z_z_speed_right
    public :: z_pressure_left, z_pressure_right

    contains

        subroutine allocate_memory()
            implicit none
            call alloc(x_qp_left, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_left.')
            call alloc(x_qp_right, 0, imx+1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'x_qp_right.')
            call alloc(y_qp_left, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_left.')
            call alloc(y_qp_right, 1, imx-1, 0, jmx+1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'y_qp_right.')
            call alloc(z_qp_left, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_left.')
            call alloc(z_qp_right, 1, imx-1, 1, jmx-1, 0, kmx+1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for ' // &
                        'z_qp_right.')
        end subroutine allocate_memory

        subroutine link_aliases()
            implicit none

            ! Link xi faces left pointers
            x_density_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 1)
            x_x_speed_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 2)
            x_y_speed_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 3)
            x_z_speed_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 4)
            x_pressure_left(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_left(:, :, :, 5)

            ! Link xi faces right pointers
            x_density_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 1)
            x_x_speed_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 2)
            x_y_speed_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 3)
            x_z_speed_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 4)
            x_pressure_right(0:imx+1, 1:jmx-1, 1:kmx-1) => x_qp_right(:, :, :, 5)

            ! Link eta faces left pointers
            y_density_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 1)
            y_x_speed_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 2)
            y_y_speed_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 3)
            y_z_speed_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 4)
            y_pressure_left(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_left(:, :, :, 5)

            ! Link eta faces right pointers
            y_density_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 1)
            y_x_speed_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 2)
            y_y_speed_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 3)
            y_z_speed_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 4)
            y_pressure_right(1:imx-1, 0:jmx+1, 1:kmx-1) => y_qp_right(:, :, :, 5)
            
            ! Link zeta faces left pointers
            z_density_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 1)
            z_x_speed_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 2)
            z_y_speed_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 3)
            z_z_speed_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 4)
            z_pressure_left(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_left(:, :, :, 5)

            ! Link zeta faces right pointers
            z_density_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 1)
            z_x_speed_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 2)
            z_y_speed_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 3)
            z_z_speed_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 4)
            z_pressure_right(1:imx-1, 1:jmx-1, 0:kmx+1) => z_qp_right(:, :, :, 5)

            !turbulent variable tk and tw linking _qp_ (6:7)
            include "turbulence_models/include/face_interpolant/link_aliases.inc"
        end subroutine link_aliases

        subroutine setup_interpolant_scheme()
            implicit none
            select case (interpolant)
                case ("none")
                    ! Do nothing
                    continue
                case ("ppm")
                    call setup_scheme_ppm()
                case ("muscl")
                    call setup_scheme_muscl()
                case ("weno")
                    call setup_scheme_weno()
                case ("weno_NM")
                    call setup_scheme_weno_NM()
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
            call dealloc(z_qp_left)
            call dealloc(z_qp_right)
        end subroutine deallocate_memory

        subroutine unlink_aliases()
            implicit none

            ! Unlink xi faces left pointers
            nullify(x_density_left)
            nullify(x_x_speed_left)
            nullify(x_y_speed_left)
            nullify(x_z_speed_left)
            nullify(x_pressure_left)

            ! Unlink xi faces right pointers
            nullify(x_density_right)
            nullify(x_x_speed_right)
            nullify(x_y_speed_right)
            nullify(x_z_speed_right)
            nullify(x_pressure_right)

            ! Unlink eta faces left pointers
            nullify(y_density_left)
            nullify(y_x_speed_left)
            nullify(y_y_speed_left)
            nullify(y_z_speed_left)
            nullify(y_pressure_left)

            ! Unlink eta faces right pointers
            nullify(y_density_right)
            nullify(y_x_speed_right)
            nullify(y_y_speed_right)
            nullify(y_z_speed_right)
            nullify(y_pressure_right)
            
            ! Unlink tau faces left pointers
            nullify(z_density_left)
            nullify(z_x_speed_left)
            nullify(z_y_speed_left)
            nullify(z_z_speed_left)
            nullify(z_pressure_left)

            ! Unlink tau faces right pointers
            nullify(z_density_right)
            nullify(z_x_speed_right)
            nullify(z_y_speed_right)
            nullify(z_z_speed_right)
            nullify(z_pressure_right)

            include "turbulence_models/include/face_interpolant/unlink_aliases.inc" 
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
                case ("muscl")
                    call destroy_scheme_muscl()
                case ("weno")
                    call destroy_scheme_weno()
                case ("weno_NM")
                    call destroy_scheme_weno_NM()
                case default
                    call dmsg(5, 'state_interpolant', &
                            'destroy_interpolant_scheme', &
                            'Interpolant not recognized.')
                    stop
            end select
        end subroutine destroy_interpolant_scheme

        subroutine extrapolate_cell_averages_to_faces()
            implicit none

            call dmsg(1, 'face_interpolant', 'extrapolate_cell_averages_to_faces')

            x_qp_left(:, :, :, :) = qp(-1:imx, 1:jmx-1, 1:kmx-1, 1:n_var)
            x_qp_right(:, :, :, :) = qp(0:imx+1, 1:jmx-1, 1:kmx-1, 1:n_var)
            y_qp_left(:, :, :, :) = qp(1:imx-1, -1:jmx, 1:kmx-1, 1:n_var)
            y_qp_right(:, :, :, :) = qp(1:imx-1, 0:jmx+1, 1:kmx-1, 1:n_var)
            z_qp_left(:, :, :, :) = qp(1:imx-1, 1:jmx-1, -1:kmx, 1:n_var)
            z_qp_right(:, :, :, :) = qp(1:imx-1, 1:jmx-1, 0:kmx+1, 1:n_var)
        end subroutine extrapolate_cell_averages_to_faces

        subroutine compute_face_interpolant()
            implicit none
            select case (interpolant)
                case ("none")
                    call extrapolate_cell_averages_to_faces()
                case ("ppm")
                    call compute_ppm_states()
                    x_qp_left(:, :, :, :) = x_qp_left_ppm(:, :, :, :)
                    x_qp_right(:, :, :, :) = x_qp_right_ppm(:, :, :, :)
                    y_qp_left(:, :, :, :) = y_qp_left_ppm(:, :, :, :)
                    y_qp_right(:, :, :, :) = y_qp_right_ppm(:, :, :, :)
                    z_qp_left(:, :, :, :) = z_qp_left_ppm(:, :, :, :)
                    z_qp_right(:, :, :, :) = z_qp_right_ppm(:, :, :, :)
                case ("muscl")
                    call compute_muscl_states()
                    x_qp_left = x_qp_left_muscl
                    x_qp_right = x_qp_right_muscl
                    y_qp_left = y_qp_left_muscl
                    y_qp_right = y_qp_right_muscl
                    z_qp_left = z_qp_left_muscl
                    z_qp_right = z_qp_right_muscl
                    x_qp_left( :, :, :, 6:) = qp(-1:imx, 1:jmx-1, 1:kmx-1, 6:n_var)
                    x_qp_right(:, :, :, 6:) = qp(0:imx+1, 1:jmx-1, 1:kmx-1, 6:n_var)
                    y_qp_left( :, :, :, 6:) = qp(1:imx-1, -1:jmx, 1:kmx-1, 6:n_var)
                    y_qp_right(:, :, :, 6:) = qp(1:imx-1, 0:jmx+1, 1:kmx-1, 6:n_var)
                    z_qp_left( :, :, :, 6:) = qp(1:imx-1, 1:jmx-1, -1:kmx, 6:n_var)
                    z_qp_right(:, :, :, 6:) = qp(1:imx-1, 1:jmx-1, 0:kmx+1, 6:n_var)
                case ("weno")
                    call compute_weno_states()
                    x_qp_left  =  x_qp_left_weno
                    x_qp_right = x_qp_right_weno
                    y_qp_left  =  y_qp_left_weno
                    y_qp_right = y_qp_right_weno
                    z_qp_left  =  z_qp_left_weno
                    z_qp_right = z_qp_right_weno
                    x_qp_left( :, :, :, 6:) = qp(-1:imx, 1:jmx-1, 1:kmx-1, 6:n_var)
                    x_qp_right(:, :, :, 6:) = qp(0:imx+1, 1:jmx-1, 1:kmx-1, 6:n_var)
                    y_qp_left( :, :, :, 6:) = qp(1:imx-1, -1:jmx, 1:kmx-1, 6:n_var)
                    y_qp_right(:, :, :, 6:) = qp(1:imx-1, 0:jmx+1, 1:kmx-1, 6:n_var)
                    z_qp_left( :, :, :, 6:) = qp(1:imx-1, 1:jmx-1, -1:kmx, 6:n_var)
                    z_qp_right(:, :, :, 6:) = qp(1:imx-1, 1:jmx-1, 0:kmx+1, 6:n_var)
                case ("weno_NM")
                    call compute_weno_NM_states()
                    x_qp_left  =  x_qp_left_weno_NM
                    x_qp_right = x_qp_right_weno_NM
                    y_qp_left  =  y_qp_left_weno_NM
                    y_qp_right = y_qp_right_weno_NM
                    z_qp_left  =  z_qp_left_weno_NM
                    z_qp_right = z_qp_right_weno_NM
                case default
                    call dmsg(5, 'state_interpolant', &
                            'compute_face_interpolant', &
                            'Interpolant not recognized.')
                    stop
            end select
        end subroutine compute_face_interpolant

        function x_sound_speed_left()
            implicit none
            real, dimension(1:imx, 1:jmx-1, 1:kmx-1) :: x_sound_speed_left
            x_sound_speed_left = sqrt(gm * x_pressure_left / x_density_left)
        end function x_sound_speed_left

        function x_sound_speed_right()
            implicit none
            real, dimension(1:imx, 1:jmx-1, 1:kmx-1) :: x_sound_speed_right
            x_sound_speed_right = sqrt(gm * x_pressure_right / x_density_right)
        end function x_sound_speed_right

        function y_sound_speed_left()
            implicit none
            real, dimension(1:imx-1, 1:jmx, 1:kmx-1) :: y_sound_speed_left
            y_sound_speed_left = sqrt(gm * y_pressure_left / y_density_left)
        end function y_sound_speed_left

        function y_sound_speed_right()
            implicit none
            real, dimension(1:imx-1, 1:jmx, 1:kmx-1) :: y_sound_speed_right
            y_sound_speed_right = sqrt(gm * y_pressure_right / y_density_right)
        end function y_sound_speed_right

        function z_sound_speed_left()
            implicit none
            real, dimension(1:imx-1, 1:jmx-1, 1:kmx) :: z_sound_speed_left
            z_sound_speed_left = sqrt(gm * z_pressure_left / z_density_left)
        end function z_sound_speed_left

        function z_sound_speed_right()
            implicit none
            real, dimension(1:imx-1, 1:jmx-1, 1:kmx) :: z_sound_speed_right
            z_sound_speed_right = sqrt(gm * z_pressure_right / z_density_right)
        end function z_sound_speed_right

end module face_interpolant
