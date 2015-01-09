module state
    !-------------------------------------------------------------------
    ! The state module contains the state variables and the methods that
    ! act on them. 
    !
    ! The state of the system is defined using the density, velocity and
    ! pressure (primitive variables qp) at the grid points. This current
    ! version assumes the grid is in atmost two dimensions. 
    !-------------------------------------------------------------------
    
    use global, only: FILE_NAME_LENGTH, STATE_FILE_UNIT, OUT_FILE_UNIT
    use utils, only: alloc, dealloc, dmsg
    use string
    use grid, only: imx, jmx
    use geometry, only: xnx, xny, ynx, yny

    implicit none
    private

    ! State variables
    real, public, dimension(:, :, :), allocatable, target :: qp
    ! Infinity variables (free stream conditions)
    real, public, dimension(4), target :: qp_inf
    ! State variable component aliases
    real, public, dimension(:, :), pointer :: density
    real, public, dimension(:, :), pointer :: x_speed
    real, public, dimension(:, :), pointer :: y_speed
    real, public, dimension(:, :), pointer :: pressure
    real, public, pointer :: density_inf
    real, public, pointer :: x_speed_inf
    real, public, pointer :: y_speed_inf
    real, public, pointer :: pressure_inf
    ! Supersonic flag
    logical :: supersonic_flag
    ! Speed of sound at xi faces
    real, public, dimension(:, :), allocatable :: x_a
    ! Speed of sound at eta faces
    real, public, dimension(:, :), allocatable :: y_a
    ! Ratio of specific heats (gamma)
    real, public :: gm
    ! Specific gas constant
    real, public :: R_gas

    ! Public methods
    public :: setup_state
    public :: destroy_state
    public :: set_ghost_cell_data
    public :: compute_sound_speeds
    public :: xi_face_normal_speeds
    public :: eta_face_normal_speeds
    public :: writestate

    contains

        subroutine link_aliases()
            implicit none
            call dmsg(1, 'state', 'link_aliases')
            density(0:imx, 0:jmx) => qp(:, :, 1)
            x_speed(0:imx, 0:jmx) => qp(:, :, 2)
            y_speed(0:imx, 0:jmx) => qp(:, :, 3)
            pressure(0:imx, 0:jmx) => qp(:, :, 4)
            density_inf => qp_inf(1)
            x_speed_inf => qp_inf(2)
            y_speed_inf => qp_inf(3)
            pressure_inf => qp_inf(4)
        end subroutine link_aliases

        subroutine unlink_aliases()
            implicit none
            call dmsg(1, 'state', 'unlink_aliases')
            nullify(density)
            nullify(x_speed)
            nullify(y_speed)
            nullify(pressure)
            nullify(density_inf)
            nullify(x_speed_inf)
            nullify(y_speed_inf)
            nullify(pressure_inf)
        end subroutine unlink_aliases

        subroutine allocate_memory()
            !-----------------------------------------------------------
            ! Allocate memory for the state variables
            !
            ! This assumes that imx and jmx (the grid size) has been set
            ! within the state module.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'state', 'allocate_memory')

            ! The state of the system is defined by the primitive 
            ! variables (density, velocity and pressure) at the grid
            ! cell centers. 
            ! There are (imx - 1) x (jmx - 1) grid cells within the 
            ! domain. We require a row of ghost cells on each boundary.
            ! This current implementation is for a 2D/1D case. 
            call alloc(qp, 0, imx, 0, jmx, 1, 4, &
                    errmsg='Error: Unable to allocate memory for state ' // &
                        'variable qp.')
            call alloc(x_a, 1, imx, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for x_a.')
            call alloc(y_a, 1, imx-1, 1, jmx, &
                    errmsg='Error: Unable to allocate memory for y_a.')
        end subroutine allocate_memory

        subroutine deallocate_memory()

            implicit none

            call dmsg(1, 'state', 'deallocate_memory')

            call dealloc(qp)
            call dealloc(x_a)
            call dealloc(y_a)

        end subroutine deallocate_memory

        subroutine setup_state(free_stream_density, free_stream_x_speed, &
                free_stream_y_speed, free_stream_pressure, state_file)
            !-----------------------------------------------------------
            ! Setup the state module.
            !
            ! This subroutine should be run before the state variables
            ! are initilized. This subroutine allocates the memory for 
            ! state variables and sets up the aliases to refer to the 
            ! components of the state.
            !-----------------------------------------------------------

            implicit none
            real, intent(in) :: free_stream_density
            real, intent(in) :: free_stream_x_speed, free_stream_y_speed
            real, intent(in) :: free_stream_pressure
            character(len=FILE_NAME_LENGTH), intent(in) :: state_file

            call dmsg(1, 'state', 'setup_state')

            call allocate_memory()
            call link_aliases()
            call init_infinity_values(free_stream_density, &
                    free_stream_x_speed, free_stream_y_speed, &
                    free_stream_pressure)
            call set_supersonic_flag()
            call initstate(state_file)

        end subroutine setup_state

        subroutine destroy_state()
            !-----------------------------------------------------------
            ! Destroy the state module.
            !
            ! This subroutine destroys the state module which includes
            ! unlinking the aliases for the state components and 
            ! deallocating the memory held by the state variables. 
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'state', 'destroy_state')

            call unlink_aliases()
            call deallocate_memory()

        end subroutine destroy_state

        subroutine set_inlet_and_exit_state_variables()
            !-----------------------------------------------------------
            ! Set / extrapolate inlet and exit state variables
            !
            ! The inlet and exit variables are set based on the 
            ! prescribed infinity (free stream) values and extrapolation
            ! of the neighbouring cells. The decision to impose or 
            ! extrapolate is based on the wave speeds and whether or not
            ! the flow is supersonic.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'state', 'set_inlet_and_exit_state_variables')

            ! Impose the density, x_speed and y_speed at the inlet
            density(0, :) = density_inf
            x_speed(0, :) = x_speed_inf
            y_speed(0, :) = y_speed_inf
            ! Extrapolate these quantities at the exit
            density(imx, :) = density(imx - 1, :)
            x_speed(imx, :) = x_speed(imx - 1, :)
            y_speed(imx, :) = y_speed(imx - 1, :)
            ! If the flow is subsonic, impose the back pressure
            ! Else, impose the inlet pressure
            ! Extrapolate at the other end
            if (supersonic_flag .eqv. .TRUE.) then
                pressure(0, :) = pressure_inf
                pressure(imx, :) = pressure(imx - 1, :)
            else
                pressure(imx, :) = pressure_inf
                pressure(0, :) = pressure(1, :)
            end if

        end subroutine set_inlet_and_exit_state_variables

        subroutine set_top_and_bottom_ghost_cell_data()
            !-----------------------------------------------------------
            ! Set the state variables for the top and bottom ghosh cells
            !
            ! The pressure and density for the top and bottom ghost 
            ! cells is extrapolated from the neighboring interior cells.
            ! The velocity components are computed by applying the 
            ! flow tangency conditions.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'state', 'set_top_and_bottom_ghost_cell_data')

            pressure(1:imx-1, 0) = pressure(1:imx-1, 1)
            pressure(1:imx-1, jmx) = pressure(1:imx-1, jmx-1)
            density(1:imx-1, 0) = density(1:imx-1, 1)
            density(1:imx-1, jmx) = density(1:imx-1, jmx-1)
            call apply_flow_tangency_conditions()

        end subroutine set_top_and_bottom_ghost_cell_data

        subroutine apply_flow_tangency_conditions()
            !-----------------------------------------------------------
            ! Apply the flow tangency conditions
            !
            ! The flow tangency conditions ensure that there is no flow
            ! across the boundaries. This is done by ensuring that the
            ! flow is parallel to the boundary.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'state', 'apply_flow_tangency_conditions')

            ! For the top cells
            x_speed(1:imx-1, jmx) = x_speed(1:imx-1, jmx-1) - &
                    (2. * &
                        ((x_speed(1:imx-1, jmx-1) * &
                            ynx(1:imx-1, jmx) ** 2.) &
                        + (y_speed(1:imx-1, jmx-1) * &
                            yny(1:imx-1, jmx) * ynx(1:imx-1, jmx)) &
                        ) &
                    )
            y_speed(1:imx-1, jmx) = y_speed(1:imx-1, jmx-1) - &
                    (2. * &
                        ((x_speed(1:imx-1, jmx-1) * &
                            ynx(1:imx-1, jmx) * yny(1:imx-1, jmx)) &
                        + (y_speed(1:imx-1, jmx-1) * &
                            yny(1:imx-1, jmx) ** 2.) &
                        ) &
                    )
            ! For the bottom cells
            x_speed(1:imx-1, 0) = x_speed(1:imx-1, 1) - &
                    (2. * &
                        ((x_speed(1:imx-1, 1) * &
                            ynx(1:imx-1, 1) ** 2.) &
                        + (y_speed(1:imx-1, 1) * &
                            yny(1:imx-1, 1) * ynx(1:imx-1, 1)) &
                        ) &
                    )
            y_speed(1:imx-1, 0) = y_speed(1:imx-1, 1) - &
                    (2. * &
                        ((x_speed(1:imx-1, 1) * &
                            ynx(1:imx-1, 1) * yny(1:imx-1, 1)) &
                        + (y_speed(1:imx-1, 1) * &
                            yny(1:imx-1, 1) ** 2.) &
                        ) &
                    )

        end subroutine apply_flow_tangency_conditions

        subroutine set_ghost_cell_data()
            !-----------------------------------------------------------
            ! Set the data in the ghost cells
            !
            ! The ghost cell data is either imposed or extrapolated from
            ! the neighboring cells inside the domain. The decision to
            ! impose or extrapolate is taken based on certain 
            ! conditions.
            !-----------------------------------------------------------

            implicit none

            call dmsg(1, 'state', 'set_ghost_cell_data')

            call set_inlet_and_exit_state_variables()
            call set_top_and_bottom_ghost_cell_data()

        end subroutine set_ghost_cell_data

        subroutine init_infinity_values(free_stream_density, &
                free_stream_x_speed, free_stream_y_speed, free_stream_pressure)
            !-----------------------------------------------------------
            ! Set the values of the infinity variables
            !-----------------------------------------------------------

            implicit none
            real, intent(in) :: free_stream_density
            real, intent(in) :: free_stream_x_speed
            real, intent(in) :: free_stream_y_speed
            real, intent(in) :: free_stream_pressure
            
            call dmsg(1, 'state', 'init_infinity_values')

            density_inf = free_stream_density
            x_speed_inf = free_stream_x_speed
            y_speed_inf = free_stream_y_speed
            pressure_inf = free_stream_pressure

        end subroutine init_infinity_values
        
        function sound_speed_inf() result(a)
            !-----------------------------------------------------------
            ! Return the free stream speed of sound.
            !-----------------------------------------------------------

            implicit none
            real :: a

            a = sqrt(gm * pressure_inf / density_inf)

        end function sound_speed_inf

        subroutine compute_sound_speeds()
            !-----------------------------------------------------------
            ! Compute the speed of sound at each cell face
            !-----------------------------------------------------------

            implicit none
            real, dimension(0:imx, 0:jmx) :: sound_speed
            
            call dmsg(1, 'state', 'compute_sound_speeds')

            ! Compute the sound speed at each cell center
            sound_speed(:, :) = sqrt(gm * pressure(:, :) / density(:, :))
            if (any(sound_speed /= sound_speed)) then
                call dmsg(5, 'state', 'compute_sound_speeds', &
                        msg='Error: NaN detected in sound_speed.')
            end if
            call compute_sound_speed_at_xi_faces(sound_speed)
            call compute_sound_speed_at_eta_faces(sound_speed)

        end subroutine compute_sound_speeds

        subroutine compute_sound_speed_at_xi_faces(sound_speed)
            implicit none
            real, dimension(0:imx, 0:jmx), intent(in) :: sound_speed
            
            call dmsg(1, 'state', 'compute_sound_speed_at_xi_faces')

            x_a(:, :) = 0.5 * &
                    (sound_speed(0:imx-1, 1:jmx-1) + &
                        sound_speed(1:imx, 1:jmx-1))
            if (any(isnan(x_a)) .or. any(x_a < 0)) stop
        end subroutine compute_sound_speed_at_xi_faces

        subroutine compute_sound_speed_at_eta_faces(sound_speed)
            implicit none
            real, dimension(0:imx, 0:jmx), intent(in) :: sound_speed
            
            call dmsg(1, 'state', 'compute_sound_speed_at_eta_faces')

            y_a(:, :) = 0.5 * &
                    (sound_speed(1:imx-1, 0:jmx-1) + &
                        sound_speed(1:imx-1, 1:jmx))
            if (any(isnan(y_a)) .or. any(y_a < 0)) stop
        end subroutine compute_sound_speed_at_eta_faces

        subroutine set_supersonic_flag()
            !-----------------------------------------------------------
            ! Set the supersonic flag based on the infinity conditions.
            !
            ! In the ghost cells, the values of the primitive variables
            ! are set either based on the neighbouring cells' values or 
            ! the infinity values. This splitting is based on the wave
            ! speeds. The supersonic flag is used as an indication as 
            ! to how these cells should get their values. 
            !-----------------------------------------------------------

            implicit none
            real :: avg_inlet_mach

            call dmsg(1, 'state', 'set_supersonic_flag')

            avg_inlet_mach = sqrt(x_speed_inf ** 2. + y_speed_inf ** 2.) / &
                    sound_speed_inf()

            if (avg_inlet_mach >= 1) then
                supersonic_flag = .TRUE.
            else
                supersonic_flag = .FALSE.
            end if

            call dmsg(5, 'state', 'set_supersonic_flag', &
                    'Supersonic flag set to ' + supersonic_flag)

        end subroutine set_supersonic_flag

        function xi_face_normal_speeds(dir) result(face_normal_speeds)
            !-----------------------------------------------------------
            ! Return the xi face normal speed
            !
            ! The face normal speed is the component of the flow speed 
            ! at a face perpendicular to it. This function returns the 
            ! face normal speeds for xi faces.
            !
            ! Since the state variables are specified at cell 
            ! centers, the direction of extrapolation must be specified 
            ! using the following parameter:
            ! 
            !   dir: character
            !       Extrapolate from previous or following cell
            !       Use '+' to extrapolate from previous cell (in the 
            !       positive direction) and '-' to extrapolate from the
            !       following cell (in the negative direction).
            !-----------------------------------------------------------

            implicit none
            character, intent(in) :: dir
            real, dimension(imx, jmx-1) :: face_normal_speeds
            integer :: i
            integer, dimension(imx) :: indices
            
            call dmsg(1, 'state', 'xi_face_normal_speeds')

            ! Find the indices of the cell centers whose speeds will be
            ! used
            indices = (/ (i, i = 0, imx-1) /)
            if (dir .eq. '-') then
                ! We have to compute the face normal speeds in the 
                ! negative direction, i.e., using the following cells.
                ! Adding 1 to the above indices gives us the correct 
                ! ones.
                indices = indices + 1
            ! else if the direction is plus, the indices are correct
            else if (dir .ne. '+') then
                ! An incorrect direction is specified
                print *, 'Error: An incorrect direction was specified to ', &
                        'xi_face_normal_speeds().'
                stop
            end if

            face_normal_speeds = x_speed(indices, 1:jmx-1) * xnx + &
                    y_speed(indices, 1:jmx-1) * xny

        end function xi_face_normal_speeds

        function eta_face_normal_speeds(dir) result(face_normal_speeds)
            !-----------------------------------------------------------
            ! Return the eta face normal speed
            !
            ! The face normal speed is the component of the flow speed 
            ! at a face perpendicular to it. This function returns the 
            ! face normal speeds for eta faces.
            !
            ! Since the state variables are specified at cell 
            ! centers, the direction of extrapolation must be specified 
            ! using the following parameter:
            ! 
            !   dir: character
            !       Extrapolate from previous or following cell
            !       Use '+' to extrapolate from previous cell (in the 
            !       positive direction) and '-' to extrapolate from the
            !       following cell (in the negative direction).
            !-----------------------------------------------------------

            implicit none
            character, intent(in) :: dir
            real, dimension(imx-1, jmx) :: face_normal_speeds
            integer :: j
            ! Find the indices of the cell centers whose speeds will be
            ! used
            integer, dimension(jmx) :: indices
            
            call dmsg(1, 'state', 'eta_face_normal_speeds')

            indices = (/ (j, j = 0, jmx-1) /)
            if (dir .eq. '-') then
                ! We have to compute the face normal speeds in the 
                ! negative direction, i.e., using the following cells.
                ! Adding 1 to the above indices gives us the correct 
                ! ones.
                indices = indices + 1
            ! else if the direction is plus, the indices are correct
            else if (dir .ne. '+') then
                ! An incorrect direction is specified
                print *, 'Error: An incorrect direction was specified to ', &
                        'eta_face_normal_speeds().'
                stop
            end if

            face_normal_speeds = x_speed(1:imx-1, indices) * ynx + &
                    y_speed(1:imx-1, indices) * yny

        end function eta_face_normal_speeds

        subroutine initstate(state_file)
            !-----------------------------------------------------------
            ! Initialize the state
            !
            ! If state_file is a tilde (~), then the state should be 
            ! set to the infinity values. Otherwise, read the state_file
            ! to get the state values.
            !-----------------------------------------------------------

            implicit none
            character(len=FILE_NAME_LENGTH), intent(in) :: state_file
            
            call dmsg(1, 'state', 'initstate')

            if (state_file .eq. '~') then
                ! Set the state to the infinity values
                call init_state_with_infinity_values()
            else
                call readstate(state_file)
            end if

        end subroutine initstate

        subroutine init_state_with_infinity_values()
            !-----------------------------------------------------------
            ! Initialize the state based on the infinity values.
            !-----------------------------------------------------------
            
            implicit none
            
            call dmsg(1, 'state', 'init_state_with_infinity_values')

            qp(:, :, 1) = qp_inf(1)
            qp(:, :, 2) = qp_inf(2)
            qp(:, :, 3) = qp_inf(3)
            qp(:, :, 4) = qp_inf(4)

        end subroutine init_state_with_infinity_values

        subroutine writestate(outfile)
            !-----------------------------------------------------------
            ! Write the state of the system to a file
            !-----------------------------------------------------------

            implicit none
            character(len=FILE_NAME_LENGTH), intent(in) :: outfile
            integer :: i, j
            
            call dmsg(1, 'state', 'writestate')

            open(OUT_FILE_UNIT, file=outfile)

            write(OUT_FILE_UNIT, *) 'CELLDATA'
            write(OUT_FILE_UNIT, *) 'Density'
            do j = 1, jmx - 1
                do i = 1, imx - 1
                    write(OUT_FILE_UNIT, *) density(i, j)
                end do
            end do

            write(OUT_FILE_UNIT, *) 'CELLDATA'
            write(OUT_FILE_UNIT, *) 'Velocity'
            do j = 1, jmx - 1
                do i = 1, imx - 1
                    write(OUT_FILE_UNIT, *) x_speed(i, j), y_speed(i, j)
                end do
            end do

            write(OUT_FILE_UNIT, *) 'CELLDATA'
            write(OUT_FILE_UNIT, *) 'Pressure'
            do j = 1, jmx - 1
                do i = 1, imx - 1
                    write(OUT_FILE_UNIT, *) pressure(i, j)
                end do
            end do
            
            close(OUT_FILE_UNIT)

        end subroutine writestate

        subroutine readstate(state_file)
            !-----------------------------------------------------------
            ! Initialize the state using a state file.
            !
            ! Prior to running this subroutine, the memory for the 
            ! state variables should have been allocated and the 
            ! pointers to alias the components of the state should have 
            ! been associated.
            !-----------------------------------------------------------
            
            implicit none
            character(len=FILE_NAME_LENGTH), intent(in) :: state_file
            integer :: i, j
            
            call dmsg(1, 'state', 'readstate')

            open(STATE_FILE_UNIT, file=state_file)

            ! Skip the section header
            read(STATE_FILE_UNIT, *)
            read(STATE_FILE_UNIT, *)
            do j = 1, jmx - 1
                do i = 1, imx - 1
                    read(STATE_FILE_UNIT, *) density(i, j)
                end do
            end do

            ! Skip the section header
            read(STATE_FILE_UNIT, *)
            read(STATE_FILE_UNIT, *)
            do j = 1, jmx - 1
                do i = 1, imx - 1
                    read(STATE_FILE_UNIT, *) x_speed(i, j), y_speed(i, j)
                end do
            end do

            ! Skip the section header
            read(STATE_FILE_UNIT, *)
            read(STATE_FILE_UNIT, *)
            do j = 1, jmx - 1
                do i = 1, imx - 1
                    read(STATE_FILE_UNIT, *) pressure(i, j)
                end do
            end do
            
            close(STATE_FILE_UNIT)

        end subroutine readstate

end module state
