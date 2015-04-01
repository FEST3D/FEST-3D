module solver

    use global, only: CONFIG_FILE_UNIT, RESNORM_FILE_UNIT, FILE_NAME_LENGTH, &
            STRING_BUFFER_LENGTH
    use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
    use string
    use grid, only: imx, jmx, setup_grid, destroy_grid
    use geometry, only: xnx, xny, ynx, yny, xA, yA, volume, setup_geometry, &
            destroy_geometry
    use state, only: qp, qp_inf, density, x_speed, y_speed, pressure, &
            density_inf, x_speed_inf, y_speed_inf, pressure_inf, gm, R_gas, &
            setup_state, destroy_state, set_ghost_cell_data
    use face_interpolant, only: interpolant, &
            x_sound_speed_left, x_sound_speed_right, &
            y_sound_speed_left, y_sound_speed_right
    use scheme, only: scheme_name, residue, setup_scheme, destroy_scheme, &
            compute_residue
    
    implicit none
    private

    real, public :: CFL
    character, public :: time_stepping_method
    real :: tolerance
    integer, public :: max_iters
    real, public :: resnorm, resnorm_0
    real, public, dimension(:, :), allocatable :: delta_t
    integer, public :: iter

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: step
    public :: converged

    contains

        subroutine get_next_token(buf)
            !-----------------------------------------------------------
            ! Extract the next token from the config file
            !
            ! Each token is on a separate line.
            ! There may be multiple comments (lines beginning with #) 
            ! and blank lines in between.
            ! The purpose of this subroutine is to ignore all these 
            ! lines and return the next "useful" line.
            !-----------------------------------------------------------

            implicit none
            character(len=STRING_BUFFER_LENGTH), intent(out) :: buf
            integer :: ios

            do
                read(CONFIG_FILE_UNIT, *, iostat=ios) buf
                if (ios /= 0) then
                    print *, 'Error while reading config file.'
                    print *, 'Current buffer length is set to: ', &
                            STRING_BUFFER_LENGTH
                    stop
                end if
                if (index(buf, '#') == 1) then
                    ! The current line begins with a hash
                    ! Ignore it
                    continue
                else if (len_trim(buf) == 0) then
                    ! The current line is empty
                    ! Ignore it
                    continue
                else
                    ! A new token has been found
                    ! Break out
                    exit
                end if
            end do
            call dmsg(0, 'solver', 'get_next_token', 'Returning: ' // trim(buf))

        end subroutine get_next_token

        subroutine read_config_file(free_stream_density, &
                free_stream_x_speed, free_stream_y_speed, &
                free_stream_pressure, grid_file, state_load_file)

            implicit none
            real, intent(out) :: free_stream_density, free_stream_x_speed, &
                    free_stream_y_speed, free_stream_pressure
            character(len=FILE_NAME_LENGTH), intent(out) :: grid_file
            character(len=FILE_NAME_LENGTH), intent(out) :: state_load_file
            character(len=FILE_NAME_LENGTH) :: config_file = "config.md"
            integer, parameter :: CONFIG_FILE_UNIT = 1
            character(len=STRING_BUFFER_LENGTH) :: buf

            call dmsg(1, 'solver', 'read_config_file')
            
            open(CONFIG_FILE_UNIT, file=config_file)

            ! Ignore the config file header
            read(CONFIG_FILE_UNIT, *)
            read(CONFIG_FILE_UNIT, *)
            
            ! Read the parameters from the file

            call get_next_token(buf)
            read(buf, *) scheme_name
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='scheme_name = ' + scheme_name)

            call get_next_token(buf)
            read(buf, *) interpolant
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='interpolant = ' + interpolant)

            call get_next_token(buf)
            read(buf, *) CFL
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='CFL = ' + CFL)

            call get_next_token(buf)
            read(buf, *) time_stepping_method
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='time_stepping_method = ' + time_stepping_method)

            call get_next_token(buf)
            read(buf, *) tolerance
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='tolerance  = ' + tolerance)

            call get_next_token(buf)
            read(buf, *) grid_file
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='grid_file = ' + grid_file)

            call get_next_token(buf)
            read(buf, *) state_load_file
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='state_load_file = ' + state_load_file)

            call get_next_token(buf)
            read(buf, *) max_iters
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='max_iters = ' + max_iters)

            call get_next_token(buf)
            read(buf, *) DEBUG_LEVEL
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='DEBUG_LEVEL = ' + DEBUG_LEVEL)

            call get_next_token(buf)
            read(buf, *) gm
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='gamma = ' + gm)

            call get_next_token(buf)
            read(buf, *) R_gas
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='R_gas = ' + R_gas)

            call get_next_token(buf)
            read(buf, *) free_stream_density
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_density = ' + free_stream_density)

            call get_next_token(buf)
            read(buf, *) free_stream_x_speed
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_x_speed = ' + free_stream_x_speed)

            call get_next_token(buf)
            read(buf, *) free_stream_y_speed
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_y_speed = ' + free_stream_y_speed)

            call get_next_token(buf)
            read(buf, *) free_stream_pressure
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='free_stream_pressure = ' + free_stream_pressure)

            close(CONFIG_FILE_UNIT)

        end subroutine read_config_file

        subroutine setup_solver()
            
            implicit none
            real :: free_stream_density
            real :: free_stream_x_speed, free_stream_y_speed
            real :: free_stream_pressure
            character(len=FILE_NAME_LENGTH) :: grid_file
            character(len=FILE_NAME_LENGTH) :: state_load_file
            
            call dmsg(1, 'solver', 'setup_solver')

            call read_config_file(free_stream_density, free_stream_x_speed, &
                    free_stream_y_speed, free_stream_pressure, grid_file, &
                    state_load_file)
            call setup_grid(grid_file)
            call setup_geometry()
            call setup_state(free_stream_density, free_stream_x_speed, &
                    free_stream_y_speed, free_stream_pressure, state_load_file)
            call allocate_memory()
            call setup_scheme()
            call initmisc()
            open(RESNORM_FILE_UNIT, file='resnorms')

        end subroutine setup_solver

        subroutine destroy_solver()

            implicit none
            
            call dmsg(1, 'solver', 'destroy_solver')

            call destroy_scheme()
            call deallocate_misc()
            call destroy_state()
            call destroy_geometry()
            call destroy_grid()
            close(RESNORM_FILE_UNIT)

        end subroutine destroy_solver

        subroutine initmisc()
            
            implicit none
            
            call dmsg(1, 'solver', 'initmisc')

            iter = 0
            resnorm = 1.
            resnorm_0 = 1.

        end subroutine initmisc

        subroutine deallocate_misc()

            implicit none
            
            call dmsg(1, 'solver', 'deallocate_misc')

            call dealloc(residue)
            call dealloc(delta_t)

        end subroutine deallocate_misc

        subroutine allocate_memory()

            implicit none
            
            call dmsg(1, 'solver', 'allocate_memory')

            call alloc(delta_t, 1, imx-1, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')

        end subroutine allocate_memory

        subroutine compute_local_time_step()
            !-----------------------------------------------------------
            ! Compute the time step to be used at each cell center
            !
            ! Local time stepping can be used to get the solution 
            ! advance towards steady state faster. If only the steady
            ! state solution is required, i.e., transients are 
            ! irrelevant, use local time stepping. 
            !-----------------------------------------------------------

            implicit none
            real, dimension(imx-1, jmx-1) :: lmx1, lmx2, lmx3, lmx4, lmxsum
            real, dimension(imx, jmx-1) :: x_sound_speed_avg
            real, dimension(imx-1, jmx) :: y_sound_speed_avg
            
            call dmsg(1, 'solver', 'compute_local_time_step')

            x_sound_speed_avg = 0.5 * &
                    (x_sound_speed_left() + x_sound_speed_right())
            y_sound_speed_avg = 0.5 * &
                    (y_sound_speed_left() + y_sound_speed_right())

            ! For left face
            lmx1(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * xnx(1:imx-1, 1:jmx-1)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * xny(1:imx-1, 1:jmx-1))) + &
                    x_sound_speed_avg(1:imx-1, 1:jmx-1)
            ! For bottom face
            lmx2(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * ynx(1:imx-1, 1:jmx-1)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * yny(1:imx-1, 1:jmx-1))) + &
                    y_sound_speed_avg(1:imx-1, 1:jmx-1)
            ! For right face
            lmx3(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * xnx(2:imx, 1:jmx-1)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * xny(2:imx, 1:jmx-1))) + &
                    x_sound_speed_avg(2:imx, 1:jmx-1)
            ! For top face
            lmx4(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * ynx(1:imx-1, 2:jmx)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * yny(1:imx-1, 2:jmx))) + &
                    y_sound_speed_avg(1:imx-1, 2:jmx)
            
            lmxsum(:, :) = (xA(1:imx-1, 1:jmx-1) * lmx1) + &
                    (yA(1:imx-1, 1:jmx-1) * lmx2) + &
                    (xA(2:imx, 1:jmx-1) * lmx3) + &
                    (yA(1:imx-1, 2:jmx) * lmx4)
            
            delta_t = 2. / lmxsum
            delta_t = delta_t * volume * CFL

        end subroutine compute_local_time_step

        subroutine compute_global_time_step()
            !-----------------------------------------------------------
            ! Compute a common time step to be used at all cell centers
            !
            ! Global time stepping is generally used to get time 
            ! accurate solutions; transients can be studied by 
            ! employing this strategy.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'compute_global_time_step')

            call compute_local_time_step()
            ! The global time step is the minimum of all the local time
            ! steps.
            delta_t = minval(delta_t)

        end subroutine compute_global_time_step

        subroutine compute_time_step()
            !-----------------------------------------------------------
            ! Compute the time step to be used
            !
            ! This calls either compute_global_time_step() or 
            ! compute_local_time_step() based on what 
            ! time_stepping_method is set to.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'compute_time_step')

            if (time_stepping_method .eq. 'g') then
                call compute_global_time_step()
            else if (time_stepping_method .eq. 'l') then
                call compute_local_time_step()
            else
                call dmsg(5, 'solver', 'compute_time_step', &
                        msg='Value for time_stepping_method (' // &
                            time_stepping_method // ') not recognized.')
                stop
            end if

        end subroutine compute_time_step

        subroutine update_solution()
            !-----------------------------------------------------------
            ! Update the solution using the residue and time step
            !-----------------------------------------------------------

            implicit none
            !real, dimension(4, 4) :: psi_inv
            !integer :: i, j
            real, dimension(0:imx, 0:jmx, 4) :: qp_temp
            
            call dmsg(1, 'solver', 'update_solution')

!            do j = 1, jmx - 1
!                do i = 1, imx - 1
!                    print *, 'Location: ', i, j
!                    psi_inv = transpose(reshape((/ &
!                            1., 0., 0., 0., &
!                            - x_speed(i, j) / density(i, j), &
!                                1. / density(i, j), 0., 0., &
!                            - y_speed(i, j) / density(i, j), &
!                                0., 1. / density(i, j), 0., &
!                            0., 0., 0., gm - 1 &
!                            /), shape(psi_inv)))
!                    if (any(abs(residue(i, j, :)) > 1e-6)) then
!                        print *, 'old qp: ', qp(i, j, :)
!                        print *, 'psi_inv: ', psi_inv
!                        print *, 'residue: ', residue(i, j, :)
!                        print *, 'delta_t: ', delta_t(i, j)
!                        print *, 'volume: ', volume(i, j)
!                    end if
!                    qp(i, j, :) = qp(i, j, :) - &
!                            (matmul(psi_inv, residue(i, j, :)) * &
!                                delta_t(i, j) / volume(i, j))
!                    if (any(abs(residue(i, j, :)) > 1e-6)) then
!                        print *, 'new qp: ', qp(i, j, :)
!                    end if
!                    if (qp(i, j, 1) < 0 .or. qp(i, j, 4) < 0) stop
!                end do
!            end do
            qp_temp(1:imx-1, 1:jmx-1, 1) = qp(1:imx-1, 1:jmx-1, 1) - &
                    (residue(:, :, 1) * &
                    delta_t(:, :) / volume(:, :))
            qp_temp(1:imx-1, 1:jmx-1, 2) = qp(1:imx-1, 1:jmx-1, 2) - &
                    (( (-1. * qp(1:imx-1, 1:jmx-1, 2) / qp(1:imx-1, 1:jmx-1, 1) * residue(:, :, 1)) + &
                       ( residue(:, :, 2) / qp(1:imx-1, 1:jmx-1, 1) )) * &
                    delta_t(:, :) / volume(:, :))
            qp_temp(1:imx-1, 1:jmx-1, 3) = qp(1:imx-1, 1:jmx-1, 3) - &
                    (( (-1. * qp(1:imx-1, 1:jmx-1, 3) / qp(1:imx-1, 1:jmx-1, 1) * residue(:, :, 1)) + &
                       ( residue(:, :, 3) / qp(1:imx-1, 1:jmx-1, 1) )) * &
                    delta_t(:, :) / volume(:, :))
            qp_temp(1:imx-1, 1:jmx-1, 4) = qp(1:imx-1, 1:jmx-1, 4) - &
                    ( ( (0.5 * (1.4 - 1.) * ( qp(1:imx-1, 1:jmx-1, 2)**2. + qp(1:imx-1, 1:jmx-1, 3)**2. ) * residue(:, :, 1)) + &
                       (- (1.4 - 1.) * qp(1:imx-1, 1:jmx-1, 2) * residue(:, :, 2)) + &
                       (- (1.4 - 1.) * qp(1:imx-1, 1:jmx-1, 3) * residue(:, :, 3)) + &
                       ((1.4 - 1.) * residue(:, :, 4)) ) * &
                    delta_t(:, :) / volume(:, :) )

            qp(1:imx-1, 1:jmx-1, :) = qp_temp(1:imx-1, 1:jmx-1, :)
            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
                print *, 'new qp: ', qp
                print *, 'residue: ', residue
                print *, 'delta_t: ', delta_t
                print *, 'volume: ', volume
            end if

        end subroutine update_solution

        subroutine step()
            !-----------------------------------------------------------
            ! Perform one time step iteration
            !
            ! This subroutine performs one iteration by stepping through
            ! time once.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'step')

            iter = iter + 1
            call set_ghost_cell_data()
            call compute_residue()
            call dmsg(1, 'solver', 'step', 'Residue computed.')
            call compute_time_step()
            call update_solution()
            call compute_residue_norm()
            if (iter .eq. 1) then
                resnorm_0 = resnorm
            end if
            write(RESNORM_FILE_UNIT, *) resnorm

        end subroutine step

        subroutine compute_residue_norm()

            implicit none
            
            call dmsg(1, 'solver', 'compute_residue_norm')

            resnorm = sum(sqrt( &
                    (residue(:, :, 1) / &
                        (density_inf * x_speed_inf)) ** 2. + &
                    (residue(:, :, 2) / &
                        (density_inf * x_speed_inf ** 2.)) ** 2. + &
                    (residue(:, :, 3) / &
                        (density_inf * x_speed_inf ** 2.)) ** 2. + &
                    (residue(:, :, 4) / &
                        (density_inf ** 2. * x_speed_inf ** 3.)) ** 2. &
                    ))

        end subroutine compute_residue_norm

        function converged() result(c)
            !-----------------------------------------------------------
            ! Check if the solution seems to have converged
            !
            ! The solution is said to have converged if the change in 
            ! the residue norm is "negligible".
            !-----------------------------------------------------------

            implicit none
            logical :: c
            
            call dmsg(1, 'solver', 'converged')

            if (resnorm / resnorm_0 < tolerance) then
                c = .TRUE.
            end if
            c = .FALSE.

        end function converged

end module solver

