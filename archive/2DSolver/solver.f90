module solver

    use global, only: CONFIG_FILE_UNIT, RESNORM_FILE_UNIT, FILE_NAME_LENGTH, &
            STRING_BUFFER_LENGTH, INTERPOLANT_NAME_LENGTH
    use utils, only: alloc, dealloc, dmsg, DEBUG_LEVEL
    use string
    use grid, only: imx, jmx, setup_grid, destroy_grid, grid_x
    use geometry, only: xnx, xny, ynx, yny, xA, yA, volume, setup_geometry, &
            destroy_geometry
    use state, only: qp, qp_inf, density, x_speed, y_speed, pressure, &
            density_inf, x_speed_inf, y_speed_inf, pressure_inf, gm, R_gas, &
            setup_state, destroy_state, set_ghost_cell_data, writestate, &
            mu_ref, T_ref, Sutherland_temp, Pr
    use face_interpolant, only: interpolant, &
            x_sound_speed_left, x_sound_speed_right, &
            y_sound_speed_left, y_sound_speed_right, &
            y_pressure_left, y_pressure_right
    use scheme, only: scheme_name, residue, setup_scheme, destroy_scheme, &
            compute_residue
    
    use van_leer, only: mass_flux_sum

    implicit none
    private

    real, public :: CFL
    character, public :: time_stepping_method
    real, public :: global_time_step
    character(len=INTERPOLANT_NAME_LENGTH) :: time_step_accuracy
    real, dimension(:, :, :), allocatable :: qp_n, dEdx_1, dEdx_2, dEdx_3
    real :: tolerance
    integer, public :: max_iters
    integer, public :: checkpoint_iter, checkpoint_iter_count
    real, public :: resnorm, resnorm_0
    real, public, dimension(:, :), allocatable :: delta_t
    integer, public :: iter
    real :: sim_clock

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
                read(CONFIG_FILE_UNIT, '(A)', iostat=ios) buf
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
            character(len=STRING_BUFFER_LENGTH) :: buf
            integer :: ios

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
            read(buf, *, iostat=ios) time_stepping_method, global_time_step
            if (ios /= 0) then
                read(buf, *) time_stepping_method
                global_time_step = -1
            end if
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='time_stepping_method = ' + time_stepping_method)
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='global_time_step = ' + global_time_step)

            call get_next_token(buf)
            read(buf, *) time_step_accuracy
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='time_step_accuracy  = ' + time_step_accuracy)

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
            read(buf, *) checkpoint_iter
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='checkpoint_iter = ' + checkpoint_iter)

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

            call get_next_token(buf)
            read(buf, *) mu_ref
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='mu_reference = ' + mu_ref)

            call get_next_token(buf)
            read(buf, *) T_ref
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='T_reference = ' + T_ref)

            call get_next_token(buf)
            read(buf, *) Sutherland_temp
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='Sutherland temperature  = ' + Sutherland_temp)

            call get_next_token(buf)
            read(buf, *) Pr
            call dmsg(5, 'solver', 'read_config_file', &
                    msg='Prandtl Number = ' + Pr)

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
            checkpoint_iter_count = 0
            call checkpoint()  ! Create an initial dump fil

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

            sim_clock = 0.
            iter = 0
            resnorm = 1.
            resnorm_0 = 1.

        end subroutine initmisc

        subroutine deallocate_misc()

            implicit none
            
            call dmsg(1, 'solver', 'deallocate_misc')

            call dealloc(delta_t)
            
            select case (time_step_accuracy)
                case ("none")
                    ! Do nothing
                    continue
                case ("RK4")
                    call destroy_RK4_time_step()
                case default
                    call dmsg(5, 'solver', 'time_setup_deallocate_memory', &
                                'time step accuracy not recognized.')
                    stop
            end select

        end subroutine deallocate_misc

        subroutine destroy_RK4_time_step()
    
            implicit none

            call dealloc(qp_n)
            call dealloc(dEdx_1)
            call dealloc(dEdx_2)
            call dealloc(dEdx_3)

        end subroutine destroy_RK4_time_step

        subroutine setup_RK4_time_step()
    
            implicit none

            call alloc(qp_n, 0, imx, 0, jmx, 1, 4, &
                    errmsg='Error: Unable to allocate memory for qp_n.')
            call alloc(dEdx_1, 1, imx-1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for dEdx_1.')
            call alloc(dEdx_2, 1, imx-1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for dEdx_2.')
            call alloc(dEdx_3, 1, imx-1, 1, jmx-1, 1, 4, &
                    errmsg='Error: Unable to allocate memory for dEdx_3.')

        end subroutine setup_RK4_time_step

        subroutine allocate_memory()

            implicit none
            
            call dmsg(1, 'solver', 'allocate_memory')

            call alloc(delta_t, 1, imx-1, 1, jmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')

            select case (time_step_accuracy)
                case ("none")
                    ! Do nothing
                    continue
                case ("RK4")
                    call setup_RK4_time_step()
                case default
                    call dmsg(5, 'solver', 'time_setup_allocate_memory', &
                                'time step accuracy not recognized.')
                    stop
            end select

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
            real, dimension(imx-1, jmx-1) :: cell_sound_speed
            
            call dmsg(1, 'solver', 'compute_local_time_step')

            cell_sound_speed = sqrt(gm*pressure(1:imx-1, 1:jmx-1) / &
                                       density(1:imx-1, 1:jmx-1))

            ! For left face
            lmx1(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * xnx(1:imx-1, 1:jmx-1)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * xny(1:imx-1, 1:jmx-1))) + &
                    cell_sound_speed
            ! For bottom face
            lmx2(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * ynx(1:imx-1, 1:jmx-1)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * yny(1:imx-1, 1:jmx-1))) + &
                    cell_sound_speed
            ! For right face
            lmx3(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * xnx(2:imx, 1:jmx-1)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * xny(2:imx, 1:jmx-1))) + &
                    cell_sound_speed
            ! For top face
            lmx4(:, :) = abs( &
                    (x_speed(1:imx-1, 1:jmx-1) * ynx(1:imx-1, 2:jmx)) + &
                    (y_speed(1:imx-1, 1:jmx-1) * yny(1:imx-1, 2:jmx))) + &
                    cell_sound_speed
            lmxsum(:, :) = (xA(1:imx-1, 1:jmx-1) * lmx1) + &
                    (yA(1:imx-1, 1:jmx-1) * lmx2) + &
                    (xA(2:imx, 1:jmx-1) * lmx3) + &
                    (yA(1:imx-1, 2:jmx) * lmx4)
            
            delta_t = 1. / lmxsum
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

            if (global_time_step > 0) then
                delta_t = global_time_step
            else
                call compute_local_time_step()
                ! The global time step is the minimum of all the local time
                ! steps.
                delta_t = minval(delta_t)
            end if

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

        subroutine update_simulation_clock
            !-----------------------------------------------------------
            ! Update the simulation clock
            !
            ! It is sometimes useful to know what the simulation time is
            ! at every iteration so that a comparison with an analytical
            ! solution is possible. Since, the global timesteps used may
            ! not be uniform, we need to track this explicitly.
            !
            ! Of course, it makes sense to track this only if the time 
            ! stepping is global and not local. If the time stepping is
            ! local, the simulation clock is set to -1. If it is global
            ! it is incremented according to the time step found.
            !-----------------------------------------------------------

            implicit none
            if (time_stepping_method .eq. 'g' .and. sim_clock >= 0.) then
                sim_clock = sim_clock + minval(delta_t)
            else if (time_stepping_method .eq. 'l') then
                sim_clock = -1
            end if

        end subroutine update_simulation_clock

        subroutine get_next_solution()

            implicit none

            select case (time_step_accuracy)
                case ("none")
                    call update_solution
                case ("RK4")
                    call RK4_update_solution()
                case default
                    call dmsg(5, 'solver', 'get_next solution', &
                                'time step accuracy not recognized.')
                    stop
            end select

        end subroutine get_next_solution

        subroutine RK4_update_solution()

            implicit none
            real, dimension(1:imx-1, 1:jmx-1, 4) :: qp_temp

            ! qp at various stages is not stored but over written
            ! The residue multiplied by the inverse of the jacobian
            ! is stored for the final update equation

            ! Stage 1 is identical to stage (n)
            ! Store qp(n)
            qp_n = qp
            dEdx_1 = get_residue_primitive()
            
            ! Stage 2
            ! Not computing delta_t since qp(1) = qp(n)
            ! Update solution will over write qp
            delta_t = 0.5 * delta_t  ! delta_t(1)
            call update_solution()

            ! Stage 3
            call set_ghost_cell_data()
            call compute_residue()
            call compute_time_step() ! delta_t(2)
            dEdx_2 = get_residue_primitive()
            delta_t = 0.5 * delta_t
            call update_solution()

            ! Stage 4
            call set_ghost_cell_data()
            call compute_residue()
            call compute_time_step() ! delta_t(3)
            dEdx_3 = get_residue_primitive()
            call update_solution()
            ! qp now is qp_4

            ! Use qp(4)
            call set_ghost_cell_data()
            call compute_residue()
            call compute_time_step() ! delta_t(4)

            ! Calculating dEdx_4 in-situ 
            qp_temp(:, :, 1) = qp_n(1:imx-1, 1:jmx-1, 1) - &
                               (((dEdx_1(:, :, 1) / 6.0) + &
                                 (dEdx_2(:, :, 1) / 3.0) + &
                                 (dEdx_3(:, :, 1) / 3.0) + &
                                 (residue(:, :, 1) / 6.0)) * &
                                delta_t / volume)
            qp_temp(:, :, 2) = qp_n(1:imx-1, 1:jmx-1, 2) - &
                               (((dEdx_1(:, :, 2) / 6.0) + &
                                 (dEdx_2(:, :, 2) / 3.0) + &
                                 (dEdx_3(:, :, 2) / 3.0) + &
                                 (( (-1. * qp(1:imx-1, 1:jmx-1, 2) / &
                                     qp(1:imx-1, 1:jmx-1, 1) * residue(:, :, 1)) + &
                                     (residue(:, :, 2) / qp(1:imx-1, 1:jmx-1, 1)) ) / 6.0) &
                                ) * delta_t / volume)
            qp_temp(:, :, 3) = qp_n(1:imx-1, 1:jmx-1, 3) - &
                               (((dEdx_1(:, :, 3) / 6.0) + &
                                 (dEdx_2(:, :, 3) / 3.0) + &
                                 (dEdx_3(:, :, 3) / 3.0) + &
                                 (( (-1. * qp(1:imx-1, 1:jmx-1, 3) / &
                                     qp(1:imx-1, 1:jmx-1, 1) * residue(:, :, 1)) + &
                                     (residue(:, :, 3) / qp(1:imx-1, 1:jmx-1, 1)) ) / 6.0) &
                                ) * delta_t / volume)
            qp_temp(:, :, 4) = qp_n(1:imx-1, 1:jmx-1, 4) - &
                               (((dEdx_1(:, :, 4) / 6.0) + &
                                 (dEdx_2(:, :, 4) / 3.0) + &
                                 (dEdx_3(:, :, 4) / 3.0) + &
                                 (( (0.5 * (gm - 1.) * (qp(1:imx-1, 1:jmx-1, 2)**2. + &
                                     qp(1:imx-1, 1:jmx-1, 3)**2.) * residue(:, :, 1)) + &
                                    ( -(gm - 1.) * qp(1:imx-1, 1:jmx-1, 2) * residue(:, :, 2)) + &
                                    (- (gm - 1.) * qp(1:imx-1, 1:jmx-1, 3) * residue(:, :, 3)) + &
                                    ((gm - 1.) * residue(:, :, 4)) ) / 6.0) &
                                ) * delta_t / volume)
            
            qp(1:imx-1, 1:jmx-1, :) = qp_temp(1:imx-1, 1:jmx-1, :)
            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
            end if

        end subroutine RK4_update_solution

        function get_residue_primitive() result(dEdx)

            implicit none

            real, dimension(1:imx-1, 1:jmx-1, 4) :: dEdx
            dEdx(:, :, 1) = residue(:, :, 1)
            dEdx(:, :, 2) = ( (-1. * qp(1:imx-1, 1:jmx-1, 2) / &
                                qp(1:imx-1, 1:jmx-1, 1) * residue(:, :, 1)) + &
                              (residue(:, :, 2) / qp(1:imx-1, 1:jmx-1, 1)) )
            dEdx(:, :, 3) = ( (-1. * qp(1:imx-1, 1:jmx-1, 3) / &
                                qp(1:imx-1, 1:jmx-1, 1) * residue(:, :, 1)) + &
                              (residue(:, :, 3) / qp(1:imx-1, 1:jmx-1, 1)) )
            dEdx(:, :, 4) = ( (0.5 * (gm - 1.) * (qp(1:imx-1, 1:jmx-1, 2)**2. + &
                              qp(1:imx-1, 1:jmx-1, 3)**2.) * residue(:, :, 1)) + &
                            ( -(gm - 1.) * qp(1:imx-1, 1:jmx-1, 2) * residue(:, :, 2)) + &
                            (- (gm - 1.) * qp(1:imx-1, 1:jmx-1, 3) * residue(:, :, 3)) + &
                            ((gm - 1.) * residue(:, :, 4)) )

        end function get_residue_primitive

        subroutine update_solution()
            !-----------------------------------------------------------
            ! Update the solution using the residue and time step
            !-----------------------------------------------------------

            implicit none
            !real, dimension(4, 4) :: psi_inv
            !integer :: i, j
            real, dimension(0:imx, 0:jmx, 4) :: qp_temp
            
            call dmsg(1, 'solver', 'update_solution')

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
            end if

        end subroutine update_solution

        subroutine checkpoint()
            !-----------------------------------------------------------
            ! Create a checkpoint dump file if the time has come
            !-----------------------------------------------------------

            implicit none

            character(len=FILE_NAME_LENGTH) :: filename

            if (checkpoint_iter .ne. 0) then
                if (mod(iter, checkpoint_iter) == 0) then
                    write(filename, '(A,I5.5,A)') 'output', checkpoint_iter_count, '.fvtk'
                    checkpoint_iter_count = checkpoint_iter_count + 1
                    call writestate(filename, 'Simulation clock: ' + sim_clock)
                    call dmsg(3, 'solver', 'checkpoint', &
                            'Checkpoint created at iteration: ' + iter)
                    call write_surface_pressure()
                end if
            end if

        end subroutine checkpoint

        subroutine write_surface_pressure()

            implicit none

            integer :: i
            open(71, file='pressure-'//interpolant)
            do i = 1, imx-1                
                write(71, *) grid_x(i, 1), y_pressure_right(i, 1)
            end do

            close(71)

        end subroutine write_surface_pressure

        subroutine step()
            !-----------------------------------------------------------
            ! Perform one time step iteration
            !
            ! This subroutine performs one iteration by stepping through
            ! time once.
            !-----------------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'step')

            call set_ghost_cell_data()
            call compute_residue()
            call dmsg(1, 'solver', 'step', 'Residue computed.')
            call compute_time_step()
         !  call update_solution()
            call get_next_solution()
            call update_simulation_clock()
            iter = iter + 1
            call compute_residue_norm()
            if (iter .eq. 1) then
                resnorm_0 = resnorm
            end if
          ! resnorm = resnorm / resnorm_0
            write(RESNORM_FILE_UNIT, *) resnorm
            call checkpoint()
            if (iter .eq. max_iters) then
                call write_surface_pressure()
            end if

        end subroutine step
        
        subroutine compute_residue_norm()

            implicit none
            
            call dmsg(1, 'solver', 'compute_residue_norm')

            resnorm = sum(sqrt( &
                    (residue(:, :, 1) / &
                        (density_inf * sqrt(x_speed_inf*x_speed_inf + &
                         y_speed_inf*y_speed_inf) )) ** 2. + &
                    (residue(:, :, 2) / &
                        (density_inf * sqrt(x_speed_inf*x_speed_inf + &
                        y_speed_inf*y_speed_inf) ** 2.)) ** 2. + &
                    (residue(:, :, 3) / &
                        (density_inf * sqrt(x_speed_inf*x_speed_inf + &
                        y_speed_inf*y_speed_inf) ** 2.)) ** 2. + &
                    (residue(:, :, 4) / &
                        (density_inf ** 2. * sqrt(x_speed_inf*x_speed_inf + &
                        y_speed_inf*y_speed_inf) ** 3.)) ** 2. &
                    ))

         !   resnorm = sum(residue(:, :, 1))


    !       resnorm = mass_flux_sum / resnorm_0

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

