module solver

  use global, only: CONFIG_FILE_UNIT 
  use global, only: RESNORM_FILE_UNIT 
  use global, only: FILE_NAME_LENGTH
  use global, only: STRING_BUFFER_LENGTH 
  use global, only: INTERPOLANT_NAME_LENGTH

  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx

  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
  use global_vars, only : volume
    
  use global_vars, only : n_var
  use global_vars, only : sst_n_var
  use global_vars, only : qp
  use global_vars, only : qp_inf
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : pressure
  use global_vars, only : tk
  use global_vars, only : tw
  use global_vars, only : tk_inf
  use global_vars, only : tw_inf
  use global_vars, only : gm
  use global_vars, only : R_gas
  use global_vars, only : mu_ref
  use global_vars, only : T_ref
  use global_vars, only : Sutherland_temp
  use global_vars, only : Pr

  use global_vars, only : qp_n
  use global_vars, only : dEdx_1
  use global_vars, only : dEdx_2
  use global_vars, only : dEdx_3
  use global_vars, only : resnorm, resnorm_0
  use global_vars, only : cont_resnorm, cont_resnorm_0
  use global_vars, only : x_mom_resnorm, x_mom_resnorm_0
  use global_vars, only : y_mom_resnorm, y_mom_resnorm_0
  use global_vars, only : z_mom_resnorm, z_mom_resnorm_0
  use global_vars, only : energy_resnorm, energy_resnorm_0
  use global_vars, only : write_percision
  use global_vars, only : CFL
  use global_vars, only : tolerance
  use global_vars, only : min_iter
  use global_vars, only : max_iters
  use global_vars, only : current_iter
  use global_vars, only : checkpoint_iter
  use global_vars, only : checkpoint_iter_count
  use global_vars, only : time_stepping_method
  use global_vars, only : time_step_accuracy
  use global_vars, only : global_time_step
  use global_vars, only : delta_t
  use global_vars, only : sim_clock
  use global_vars, only : turbulence
  use global_vars, only : supersonic_flag

  use global_vars, only: F_p
  use global_vars, only: G_p
  use global_vars, only: H_p
  use global_vars, only: mass_residue
  use global_vars, only: x_mom_residue
  use global_vars, only: y_mom_residue
  use global_vars, only: z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only: TKE_residue
  use global_vars, only: omega_residue
  use global_vars, only: res_write_interval

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL

  use string
  use read, only : read_controls
  use read, only : read_scheme
  use read, only : read_flow

  use grid, only: setup_grid, destroy_grid
  use geometry, only: setup_geometry, destroy_geometry
  use state, only:  setup_state, destroy_state

  use face_interpolant, only: interpolant, &
          x_qp_left, x_qp_right, &
          y_qp_left, y_qp_right, &
          z_qp_left, z_qp_right, compute_face_interpolant, &
          extrapolate_cell_averages_to_faces
  use scheme, only: scheme_name, setup_scheme, destroy_scheme, &
          compute_fluxes, compute_residue
  use source, only: setup_source, destroy_source, add_source_term_residue, &
                    compute_gradients_cell_centre
  use boundary_conditions, only: setup_boundary_conditions, &
          apply_boundary_conditions, set_wall_bc_at_faces, &
          destroy_boundary_conditions
  use wall_dist, only: setup_wall_dist, destroy_wall_dist, find_wall_dist
  use viscous, only: compute_viscous_fluxes
  use boundary_state_reconstruction, only: reconstruct_boundary_state
  use layout, only: process_id, grid_file_buf, bc_file, &
  get_process_data, read_layout_file, total_process
  use parallel, only: allocate_buffer_cells,send_recv
!  use state, only: turbulence
  use resnorm_, only : write_resnorm, setup_resnorm, destroy_resnorm
  use dump_solution, only : checkpoint
  include "turbulence_models/include/solver/import_module.inc"

    use mpi
    implicit none
!    include "mpif.h"
    private

    real, dimension(:), allocatable, target :: qp_temp
    real, pointer :: density_temp, x_speed_temp, &
                                         y_speed_temp, z_speed_temp, &
                                         pressure_temp
    include "turbulence_models/include/solver/variables_deceleration.inc"

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: step
    public :: converged

    contains


        subroutine setup_solver()
            
            implicit none

            call dmsg(1, 'solver', 'setup_solver')
            call get_process_data() ! parallel calls
            call read_layout_file(process_id) ! reads layout file calls
            
            call read_controls()
            call read_scheme()
            call read_flow()
                  !todo make it general for all turbulence model
                  if(turbulence=="sst")then
                    n_var=n_var+sst_n_var
                  end if
            call setup_grid(grid_file_buf)
            call setup_geometry()
            call setup_state()
            call setup_boundary_conditions(bc_file)
            call allocate_memory()
            call allocate_buffer_cells(3) !parallel buffers
            call setup_scheme()
            if(turbulence /= 'none') then
              call setup_wall_dist
              call find_wall_dist()
            end if
            if(mu_ref /= 0. .or. turbulence /= 'none') then
              call setup_source()
            end if
            call link_aliases_solver()
            call setup_resnorm()
            call initmisc()
            checkpoint_iter_count = 0
            call checkpoint()  ! Create an initial dump file
            call dmsg(1, 'solver', 'setup_solver', 'Setup solver complete')

        end subroutine setup_solver

        subroutine destroy_solver()

            implicit none
            
            call dmsg(1, 'solver', 'destroy_solver')

            if(mu_ref /= 0. .or. turbulence /= 'none')  then 
              call destroy_source()
            end if
            if(turbulence /= 'none') then
              call destroy_wall_dist()
            end if
            call destroy_scheme()
            call deallocate_misc()
            call unlink_aliases_solver()
            call destroy_state()
            call destroy_geometry()
            call destroy_grid()
            call destroy_resnorm()

        end subroutine destroy_solver

        subroutine initmisc()
            
            implicit none
            
            call dmsg(1, 'solver', 'initmisc')

            sim_clock = 0.
            current_iter = 0
!            resnorm = 1.
!            resnorm_0 = 1.

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

            call alloc(qp_n, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for qp_n.')
            call alloc(dEdx_1, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for dEdx_1.')
            call alloc(dEdx_2, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for dEdx_2.')
            call alloc(dEdx_3, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for dEdx_3.')

        end subroutine setup_RK4_time_step

        subroutine allocate_memory()

            implicit none
            
            call dmsg(1, 'solver', 'allocate_memory')

            call alloc(delta_t, 1, imx-1, 1, jmx-1, 1, kmx-1, &
                    errmsg='Error: Unable to allocate memory for delta_t.')
            call alloc(qp_temp, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for qp_temp.')

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

        subroutine unlink_aliases_solver()

            implicit none

            nullify(density_temp)
            nullify(x_speed_temp)
            nullify(y_speed_temp)
            nullify(z_speed_temp)
            nullify(pressure_temp)
            include "turbulence_models/include/solver/unlink_aliases_solver.inc"

        end subroutine unlink_aliases_solver

        subroutine link_aliases_solver()

            implicit none

            call dmsg(1, 'solver', 'link_aliases_solver')

            density_temp => qp_temp(1)
            x_speed_temp => qp_temp(2)
            y_speed_temp => qp_temp(3)
            z_speed_temp => qp_temp(4)
            pressure_temp => qp_temp(5)
            include "turbulence_models/include/solver/link_aliases_solver.inc"
        end subroutine link_aliases_solver

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

            real :: lmx1, lmx2, lmx3, lmx4, lmx5, lmx6, lmxsum
            real :: x_sound_speed_avg, y_sound_speed_avg, z_sound_speed_avg
            integer :: i, j, k

            call dmsg(1, 'solver', 'compute_local_time_step')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
               ! For orientation, refer to the report. The standard i,j,k 
               ! direction are marked. All orientation notations are w.r.t 
               ! to the perspective shown in the image.

               ! Faces with lower index
               x_sound_speed_avg = 0.5 * (sqrt(gm * x_qp_left(i, j, k, 5) / &
                                                    x_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * x_qp_right(i, j, k, 5) / &
                                                    x_qp_right(i, j, k, 1)) )
               y_sound_speed_avg = 0.5 * (sqrt(gm * y_qp_left(i, j, k, 5) / &
                                                    y_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * y_qp_right(i, j, k, 5) / &
                                                    y_qp_right(i, j, k, 1)) )
               z_sound_speed_avg = 0.5 * (sqrt(gm * z_qp_left(i, j, k, 5) / &
                                                    z_qp_left(i, j, k, 1)) + &
                                          sqrt(gm * z_qp_right(i, j, k, 5) / &
                                                    z_qp_right(i, j, k, 1)) )
               
               ! For left face: i.e., lower index face along xi direction
               lmx1 = abs( &
                    (x_speed(i, j, k) * xnx(i, j, k)) + &
                    (y_speed(i, j, k) * xny(i, j, k)) + &
                    (z_speed(i, j, k) * xnz(i, j, k))) + &
                    x_sound_speed_avg
               ! For front face, i.e., lower index face along eta direction
               lmx2 = abs( &
                    (x_speed(i, j, k) * ynx(i, j, k)) + &
                    (y_speed(i, j, k) * yny(i, j, k)) + &
                    (z_speed(i, j, k) * ynz(i, j, k))) + &
                    y_sound_speed_avg
               ! For bottom face, i.e., lower index face along zeta direction
               lmx3 = abs( &
                    (x_speed(i, j, k) * znx(i, j, k)) + &
                    (y_speed(i, j, k) * zny(i, j, k)) + &
                    (z_speed(i, j, k) * znz(i, j, k))) + &
                    z_sound_speed_avg

               ! Faces with higher index
               x_sound_speed_avg = 0.5 * (sqrt(gm * x_qp_left(i+1,j,k,5) / x_qp_left(i+1,j,k,1)) + &
                                          sqrt(gm * x_qp_right(i+1,j,k,5) / x_qp_right(i+1,j,k,1)) )
               y_sound_speed_avg = 0.5 * (sqrt(gm * y_qp_left(i,j+1,k,5) / y_qp_left(i,j+1,k,1)) + &
                                          sqrt(gm * y_qp_right(i,j+1,k,5) / y_qp_right(i,j+1,k,1)) )
               z_sound_speed_avg = 0.5 * (sqrt(gm * z_qp_left(i,j,k+1,5) / z_qp_left(i,j,k+1,1)) + &
                                          sqrt(gm * z_qp_right(i,j,k+1,5) / z_qp_right(i,j,k+1,1)) )
               
               ! For right face, i.e., higher index face along xi direction
               lmx4 = abs( &
                    (x_speed(i+1, j, k) * xnx(i+1, j, k)) + &
                    (y_speed(i+1, j, k) * xny(i+1, j, k)) + &
                    (z_speed(i+1, j, k) * xnz(i+1, j, k))) + &
                    x_sound_speed_avg
               ! For back face, i.e., higher index face along eta direction
               lmx5 = abs( &
                    (x_speed(i, j+1, k) * ynx(i, j+1, k)) + &
                    (y_speed(i, j+1, k) * yny(i, j+1, k)) + &
                    (z_speed(i, j+1, k) * ynz(i, j+1, k))) + &
                    y_sound_speed_avg
               ! For top face, i.e., higher index face along zeta direction
               lmx6 = abs( &
                    (x_speed(i, j, k+1) * znx(i, j, k+1)) + &
                    (y_speed(i, j, k+1) * zny(i, j, k+1)) + &
                    (z_speed(i, j, k+1) * znz(i, j, k+1))) + &
                    z_sound_speed_avg

               lmxsum = (xA(i, j, k) * lmx1) + &
                        (yA(i, j, k) * lmx2) + &
                        (zA(i, j, k) * lmx3) + &
                        (xA(i+1, j, k) * lmx4) + &
                        (yA(i, j+1, k) * lmx5) + &
                        (zA(i, j, k+1) * lmx6)
            
               delta_t(i, j, k) = 1. / lmxsum
               delta_t(i, j, k) = delta_t(i, j, k) * volume(i, j, k) * CFL
              end do
             end do
            end do

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
                    call update_solution()
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
            integer :: i, j, k

            ! qp at various stages is not stored but over written
            ! The residue multiplied by the inverse of the jacobian
            ! is stored for the final update equation

            ! Stage 1 is identical to stage (n)
            ! Store qp(n)
            qp_n = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1:n_var)
            dEdx_1 = get_residue_primitive()
            
            ! Stage 2
            ! Not computing delta_t since qp(1) = qp(n)
            ! Update solution will over write qp
            delta_t = 0.5 * delta_t  ! delta_t(1)
            call update_solution()

            ! Stage 3
            call sub_step()
            dEdx_2 = get_residue_primitive()
            delta_t = 0.5 * delta_t
            call update_solution()

            ! Stage 4
            call sub_step()
            dEdx_3 = get_residue_primitive()
            call update_solution()

            ! qp now is qp_4
            ! Use qp(4)
            call sub_step()

            ! Calculating dEdx_4 in-situ and updating the solution
            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                density_temp  = qp_n(i, j, k, 1) - &
                               (((dEdx_1(i, j, k, 1) / 6.0) + &
                                 (dEdx_2(i, j, k, 1) / 3.0) + &
                                 (dEdx_3(i, j, k, 1) / 3.0) + &
                                 (mass_residue(i, j, k) / 6.0)) * &
                                delta_t(i, j, k) / volume(i, j, k))
                x_speed_temp = qp_n(i, j, k, 2) - &
                               (((dEdx_1(i, j, k, 2) / 6.0) + &
                                 (dEdx_2(i, j, k, 2) / 3.0) + &
                                 (dEdx_3(i, j, k, 2) / 3.0) + &
                                 (( (-1 * x_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( x_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                y_speed_temp = qp_n(i, j, k, 3) - &
                               (((dEdx_1(i, j, k, 3) / 6.0) + &
                                 (dEdx_2(i, j, k, 3) / 3.0) + &
                                 (dEdx_3(i, j, k, 3) / 3.0) + &
                                 (( (-1 * y_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( y_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                z_speed_temp = qp_n(i, j, k, 4) - &
                               (((dEdx_1(i, j, k, 4) / 6.0) + &
                                 (dEdx_2(i, j, k, 4) / 3.0) + &
                                 (dEdx_3(i, j, k, 4) / 3.0) + &
                                 (( (-1 * z_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( z_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
                pressure_temp = qp_n(i, j, k, 5) - &
                               (((dEdx_1(i, j, k, 5) / 6.0) + &
                                 (dEdx_2(i, j, k, 5) / 3.0) + &
                                 (dEdx_3(i, j, k, 5) / 3.0) + &
                                 (( (0.5 * (gm - 1.) * ( x_speed(i, j, k) ** 2. + &
                                                         y_speed(i, j, k) ** 2. + &
                                                         z_speed(i, j, k) ** 2.) * &
                                                        mass_residue(i, j, k)) + &
                       (- (gm - 1.) * x_speed(i, j, k) * x_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * y_speed(i, j, k) * y_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * z_speed(i, j, k) * z_mom_residue(i, j, k)) + &
                       ((gm - 1.) * energy_residue(i, j, k)) ) / 6.0) &
                                ) * delta_t(i, j, k) / volume(i, j, k))
            
                density(i, j, k) = density_temp
                x_speed(i, j, k) = x_speed_temp
                y_speed(i, j, k) = y_speed_temp
                z_speed(i, j, k) = z_speed_temp
                pressure(i, j, k) = pressure_temp
                include "turbulence_models/include/solver/RK4_update_solution.inc"
              end do
             end do
            end do

            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
            end if

        end subroutine RK4_update_solution

        function get_residue_primitive() result(dEdx)

            implicit none

            real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, n_var) :: dEdx
            dEdx(:, :, :, 1) = mass_residue
            dEdx(:, :, :, 2) = ( (-1 * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( x_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 3) = ( (-1 * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( y_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 4) = ( (-1 * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                                     mass_residue) + &
                             ( z_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
            dEdx(:, :, :, 5) = ( (0.5 * (gm - 1.) * ( x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
                                                      y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
                                                      z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2.) * &
                                                    mass_residue) + &
                       (- (gm - 1.) * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * x_mom_residue) + &
                       (- (gm - 1.) * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * y_mom_residue) + &
                       (- (gm - 1.) * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * z_mom_residue) + &
                       ((gm - 1.) * energy_residue) )

            include "turbulence_models/include/solver/get_residue_primitive.inc"

        end function get_residue_primitive

        subroutine update_solution()
            !-----------------------------------------------------------
            ! Update the solution using the residue and time step
            !-----------------------------------------------------------

            implicit none
            integer :: i, j, k
            
            call dmsg(1, 'solver', 'update_solution')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
               density_temp = density(i, j, k) - &
                            (mass_residue(i, j, k) * &
                            delta_t(i, j, k) / volume(i, j, k))

               x_speed_temp = x_speed(i, j, k) - &
                            (( (-1 * x_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( x_mom_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               y_speed_temp = y_speed(i, j, k) - &
                            (( (-1 * y_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( y_mom_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               z_speed_temp = z_speed(i, j, k) - &
                            (( (-1 * z_speed(i, j, k) / density(i, j, k) * &
                                     mass_residue(i, j, k)) + &
                             ( z_mom_residue(i, j, k) / density(i, j, k)) ) * &
                            delta_t(i, j, k) / volume(i, j, k))

               pressure_temp = pressure(i, j, k) - &
                   ( ( (0.5 * (gm - 1.) * ( x_speed(i, j, k) ** 2. + &
                                            y_speed(i, j, k) ** 2. + &
                                            z_speed(i, j, k) ** 2.) * &
                                          mass_residue(i, j, k)) + &
                       (- (gm - 1.) * x_speed(i, j, k) * x_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * y_speed(i, j, k) * y_mom_residue(i, j, k)) + &
                       (- (gm - 1.) * z_speed(i, j, k) * z_mom_residue(i, j, k)) + &
                       ((gm - 1.) * energy_residue(i, j, k)) ) * &
                       delta_t(i, j, k) / volume(i, j, k) ) 

               density(i, j, k) = density_temp
               x_speed(i, j, k) = x_speed_temp
               y_speed(i, j, k) = y_speed_temp
               z_speed(i, j, k) = z_speed_temp
               pressure(i, j, k) = pressure_temp
               include "turbulence_models/include/solver/update_solution.inc"
              end do
             end do
            end do

            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
                stop
            end if

        end subroutine update_solution


        subroutine sub_step()

            implicit none

            !TODO: Better name for this??

            call dmsg(1, 'solver', 'sub_step')
            call send_recv(3) ! parallel call-argument:no of layers 
            call apply_boundary_conditions()
            call compute_face_interpolant()
            call reconstruct_boundary_state(interpolant)
            call set_wall_bc_at_faces()
            call compute_fluxes()
            if (mu_ref /= 0.0) then
                if (interpolant /= "none") then
                    call extrapolate_cell_averages_to_faces()
                    call set_wall_bc_at_faces()
                end if
                call compute_gradients_cell_centre()
                call add_source_term_residue()
                call compute_viscous_fluxes(F_p, G_p, H_p)
            end if
            call compute_residue()
            call dmsg(1, 'solver', 'step', 'Residue computed.')
            call compute_time_step()

        end subroutine sub_step
        
        subroutine step()
            !-----------------------------------------------------------
            ! Perform one time step iteration
            !
            ! This subroutine performs one iteration by stepping through
            ! time once.
            !-----------------------------------------------------------

            implicit none
            call dmsg(1, 'solver', 'step')

            if (process_id==0) then
              print*, current_iter
            end if
            call sub_step()

            call get_next_solution()
            call update_simulation_clock()
            current_iter = current_iter + 1

            if ( mod(current_iter,res_write_interval) == 0 .or. current_iter == 1) then
              call write_resnorm()
            end if

            call checkpoint()

        end subroutine step

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

