module solver

  use global, only: CONFIG_FILE_UNIT 
  use global, only: RESNORM_FILE_UNIT 
  use global, only: FILE_NAME_LENGTH
  use global, only: STRING_BUFFER_LENGTH 
  use global, only: INTERPOLANT_NAME_LENGTH
  use global, only: STOP_FILE_UNIT
  use global, only: stop_file
  use global_vars, only : want_to_stop
  use global_vars, only : Halt

  use global_kkl , only : cmu
  use global_kkl , only : cd1
  use global_kkl , only : eta
  use global_kkl , only : fphi
  use global_sst , only : beta1
  use global_sst , only : beta2
  use global_sst , only : bstar
!  use global_sst , only : sst_F1
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx

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
  use global_vars, only : tkl
  use global_vars, only : tk_inf
  use global_vars, only : tw_inf
  use global_vars, only : tkl_inf
  use global_vars, only : gm
  use global_vars, only : R_gas
  use global_vars, only : mu_ref
  use global_vars, only : mu

  use global_vars, only : qp_n
  use global_vars, only : dEdx_1
  use global_vars, only : dEdx_2
  use global_vars, only : dEdx_3
  use global_vars, only : CFL
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
  use global_vars, only: KL_residue
  use global_vars, only: dissipation_residue
  use global_vars, only: tv_residue
  use global_vars, only: res_write_interval
  use global_vars, only: r_list
  use global_vars, only: w_list
  use global_vars, only: Res_itr

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL

  use CC,    only: setupCC
  use CC,    only: destroyCC

  use string
  use read, only : read_input_and_controls

  use grid, only: setup_grid, destroy_grid
  use geometry, only: setup_geometry, destroy_geometry
  use state, only:  setup_state, destroy_state
  use gradients, only : setup_gradients
  use gradients, only : destroy_gradients

  use face_interpolant, only: interpolant, &
          x_qp_left, x_qp_right, &
          y_qp_left, y_qp_right, &
          z_qp_left, z_qp_right, compute_face_interpolant, &
          extrapolate_cell_averages_to_faces
  use scheme, only: scheme_name, setup_scheme, destroy_scheme, &
          compute_fluxes, compute_residue
  use source, only: add_source_term_residue, setup_source, destroy_source
  use wall_dist, only: setup_wall_dist, destroy_wall_dist, find_wall_dist
  use viscous, only: compute_viscous_fluxes
!  use turbulent_fluxes, only: compute_turbulent_fluxes
  use boundary_state_reconstruction, only: reconstruct_boundary_state
  use layout, only: process_id, grid_file_buf, bc_file, &
  get_process_data, read_layout_file, total_process
!  use parallel, only: allocate_buffer_cells,send_recv
  use interface, only: setup_interface
  use interface, only: destroy_interface
!  use state, only: turbulence
  use resnorm, only : find_resnorm, setup_resnorm, destroy_resnorm
  use dump_solution, only : checkpoint
  !use transport    , only : setup_transport
  !use transport    , only : destroy_transport
  !use transport    , only : calculate_transport
  use viscosity    , only : setup_viscosity
  use viscosity    , only : destroy_viscosity
  use viscosity    , only : calculate_viscosity
!  use blending_function , only : setup_sst_F1
!  use blending_function , only : destroy_sst_F1
!  use blending_function , only : calculate_sst_F1
  use wall        , only : write_surfnode
  include "turbulence_models/include/solver/import_module.inc"
  use bc, only: setup_bc
  use bc, only: destroy_bc
  use bc_primitive, only: populate_ghost_primitive
  use summon_grad_evaluation, only : evaluate_all_gradients
  use time , only : setup_time
  use time , only : destroy_time
  use time , only : compute_time_step
  use global_vars, only: dist
  use update, only : get_next_solution
  use update, only : setup_update
  use update, only : destroy_update
  use mapping, only : read_interface_map
!  use residual_smoothing, only: setup_implicit_residual_smoothing
!  use residual_smoothing, only: destroy_implicit_residual_smoothing
#include "error.inc"
#include "mpi.inc"
    private

    real, dimension(:), allocatable, target :: qp_temp
    real, pointer :: density_temp, x_speed_temp, &
                                         y_speed_temp, z_speed_temp, &
                                         pressure_temp
    include "turbulence_models/include/solver/variables_deceleration.inc"

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: iterate_one_more_time_step

    contains


        subroutine setup_solver()
            
            implicit none
            integer :: ierr

            call dmsg(1, 'solver', 'setup_solver')
            call get_process_data() ! parallel calls
            call read_layout_file(process_id) ! reads layout file calls
            
            call read_input_and_controls()
            call setup_grid(grid_file_buf)
            call setup_geometry()
            call setup_viscosity()
            call setup_state()
            call setup_gradients()
            call setup_source
            call setup_bc()
            call setup_time()
            call setup_update()
            !call setup_implicit_residual_smoothing()
            !call allocate_memory()
            !call allocate_buffer_cells(3) !parallel buffers
            call setup_interface()
            call setup_scheme()
            if(turbulence /= 'none') then
              call write_surfnode()
              call setup_wall_dist()
              call mpi_barrier(MPI_COMM_WORLD,ierr)
              call find_wall_dist()
            end if
            call setupCC()
!            if(mu_ref /= 0. .or. turbulence /= 'none') then
!              call setup_source()
!            end if
!            call setup_sst_F1()
            !call link_aliases_solver()
            call setup_resnorm()
            call initmisc()
            checkpoint_iter_count = 0
            call checkpoint()  ! Create an initial dump file
            current_iter=1
            call dmsg(1, 'solver', 'setup_solver', 'checkpoint')
            if(process_id==0)then
              open(STOP_FILE_UNIT, file=stop_file)
            end if
            call dmsg(1, 'solver', 'setup_solver', 'Setup solver complete')

        end subroutine setup_solver

        subroutine destroy_solver()

            implicit none
            
            call dmsg(1, 'solver', 'destroy_solver')

            if(process_id==0)then
              close(STOP_FILE_UNIT)
            end if
            call destroy_update()
            !call destroy_implicit_residual_smoothing()
            call destroy_viscosity()
!            if(mu_ref /= 0. .or. turbulence /= 'none')  then 
!              call destroy_source()
!            end if
            call destroy_gradients()
            call destroyCC()
            if(turbulence /= 'none') then
              call destroy_wall_dist()
            end if
            call destroy_scheme()
            !call deallocate_misc()
           ! call unlink_aliases_solver()
            call destroy_source()
            call destroy_state()
            call destroy_geometry()
            call destroy_grid()
            call destroy_resnorm()
            call destroy_interface()
!            if(turbulence=='sst')then
!              call destroy_sst_F1()
!            end if
            call destroy_time()
            call destroy_bc()

            if(allocated(r_list)) deallocate(r_list)
            if(allocated(w_list)) deallocate(w_list)

        end subroutine destroy_solver

        subroutine initmisc()
            
            implicit none
            
            call dmsg(1, 'solver', 'initmisc')

            sim_clock = 0.
            current_iter = 0
!            resnorm = 1.
!            resnorm_0 = 1.

        end subroutine initmisc

!        subroutine deallocate_misc()
!
!            implicit none
!            
!            call dmsg(1, 'solver', 'deallocate_misc')
!
!            call dealloc(delta_t)
!
!            select case (time_step_accuracy)
!                case ("none")
!                    ! Do nothing
!                    continue
!                case ("RK4")
!                    call destroy_RK4_time_step()
!                case default
!                    call dmsg(5, 'solver', 'time_setup_deallocate_memory', &
!                                'time step accuracy not recognized.')
!                    stop
!            end select
!
!        end subroutine deallocate_misc
!
!        subroutine destroy_RK4_time_step()
!    
!            implicit none
!
!            call dealloc(qp_n)
!            call dealloc(dEdx_1)
!            call dealloc(dEdx_2)
!            call dealloc(dEdx_3)
!
!        end subroutine destroy_RK4_time_step
!
!        subroutine setup_RK4_time_step()
!    
!            implicit none
!
!            call alloc(qp_n, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for qp_n.')
!            call alloc(dEdx_1, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for dEdx_1.')
!            call alloc(dEdx_2, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for dEdx_2.')
!            call alloc(dEdx_3, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for dEdx_3.')
!
!        end subroutine setup_RK4_time_step
!
!        subroutine allocate_memory()
!
!            implicit none
!            
!            call dmsg(1, 'solver', 'allocate_memory')
!
!!            call alloc(delta_t, 1, imx-1, 1, jmx-1, 1, kmx-1, &
!!                    errmsg='Error: Unable to allocate memory for delta_t.')
!            call alloc(qp_temp, 1, n_var, &
!                    errmsg='Error: Unable to allocate memory for qp_temp.')
!
!            select case (time_step_accuracy)
!                case ("none")
!                    ! Do nothing
!                    continue
!                case ("RK4")
!                    call setup_RK4_time_step()
!                case default
!                    call dmsg(5, 'solver', 'time_setup_allocate_memory', &
!                                'time step accuracy not recognized.')
!                    stop
!            end select
!
!        end subroutine allocate_memory
!
!        subroutine unlink_aliases_solver()
!
!            implicit none
!
!            nullify(density_temp)
!            nullify(x_speed_temp)
!            nullify(y_speed_temp)
!            nullify(z_speed_temp)
!            nullify(pressure_temp)
!            include "turbulence_models/include/solver/unlink_aliases_solver.inc"
!
!        end subroutine unlink_aliases_solver
!
!        subroutine link_aliases_solver()
!
!            implicit none
!
!            call dmsg(1, 'solver', 'link_aliases_solver')
!
!            density_temp => qp_temp(1)
!            x_speed_temp => qp_temp(2)
!            y_speed_temp => qp_temp(3)
!            z_speed_temp => qp_temp(4)
!            pressure_temp => qp_temp(5)
!            include "turbulence_models/include/solver/link_aliases_solver.inc"
!        end subroutine link_aliases_solver
!
!
!        subroutine get_next_solution()
!
!            implicit none
!
!            select case (time_step_accuracy)
!                case ("none")
!                    call update_solution()
!                case ("RK4")
!                    call RK4_update_solution()
!                case default
!                    call dmsg(5, 'solver', 'get_next solution', &
!                                'time step accuracy not recognized.')
!                    stop
!            end select
!
!        end subroutine get_next_solution
!
!        subroutine RK4_update_solution()
!
!            implicit none
!            integer :: i, j, k
!            real, dimension(1:imx-1,1:jmx-1,1:kmx-1) :: delta_t_0
!
!            delta_t_0 = delta_t
!            ! qp at various stages is not stored but over written
!            ! The residue multiplied by the inverse of the jacobian
!            ! is stored for the final update equation
!
!            ! Stage 1 is identical to stage (n)
!            ! Store qp(n)
!            qp_n = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1:n_var)
!            dEdx_1 = get_residue_primitive()
!            
!            ! Stage 2
!            ! Not computing delta_t since qp(1) = qp(n)
!            ! Update solution will over write qp
!            delta_t = 0.5 * delta_t_0  ! delta_t(1)
!            call update_solution()
!
!            ! Stage 3
!            call get_total_conservative_Residue()
!            dEdx_2 = get_residue_primitive()
!            delta_t = 0.5 * delta_t_0
!            call update_solution()
!
!            ! Stage 4
!            call get_total_conservative_Residue()
!            dEdx_3 = get_residue_primitive()
!            delta_t = delta_t_0
!            call update_solution()
!
!            ! qp now is qp_4
!            ! Use qp(4)
!            call get_total_conservative_Residue()
!
!            ! Calculating dEdx_4 in-situ and updating the solution
!            do k = 1, kmx - 1
!             do j = 1, jmx - 1
!              do i = 1, imx - 1
!                density_temp  = qp_n(i, j, k, 1) - &
!                               (((dEdx_1(i, j, k, 1) / 6.0) + &
!                                 (dEdx_2(i, j, k, 1) / 3.0) + &
!                                 (dEdx_3(i, j, k, 1) / 3.0) + &
!                                 (mass_residue(i, j, k) / 6.0)) * &
!                                delta_t(i, j, k) / volume(i, j, k))
!                x_speed_temp = qp_n(i, j, k, 2) - &
!                               (((dEdx_1(i, j, k, 2) / 6.0) + &
!                                 (dEdx_2(i, j, k, 2) / 3.0) + &
!                                 (dEdx_3(i, j, k, 2) / 3.0) + &
!                                 (( (-1 * x_speed(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( x_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
!                                ) * delta_t(i, j, k) / volume(i, j, k))
!                y_speed_temp = qp_n(i, j, k, 3) - &
!                               (((dEdx_1(i, j, k, 3) / 6.0) + &
!                                 (dEdx_2(i, j, k, 3) / 3.0) + &
!                                 (dEdx_3(i, j, k, 3) / 3.0) + &
!                                 (( (-1 * y_speed(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( y_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
!                                ) * delta_t(i, j, k) / volume(i, j, k))
!                z_speed_temp = qp_n(i, j, k, 4) - &
!                               (((dEdx_1(i, j, k, 4) / 6.0) + &
!                                 (dEdx_2(i, j, k, 4) / 3.0) + &
!                                 (dEdx_3(i, j, k, 4) / 3.0) + &
!                                 (( (-1 * z_speed(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( z_mom_residue(i, j, k) / density(i, j, k)) ) / 6.0) &
!                                ) * delta_t(i, j, k) / volume(i, j, k))
!                pressure_temp = qp_n(i, j, k, 5) - &
!                               (((dEdx_1(i, j, k, 5) / 6.0) + &
!                                 (dEdx_2(i, j, k, 5) / 3.0) + &
!                                 (dEdx_3(i, j, k, 5) / 3.0) + &
!                                 (( (0.5 * (gm - 1.) * ( x_speed(i, j, k) ** 2. + &
!                                                         y_speed(i, j, k) ** 2. + &
!                                                         z_speed(i, j, k) ** 2.) * &
!                                                        mass_residue(i, j, k)) + &
!                       (- (gm - 1.) * x_speed(i, j, k) * x_mom_residue(i, j, k)) + &
!                       (- (gm - 1.) * y_speed(i, j, k) * y_mom_residue(i, j, k)) + &
!                       (- (gm - 1.) * z_speed(i, j, k) * z_mom_residue(i, j, k)) + &
!                       ((gm - 1.) * energy_residue(i, j, k)) ) / 6.0) &
!                                ) * delta_t(i, j, k) / volume(i, j, k))
!            
!                density(i, j, k) = density_temp
!                x_speed(i, j, k) = x_speed_temp
!                y_speed(i, j, k) = y_speed_temp
!                z_speed(i, j, k) = z_speed_temp
!                pressure(i, j, k) = pressure_temp
!                include "turbulence_models/include/solver/RK4_update_solution.inc"
!              end do
!             end do
!            end do
!
!            if (any(density < 0) .or. any(pressure < 0)) then
!                call dmsg(5, 'solver', 'update_solution', &
!                        'ERROR: Some density or pressure is negative.')
!            end if
!
!        end subroutine RK4_update_solution
!
!        function get_residue_primitive() result(dEdx)
!
!            implicit none
!
!            real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, n_var) :: dEdx
!            real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1) :: beta
!            dEdx(:, :, :, 1) = mass_residue
!            dEdx(:, :, :, 2) = ( (-1 * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
!                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
!                                     mass_residue) + &
!                             ( x_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
!            dEdx(:, :, :, 3) = ( (-1 * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
!                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
!                                     mass_residue) + &
!                             ( y_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
!            dEdx(:, :, :, 4) = ( (-1 * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
!                                       density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
!                                     mass_residue) + &
!                             ( z_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )
!            dEdx(:, :, :, 5) = ( (0.5 * (gm - 1.) * ( x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
!                                                      y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
!                                                      z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2.) * &
!                                                    mass_residue) + &
!                       (- (gm - 1.) * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * x_mom_residue) + &
!                       (- (gm - 1.) * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * y_mom_residue) + &
!                       (- (gm - 1.) * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * z_mom_residue) + &
!                       ((gm - 1.) * energy_residue) )
!
!            include "turbulence_models/include/solver/get_residue_primitive.inc"
!
!        end function get_residue_primitive
!
!        subroutine update_solution()
!            !-----------------------------------------------------------
!            ! Update the solution using the residue and time step
!            !-----------------------------------------------------------
!
!            implicit none
!            integer :: i, j, k
!            real :: beta
!            
!            call dmsg(1, 'solver', 'update_solution')
!
!            do k = 1, kmx - 1
!             do j = 1, jmx - 1
!              do i = 1, imx - 1
!               density_temp = density(i, j, k) - &
!                            (mass_residue(i, j, k) * &
!                            delta_t(i, j, k) / volume(i, j, k))
!
!               x_speed_temp = x_speed(i, j, k) - &
!                            (( (-1 * x_speed(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( x_mom_residue(i, j, k) / density(i, j, k)) ) * &
!                            delta_t(i, j, k) / volume(i, j, k))
!
!               y_speed_temp = y_speed(i, j, k) - &
!                            (( (-1 * y_speed(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( y_mom_residue(i, j, k) / density(i, j, k)) ) * &
!                            delta_t(i, j, k) / volume(i, j, k))
!
!               z_speed_temp = z_speed(i, j, k) - &
!                            (( (-1 * z_speed(i, j, k) / density(i, j, k) * &
!                                     mass_residue(i, j, k)) + &
!                             ( z_mom_residue(i, j, k) / density(i, j, k)) ) * &
!                            delta_t(i, j, k) / volume(i, j, k))
!
!               pressure_temp = pressure(i, j, k) - &
!                   ( ( (0.5 * (gm - 1.) * ( x_speed(i, j, k) ** 2. + &
!                                            y_speed(i, j, k) ** 2. + &
!                                            z_speed(i, j, k) ** 2.) * &
!                                          mass_residue(i, j, k)) + &
!                       (- (gm - 1.) * x_speed(i, j, k) * x_mom_residue(i, j, k)) + &
!                       (- (gm - 1.) * y_speed(i, j, k) * y_mom_residue(i, j, k)) + &
!                       (- (gm - 1.) * z_speed(i, j, k) * z_mom_residue(i, j, k)) + &
!                       ((gm - 1.) * energy_residue(i, j, k)) ) * &
!                       delta_t(i, j, k) / volume(i, j, k) ) 
!
!               include "turbulence_models/include/solver/update_solution.inc"
!               density(i, j, k) = density_temp
!               x_speed(i, j, k) = x_speed_temp
!               y_speed(i, j, k) = y_speed_temp
!               z_speed(i, j, k) = z_speed_temp
!               pressure(i, j, k) = pressure_temp
!              end do
!             end do
!            end do
!
!            if (any(density < 0.) .or. any(pressure < 0.)) then
!                call dmsg(5, 'solver', 'update_solution', &
!                        'ERROR: Some density or pressure is negative.')
!                !stop
!            end if
!
!            do k = -2,kmx+2
!              do j = -2,jmx+2
!                do i = -2,imx+2
!                  if (density(i,j,k)<0.) then
!                    print*, process_id, i,j,k, "density: ", density(i,j,k)
!                  end if
!                  if (pressure(i,j,k)<0.) then
!                    print*, process_id, i,j,k, "pressure: ", pressure(i,j,k)
!                  end if
!                end do
!              end do
!            end do
!
!        end subroutine update_solution
!
!
!        subroutine get_total_conservative_Residue()
!
!            implicit none
!
!            call dmsg(1, 'solver', 'get_total_conservative_Residue')
!            call send_recv(3) ! parallel call-argument:no of layers 
!            call populate_ghost_primitive()
!            call compute_face_interpolant()
!            call reconstruct_boundary_state(interpolant)
!            call compute_fluxes()
!            if (mu_ref /= 0.0) then
!              call evaluate_all_gradients()
!              call calculate_transport()
!              call calculate_sst_F1()
!              call compute_viscous_fluxes(F_p, G_p, H_p)
!              call compute_turbulent_fluxes(F_p, G_p, H_p)
!            end if
!            call compute_residue()
!            call add_source_term_residue()
!            call dmsg(1, 'solver', 'step', 'Residue computed.')
!
!        end subroutine get_total_conservative_Residue
        
        subroutine iterate_one_more_time_step()
            !-----------------------------------------------------------
            ! Perform one time step iteration
            !
            ! This subroutine performs one iteration by stepping through
            ! time once.
            !-----------------------------------------------------------

            implicit none
            integer :: ierr
            call dmsg(1, 'solver', 'iterate_one_more_time_step')

            if (process_id==0) then
              print*, current_iter
            end if
!            call get_total_conservative_Residue()
!            call compute_time_step()
            call get_next_solution()
            if((mod(current_iter,res_write_interval)==0 .or. &
                    current_iter==Res_itr .or.               &
                    current_iter==1))      then
              call find_resnorm()
            end if
            call checkpoint()
            current_iter = current_iter + 1
            if(process_id==0)then
              REWIND(STOP_FILE_UNIT)
              read(STOP_FILE_UNIT,*) want_to_stop
            end if
            call MPI_BCAST(want_to_stop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            !if (want_to_stop==1) max_iters=current_iter-1
            if (want_to_stop==1) Halt = .TRUE.

        end subroutine iterate_one_more_time_step


end module solver

