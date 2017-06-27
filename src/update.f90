module update
  use global, only: CONFIG_FILE_UNIT 
  use global, only: RESNORM_FILE_UNIT 
  use global, only: FILE_NAME_LENGTH
  use global, only: STRING_BUFFER_LENGTH 
  use global, only: INTERPOLANT_NAME_LENGTH
  use global, only: STOP_FILE_UNIT
  use global, only: stop_file
  use global_vars, only : want_to_stop

  use global_sst , only : beta1
  use global_sst , only : beta2
  use global_sst , only : bstar
  use global_sst , only : sst_F1
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
  use global_vars, only: rw_list
  use global_vars, only: merror

  use utils, only: alloc
  use utils, only:  dealloc 
  use utils, only:  dmsg
  use utils, only:  DEBUG_LEVEL

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
  use source, only: add_source_term_residue
  use wall_dist, only: setup_wall_dist, destroy_wall_dist, find_wall_dist
  use viscous, only: compute_viscous_fluxes
  use turbulent_flux, only: compute_turbulent_fluxes
  use boundary_state_reconstruction, only: reconstruct_boundary_state
  use layout, only: process_id, grid_file_buf, bc_file, &
  get_process_data, read_layout_file, total_process
  use parallel, only: allocate_buffer_cells,send_recv
!  use state, only: turbulence
  use resnorm_, only : write_resnorm, setup_resnorm, destroy_resnorm
  use dump_solution, only : checkpoint
  use transport    , only : setup_transport
  use transport    , only : destroy_transport
  use transport    , only : calculate_transport
  use blending_function , only : setup_sst_F1
  use blending_function , only : destroy_sst_F1
  use blending_function , only : calculate_sst_F1
  use wall        , only : write_surfnode
  include "turbulence_models/include/solver/import_module.inc"
  use bc, only: setup_bc
  use bc_primitive, only: populate_ghost_primitive
  use summon_grad_evaluation, only : evaluate_all_gradients

#ifdef __GFORTRAN__
    use mpi
#endif    
    implicit none
#ifdef __INTEL_COMPILER
    include "mpif.h"
#endif
    private

    real, dimension(:), allocatable, target :: qp_temp
    real, pointer :: density_temp, x_speed_temp, &
                                         y_speed_temp, z_speed_temp, &
                                         pressure_temp
    real, dimension(:,:,:,:), allocatable  :: pR !primitve residue 
    include "turbulence_models/include/solver/variables_deceleration.inc"

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: step
    public :: converged

    contains


      subroutine setup_update()
        implicit none

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
      end subroutine setup_update


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

            call dealloc(qc_n)
            call dealloc(R1)

        end subroutine destroy_RK4_time_step

        subroutine setup_RK4_time_step()
    
            implicit none

            call alloc(qc_n, -2, imx+2, -2, jmx+2, -2, kmx+2, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for qc_n.')
            call alloc(R1, 1, imx-1, 1, jmx-1, 1, kmx-1, 1, n_var, &
                    errmsg='Error: Unable to allocate memory for R1.')

        end subroutine setup_RK4_time_step


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

        function pR()
          implicit none
          real, dimension(1:imx-1, 1:jmx-1, 1:kmx-1, n_var) :: pR

          pR(:, :, :, 1) = mass_residue

          pR(:, :, :, 2) = ((-1 * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                               mass_residue) + &
                       (x_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)))

          pR(:, :, :, 3) = ( (-1 * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                 density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                               mass_residue) + &
                       ( y_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )

          pR(:, :, :, 4) = ( (-1 * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) / &
                                 density(1:imx-1, 1:jmx-1, 1:kmx-1) * &
                               mass_residue) + &
                       ( z_mom_residue / density(1:imx-1, 1:jmx-1, 1:kmx-1)) )

          pR(:, :, :, 5) = ((0.5*(gm - 1.)*(x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
                            y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2. + &
                            z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) ** 2.) * &
                                                  mass_residue) + &
                     (- (gm - 1.) * x_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * x_mom_residue) + &
                     (- (gm - 1.) * y_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * y_mom_residue) + &
                     (- (gm - 1.) * z_speed(1:imx-1, 1:jmx-1, 1:kmx-1) * z_mom_residue) + &
                     ((gm - 1.) * energy_residue) )
          
        end function pR

        subroutine RK3()
          implicit none

            qp_n = qp(1:imx-1, 1:jmx-1, 1:kmx-1,:)
            delta_t = delta_t_0
            call update_solution()
            !qp(1:imx-1,1:jmx-1,1:kmx-1,:) = qp(1:imx-1, 1:jmx-1, 1:kmx-1,:)
            call sub_step()
            delta_t = delta_t_0
            call update_solution()
            qp(1:imx-1,1:jmx-1,1:kmx-1,:) = 0.75*qp_n + 0.25*(qp(1:imx-1,1:jmx-1,1:kmx-1,:))
            call sub_step()
            delta_t = delta_t_0
            call update_solution()
            qp(1:imx-1,1:jmx-1,1:kmx-1,:) = (0.33334)*qp_n + (0.66666)*(qp(1:imx-1,1:jmx-1,1:kmx-1,:))
        end subroutine RK3


        subroutine update_solution(time_factor)
            !-----------------------------------------------------------
            ! Update the solution using the residue and time step
            !-----------------------------------------------------------

            implicit none
            integer, intent(in) :: time_factor
            integer :: i, j, k
            real :: beta
            real :: LS    ! length scale
            
            call dmsg(1, 'solver', 'update_solution')

            do k = 1, kmx - 1
             do j = 1, jmx - 1
              do i = 1, imx - 1
                LS  = time_factor*delta(i,j,k)/volume(i,j,k)
                density_temp = density(i, j, k) - LS*mass_residue(i, j, k)

               x_speed_temp = x_speed(i, j, k) - LS*&
                             ((-1 * x_speed(i, j, k)/density(i, j, k)*mass_residue(i, j, k)) + &
                             ( x_mom_residue(i, j, k) / density(i, j, k)) )

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

               include "turbulence_models/include/solver/update_solution.inc"
               density(i, j, k) = density_temp
               x_speed(i, j, k) = x_speed_temp
               y_speed(i, j, k) = y_speed_temp
               z_speed(i, j, k) = z_speed_temp
               pressure(i, j, k) = pressure_temp
              end do
             end do
            end do

            if (any(density < 0.) .or. any(pressure < 0.)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
                !stop
            end if

            do k = -2,kmx+2
              do j = -2,jmx+2
                do i = -2,imx+2
                  if (density(i,j,k)<0.) then
                    print*, process_id, i,j,k, "density: ", density(i,j,k)
                  end if
                  if (pressure(i,j,k)<0.) then
                    print*, process_id, i,j,k, "pressure: ", pressure(i,j,k)
                  end if
                end do
              end do
            end do

        end subroutine update_solution


        subroutine RK4_update_solution()

            implicit none
            integer :: i, j, k

            ! qp at various stages is not stored but over written
            ! The residue multiplied by the inverse of the jacobian
            ! is stored for the final update equation

            ! Stage 1 is identical to stage (n)
            ! Store qp(n)
            call compute_time_step()
            qc_n = qp(1:imx-1, 1:jmx-1, 1:kmx-1, 1:n_var)
            dEdx_1 = get_residue_primitive()
            
            ! Stage 2
            ! Not computing delta_t since qp(1) = qp(n)
            ! Update solution will over write qp
            delta_t = 0.5 * delta_t_0  ! delta_t(1)
            call update_solution()

            ! Stage 3
            call sub_step()
            dEdx_2 = get_residue_primitive()
            delta_t = 0.5 * delta_t_0
            call update_solution()

            ! Stage 4
            call sub_step()
            dEdx_3 = get_residue_primitive()
            delta_t = delta_t_0
            call update_solution()

            ! qp now is qp_4
            ! Use qp(4)
            call sub_step()

              end do
             end do
            end do

            if (any(density < 0) .or. any(pressure < 0)) then
                call dmsg(5, 'solver', 'update_solution', &
                        'ERROR: Some density or pressure is negative.')
            end if

        end subroutine RK4_update_solution


        subroutine get_primitive_from_conservative()
          implicit none
          integer :: i,j,k
          real    :: r,u,v,w,k,p

          qp(:,:,:,1)  = qc(:,:,:,1)               ! density
          qp(:,:,:,2:) = qc(:,:,:,2:)/qc(:,:,:,1)  ! all variables except density

          ! get turbulent kinetic energy
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = 1,imx-1

                select case (trim(turbulence))
                  case('sst', 'kkl')
                    k = qp(i,j,k,6)
                  case DEFAULT
                    k = 0.
                end select

                r = qp(i,j,k,1) !density
                u = qp(i,j,k,2) !x_speed
                v = qp(i,j,k,3) !y_speed
                w = qp(i,j,k,4) !z_speed

                p = (gm-1.)*(qc(i,j,k,5) - (0.5*r*(u**2+v**2+w**2)) - r*k )
                qp(i,j,k,5) = p !pressure

              end do
            end do
          end do

        end subroutine get_primitive_from_conservative


        subroutine get_conservative_from_primitive()
          implicit none
          integer :: i,j,k
          real    :: r,u,v,w,k,p

          qc(:,:,:,1)  = qp(:,:,:,1)               ! density
          qc(:,:,:,2:) = qp(:,:,:,2:)*qp(:,:,:,1)  ! all variables except density

          ! get turbulent kinetic energy
          do k = 1,kmx-1
            do j = 1,jmx-1
              do i = 1,imx-1

                select case (trim(turbulence))
                  case('sst', 'kkl')
                    k = qp(i,j,k,6)
                  case DEFAULT
                    k = 0.
                end select

                r = qp(i,j,k,1) !density
                u = qp(i,j,k,2) !x_speed
                v = qp(i,j,k,3) !y_speed
                w = qp(i,j,k,4) !z_speed
                p = qp(i,j,k,5) !pressure

                qc(i,j,k,5) = p/(gm-1.) + (0.5*r*(u**2+v**2+w**2)) + r*k

              end do
            end do
          end do
        end subroutine get_conservative_from_primitive




end module update
