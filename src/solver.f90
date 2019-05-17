module solver
  !< Setup, run, and destroy the solver
  !< allocate/deallcoate memory, initialize, iterate
  !-------------------------------------------------

  use global, only: STOP_FILE_UNIT
  use global, only: stop_file
  use global_vars, only : want_to_stop
  use global_vars, only : Halt

  use global_vars, only : max_iters
  use global_vars, only : current_iter
  use global_vars, only : checkpoint_iter
  use global_vars, only : checkpoint_iter_count
  use global_vars, only : sim_clock
  use global_vars, only : turbulence
  use global_vars, only : supersonic_flag

  use global_vars, only: res_write_interval
  use global_vars, only: r_list
  use global_vars, only: w_list
  use global_vars, only: Res_itr

  use utils, only:  dmsg

  use CC,    only: setupCC
  use CC,    only: destroyCC

  use string
  use read, only : read_input_and_controls

  use grid,      only : setup_grid, destroy_grid
  use geometry,  only : setup_geometry, destroy_geometry
  use state,     only : setup_state, destroy_state
  use gradients, only : setup_gradients
  use gradients, only : destroy_gradients
  use Scheme,    only : setup_scheme, destroy_scheme

  use source,        only: add_source_term_residue, setup_source, destroy_source
  use wall_dist,     only: setup_wall_dist, destroy_wall_dist, find_wall_dist
  use viscous,       only: compute_viscous_fluxes
  use layout,        only: process_id, grid_file_buf, bc_file, &
                           get_process_data, read_layout_file, total_process
  use interface1,    only : setup_interface
  use interface1,    only : destroy_interface
  use resnorm,       only : find_resnorm, setup_resnorm, destroy_resnorm
  use dump_solution, only : checkpoint
  use viscosity    , only : setup_viscosity
  use viscosity    , only : destroy_viscosity
  use viscosity    , only : calculate_viscosity
  use wall         , only : write_surfnode
  use bc,            only : setup_bc
  use bc,            only : destroy_bc
  use time ,         only : setup_time
  use time ,         only : destroy_time
  use time ,         only : compute_time_step
  use update,        only : get_next_solution
  use update,        only : setup_update
  use update,        only : destroy_update
  use mapping,       only : read_interface_map
  use bc_primitive,  only : populate_ghost_primitive
  use summon_grad_evaluation, only : evaluate_all_gradients
  use boundary_state_reconstruction, only: reconstruct_boundary_state
#include "error.inc"
#include "mpi.inc"
    private

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: iterate_one_more_time_step

    contains


        subroutine setup_solver()
          !< Call to allocate memoery and initialize domain
          !--------------------------------------------------
            
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
            call setup_interface()
            call setup_scheme()
            if(turbulence /= 'none') then
              call write_surfnode()
              call setup_wall_dist()
              call mpi_barrier(MPI_COMM_WORLD,ierr)
              call find_wall_dist()
            end if
            call setupCC()
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
          !< Call to different modules to deallocate memory
          !--------------------------------------------------

            implicit none
            
            call dmsg(1, 'solver', 'destroy_solver')

            if(process_id==0)then
              close(STOP_FILE_UNIT)
            end if
            call destroy_update()
            call destroy_viscosity()
            call destroy_gradients()
            call destroyCC()
            if(turbulence /= 'none') then
              call destroy_wall_dist()
            end if
            call destroy_scheme()
            call destroy_source()
            call destroy_state()
            call destroy_geometry()
            call destroy_grid()
            call destroy_resnorm()
            call destroy_interface()
            call destroy_time()
            call destroy_bc()

            if(allocated(r_list)) deallocate(r_list)
            if(allocated(w_list)) deallocate(w_list)

        end subroutine destroy_solver

        subroutine initmisc()
          !< Initilize miscellaneous variables
          !----------------------------------
            
            implicit none
            
            call dmsg(1, 'solver', 'initmisc')

            sim_clock = 0.
            current_iter = 0

        end subroutine initmisc

        
        subroutine iterate_one_more_time_step()
            !< Perform one time step iteration
            !  This subroutine performs one iteration by stepping through
            !  time once.
            !-----------------------------------------------------------

            implicit none
            integer :: ierr
            call dmsg(1, 'solver', 'iterate_one_more_time_step')

            if (process_id==0) then
              print*, current_iter
            end if
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

