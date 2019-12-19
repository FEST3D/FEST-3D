module solver
  !< Setup, run, and destroy the solver
  !< allocate/deallcoate memory, initialize, iterate
  !-------------------------------------------------

  use vartypes
!  use global, only: STOP_FILE_UNIT
!  use global, only: stop_file
!  use global, only: mapfile
!  use global, only: periodicfile
  use global_vars, only : want_to_stop
  use global_vars, only : Halt
  use global_vars, only : sim_clock
  use global_vars, only: mu_t
  use global_vars, only: mu

  use CC,    only: setupCC

  use read, only : read_input_and_controls

  use grid,      only : setup_grid, destroy_grid
  use geometry,  only : setup_geometry!, destroy_geometry
  use state,     only : setup_state!, destroy_state
  use gradients, only : setup_gradients
  use gradients, only : evaluate_all_gradients
  use Scheme,    only : setup_scheme!, destroy_scheme

  use wall_dist,     only: setup_wall_dist, destroy_wall_dist, find_wall_dist
  use viscous,       only: compute_viscous_fluxes
  use layout,        only: process_id, get_process_data, read_layout_file, total_process
  use interface1,    only : setup_interface
  use interface1,    only : destroy_interface
  use resnorm,       only : find_resnorm, setup_resnorm!, destroy_resnorm
  use dump_solution, only : checkpoint
  use viscosity    , only : setup_viscosity
!  use viscosity    , only : destroy_viscosity
  use viscosity    , only : calculate_viscosity
  use wall         , only : write_surfnode
  use bc,            only : setup_bc
  use bc,            only : destroy_bc
  use time ,         only : setup_time
  use time ,         only : destroy_time
  use update,        only : get_next_solution
  use update,        only : setup_update
!  use update,        only : destroy_update
  use mapping,       only : read_interface_map
  use bc_primitive,  only : populate_ghost_primitive
  use boundary_state_reconstruction, only: reconstruct_boundary_state
#include "debug.h"
#include "error.h"
#include "mpi.inc"
    private

    type(extent) :: dims
    type(nodetype), dimension(:,:,:), allocatable :: nodes
    type(celltype), dimension(:,:,:), allocatable :: cells
    type(facetype), dimension(:,:,:), allocatable :: Ifaces, Jfaces, Kfaces
    type(controltype), public :: control
    type(schemetype), public :: scheme
    type(flowtype) :: flow
    type(filetype) :: files

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

            DebugCall('setup_solver: Start')
            call get_process_data() ! parallel calls
            call read_layout_file(files, process_id) ! reads layout file calls
            
            call read_input_and_controls(files, control, scheme, flow)
            call setup_grid(files, nodes, dims)
            call setup_geometry(cells, Ifaces, Jfaces, Kfaces, nodes, dims)
            call setup_viscosity(mu, mu_t, scheme, flow, dims)
            call setup_state(files, control, scheme, flow, dims)
            call setup_gradients(control,scheme,flow,dims)
            !call setup_source
            call setup_bc(files, scheme, flow, dims)
            call setup_time(control,dims)
            call setup_update(control,scheme,flow, dims)
            call setup_interface(control,dims)
            call setup_scheme(control, scheme, dims)
            if(scheme%turbulence /= 'none') then
              call write_surfnode(files, nodes, dims)
              call setup_wall_dist(files, dims)
              call mpi_barrier(MPI_COMM_WORLD,ierr)
              call find_wall_dist(nodes, dims)
            end if
            call setupCC(scheme, dims)
            call setup_resnorm(files, control, scheme, flow)
            call initmisc()
            control%checkpoint_iter_count = 0
            call checkpoint(files, nodes, control, scheme, dims)  ! Create an initial dump file
            control%current_iter=1
            DebugCall('setup_solver: checkpoint')
            if(process_id==0)then
              open(files%STOP_FILE_UNIT, file=files%stop_file)
            end if
            DebugCall('Setup solver complete')

        end subroutine setup_solver

        subroutine destroy_solver()
          !< Call to different modules to deallocate memory
          !--------------------------------------------------

            implicit none
            
            DebugCall('destroy_solver')
!
!            if(process_id==0)then
!              close(STOP_FILE_UNIT)
!            end if
!            call destroy_update()
!            call destroy_viscosity()
!            !call destroy_gradients()
!            call destroyCC()
!            if(turbulence /= 'none') then
!              call destroy_wall_dist()
!            end if
!            call destroy_scheme()
!            !call destroy_source()
!            call destroy_state()
!            !call destroy_geometry()
!            !call destroy_grid()
!            call destroy_resnorm()
!            call destroy_interface()
!            call destroy_time()
!            call destroy_bc()
!
!            if(allocated(r_list)) deallocate(r_list)
!            if(allocated(w_list)) deallocate(w_list)

        end subroutine destroy_solver

        subroutine initmisc()
          !< Initilize miscellaneous variables
          !----------------------------------
            
            implicit none
            
            DebugCall('initmisc')

            sim_clock = 0.
            control%current_iter = 0

        end subroutine initmisc

        
        subroutine iterate_one_more_time_step()
            !< Perform one time step iteration
            !  This subroutine performs one iteration by stepping through
            !  time once.
            !-----------------------------------------------------------

            implicit none
            integer :: ierr
            DebugCall('iterate_one_more_time_step')

            if (process_id==0) then
              print*, control%current_iter
            end if
            call get_next_solution(control, scheme, flow, dims)
            !if((mod(control%current_iter,control%res_write_interval)==0 .or. &
            !        control%current_iter==Res_itr .or.               &
            !        control%current_iter==1))      then
              call find_resnorm(files%RESNORM_FILE_UNIT, control, scheme, dims)
            !end if
            call checkpoint(files, nodes, control, scheme, dims)
            control%current_iter = control%current_iter + 1
            if(process_id==0)then
              REWIND(files%STOP_FILE_UNIT)
              read(files%STOP_FILE_UNIT,*) want_to_stop
            end if
            call MPI_BCAST(want_to_stop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            !if (want_to_stop==1) max_iters=current_iter-1
            if (want_to_stop==1) Halt = .TRUE.

        end subroutine iterate_one_more_time_step


end module solver

