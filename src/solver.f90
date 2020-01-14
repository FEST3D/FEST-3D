module solver
  !< Setup, run, and destroy the solver
  !< allocate/deallcoate memory, initialize, iterate
  !-------------------------------------------------
  use vartypes
  use mpi
  use CC,            only : setupCC
  use read,          only : read_input_and_controls
  use grid,          only : setup_grid
  use geometry,      only : setup_geometry
  use state,         only : setup_state
  use gradients,     only : setup_gradients
  use Scheme,        only : setup_scheme
  use wall_dist,     only : setup_wall_dist, find_wall_dist
  use viscous,       only : compute_viscous_fluxes
  use layout,        only : get_process_data, read_layout_file
  use interface1,    only : setup_interface
  use resnorm,       only : find_resnorm, setup_resnorm!, destroy_resnorm
  use dump_solution, only : checkpoint
  use viscosity    , only : setup_viscosity
  use viscosity    , only : calculate_viscosity
  use wall         , only : write_surfnode
  use bc,            only : setup_bc
  use time ,         only : setup_time
  use time ,         only : destroy_time
  use update,        only : get_next_solution
  use update,        only : setup_update
#include "debug.h"
#include "error.h"
    private

    type(extent) :: dims
    !< Extent of the domain:imx,jmx,kmx
    type(nodetype), dimension(:,:,:), allocatable :: nodes
    !< Grid points 
    type(celltype), dimension(:,:,:), allocatable :: cells
    !< Cell center quantities: volume, cellCenter
    type(facetype), dimension(:,:,:), allocatable :: Ifaces, Jfaces, Kfaces
    !< Face quantities: area and unit normal
    type(controltype), public :: control
    !< Control parameters
    type(schemetype), public :: schemes
    !< finite-volume Schemes
    type(flowtype) :: flow
    !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
    type(filetype) :: files
    !< Files' name and handler
    type(boundarytype) :: boundary
    !< boundary conditions and fixed values
    real(wp), dimension(:, :, :, :), allocatable :: qp           
    !< Store primitive variable at cell center
    real(wp), dimension(:, :, :   ), allocatable :: Temp
    !< Store Temperature variable at cell center
    real(wp), public, dimension(:, :, :, :), allocatable, target :: F
    !< Store fluxes throught the I faces
    real(wp), public, dimension(:, :, :, :), allocatable, target :: G
    !< Store fluxes throught the J faces
    real(wp), public, dimension(:, :, :, :), allocatable, target :: H
    !< Store fluxes throught the K faces
    real(wp), public, dimension(:, :, :, :), allocatable, target :: residue
    !< Store residue at each cell-center
    real(wp), dimension(:, :, :), allocatable     :: delta_t  
    !< Local time increment value at each cell center

    ! Public methods
    public :: setup_solver
    public :: destroy_solver
    public :: iterate_one_more_time_step
    public :: abort_run
    public :: finish_run
    public :: start_run

    contains

        subroutine abort_run()
          !< Aborting the solver
          implicit none

          call finish_run()
          stop

        end subroutine abort_run

        subroutine finish_run()
          !< Finishing the solution computation
          implicit none
          integer :: ierr

          call destroy_solver()
          call MPI_FINALIZE(ierr)

        end subroutine finish_run

        subroutine start_run()
          !< Starting the solver setup
          implicit none
          integer :: ierr

          call MPI_INIT(ierr)
          call setup_solver()

        end subroutine start_run

        subroutine setup_solver()
          !< Call to allocate memoery and initialize domain
          !--------------------------------------------------
            
            implicit none
            integer :: ierr

            DebugCall('setup_solver: Start')
            call get_process_data(control) ! parallel calls
            call read_layout_file(files, control, boundary) ! reads layout file calls
            
            call read_input_and_controls(files, control, schemes, flow)
            call setup_grid(files, nodes, control, boundary, dims)
            call setup_geometry(cells, Ifaces, Jfaces, Kfaces, nodes, boundary, dims)
            !call setup_viscosity(mu, mu_t, schemes, flow, dims)
            call setup_viscosity(schemes, flow, dims)
            call setup_state(files, qp, control, schemes, flow, dims)
            allocate(Temp(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2))
            call setup_gradients(control,schemes,flow,dims)
            !call setup_source
            call setup_bc(files, schemes, flow, boundary, dims)
            call setup_time(delta_t, control,dims)
            call setup_update(control,schemes,flow, dims)
            call setup_interface(control,dims)
            call setup_scheme(residue, F,G,H, control, dims)
            if(schemes%turbulence /= 'none') then
              call write_surfnode(files, nodes, control, boundary, dims)
              call setup_wall_dist(files, dims)
              call mpi_barrier(MPI_COMM_WORLD,ierr)
              call find_wall_dist(nodes, dims)
            end if
            call setupCC(schemes, cells, Ifaces,Jfaces,Kfaces, dims)
            call setup_resnorm(files, control, schemes, flow)
            call initmisc()
            control%checkpoint_iter_count = 0
            call checkpoint(files, qp, nodes, control, schemes, dims)  ! Create an initial dump file
            control%current_iter=1
            DebugCall('setup_solver: checkpoint')
            DebugCall('Setup solver complete')

        end subroutine setup_solver

        subroutine destroy_solver()
          !< Call to different modules to deallocate memory
          !--------------------------------------------------

            implicit none
            
            DebugCall('destroy_solver')
            call destroy_time(control)

        end subroutine destroy_solver

        subroutine initmisc()
          !< Initilize miscellaneous variables
          !----------------------------------
            
            implicit none
            
            DebugCall('initmisc')

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

            call get_next_solution(qp, Temp, residue, delta_t, cells, F,G,H, Ifaces,Jfaces,Kfaces,control, schemes, flow, boundary, dims)
            call find_resnorm(files%RESNORM_FILE_UNIT, residue, F,G,H, control, schemes, dims)
            call checkpoint(files, qp, nodes, control, schemes, dims)
            control%current_iter = control%current_iter + 1
            if(process_id==0)then
              open(files%STOP_FILE_UNIT, file=files%stop_file)
              read(files%STOP_FILE_UNIT,*) control%want_to_stop
              close(files%STOP_FILE_UNIT)
            end if
            call MPI_BCAST(control%want_to_stop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (control%want_to_stop==1) control%Halt = .TRUE.

        end subroutine iterate_one_more_time_step


end module solver

