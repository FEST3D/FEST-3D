  !< Create and destroy the solver setup
module create_and_erase
  !< Create and destroy the solver setup
  !-----------------------------------------------------------
  ! 170609  -Jatinder Pal Singh Sandhu
  ! AIM : 1)to setup, create, allocate memory, link pointer
  !         (everything that is required before first iteration
  !       2)free memory and free pointers
  !------------------------------------------------------------
  use global, only: STOP_FILE_UNIT
  use global, only: stop_file
    
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
  use global_vars, only : r_list
  use global_vars, only : w_list

  use utils      , only :   alloc
  use utils      , only : dealloc 
  use utils      , only : dmsg

  use string
  use read       , only : read_input_and_controls

  use grid       , only :   setup_grid
  use grid       , only : destroy_grid
  use geometry   , only :   setup_geometry
  use geometry   , only : destroy_geometry
  use state      , only :    setup_state
  use state      , only :  destroy_state
  use gradients  , only :   setup_gradients
  use gradients  , only : destroy_gradients
  use scheme     , only :   setup_scheme
  use scheme     , only : destroy_scheme 
  use source     , only : add_source_term_residue
  use wall_dist  , only :   setup_wall_dist
  use wall_dist  , only : destroy_wall_dist
  use wall_dist  , only :    find_wall_dist
  use layout     , only : process_id
  use layout     , only : grid_file_buf
  use layout     , only : bc_file
  use layout     , only : get_process_data
  use layout     , only : read_layout_file
  use layout     , only : total_process
  use parallel   , only : allocate_buffer_cells
  use resnorm_   , only : destroy_resnorm
  use resnorm_   , only :   setup_resnorm
  use transport  , only : setup_transport
  use transport  , only : destroy_transport
  use bc         , only : setup_bc
  use blending_function , only : setup_sst_F1
  use blending_function , only : destroy_sst_F1
  use wall       , only : write_surfnode
  use time       , only :   setup_time
  use time       , only : destroy_time

  private

  public :: setup_all
  public :: destroy_all

  contains

      subroutine setup_all()
        !< To setup, create, allocate memory, link pointer
        !<    (everything that is required before first iteration
          implicit none

          call dmsg(1, 'create_erase', 'setup_all')
          call get_process_data()           ! parallel calls
          call read_layout_file(process_id) ! reads layout file calls
          call read_input_and_controls()    ! all input config file are read
          call setup_grid(grid_file_buf)    ! read grid 
          call setup_geometry()             ! calculate geometric quantities (area, normal and volume)
          call setup_state()                ! allocate memroy and initialize state variable
          call setup_transport()            ! allocate memroy to viscosity
          call setup_gradients()            ! allocate memroy to gradients
          call setup_bc()                   ! set id and face_names array
          call allocate_memory()            
          call allocate_buffer_cells(3)     ! parallel buffers (MPI interafce communication)
          call setup_scheme()               ! face convective flux: memory and scheme
          if(turbulence /= 'none') then
            call write_surfnode()
            call setup_wall_dist()
            call find_wall_dist()
          end if
          call setup_sst_F1()
          call link_aliases_solver()
          call setup_resnorm()
          call initmisc()
          checkpoint_iter_count = 0
          call checkpoint()  ! Create an initial dump file
          call setup_time()
          call dmsg(5, 'create_erase', 'setup_all', 'Setup complete')

      end subroutine setup_all



      subroutine destroy_all()
          !< free memory and free pointers
          implicit none
          
          call dmsg(1, 'create_erase', 'destroy_all')
          call destroy_time()
          call destroy_transport()
          call destroy_gradients()
          call destroy_wall_dist()
          call destroy_scheme()
          call deallocate_misc()
          call unlink_aliases_solver()
          call destroy_state()
          call destroy_geometry()
          call destroy_grid()
          call destroy_resnorm()
          call destroy_sst_F1()

          if(allocated(r_list)) deallocate(r_list)
          if(allocated(w_list)) deallocate(w_list)
          call dmsg(5, 'create_erase', 'destroy_all', 'Memory_freed')

      end subroutine destroy_all

end module create_and_erase
