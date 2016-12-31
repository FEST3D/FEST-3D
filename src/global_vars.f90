module global_vars
  !----------------------------------------------
  ! contains all the public/global variables used 
  ! by more than one module
  !----------------------------------------------

  use global, only : INTERPOLANT_NAME_LENGTH
  use global, only : FORMAT_LENGTH
  use global, only : SCHEME_NAME_LENGTH
  use global, only : FILE_NAME_LENGTH

  implicit none
  public

  ! Parallel processig variables
  integer :: total_process
  integer :: total_entries
  integer :: process_id
  integer :: imin_id
  integer :: imax_id
  integer :: jmin_id
  integer :: jmax_id
  integer :: kmin_id
  integer :: kmax_id

  ! Time controls
  integer :: min_iter=1     !Minimum iteration value, starting iteration value
  integer :: max_iters       !Maximum iteration value, stop at
  integer :: start_from=0   ! folder load level from time_directories eg 1 for folder 0001
  integer :: checkpoint_iter
  integer :: checkpoint_iter_count
  integer :: current_iter

  !write controls
!  integer :: write_interval
  integer :: res_write_interval
  integer :: purge_write
  integer :: write_percision
  character(len=FORMAT_LENGTH):: write_data_format   ! either ascii or binary
  character(len=FORMAT_LENGTH):: write_file_format   ! either vtk or tecplot
  character(len=FILE_NAME_LENGTH):: outfile   ! either vtk or tecplot

  ! solver specific
  real :: CFL
  real :: tolerance

  !solver time secific
  character                                 :: time_stepping_method
  character(len=INTERPOLANT_NAME_LENGTH)    :: time_step_accuracy
  real                                      :: global_time_step
  real, dimension(:, :, :), allocatable     :: delta_t
  real                                      :: sim_clock

  !scheme
  character(len=SCHEME_NAME_LENGTH) :: scheme_name
  character(len=INTERPOLANT_NAME_LENGTH) :: interpolant

  !solution specific
  real, dimension(:, :, :, :), allocatable  :: qp_n
  real, dimension(:, :, :, :), allocatable  :: dEdx_1
  real, dimension(:, :, :, :), allocatable  :: dEdx_2
  real, dimension(:, :, :, :), allocatable  :: dEdx_3

  ! state variables viscous
  integer                                           :: n_var
  real, dimension(:, :, :, :), allocatable, target  :: qp
  real, dimension(:)         , allocatable, target  :: qp_inf
  real, dimension(:, :, :)                , pointer :: density
  real, dimension(:, :, :)                , pointer :: x_speed
  real, dimension(:, :, :)                , pointer :: y_speed
  real, dimension(:, :, :)                , pointer :: z_speed
  real, dimension(:, :, :)                , pointer :: pressure
  real                                    , pointer :: density_inf
  real                                    , pointer :: x_speed_inf
  real                                    , pointer :: y_speed_inf
  real                                    , pointer :: z_speed_inf
  real                                    , pointer :: pressure_inf

  ! Freestram variable used to read file before inf pointer are linked and allocated
  real                                              :: free_stream_density
  real                                              :: free_stream_x_speed
  real                                              :: free_stream_y_speed
  real                                              :: free_stream_z_speed
  real                                              :: free_stream_pressure
  real                                              :: free_stream_tk
  real                                              :: free_stream_tw
  real, dimension(:, :, :), allocatable             :: dist

  ! state variable turbulent
  integer                                           :: sst_n_var=2
!  real, dimension(:, :, :, :), allocatable, target  :: tqp       ! turbulent primitive
!  real, dimension(:)         , allocatable, target  :: tqp_inf   ! turbulent primitive at inf
  real, dimension(:, :, :)                , pointer :: tk        ! TKE/mass
  real, dimension(:, :, :)                , pointer :: tw        ! omega
  real                                    , pointer :: tk_inf
  real                                    , pointer :: tw_inf

  ! residue variables
  real, dimension(:, :, :, :)             , pointer :: F_p
  real, dimension(:, :, :, :)             , pointer :: G_p
  real, dimension(:, :, :, :)             , pointer :: H_p
  real, dimension(:, :, :)                , pointer :: mass_residue
  real, dimension(:, :, :)                , pointer :: x_mom_residue
  real, dimension(:, :, :)                , pointer :: y_mom_residue
  real, dimension(:, :, :)                , pointer :: z_mom_residue
  real, dimension(:, :, :)                , pointer :: energy_residue
  real, dimension(:, :, :)                , pointer :: TKE_residue
  real, dimension(:, :, :)                , pointer :: omega_residue

  ! thermal properties
  real                                              :: gm    !gamma
  real                                              :: R_gas !univarsal gas constant
  
  ! Transport properties
  real                                              :: mu_ref !viscoity

  ! sutherland law variable
  real                                              :: T_ref
  real                                              :: Sutherland_temp

  ! nondimensional numbers
  real                                              :: Pr !prandtl number

  ! switches
  logical                                           :: supersonic_flag
  integer                                           :: ilimiter_switch
  integer                                           :: PB_switch
  character(len=5)                                  :: turbulence ! todo character length
  

  !residual specific
  real, pointer ::        resnorm
  real, pointer ::    vis_resnorm
  real, pointer ::   turb_resnorm
  real, pointer ::   cont_resnorm     
  real, pointer ::  x_mom_resnorm    
  real, pointer ::  y_mom_resnorm    
  real, pointer ::  z_mom_resnorm    
  real, pointer :: energy_resnorm   
  real, pointer ::    TKE_resnorm    
  real, pointer ::  omega_resnorm    
  real, pointer ::        resnorm_d1
  real, pointer ::    vis_resnorm_d1
  real, pointer ::   turb_resnorm_d1 
  real, pointer ::   cont_resnorm_d1
  real, pointer ::  x_mom_resnorm_d1
  real, pointer ::  y_mom_resnorm_d1
  real, pointer ::  z_mom_resnorm_d1
  real, pointer :: energy_resnorm_d1
  real, pointer ::    TKE_resnorm_d1
  real, pointer ::  omega_resnorm_d1
  real          ::        resnorm_0
  real          ::    vis_resnorm_0
  real          ::   turb_resnorm_0 
  real          ::   cont_resnorm_0
  real          ::  x_mom_resnorm_0
  real          ::  y_mom_resnorm_0
  real          ::  z_mom_resnorm_0
  real          :: energy_resnorm_0
  real          ::    TKE_resnorm_0
  real          ::  omega_resnorm_0
  real          ::        resnorm_0s
  real          ::    vis_resnorm_0s
  real          ::   turb_resnorm_0s 
  real          ::   cont_resnorm_0s
  real          ::  x_mom_resnorm_0s
  real          ::  y_mom_resnorm_0s
  real          ::  z_mom_resnorm_0s
  real          :: energy_resnorm_0s
  real          ::    TKE_resnorm_0s
  real          ::  omega_resnorm_0s

  ! grid variables
  integer                                 :: imx, jmx, kmx        ! no. of points
  real, dimension(:, :, :), allocatable  :: grid_x, grid_y, grid_z! point coordinates

  ! geometry variables
  real, dimension(:, :, :), allocatable, target :: xnx, xny, xnz !face unit norm x
  real, dimension(:, :, :), allocatable, target :: ynx, yny, ynz !face unit norm y
  real, dimension(:, :, :), allocatable, target :: znx, zny, znz !face unit norm z
  real, dimension(:, :, :), allocatable, target :: xA, yA, zA    !face area
  real, dimension(:, :, :), allocatable, target :: volume
  real, dimension(:, :, :), allocatable, target ::   left_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::  right_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::  front_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::   back_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::    top_ghost_centroid
  real, dimension(:, :, :), allocatable, target :: bottom_ghost_centroid

end module global_vars

