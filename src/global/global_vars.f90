module global_vars
  !----------------------------------------------
  ! contains all the public/global variables used 
  ! by more than one module
  !----------------------------------------------

  use global, only : INTERPOLANT_NAME_LENGTH
  use global, only : FORMAT_LENGTH
  use global, only : SCHEME_NAME_LENGTH
  use global, only : FILE_NAME_LENGTH
  use global, only : STATE_NAME_LENGTH
  use global, only : FLOW_TYPE_LENGTH

  implicit none
  public

  ! Parallel processig variables
  integer :: total_process      ! total no. of process to be used for computation
  integer :: total_entries      ! total enteries in layout.md for each processor
  integer :: process_id         ! id no. of each processor
  integer :: imin_id            ! bc_list at imin for particulat processor
  integer :: imax_id            ! bc_list at imax for particulat processor
  integer :: jmin_id            ! bc_list at jmin for particulat processor
  integer :: jmax_id            ! bc_list at jmax for particulat processor
  integer :: kmin_id            ! bc_list at kmin for particulat processor
  integer :: kmax_id            ! bc_list at kmax for particulat processor

  ! Time controls
  integer :: min_iter=1         !Minimum iteration value, starting iteration value
  integer :: max_iters          !Maximum iteration value, stop at
  integer :: start_from=0       ! folder load level from time_directories eg 1 for folder 0001
  integer :: checkpoint_iter    ! Write interval for output file
  integer :: checkpoint_iter_count ! write file counter
  integer :: current_iter       ! current iteration number

  !write controls
  integer :: r_count=0
  integer :: w_count=0
  integer :: res_write_interval                      ! resnorm write interval
  integer :: purge_write                             ! number of output files per process to keep
  integer :: write_percision                         ! number of place after decimal 
  character(len=FORMAT_LENGTH):: write_data_format   ! either ascii or binary
  character(len=FORMAT_LENGTH):: write_file_format   ! either vtk or tecplot
  character(len=FORMAT_LENGTH)::  read_data_format='ASCII'   ! either ascii or binary
  character(len=FORMAT_LENGTH)::  read_file_format="vtk"   ! either vtk or tecplot
  character(len=FILE_NAME_LENGTH)::     outfile          ! name of output file
  character(len=FILE_NAME_LENGTH)::      infile          ! name of load file
  character(len=FILE_NAME_LENGTH):: restartfile          ! name of restart log file
  character(len=STATE_NAME_LENGTH), dimension(:), allocatable ::  r_list ! read  control list
  character(len=STATE_NAME_LENGTH), dimension(:), allocatable ::  w_list ! write control list
  character(len=FLOW_TYPE_LENGTH):: previous_flow_type="none"

  ! solver specific
  real :: CFL                  ! courrent number
  real :: tolerance            ! minimum value of resnorm after which simulation stop
  integer :: want_to_stop=0

  !solver time secific
  character                                 :: time_stepping_method  ! local or global
  character(len=INTERPOLANT_NAME_LENGTH)    :: time_step_accuracy    ! Ranga-Kutta 4th order or first order
  real                                      :: global_time_step      ! value of global time step
  real, dimension(:, :, :), allocatable     :: delta_t               ! time increment value
  real                                      :: sim_clock             

  !scheme
  character(len=SCHEME_NAME_LENGTH) :: scheme_name         !flux calculation -> ausm, ldfss0, vanleer
  character(len=INTERPOLANT_NAME_LENGTH) :: interpolant    !face state reconstruction -> muscl, ppm, none

  !solution specific (used for Ranga_kutta_4th order)
  real, dimension(:, :, :, :), allocatable  :: qp_n
  real, dimension(:, :, :, :), allocatable  :: dEdx_1
  real, dimension(:, :, :, :), allocatable  :: dEdx_2
  real, dimension(:, :, :, :), allocatable  :: dEdx_3

  ! state variables viscous
   ! number of variable to solve for
   ! primitive variable at cell center
   ! primitive variable at infinity  
   ! rho pointer, point to slice of qp
   ! u pointer, point to slice of qp 
   ! v pointer, point to slice of qp 
   ! w pointer, point to slice of qp 
   ! P pointer, point to slice of qp 
   ! rho pointer, point to slice of qp_inf
   ! u pointer, point to slice of qp_inf
   ! v pointer, point to slice of qp_inf
   ! w pointer, point to slice of qp_inf
   ! p pointer, point to slice of qp_inf
  integer                                           :: n_var=5
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
   ! read Rho_inf from control file
   ! read U_inf from control file
   ! read V_inf from control file
   ! read W_inf from control file
   ! read P_inf from control file
   ! read tk_inf from control file
   ! read tw_inf from control file
   ! wall distance for each cell center
  real                                              :: free_stream_density  
  real                                              :: free_stream_x_speed  
  real                                              :: free_stream_y_speed  
  real                                              :: free_stream_z_speed  
  real                                              :: free_stream_pressure 
  real                                              :: free_stream_tk       
  real                                              :: free_stream_tw       
  real                                              :: free_stream_te
  real                                              :: free_stream_tv
  real                                              :: free_stream_tkl
  real                                              :: vel_mag ! free_stream velocity magnitude
  real, dimension(:, :, :), allocatable             :: dist 

  ! state variable turbulent
  integer                                           :: sst_n_var=2
  integer                                           :: sst_n_grad=2
!  real, dimension(:, :, :, :), allocatable, target  :: tqp       ! turbulent primitive
!  real, dimension(:)         , allocatable, target  :: tqp_inf   ! turbulent primitive at inf
  real, dimension(:, :, :)                , pointer :: tk        ! TKE/mass
  real, dimension(:, :, :)                , pointer :: tw        ! omega
  real, dimension(:, :, :)                , pointer :: te        ! Dissipation
  real, dimension(:, :, :)                , pointer :: tv        ! sa visocity
  real, dimension(:, :, :)                , pointer :: tkl       ! KL K-KL method
  real                                    , pointer :: tk_inf    ! TKE/mass at inf
  real                                    , pointer :: tw_inf    ! omega at inf
  real                                    , pointer :: te_inf    ! dissipation at inf
  real                                    , pointer :: tv_inf    ! SA viscosity at inf
  real                                    , pointer :: tkl_inf   ! kl at inf

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
  real, dimension(:, :, :)                , pointer :: KL_residue
  real, dimension(:, :, :)                , pointer :: dissipation_residue
  real, dimension(:, :, :)                , pointer :: tv_residue

  ! thermal properties
  real                                              :: gm    !gamma
  real                                              :: R_gas !univarsal gas constant
  
  ! Transport properties
  real                                              :: mu_ref !viscoity
  character(len=FILE_NAME_LENGTH)                   :: mu_variation

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
  
  real, dimension(:, :, :), allocatable             :: mu
  real, dimension(:, :, :), allocatable, target     :: mu_t
  real, dimension(:, :, :)              , pointer   :: sst_mu
  real, dimension(:, :, :)              , pointer   :: kkl_mu

  !residual specific
  real, pointer ::        resnorm     !            residue normalized
  real, pointer ::    vis_resnorm     ! {rho+V+P}  residue normalized
  real, pointer ::   turb_resnorm     !  turbulent residue normalized
  real, pointer ::   cont_resnorm     !       mass residue normalized
  real, pointer ::  x_mom_resnorm     ! x momentum residue normalized
  real, pointer ::  y_mom_resnorm     ! Y momentum residue normalized
  real, pointer ::  z_mom_resnorm     ! Z momentum residue normalized
  real, pointer :: energy_resnorm     !     energy residue normalized
  real, pointer ::    TKE_resnorm     !        TKE residue normalized
  real, pointer ::  omega_resnorm     !      omega residue normalized
  real, pointer ::        resnorm_d1  !            residue normalized/same at iter 1
  real, pointer ::    vis_resnorm_d1  ! {rho+V+P}  residue normalized/same at iter 1
  real, pointer ::   turb_resnorm_d1  !  turbulent residue normalized/same at iter 1 
  real, pointer ::   cont_resnorm_d1  !       mass residue normalized/same at iter 1
  real, pointer ::  x_mom_resnorm_d1  ! x momentum residue normalized/same at iter 1
  real, pointer ::  y_mom_resnorm_d1  ! Y momentum residue normalized/same at iter 1
  real, pointer ::  z_mom_resnorm_d1  ! Z momentum residue normalized/same at iter 1
  real, pointer :: energy_resnorm_d1  !     energy residue normalized/same at iter 1
  real, pointer ::    TKE_resnorm_d1  !        TKE residue normalized/same at iter 1
  real, pointer ::  omega_resnorm_d1  !      omega residue normalized/same at iter 1
  real          ::        resnorm_0   !            residue normalized at iter 1
  real          ::    vis_resnorm_0   ! {rho+V+P}  residue normalized at iter 1
  real          ::   turb_resnorm_0   !  turbulent residue normalized at iter 1 
  real          ::   cont_resnorm_0   !       mass residue normalized at iter 1
  real          ::  x_mom_resnorm_0   ! x momentum residue normalized at iter 1
  real          ::  y_mom_resnorm_0   ! Y momentum residue normalized at iter 1
  real          ::  z_mom_resnorm_0   ! Z momentum residue normalized at iter 1
  real          :: energy_resnorm_0   !     energy residue normalized at iter 1
  real          ::    TKE_resnorm_0   !        TKE residue normalized at iter 1
  real          ::  omega_resnorm_0   !      omega residue normalized at iter 1
  !used for MPI manipulation
  real          ::        resnorm_0s  !            residue normalized at iter 1 
  real          ::    vis_resnorm_0s  ! {rho+V+P}  residue normalized at iter 1 
  real          ::   turb_resnorm_0s  !  turbulent residue normalized at iter 1  
  real          ::   cont_resnorm_0s  !       mass residue normalized at iter 1 
  real          ::  x_mom_resnorm_0s  ! x momentum residue normalized at iter 1 
  real          ::  y_mom_resnorm_0s  ! Y momentum residue normalized at iter 1 
  real          ::  z_mom_resnorm_0s  ! Z momentum residue normalized at iter 1 
  real          :: energy_resnorm_0s  !     energy residue normalized at iter 1 
  real          ::    TKE_resnorm_0s  !        TKE residue normalized at iter 1 
  real          ::  omega_resnorm_0s  !      omega residue normalized at iter 1 
  real          :: merror=0.

  ! grid variables
  integer                                 :: imx, jmx, kmx        ! no. of points
  integer                                 :: imn, jmn, kmn
  real, dimension(:, :, :), allocatable  :: grid_x, grid_y, grid_z! point coordinates

  ! geometry variables
  real, dimension(:, :, :,:), allocatable, target :: xn           !face unit norm x
  real, dimension(:, :, :,:), allocatable, target :: yn           !face unit norm y
  real, dimension(:, :, :,:), allocatable, target :: zn           !face unit norm z
  real, dimension(:, :, :)             , pointer :: xnx, xny, xnz !face unit norm x
  real, dimension(:, :, :)             , pointer :: ynx, yny, ynz !face unit norm y
  real, dimension(:, :, :)             , pointer :: znx, zny, znz !face unit norm z
  real, dimension(:, :, :), allocatable, target :: xA, yA, zA    !face area
  real, dimension(:, :, :), allocatable, target :: volume
  real, dimension(:, :, :), allocatable, target ::   left_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::  right_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::  front_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::   back_ghost_centroid
  real, dimension(:, :, :), allocatable, target ::    top_ghost_centroid
  real, dimension(:, :, :), allocatable, target :: bottom_ghost_centroid
  

  ! gradients
  integer                                           :: n_grad=4
  real, dimension(:, :, :, :), allocatable, target  :: gradqp_x
  real, dimension(:, :, :, :), allocatable, target  :: gradqp_y
  real, dimension(:, :, :, :), allocatable, target  :: gradqp_z

  real, dimension(:, :, :),                 pointer :: gradu_x
  real, dimension(:, :, :),                 pointer :: gradu_y
  real, dimension(:, :, :),                 pointer :: gradu_z 
  real, dimension(:, :, :),                 pointer :: gradv_x 
  real, dimension(:, :, :),                 pointer :: gradv_y
  real, dimension(:, :, :),                 pointer :: gradv_z
  real, dimension(:, :, :),                 pointer :: gradw_x
  real, dimension(:, :, :),                 pointer :: gradw_y
  real, dimension(:, :, :),                 pointer :: gradw_z
  real, dimension(:, :, :),                 pointer :: gradT_x
  real, dimension(:, :, :),                 pointer :: gradT_y
  real, dimension(:, :, :),                 pointer :: gradT_z
  real, dimension(:, :, :),                 pointer :: gradtk_x
  real, dimension(:, :, :),                 pointer :: gradtk_y
  real, dimension(:, :, :),                 pointer :: gradtk_z
  real, dimension(:, :, :),                 pointer :: gradtw_x
  real, dimension(:, :, :),                 pointer :: gradtw_y
  real, dimension(:, :, :),                 pointer :: gradtw_z
  real, dimension(:, :, :),                 pointer :: gradtkl_x
  real, dimension(:, :, :),                 pointer :: gradtkl_y
  real, dimension(:, :, :),                 pointer :: gradtkl_z
  real, dimension(:, :, :),                 pointer :: gradte_x
  real, dimension(:, :, :),                 pointer :: gradte_y
  real, dimension(:, :, :),                 pointer :: gradte_z
  real, dimension(:, :, :),                 pointer :: gradtv_x
  real, dimension(:, :, :),                 pointer :: gradtv_y
  real, dimension(:, :, :),                 pointer :: gradtv_z

  ! higher order boundary condtioion
  integer  :: accur=1
  character(len=4), dimension(6) :: face_names
  integer,          dimension(6) :: id
  real                           :: c1
  real                           :: c2
  real                           :: c3

  ! store fix values for 6 faces of domain
  real, dimension(6) :: fixed_density  = 0.
  real, dimension(6) :: fixed_pressure = 0.
  real, dimension(6) :: fixed_x_speed  = 0.
  real, dimension(6) :: fixed_y_speed  = 0.
  real, dimension(6) :: fixed_z_speed  = 0.
  real, dimension(6) :: fixed_tk       = 0.
  real, dimension(6) :: fixed_tw       = 0.

  ! variable for post_processing
  integer :: N_blocks
  integer :: I_blocks
  integer :: J_blocks
  integer :: K_blocks
  integer, dimension(:), allocatable:: imin
  integer, dimension(:), allocatable:: imax
  integer, dimension(:), allocatable:: jmin
  integer, dimension(:), allocatable:: jmax
  integer, dimension(:), allocatable:: kmin
  integer, dimension(:), allocatable:: kmax

end module global_vars

