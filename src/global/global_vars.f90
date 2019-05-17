  !< Common variables (can be re-assigned) used by other modules
module global_vars
  !< Contains all the public/global variables used by more than one module
  !----------------------------------------------

  use global, only : INTERPOLANT_NAME_LENGTH
  use global, only : FORMAT_LENGTH
  use global, only : SCHEME_NAME_LENGTH
  use global, only : FILE_NAME_LENGTH
  use global, only : STATE_NAME_LENGTH
  use global, only : FLOW_TYPE_LENGTH
  use global, only : TOLERANCE_LENGTH

  implicit none
  public

  ! Parallel processig variables
  integer :: total_process      !< Total number of process to be used for computation
  integer :: total_entries      !< total enteries in layout.md for each processor
  integer :: process_id         !< id no. of each processor assinged by MPICH library
  integer :: imin_id            !< Boundary condition number/ID at imin for particulat processor
  integer :: imax_id            !< Boundary condition number/ID at imax for particulat processor
  integer :: jmin_id            !< Boundary condition number/ID at jmin for particulat processor
  integer :: jmax_id            !< Boundary condition number/ID at jmax for particulat processor
  integer :: kmin_id            !< Boundary condition number/ID at kmin for particulat processor
  integer :: kmax_id            !< Boundary condition number/ID at kmax for particulat processor
  integer :: layers=3           !< Number of ghost cell layers to transfer with mpi

  ! Time controls
  integer :: min_iter=1              !< Minimum iteration value, starting iteration value
  integer :: max_iters=1             !< Maximum iteration value, stop after these many iteration
  integer :: start_from=0            !< Number of the folder (in time_directories) to load stored state from to restart computation
  integer :: checkpoint_iter=0       !< Write interval for output file. Number of iteration after which solver will dump/store a state in a folder in time_directories
  integer :: checkpoint_iter_count=0 !< Counter of folder number to write in time_directories/
  integer :: current_iter=0          !< current iteration number

  !write controls
  integer :: r_count=0                               !< Number of variable to read from the restart file
  integer :: w_count=0                               !< Number of variable to write in the output file
  integer :: res_write_interval                      !< Write resnorm after every "res_write_interval" iteration
  integer :: purge_write                             !< Remove unwanted folder. If Purge_write=2, latest two folder in time_direcotires are kept and 0=no purge
  integer :: last_iter=0                             !< Last iteration that is stored in the restart file
  integer :: write_percision=6                       !< number of place after decimal. Only used for resnorm file
  character(len=FORMAT_LENGTH):: write_data_format   !< write data type. Either ASCII or BINARY
  character(len=FORMAT_LENGTH):: write_file_format   !< Write file type. Either vtk or tecplot
  character(len=FORMAT_LENGTH)::  read_data_format='ASCII'   !< Read data type in file. Either ASCII or BINARY
  character(len=FORMAT_LENGTH)::  read_file_format="vtk"     !< Read file type. Either vtk or tecplot
  character(len=FILE_NAME_LENGTH)::     outfile          !< String to store name of output file
  character(len=FILE_NAME_LENGTH)::      infile          !< String to store the name of restart/load file
  character(len=FILE_NAME_LENGTH):: restartfile          !< Sting to store the name of restart log file
  character(len=STATE_NAME_LENGTH), dimension(:), allocatable ::  r_list !< Read variable list
  character(len=STATE_NAME_LENGTH), dimension(:), allocatable ::  w_list !< Write variable list
  character(len=FLOW_TYPE_LENGTH):: previous_flow_type="none"   !< Type of flow:inviscid, laminar, etc, stored in the load file 

  ! solver specific
  real :: CFL                  !< Courant–Friedrichs–Lewy (CFL) (Read from input)
  real :: tolerance            !< Minimum value of resnorm after which simulation stop
  character(len=TOLERANCE_LENGTH) :: tolerance_type="abs" !< Type of tolerance to check:absolute or relative
  integer :: want_to_stop=0    !< 0: continue the solver; 1=Stop the solver
  logical :: Halt = .FALSE.    !< Logical value used to stop the solver in main program file.

  !solver time secific
  character                                 :: time_stepping_method  !< Either local time stepping or global time stepping
  character(len=INTERPOLANT_NAME_LENGTH)    :: time_step_accuracy    !< Type of time_integration scheme: RK4, none(firt order explicit) implicit, 
  real                                      :: global_time_step      !< value of global time step to march the solution with
  real, dimension(:, :, :), allocatable     :: delta_t               !< local time increment value at each cell center.
  real                                      :: sim_clock             !< Simluation clock time.

  !scheme
  character(len=SCHEME_NAME_LENGTH) :: scheme_name         !< flux Scheme to use: ausm, ldfss0, vanleer, ausmup, ausmp, slau
  character(len=INTERPOLANT_NAME_LENGTH) :: interpolant    !< face state reconstruction  method to user: muscl, ppm, none, weno, 

  !solution specific (used for Ranga_kutta_4th order)
  real, dimension(:, :, :, :), allocatable  :: qp_n
  real, dimension(:, :, :, :), allocatable  :: dEdx_1
  real, dimension(:, :, :, :), allocatable  :: dEdx_2
  real, dimension(:, :, :, :), allocatable  :: dEdx_3

  ! state variables viscous
  integer                                           :: n_var=5
   !< number of variable to solve for
  real, dimension(:, :, :, :), allocatable, target  :: qp           
   !< Store primitive variable at cell center
  real, dimension(:)         , allocatable, target  :: qp_inf       
   !< Store primitive variable at infinity  
  real, dimension(:, :, :)                , pointer :: density      
   !< rho pointer, point to slice of qp (:,:,:,1)
  real, dimension(:, :, :)                , pointer :: x_speed      
   !< u pointer, point to slice of qp (:,:,:,2) 
  real, dimension(:, :, :)                , pointer :: y_speed      
   !< v pointer, point to slice of qp (:,:,:,3) 
  real, dimension(:, :, :)                , pointer :: z_speed      
   !< w pointer, point to slice of qp (:,:,:,4)
  real, dimension(:, :, :)                , pointer :: pressure     
   !< P pointer, point to slice of qp (:,:,:,5)
  real                                    , pointer :: density_inf  
   !< rho pointer, point to slice of qp_inf (1)
  real                                    , pointer :: x_speed_inf  
   !< u pointer, point to slice of qp_inf (2)
  real                                    , pointer :: y_speed_inf  
   !< v pointer, point to slice of qp_inf (3)
  real                                    , pointer :: z_speed_inf  
   !< w pointer, point to slice of qp_inf (4)
  real                                    , pointer :: pressure_inf 
   !< p pointer, point to slice of qp_inf (5)
  real                                              :: MInf
   !< Free-stream Mach number
  real, dimension(:, :, :), allocatable, target  :: intermittency
   !< intermiitency pointer
  real, dimension(:, :, :), allocatable, target  :: ExtraVar1
   !< Extravar1 used only for debuging or store some sepcial kind of compination of other varialbes
  real, dimension(:, :, :), allocatable, target  :: ExtraVar2
   !< Extravar2 used only for debuging or store some sepcial kind of compination of other varialbes
  real, dimension(:, :, :), allocatable, target  :: ExtraVar3
   !< Extravar3 used only for debuging or store some sepcial kind of compination of other varialbes
  real, dimension(:, :, :), allocatable, target  :: ExtraVar4
   !< Extravar4 used only for debuging or store some sepcial kind of compination of other varialbes
  real, dimension(:, :, :), allocatable, target  :: ExtraVar5
   !< Extravar5 used only for debuging or store some sepcial kind of compination of other varialbes

  ! Freestram variable used to read file before inf pointer are linked and allocated
  real                                              :: free_stream_density  
   !< read freestream Density from control file
  real                                              :: free_stream_x_speed  
   !< read freestream U from control file
  real                                              :: free_stream_y_speed  
   !< read freestream V from control file
  real                                              :: free_stream_z_speed  
   !< read freestream W from control file
  real                                              :: free_stream_pressure 
   !< read freestream Pressure from control file
  real                                              :: free_stream_tk       
   !< read freestream turbulent kinetic energy rate from control file
  real                                              :: free_stream_tw       
   !< read freestream turbulent dissipation rate from control file
  real                                              :: free_stream_te
   !< read freestream turbulent dissipation from control file
  real                                              :: free_stream_tv
   !< read freestream turbulent viscosity(SA) from control file
  real                                              :: free_stream_tkl
   !< read freestream kL variable from control file
  real                                              :: free_stream_tu
   !< read freestream turbulence intensity (percentage) from control file
  real                                              :: free_stream_tgm
   !< read freestream turbulence intermittency from control file
  real                                              :: vel_mag
   !< calulated freestream Velocity Magnitude from control file
  real                                              :: Reynolds_number 
   !< calculated free_stream Reynolds_number
  real                                              :: mu_ratio_inf 
   !< read freestream turbulent viscosity to molecular viscosity ratio
  real                                              :: Turb_intensity_inf 
   !< calculate free_stream turbulence intensity 
  real, dimension(:, :, :), allocatable             :: dist 
   !< store wall distance for each cell center
  real, dimension(:, :, :), allocatable             :: CCnormalX 
   !< Cell-Center normal nx with respect to wall; used for transition model (pressure gradient calcualtion)
  real, dimension(:, :, :), allocatable             :: CCnormalY
   !< Cell-Center normal ny with respect to wall; used for transiton model (pressure gradient calculation)
  real, dimension(:, :, :), allocatable             :: CCnormalZ
   !< Cell-Center normal nz with respect to wall; used for transiton model (pressure gradient calculation)
  real, dimension(:, :, :), allocatable             :: CCVn 
  !< store value at Cell-Center of dot product between velocity vector and cell-center normal. {vec(Velocity).normal}
  real, dimension(:, :, :), allocatable             :: DCCVnX
  !< store Derivative of Cell-Center CCVn with respect to x
  real, dimension(:, :, :), allocatable             :: DCCVnY
  !< store Derivative of Cell-Center CCVn with respect to y
  real, dimension(:, :, :), allocatable             :: DCCVnZ
  !< store Derivative of Cell-Center CCVn with respect to z

  ! state variable turbulent
  integer                                           :: sst_n_var=2 
  integer                                           :: sst_n_grad=2
!  real, dimension(:, :, :, :), allocatable, target  :: tqp       !< turbulent primitive
!  real, dimension(:)         , allocatable, target  :: tqp_inf   !< turbulent primitive at inf
  real, dimension(:, :, :)                , pointer :: tk        !< TKE/mass
  real, dimension(:, :, :)                , pointer :: tw        !< omega
  real, dimension(:, :, :)                , pointer :: te        !< Dissipation
  real, dimension(:, :, :)                , pointer :: tv        !< sa visocity
  real, dimension(:, :, :)                , pointer :: tkl       !< KL K-KL method
  real, dimension(:, :, :)                , pointer :: tgm       !< intermittency of LCTM2015
  real                                    , pointer :: tk_inf    !< TKE/mass at inf
  real                                    , pointer :: tw_inf    !< omega at inf
  real                                    , pointer :: te_inf    !< dissipation at inf
  real                                    , pointer :: tv_inf    !< SA viscosity at inf
  real                                    , pointer :: tkl_inf   !< kl at inf
  real                                    , pointer :: tgm_inf   !< intermittency at inf

  ! residue variables
  real, dimension(:, :, :, :)             , pointer :: F_p
   !< Flux pointer for face in the I direction
  real, dimension(:, :, :, :)             , pointer :: G_p
   !< Flux pointer for face in the G direction
  real, dimension(:, :, :, :)             , pointer :: H_p
   !< Flux pointer for face in the K direction
  real, dimension(:, :, :, :)             , pointer :: residue
   !< Store residue at each cell-center
  real, dimension(:, :, :)                , pointer :: mass_residue
   !< Store continuity equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: x_mom_residue
   !< Store x-momentum equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: y_mom_residue
   !< Store y-momentum equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: z_mom_residue
   !< Store z-momentum equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: energy_residue
   !< Store energy equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: TKE_residue
   !< Store TKE equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: omega_residue
   !< Store Omega equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: KL_residue
   !< Store KL equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: dissipation_residue
   !< Store Disspaiton equation residual at each cell-center
  real, dimension(:, :, :)                , pointer :: tv_residue
   !< Store nut equation(SA model) residual at each cell-center

  ! thermal properties
  real                                              :: gm    !< gamma commonly 1.4
  real                                              :: R_gas !< univarsal gas constant
  
  ! Transport properties
  real                                              :: mu_ref 
  !< molecular viscoity reference
  character(len=FILE_NAME_LENGTH)                   :: mu_variation
  !< Type of viscosity variaiton: Sutherland or constant

  ! sutherland law variable
  real                                              :: T_ref
  !< Reference Temperature of flow for viscosity calculation
  real                                              :: Sutherland_temp
  !< Sutherland temperature for viscosity calculation

  ! nondimensional numbers
  real                                              :: Pr=0.7 !< prandtl number
  real                                              :: tPr=0.9 !< turbulent Prandtl number

  ! switches
  logical                                           :: supersonic_flag
   !< Switch for boundary condition. No longer in use
  integer                                           :: ilimiter_switch
   !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for I direction faces.
  integer                                           :: jlimiter_switch
   !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for J direction faces.
  integer                                           :: klimiter_switch
   !< Turn on/off application of limiter for MUSCL (higer order face state reconstiion) for K direction faces.
  integer                                           :: itlimiter_switch
   !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for I direction faces.
  integer                                           :: jtlimiter_switch
   !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for J direction faces.
  integer                                           :: ktlimiter_switch
   !< Turn on/off application of limiter for MUSCL (higer order face turbulent variable state reconstiion) for K direction faces.
  integer                                           :: iPB_switch
   !< Turn on/off application of pressure based switching for higher order methods for I direction faces.
  integer                                           :: jPB_switch
   !< Turn on/off application of pressure based switching for higher order methods for J direction faces.
  integer                                           :: kPB_switch
   !< Turn on/off application of pressure based switching for higher order methods for K direction faces.
  character(len=8)                                  :: turbulence
   !< Store Turbulence model name
  character(len=8)                                  :: transition
   !< Store Transition model name
  
  real, dimension(:, :, :), allocatable, target     :: mu
   !< cell-center molecular viscosity
  real, dimension(:, :, :), allocatable, target     :: mu_t
   !< cell-center turbulent viscosity
  real, dimension(:, :, :)              , pointer   :: sst_mu
   !< pointer to  turbulent viscosity for SST turbulence model
  real, dimension(:, :, :)              , pointer   :: kkl_mu
   !< pointer to  turbulent viscosity for KKL turbulence model
  real, dimension(:, :, :)              , pointer   :: sa_mu
   !< pointer to  turbulent viscosity for SA turbulence model

  !residual specific
  character(len=STATE_NAME_LENGTH), dimension(:), allocatable :: Res_list
  integer            :: Res_count       !< no of variable to save
  integer            :: Res_itr=3       !< iteration to save
  real, dimension(:), allocatable :: Res_abs       !< absolute value
  real, dimension(:), allocatable :: Res_rel       !< relative value
  real, dimension(:), allocatable :: Res_save      !< saved iteration for relative
  real, dimension(:), allocatable :: Res_scale     !< scaling factor
  real, pointer ::        resnorm     !<             residual normalized
  real, pointer ::    vis_resnorm     !<  {rho+V+P} equation residual normalized
  real, pointer ::   turb_resnorm     !<  turbulent residual normalized
  real, pointer ::   cont_resnorm     !<  mass residual normalized
  real, pointer ::  x_mom_resnorm     !<  x momentum residual normalized
  real, pointer ::  y_mom_resnorm     !<  Y momentum residual normalized
  real, pointer ::  z_mom_resnorm     !<  Z momentum residual normalized
  real, pointer :: energy_resnorm     !<  energy residual normalized
  real, pointer ::    TKE_resnorm     !<  TKE residual normalized
  real, pointer ::  omega_resnorm     !<  omega residual normalized
  real, pointer ::        resnorm_d1  !<  residual normalized/same at iter 1
  real, pointer ::    vis_resnorm_d1  !<  {rho+V+P}  residual normalized/same at iter 1
  real, pointer ::   turb_resnorm_d1  !<  turbulent residual normalized/same at iter 1 
  real, pointer ::   cont_resnorm_d1  !<  mass residual normalized/same at iter 1
  real, pointer ::  x_mom_resnorm_d1  !<  x momentum residual normalized/same at iter 1
  real, pointer ::  y_mom_resnorm_d1  !<  Y momentum residual normalized/same at iter 1
  real, pointer ::  z_mom_resnorm_d1  !<  Z momentum residual normalized/same at iter 1
  real, pointer :: energy_resnorm_d1  !<  energy residual normalized/same at iter 1
  real, pointer ::    TKE_resnorm_d1  !<  TKE residual normalized/same at iter 1
  real, pointer ::  omega_resnorm_d1  !<  omega residual normalized/same at iter 1
  real          ::        resnorm_0   !<  residual normalized at iter 1
  real          ::    vis_resnorm_0   !<  {rho+V+P}  residual normalized at iter 1
  real          ::   turb_resnorm_0   !<  turbulent residual normalized at iter 1 
  real          ::   cont_resnorm_0   !<  mass residual normalized at iter 1
  real          ::  x_mom_resnorm_0   !<  x momentum residual normalized at iter 1
  real          ::  y_mom_resnorm_0   !<  Y momentum residual normalized at iter 1
  real          ::  z_mom_resnorm_0   !<  Z momentum residual normalized at iter 1
  real          :: energy_resnorm_0   !<  energy residual normalized at iter 1
  real          ::    TKE_resnorm_0   !<  TKE residual normalized at iter 1
  real          ::  omega_resnorm_0   !<  omega residual normalized at iter 1
  !used for MPI manipulation
  real          ::        resnorm_0s  !<  residual normalized at iter 1 
  real          ::    vis_resnorm_0s  !<  {rho+V+P}  residual normalized at iter 1 
  real          ::   turb_resnorm_0s  !<  turbulent residual normalized at iter 1  
  real          ::   cont_resnorm_0s  !<  mass residual normalized at iter 1 
  real          ::  x_mom_resnorm_0s  !<  x momentum residual normalized at iter 1 
  real          ::  y_mom_resnorm_0s  !<  Y momentum residual normalized at iter 1 
  real          ::  z_mom_resnorm_0s  !<  Z momentum residual normalized at iter 1 
  real          :: energy_resnorm_0s  !<  energy residual normalized at iter 1 
  real          ::    TKE_resnorm_0s  !<  TKE residual normalized at iter 1 
  real          ::  omega_resnorm_0s  !<  omega residual normalized at iter 1 

  ! grid variables
  integer                                 :: imx
   !< Maximum number of grid points in the I direction
  integer                                 :: jmx
   !< Maximum number of grid points in the K direction
  integer                                 :: kmx
   !< Maximum number of grid points in the K direction
  integer                                 :: imn, jmn, kmn
  real, dimension(:, :, :), allocatable  :: grid_x
   !< X corrdinate of the grid point
  real, dimension(:, :, :), allocatable  :: grid_y
   !< Y corrdinate of the grid point
  real, dimension(:, :, :), allocatable  :: grid_z
   !< z corrdinate of the grid point

  ! geometry variables
  real, dimension(:, :, :,:), allocatable, target :: xn
   !< Store unit face normal vector for all I faces 
  real, dimension(:, :, :,:), allocatable, target :: yn
   !< Store unit face normal vector for all J faces 
  real, dimension(:, :, :,:), allocatable, target :: zn
   !< Store unit face normal vector for all K faces 
  real, dimension(:, :, :)             , pointer :: xnx
   !< Pointer to x component of face unit normal of I faces
  real, dimension(:, :, :)             , pointer :: xny
   !< Pointer to y component of face unit normal of I faces
  real, dimension(:, :, :)             , pointer :: xnz
   !< Pointer to z component of face unit normal of I faces
  real, dimension(:, :, :)             , pointer :: ynx
   !< Pointer to x component of face unit normal of J faces
  real, dimension(:, :, :)             , pointer :: yny
   !< Pointer to y component of face unit normal of J faces
  real, dimension(:, :, :)             , pointer :: ynz
   !< Pointer to z component of face unit normal of J faces
  real, dimension(:, :, :)             , pointer :: znx
   !< Pointer to x component of face unit normal of K faces
  real, dimension(:, :, :)             , pointer :: zny
   !< Pointer to y component of face unit normal of K faces
  real, dimension(:, :, :)             , pointer :: znz
   !< Pointer to z component of face unit normal of K faces
  real, dimension(:, :, :), allocatable, target :: xA
   !< Store magnitude of face area vector of I direction faces
  real, dimension(:, :, :), allocatable, target :: yA
   !< Store magnitude of face area vector of J direction faces
  real, dimension(:, :, :), allocatable, target :: zA
   !< Store magnitude of face area vector of K direction faces
  real, dimension(:, :, :), allocatable, target :: volume
   !< Store cell volume
  real, dimension(:, :, :), allocatable, target ::   left_ghost_centroid
   !< Store the cell center of the ghost cell on 1 face
  real, dimension(:, :, :), allocatable, target ::  right_ghost_centroid
   !< Store the cell center of the ghost cell on 2 face
  real, dimension(:, :, :), allocatable, target ::  front_ghost_centroid
   !< Store the cell center of the ghost cell on 3 face
  real, dimension(:, :, :), allocatable, target ::   back_ghost_centroid
   !< Store the cell center of the ghost cell on 4 face
  real, dimension(:, :, :), allocatable, target ::    top_ghost_centroid
   !< Store the cell center of the ghost cell on 5 face
  real, dimension(:, :, :), allocatable, target :: bottom_ghost_centroid
   !< Store the cell center of the ghost cell on 6 face
  

  ! gradients
  integer                                           :: n_grad=4 !< number of variable to store gradient for
  real, dimension(:, :, :, :), allocatable, target  :: gradqp_x !< Store gradient of n_grad variables with respect to direction x
  real, dimension(:, :, :, :), allocatable, target  :: gradqp_y !< Store gradient of n_grad variables with respect to direction y
  real, dimension(:, :, :, :), allocatable, target  :: gradqp_z !< Store gradient of n_grad variables with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradu_x  !< Gradient of variable U with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradu_y  !< Gradient of variable U with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradu_z  !< Gradient of variable U with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradv_x  !< Gradient of variable V with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradv_y  !< Gradient of variable V with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradv_z  !< Gradient of variable V with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradw_x  !< Gradient of variable W with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradw_y  !< Gradient of variable W with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradw_z  !< Gradient of variable W with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradT_x  !< Gradient of variable Temperature with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradT_y  !< Gradient of variable Temperature with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradT_z  !< Gradient of variable Temperature with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradtk_x !< Gradient of variable turbulent kinetic energy with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradtk_y !< Gradient of variable turbulent kinetic energy with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradtk_z !< Gradient of variable turbulent kinetic energy with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradtw_x !< Gradient of variable dissipation rate with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradtw_y !< Gradient of variable dissipation rate with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradtw_z !< Gradient of variable dissipation rate with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradtkl_x!< Gradient of variable kL with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradtkl_y!< Gradient of variable kL with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradtkl_z!< Gradient of variable kL with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradte_x !< Gradient of variable turbulent energy dissiaption with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradte_y !< Gradient of variable turbulent energy dissiaption with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradte_z !< Gradient of variable turbulent energy dissiaption with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradtv_x !< Gradient of variable turbulenct visocity(SA mode) with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradtv_y !< Gradient of variable turbulenct visocity(SA mode) with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradtv_z !< Gradient of variable turbulenct visocity(SA mode) with respect to direction z
  real, dimension(:, :, :),                 pointer :: gradtgm_x!< Gradient of variable intermittency with respect to direction x
  real, dimension(:, :, :),                 pointer :: gradtgm_y!< Gradient of variable intermittency with respect to direction y
  real, dimension(:, :, :),                 pointer :: gradtgm_z!< Gradient of variable intermittency with respect to direction z

  ! higher order boundary condtioion
  integer  :: accur=1                          !< Switch for higher order boundary condition
  character(len=4), dimension(6) :: face_names !< Store name of all six boundary faces
  integer,          dimension(6) :: id         !< Store the boundary condition ID of all six faces
  real                           :: c1         !< First coefficient user for higher order boundary condition
  real                           :: c2         !< Second coefficient user for higher order boundary condition
  real                           :: c3         !< Third coefficient user for higher order boundary condition

  ! store fix values for 6 faces of domain
  real, dimension(6) :: fixed_density  = 0.
  !< Density value to fix at particular boundary
  real, dimension(6) :: fixed_pressure = 0.
  !< Pressure value to fix at particular boundary
  real, dimension(6) :: fixed_x_speed  = 0.
  !< x component of velocity to fix at particular boundary condition
  real, dimension(6) :: fixed_y_speed  = 0.
  !< y component of velocity to fix at particular boundary condition
  real, dimension(6) :: fixed_z_speed  = 0.
  !< z component of velocity to fix at particular boundary condition
  real, dimension(6) :: fixed_tk       = 0.
  !< Turbulent kinetic energy value to fix at particular boundary condition
  real, dimension(6) :: fixed_tw       = 0.
  !< Turbulent kinetic energy dissiaption rate value to fix at particular boundary condition(k-omega and SST model)
  real, dimension(6) :: fixed_te       = 0.
  !< Turbulent kinetic energy dissiaption value to fix at particular boundary condition (K-eplision model)
  real, dimension(6) :: fixed_tv       = 0.
  !< Turbulent viscosity varialble value to fix at particular boundary condition (for SA turbulence model)
  real, dimension(6) :: fixed_tkl       = 0.
  !< (Turbulent kinetic energy x length) varialble value to fix at particular boundary condition (for k-kL turbulence model)
  real, dimension(6) :: fixed_tgm       = 0.
  !<  Fixed intermittency value to apply at particular boundary condition (for SST2003-gamma transition model)
  real, dimension(6) :: fixed_wall_temperature  = 0.
  !<  Fixed wall temperature value to apply at isothermal wall boundary condition.
  real, dimension(6) :: fixed_Tpressure         = 0.
  !<  Fixed Total Pressure value to apply at particular boundary condition
  real, dimension(6) :: fixed_Ttemperature      = 0.
  !<  Fixed Total Temperature value to apply at particular boundary condition

  ! variable for post_processing
  integer :: N_blocks
  !< Total number of blocks 
  integer :: I_blocks
  !< Total number of blocks  in I direction
  integer :: J_blocks
  !< Total number of blocks  in J direction
  integer :: K_blocks
  !< Total number of blocks  in K direction
  integer, dimension(:), allocatable:: imin
  integer, dimension(:), allocatable:: imax
  integer, dimension(:), allocatable:: jmin
  integer, dimension(:), allocatable:: jmax
  integer, dimension(:), allocatable:: kmin
  integer, dimension(:), allocatable:: kmax

  !interface mapping
  integer, dimension(6) :: ilo, ihi 
   !< Store the lower and upper bound of the indecies of I loop for the interface mapping
  integer, dimension(6) :: jlo, jhi
   !< Store the lower and upper bound of the indecies of J loop for the interface mapping
  integer, dimension(6) :: klo, khi
   !< Store the lower and upper bound of the indecies of K loop for the interface mapping
  integer, dimension(6) :: dir_switch=0
   !< Switch for each boundary face. Activated only if ( for eg I-direction in the mapping is mapped with J-direction)
  integer, dimension(6) :: otherface
   !< Store the face number with which the current interface is connected.

  !zero flux faces
  integer,dimension(:),allocatable::make_F_flux_zero
   !< Store zero to boundary face, which has wall ID, to make F flux zero
  integer,dimension(:),allocatable::make_G_flux_zero
   !< Store zero to boundary face, which has wall ID, to make G flux zero
  integer,dimension(:),allocatable::make_H_flux_zero  
   !< Store zero to boundary face, which has wall ID, to make H flux zero

  !periodic boundary condition
  integer, dimension(6) :: PbcId = -1 !< Block ID for Periodic boundary condition

end module global_vars

