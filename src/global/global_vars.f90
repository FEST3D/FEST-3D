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
  integer :: total_entries      !< Total enteries in layout.md for each processor
  integer :: process_id         !< Id no. of each processor assinged by MPICH library
  integer :: imin_id            !< Boundary condition number/ID at imin for particulat processor
  integer :: imax_id            !< Boundary condition number/ID at imax for particulat processor
  integer :: jmin_id            !< Boundary condition number/ID at jmin for particulat processor
  integer :: jmax_id            !< Boundary condition number/ID at jmax for particulat processor
  integer :: kmin_id            !< Boundary condition number/ID at kmin for particulat processor
  integer :: kmax_id            !< Boundary condition number/ID at kmax for particulat processor
  integer :: layers=3           !< Number of ghost cell layers to transfer with mpi
  character(len=FILE_NAME_LENGTH)::     outfile          !< String to store name of output file
  character(len=FILE_NAME_LENGTH)::      infile          !< String to store the name of restart/load file
  character(len=FILE_NAME_LENGTH):: restartfile          !< Sting to store the name of restart log file
  integer :: want_to_stop=0    !< 0: continue the solver; 1=Stop the solver
  logical :: Halt = .FALSE.    !< Logical value used to stop the solver in main program file.

  real, dimension(:, :, :), allocatable     :: delta_t               !< Local time increment value at each cell center
  real                                      :: sim_clock             !< Simluation clock time


  !solution specific (used for Ranga_kutta_4th order)
  real, dimension(:, :, :, :), allocatable  :: qp_n
  real, dimension(:, :, :, :), allocatable  :: dEdx_1
  real, dimension(:, :, :, :), allocatable  :: dEdx_2
  real, dimension(:, :, :, :), allocatable  :: dEdx_3
  real, dimension(:, :, :, :), allocatable, target  :: qp           
   !< Store primitive variable at cell center
  real, dimension(:, :, :)                , pointer :: density      
   !< Rho pointer, point to slice of qp (:,:,:,1)
  real, dimension(:, :, :)                , pointer :: x_speed      
   !< U pointer, point to slice of qp (:,:,:,2) 
  real, dimension(:, :, :)                , pointer :: y_speed      
   !< V pointer, point to slice of qp (:,:,:,3) 
  real, dimension(:, :, :)                , pointer :: z_speed      
   !< W pointer, point to slice of qp (:,:,:,4)
  real, dimension(:, :, :)                , pointer :: pressure     
   !< P pointer, point to slice of qp (:,:,:,5)
  real, dimension(:, :, :), allocatable, target  :: intermittency
   !< Intermiitency pointer
  real, dimension(:, :, :), allocatable             :: dist 
   !< Store wall distance for each cell center
  real, dimension(:, :, :), allocatable             :: CCnormalX 
   !< Cell-Center normal nx with respect to wall; used for transition model (pressure gradient calcualtion)
  real, dimension(:, :, :), allocatable             :: CCnormalY
   !< Cell-Center normal ny with respect to wall; used for transiton model (pressure gradient calculation)
  real, dimension(:, :, :), allocatable             :: CCnormalZ
   !< Cell-Center normal nz with respect to wall; used for transiton model (pressure gradient calculation)
  real, dimension(:, :, :), allocatable             :: CCVn 
  !< Store value at Cell-Center of dot product between velocity vector and cell-center normal. {vec(Velocity).normal}
  real, dimension(:, :, :), allocatable             :: DCCVnX
  !< Store Derivative of Cell-Center CCVn with respect to x
  real, dimension(:, :, :), allocatable             :: DCCVnY
  !< Store Derivative of Cell-Center CCVn with respect to y
  real, dimension(:, :, :), allocatable             :: DCCVnZ
  !< Store Derivative of Cell-Center CCVn with respect to z

  ! state variable turbulent
  integer                                           :: sst_n_var=2 
  integer                                           :: sst_n_grad=2
  real, dimension(:, :, :)                , pointer :: tk        !< TKE/mass
  real, dimension(:, :, :)                , pointer :: tw        !< Omega
  real, dimension(:, :, :)                , pointer :: te        !< Dissipation
  real, dimension(:, :, :)                , pointer :: tv        !< SA visocity
  real, dimension(:, :, :)                , pointer :: tkl       !< KL K-KL method
  real, dimension(:, :, :)                , pointer :: tgm       !< Intermittency of LCTM2015
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
  real, dimension(:, :, :), allocatable, target     :: mu
   !< Cell-center molecular viscosity
  real, dimension(:, :, :), allocatable, target     :: mu_t
   !< Cell-center turbulent viscosity
  real, dimension(:, :, :)              , pointer   :: sst_mu
   !< Pointer to  turbulent viscosity for SST turbulence model
  real, dimension(:, :, :)              , pointer   :: kkl_mu
   !< Pointer to  turbulent viscosity for KKL turbulence model
  real, dimension(:, :, :)              , pointer   :: sa_mu
   !< Pointer to  turbulent viscosity for SA turbulence model
  real, pointer ::        resnorm     !<             Residual normalized
  real, pointer ::    vis_resnorm     !<  {rho+V+P} equation residual normalized
  real, pointer ::   turb_resnorm     !<  Turbulent residual normalized
  real, pointer ::   cont_resnorm     !<  Mass residual normalized
  real, pointer ::  x_mom_resnorm     !<  X momentum residual normalized
  real, pointer ::  y_mom_resnorm     !<  Y momentum residual normalized
  real, pointer ::  z_mom_resnorm     !<  Z momentum residual normalized
  real, pointer :: energy_resnorm     !<  Energy residual normalized
  real, pointer ::    TKE_resnorm     !<  TKE residual normalized
  real, pointer ::  omega_resnorm     !<  Omega residual normalized
  real, pointer ::        resnorm_d1  !<  Residual normalized/same at iter 1
  real, pointer ::    vis_resnorm_d1  !<  {rho+V+P}  residual normalized/same at iter 1
  real, pointer ::   turb_resnorm_d1  !<  Turbulent residual normalized/same at iter 1 
  real, pointer ::   cont_resnorm_d1  !<  Mass residual normalized/same at iter 1
  real, pointer ::  x_mom_resnorm_d1  !<  X momentum residual normalized/same at iter 1
  real, pointer ::  y_mom_resnorm_d1  !<  Y momentum residual normalized/same at iter 1
  real, pointer ::  z_mom_resnorm_d1  !<  Z momentum residual normalized/same at iter 1
  real, pointer :: energy_resnorm_d1  !<  Energy residual normalized/same at iter 1
  real, pointer ::    TKE_resnorm_d1  !<  TKE residual normalized/same at iter 1
  real, pointer ::  omega_resnorm_d1  !<  Omega residual normalized/same at iter 1
  real          ::        resnorm_0   !<  Residual normalized at iter 1
  real          ::    vis_resnorm_0   !<  {rho+V+P}  residual normalized at iter 1
  real          ::   turb_resnorm_0   !<  Turbulent residual normalized at iter 1 
  real          ::   cont_resnorm_0   !<  Mass residual normalized at iter 1
  real          ::  x_mom_resnorm_0   !<  X momentum residual normalized at iter 1
  real          ::  y_mom_resnorm_0   !<  Y momentum residual normalized at iter 1
  real          ::  z_mom_resnorm_0   !<  Z momentum residual normalized at iter 1
  real          :: energy_resnorm_0   !<  Energy residual normalized at iter 1
  real          ::    TKE_resnorm_0   !<  TKE residual normalized at iter 1
  real          ::  omega_resnorm_0   !<  Omega residual normalized at iter 1
  !used for MPI manipulation
  real          ::        resnorm_0s  !<  Residual normalized at iter 1 
  real          ::    vis_resnorm_0s  !<  {rho+V+P}  residual normalized at iter 1 
  real          ::   turb_resnorm_0s  !<  Turbulent residual normalized at iter 1  
  real          ::   cont_resnorm_0s  !<  Mass residual normalized at iter 1 
  real          ::  x_mom_resnorm_0s  !<  X momentum residual normalized at iter 1 
  real          ::  y_mom_resnorm_0s  !<  Y momentum residual normalized at iter 1 
  real          ::  z_mom_resnorm_0s  !<  Z momentum residual normalized at iter 1 
  real          :: energy_resnorm_0s  !<  Energy residual normalized at iter 1 
  real          ::    TKE_resnorm_0s  !<  TKE residual normalized at iter 1 
  real          ::  omega_resnorm_0s  !<  Omega residual normalized at iter 1 
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
  

  ! higher order boundary condtioion
!  integer  :: accur=1                          !< Switch for higher order boundary condition
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
  !< X component of velocity to fix at particular boundary condition
  real, dimension(6) :: fixed_y_speed  = 0.
  !< Y component of velocity to fix at particular boundary condition
  real, dimension(6) :: fixed_z_speed  = 0.
  !< Z component of velocity to fix at particular boundary condition
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

