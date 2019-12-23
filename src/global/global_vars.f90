  !< Common variables (can be re-assigned) used by other modules
module global_vars
  !< Contains all the public/global variables used by more than one module
  !----------------------------------------------

  implicit none
  public

  real, dimension(:, :, :), allocatable     :: delta_t  
  !< Local time increment value at each cell center
  real, dimension(:, :, :), allocatable             :: dist 

  real, dimension(:, :, :), allocatable, target     :: mu
   !< Cell-center molecular viscosity
  real, dimension(:, :, :), allocatable, target     :: mu_t
   !< Cell-center turbulent viscosity
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


end module global_vars

