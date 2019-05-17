  !< Declare all the constants used by k-kL turbulence model
module global_kkl
  !< Declare all the constants used by k-kL turbulence model

  real, parameter  :: zeta1       = 1.2
  real, parameter  :: zeta2       = 0.97
  real, parameter  :: zeta3       = 0.13
  real, parameter  :: sigma_k    = 1.0
  real, parameter  :: sigma_phi   = 1.0
  real, parameter  :: cmu        = 0.09
  real, parameter  :: kappa      = 0.41
  real, parameter  :: c11        = 10.0
  real, parameter  :: c12        = 1.3
  real, parameter  :: cd1        = 4.7

  real             :: cphi1
  real             :: cphi2
  real             :: fphi
  real             :: eta


end module global_kkl
