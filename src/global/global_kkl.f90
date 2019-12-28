  !< Declare all the constants used by k-kL turbulence model
module global_kkl
  !< Declare all the constants used by k-kL turbulence model
  use iso_fortran_env, only : wp => real64

  real(wp), parameter  :: zeta1       = 1.2
  real(wp), parameter  :: zeta2       = 0.97
  real(wp), parameter  :: zeta3       = 0.13
  real(wp), parameter  :: sigma_k    = 1.0
  real(wp), parameter  :: sigma_phi   = 1.0
  real(wp), parameter  :: cmu        = 0.09
  real(wp), parameter  :: kappa      = 0.41
  real(wp), parameter  :: c11        = 10.0
  real(wp), parameter  :: c12        = 1.3
  real(wp), parameter  :: cd1        = 4.7

  real(wp)             :: cphi1
  real(wp)             :: cphi2
  real(wp)             :: fphi
  real(wp)             :: eta


end module global_kkl
