  !< Declare all the constants used by SA turbulence model
module global_sa
  !< Declare all the constants used by SA turbulence model
  use iso_fortran_env, only : wp => real64

  real(wp), parameter  :: cb1    = 0.1355
  real(wp), parameter  :: cb2    = 0.6220
  real(wp), parameter  :: cw2    = 0.3
  real(wp), parameter  :: cw3    = 2.0
  real(wp), parameter  :: cv1    = 7.1
  real(wp), parameter  :: ct3    = 1.2
  real(wp), parameter  :: ct4    = 0.5
  real(wp), parameter  :: sigma_sa  = 2./3.
  real(wp), parameter  :: kappa_sa  = 0.41

  real(wp), parameter  :: cw1    = (cb1/(kappa_sa**2)) + ((1+cb2)/sigma_sa)

  real(wp), parameter  :: cv1_3 = cv1**3
  real(wp), parameter  :: cw3_6 = cw3**6

end module global_sa
