  !< Declare all the constants used by Wilcox 2006 turbulence model
module global_wilcox2006
  !< Declare all the constants used by Wilcox 2006 turbulence model
  use iso_fortran_env, only : wp => real64

  real(wp), parameter  :: sigma_k = 0.6
  real(wp), parameter  :: sigma_w = 0.5
  real(wp), parameter  :: beta0   = 0.0708
  real(wp), parameter  :: bstar   = 0.09
  real(wp), parameter  :: a1      = 0.31
  real(wp), parameter  :: gama    = 13.0/25.0
  real(wp), parameter  :: clim    = 7.0/8.0

  ! to be used after blending
  real(wp)             :: w6_beta

end module global_wilcox2006
