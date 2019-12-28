  !< Declare all the constants used by SST turbulence model
module global_sst
  !< Declare all the constants used by SST turbulence model
  use iso_fortran_env, only : wp => real64

  real(wp), parameter  :: sigma_k1   = 0.85
  real(wp), parameter  :: sigma_k2   = 1.0
  real(wp), parameter  :: sigma_w1   = 0.5
  real(wp), parameter  :: sigma_w2   = 0.856
  real(wp), parameter  :: beta1      = 0.075
  real(wp), parameter  :: beta2      = 0.0828
  real(wp), parameter  :: bstar      = 0.09
  real(wp), parameter  :: kappa      = 0.41
  real(wp), parameter  :: a1         = 0.31
  real(wp)  :: gama1     = (beta1/bstar)-((sigma_w1*(kappa**2))/sqrt(bstar))
  real(wp)  :: gama2     = (beta2/bstar)-((sigma_w2*(kappa**2))/sqrt(bstar))

  ! to be used after blending with F1
  real(wp)             :: beta
  real(wp)             :: sigma_w
  real(wp)             :: sigma_k
  real(wp)             :: gama

  ! blending function
  real(wp), dimension(:,:,:), allocatable, target :: sst_F1

end module global_sst
