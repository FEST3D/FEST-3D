  !< Declare all the constants used by SST turbulence model
module global_sst
  !< Declare all the constants used by SST turbulence model

  real, parameter  :: sigma_k1   = 0.85
  real, parameter  :: sigma_k2   = 1.0
  real, parameter  :: sigma_w1   = 0.5
  real, parameter  :: sigma_w2   = 0.856
  real, parameter  :: beta1      = 0.075
  real, parameter  :: beta2      = 0.0828
  real, parameter  :: bstar      = 0.09
  real, parameter  :: kappa      = 0.41
  real, parameter  :: a1         = 0.31
  real  :: gama1     = (beta1/bstar)-((sigma_w1*(kappa**2))/sqrt(bstar))
  real  :: gama2     = (beta2/bstar)-((sigma_w2*(kappa**2))/sqrt(bstar))

  ! to be used after blending with F1
  real             :: beta
  real             :: sigma_w
  real             :: sigma_k
  real             :: gama

  ! blending function
  real, dimension(:,:,:), allocatable, target :: sst_F1

end module global_sst
