module global_sst

  real, parameter  :: sigma_k1   = 0.85
  real, parameter  :: sigma_k2   = 1.0
  real, parameter  :: sigma_w1   = 0.5
  real, parameter  :: sigma_w2   = 0.856
  real, parameter  :: beta1      = 0.075
  real, parameter  :: beta2      = 0.0828
  real, parameter  :: bstar      = 0.09
  real, parameter  :: kappa      = 0.41
  real, parameter  :: a1         = 0.31

  real, dimension(:,:,:), allocatable :: sst_F1

end module global_sst
