!should be included in viscous module (rename it as transport) 


    real :: sst_f_mu, F1_f
!    real, public, dimension(:, :, :), allocatable :: gradtk_x, gradtk_y, gradtk_z,&
!                                                  gradtw_x, gradtw_y, gradtw_z
    real, parameter :: sigma_k1 =       0.85
    real, parameter :: sigma_k2 =       1.0
    real, parameter :: sigma_omega1 =   0.5
    real, parameter :: sigma_omega2 =   0.856
    real, parameter :: beta1 =          0.075
    real, parameter :: beta2 =          0.0828
    real, parameter :: bs =             0.09 !beta*
    real, parameter :: kappa =          0.41
    real, parameter :: a1 =             0.31
    real :: g1_w, g2_w, g_w ! gamma with blend 
    real :: sigma_k, sigma_omega, beta
    !real :: F2, arg1, arg2
    !real :: vort, CD,
    !real :: KE, omega, nu, rho
    !real :: P_k, P_w !production term
    !real :: D_k, D_w !dissipation term
    !real :: S_k, S_w !source term
    !real :: lamda    !extra source due to blend

