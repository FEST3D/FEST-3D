module sst_source

  use global_sst   ,only : sigma_k1
  use global_sst   ,only : sigma_k2
  use global_sst   ,only : sigma_w1
  use global_sst   ,only : sigma_w2
  use global_sst   ,only : beta1
  use global_sst   ,only : beta2
  use global_sst   ,only : bstar
  use global_sst   ,only : kappa
  use global_sst   ,only : a1
  use global_sst   ,only : gama1
  use global_sst   ,only : gama2
  use global_sst   ,only : beta
  use global_sst   ,only : sigma_w
  use global_sst   ,only : sigma_k
  use global_sst   ,only : gama
  use global_sst   ,only : sst_F1

  use global_vars  ,only :   imx
  use global_vars  ,only :   jmx
  use global_vars  ,only :   kmx
  use global_vars  ,only :   volume
  use global_vars  ,only :   density
  use global_vars  ,only :   pressure
  use global_vars  ,only :   tk
  use global_vars  ,only :   tw
  use global_vars  ,only :   mu
  use global_vars  ,only :   sst_mu
  use global_vars  ,only :   dist
  use global_vars  ,only :   gradu_x
  use global_vars  ,only :   gradu_y
  use global_vars  ,only :   gradu_z
  use global_vars  ,only :   gradv_x
  use global_vars  ,only :   gradv_y
  use global_vars  ,only :   gradv_z
  use global_vars  ,only :   gradw_x
  use global_vars  ,only :   gradw_y
  use global_vars  ,only :   gradw_z
  use global_vars  ,only :   gradtk_x
  use global_vars  ,only :   gradtk_y
  use global_vars  ,only :   gradtk_z
  use global_vars  ,only :   gradtw_x
  use global_vars  ,only :   gradtw_y
  use global_vars  ,only :   gradtw_z
  use global_vars  ,only :   TKE_residue
  use global_vars  ,only : omega_residue

  implicit none
  private

  public :: add_sst_source
  contains

    subroutine add_sst_source()
      implicit none
      integer :: i,j,k

      real :: CD
      real :: F1
      real :: vort
      real :: S_k
      real :: S_w
      real :: D_k
      real :: D_w
      real :: P_k
      real :: P_w
      real :: lamda


      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            ! __ vorticity __
            vort = sqrt( 0.5 * ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                            + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                            + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                             )&
                       )

            CD = 2*density(i,j,k)*sigma_w2*(gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                          + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                          + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                           )/tw(i,j,k)
            CD = max(CD, 1e-20)
            F1 = sst_F1(i,j,k)


            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w1*(1. - F1)
            gama        =       gama1*F1  +       gama2*(1. - F1)
            beta        =       beta1*F1  +       beta2*(1. - F1)



            ! ____ Dissipation term ___
            D_k = bstar*density(i,j,k)*tw(i,j,k)*tk(i,j,k)
            D_w = beta*density(i,j,k)*tw(i,j,k)**2

            ! ____ PRODUCTION term____ 
            P_k = sst_mu(i,j,k)*(vort**2)
            P_k = min(P_k,20*D_k)
            P_w = min(gama*density(i,j,k)*P_k*sst_mu(i,j,k), 10.0*D_w)

            ! ____ cross diffusion term ___
            lamda = (1. - F1)*CD

            S_k = P_k - D_k           !Source term TKE
            S_w = P_w - D_w  +lamda   !source term omega

            S_k = S_k * volume(i, j, k)
            S_w = S_w * volume(i, j, k)

            TKE_residue(i, j, k)   = TKE_residue(i, j, k) - S_k
            omega_residue(i, j, k) = omega_residue(i, j, k) - S_w

          end do
        end do
      end do

    end subroutine add_sst_source

end module sst_source
