module kkl_source


  use global_kkl, only : zeta1
  use global_kkl, only : zeta2
  use global_kkl, only : zeta3
  use global_kkl, only : sigma_k
  use global_kkl, only : sigma_phi
  use global_kkl, only : cmu
  use global_kkl, only : kappa
  use global_kkl, only : c11
  use global_kkl, only : c12
  use global_kkl, only : cd1

  use global_kkl, only : cphi1
  use global_kkl, only : cphi2
  use global_kkl, only : fphi
  use global_kkl, only : eta

  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
  use global_vars, only : volume
  use global_vars  ,only : imx
  use global_vars  ,only : jmx
  use global_vars  ,only : kmx
  use global_vars  ,only : volume
  use global_vars  ,only : density
  use global_vars  ,only : pressure
  use global_vars  ,only : tk
  use global_vars  ,only : tkl
  use global_vars  ,only : mu
  use global_vars  ,only : mu_t
  use global_vars  ,only : dist
  use global_vars  ,only : gradu_x
  use global_vars  ,only : gradu_y
  use global_vars  ,only : gradu_z
  use global_vars  ,only : gradv_x
  use global_vars  ,only : gradv_y
  use global_vars  ,only : gradv_z
  use global_vars  ,only : gradw_x
  use global_vars  ,only : gradw_y
  use global_vars  ,only : gradw_z
  use global_vars  ,only : gradtk_x
  use global_vars  ,only : gradtk_y
  use global_vars  ,only : gradtk_z
  use global_vars  ,only : gradtw_x
  use global_vars  ,only : gradtw_y
  use global_vars  ,only : gradtw_z
  use global_vars  ,only : TKE_residue
  use global_vars  ,only : KL_residue

  implicit none
  private

  public :: add_kkl_source
  contains

    subroutine add_kkl_source()
      implicit none
      integer :: i,j,k

      real :: Tau11
      real :: Tau12
      real :: Tau13

      real :: Tau21
      real :: Tau22
      real :: Tau23

      real :: Tau31
      real :: Tau32
      real :: Tau33

      real :: S11
      real :: S12
      real :: S13

      real :: S21
      real :: S22
      real :: S23

      real :: S31
      real :: S32
      real :: S33

      real :: delv

      real :: d2udx2
      real :: d2udy2
      real :: d2udz2

      real :: d2vdx2
      real :: d2vdy2
      real :: d2vdz2

      real :: d2wdx2
      real :: d2wdy2
      real :: d2wdz2

      real :: Lvk
      real :: fp
      real :: ud
      real :: udd

      real :: S_k
      real :: S_kl
      real :: D_k
      real :: D_kl
      real :: P_k
      real :: P_kl


      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            S11 = 0.5*(gradu_x(i,j,k) + gradu_x(i,j,k))
            S12 = 0.5*(gradu_y(i,j,k) + gradv_x(i,j,k))
            S13 = 0.5*(gradu_z(i,j,k) + gradw_x(i,j,k))

            S21 = 0.5*(gradv_x(i,j,k) + gradu_y(i,j,k))
            S22 = 0.5*(gradv_y(i,j,k) + gradv_y(i,j,k))
            S23 = 0.5*(gradv_z(i,j,k) + gradw_y(i,j,k))

            S31 = 0.5*(gradw_x(i,j,k) + gradu_z(i,j,k))
            S32 = 0.5*(gradw_y(i,j,k) + gradv_z(i,j,k))
            S33 = 0.5*(gradw_z(i,j,k) + gradw_z(i,j,k))

            delv = gradu_x(i,j,k) + gradv_y(i,j,k) + gradw_z(i,j,k)

            Tau11 = mu_t(i,j,k)*(2*S11 - (2/3)*delv) - (2/3)*density(i,j,k)*tk(i,j,k)
            Tau12 = mu_t(i,j,k)*(2*S12)
            Tau13 = mu_t(i,j,k)*(2*S13)
            Tau21 = mu_t(i,j,k)*(2*S21)
            Tau22 = mu_t(i,j,k)*(2*S22 - (2/3)*delv) - (2/3)*density(i,j,k)*tk(i,j,k)
            Tau23 = mu_t(i,j,k)*(2*S23)
            Tau31 = mu_t(i,j,k)*(2*S31)
            Tau32 = mu_t(i,j,k)*(2*S32)
            Tau33 = mu_t(i,j,k)*(2*S33 - (2/3)*delv) - (2/3)*density(i,j,k)*tk(i,j,k)

            P_k = 0.
            P_k = P_k + Tau11*gradu_x(i,j,k) + Tau12*gradu_y(i,j,k) + Tau13*gradu_z(i,j,k)
            P_k = P_k + Tau21*gradv_x(i,j,k) + Tau22*gradv_y(i,j,k) + Tau23*gradv_z(i,j,k)
            P_k = P_k + Tau31*gradw_x(i,j,k) + Tau32*gradw_y(i,j,k) + Tau33*gradw_z(i,j,k)
            D_k = (cmu**0.75)*density(i,j,k)*(tk(i,j,k)**2.5)/tkl(i,j,k)
            P_k = min(P_k, 20*D_k)

            ! calculation of Lvk
            ! first get second order gradients 
            d2udx2 =(-(gradu_x(i-1,j  ,k  )+gradu_x(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradu_x(i  ,j-1,k  )+gradu_x(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradu_x(i  ,j  ,k-1)+gradu_x(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradu_x(i+1,j  ,k  )+gradu_x(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradu_x(i  ,j+1,k  )+gradu_x(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradu_x(i  ,j  ,k+1)+gradu_x(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2udy2 =(-(gradu_y(i-1,j  ,k  )+gradu_y(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradu_y(i  ,j-1,k  )+gradu_y(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradu_y(i  ,j  ,k-1)+gradu_y(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradu_y(i+1,j  ,k  )+gradu_y(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradu_y(i  ,j+1,k  )+gradu_y(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradu_y(i  ,j  ,k+1)+gradu_y(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2udz2 =(-(gradu_z(i-1,j  ,k  )+gradu_z(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradu_z(i  ,j-1,k  )+gradu_z(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradu_z(i  ,j  ,k-1)+gradu_z(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradu_z(i+1,j  ,k  )+gradu_z(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradu_z(i  ,j+1,k  )+gradu_z(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradu_z(i  ,j  ,k+1)+gradu_z(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            ! gradient of v component
            d2vdx2 =(-(gradv_x(i-1,j  ,k  )+gradv_x(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradv_x(i  ,j-1,k  )+gradv_x(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradv_x(i  ,j  ,k-1)+gradv_x(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradv_x(i+1,j  ,k  )+gradv_x(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradv_x(i  ,j+1,k  )+gradv_x(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradv_x(i  ,j  ,k+1)+gradv_x(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2vdy2 =(-(gradv_y(i-1,j  ,k  )+gradv_y(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradv_y(i  ,j-1,k  )+gradv_y(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradv_y(i  ,j  ,k-1)+gradv_y(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradv_y(i+1,j  ,k  )+gradv_y(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradv_y(i  ,j+1,k  )+gradv_y(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradv_y(i  ,j  ,k+1)+gradv_y(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2vdz2 =(-(gradv_z(i-1,j  ,k  )+gradv_z(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradv_z(i  ,j-1,k  )+gradv_z(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradv_z(i  ,j  ,k-1)+gradv_z(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradv_z(i+1,j  ,k  )+gradv_z(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradv_z(i  ,j+1,k  )+gradv_z(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradv_z(i  ,j  ,k+1)+gradv_z(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))


            !gradients of w components
            d2wdx2 =(-(gradw_x(i-1,j  ,k  )+gradw_x(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradw_x(i  ,j-1,k  )+gradw_x(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradw_x(i  ,j  ,k-1)+gradw_x(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradw_x(i+1,j  ,k  )+gradw_x(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradw_x(i  ,j+1,k  )+gradw_x(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradw_x(i  ,j  ,k+1)+gradw_x(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2wdy2 =(-(gradw_y(i-1,j  ,k  )+gradw_y(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradw_y(i  ,j-1,k  )+gradw_y(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradw_y(i  ,j  ,k-1)+gradw_y(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradw_y(i+1,j  ,k  )+gradw_y(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradw_y(i  ,j+1,k  )+gradw_y(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradw_y(i  ,j  ,k+1)+gradw_y(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2wdz2 =(-(gradw_z(i-1,j  ,k  )+gradw_z(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradw_z(i  ,j-1,k  )+gradw_z(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradw_z(i  ,j  ,k-1)+gradw_z(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradw_z(i+1,j  ,k  )+gradw_z(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradw_z(i  ,j+1,k  )+gradw_z(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradw_z(i  ,j  ,k+1)+gradw_z(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            udd = sqrt( (d2udx2+d2udy2+d2udz2)**2 &
                      + (d2udx2+d2udy2+d2udz2)**2 &
                      + (d2udx2+d2udy2+d2udz2)**2 )

            ud  = sqrt(2*(s11**2 + s12**2 + s13**2 &
                         +s21**2 + s22**2 + s23**2 &
                         +s31**2 + s32**2 + s33**2 ))

            Lvk = kappa*abs(ud/udd)

            fp = min(max(P_k/D_k, 0.5),1.0)
            ! Lvk limiter
            Lvk = max(Lvk, tkl(i,j,k)/(tk(i,j,k)*c11))
            Lvk = min(Lvk, c12*kappa*dist(i,j,k)*fp)

            eta = density(i,j,k)*dist(i,j,k)*sqrt(0.3*tk(i,j,k))/(20*mu(i,j,k))
            fphi = (1 + cd1*eta)/(1 + eta**4)
            cphi2 = zeta3
            cphi1 = (zeta1 - zeta2*((tkl(i,j,k)/(tk(i,j,k)*Lvk))**2))

            P_kl = cphi1*tkl(i,j,k)*P_k/tk(i,j,k)
            D_kl = cphi2*density(i,j,k)*(tk(i,j,k)**1.5)


            S_k  = P_k  - D_k  - 2*mu(i,j,k)*tk(i,j,k)/(dist(i,j,k)**2)       !Source term TKE
            S_kl = P_kl - D_kl - 6*mu(i,j,k)*tkl(i,j,k)*fphi/(dist(i,j,k)**2) !source term KL

            S_k  = S_k  * volume(i, j, k)
            S_kl = S_kl * volume(i, j, k)

            TKE_residue(i, j, k)   = TKE_residue(i, j, k) - S_k
            KL_residue(i, j, k)    = KL_residue(i, j, k) - S_kl

          end do
        end do
      end do

    end subroutine add_kkl_source

end module kkl_source
