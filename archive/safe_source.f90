module source
  !-----------------------------------------------------------------
  !170609  - jatinder Pal Singh Sandhu
  ! AIM: to add source term residue to already calculate residuals
  !      if requires
  !-----------------------------------------------------------------
#include "error.inc"
  use global_vars, only : turbulence
  use global_vars, only : process_id
  use utils      , only : dmsg
  use utils      , only : turbulence_read_error
  !--- variable required for sst source calculation ---!
  use global_sst   ,only : sigma_k1
  use global_sst   ,only : sigma_k2
  use global_sst   ,only : sigma_w1
  use global_sst   ,only : sigma_w2
  use global_sst   ,only : beta1
  use global_sst   ,only : beta2
  use global_sst   ,only : bstar
  use global_sst   ,only : a1
  use global_sst   ,only : gama1
  use global_sst   ,only : gama2
  use global_sst   ,only : beta
  use global_sst   ,only : sigma_w
  use global_sst   ,only : sigma_k
  use global_sst   ,only : gama
  use global_sst   ,only : sst_F1

  use global_vars  ,only : intermittency
  use global_vars  ,only : mu_ref
  use global_vars  ,only : density_inf
  use global_vars  ,only : vel_mag
  use global_vars  ,only :  qp
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
  use global_vars  ,only : xn
  use global_vars  ,only : yn
  use global_vars  ,only : zn

  !--- variables required for kkl source calculation ---!
  use global_kkl, only : zeta1
  use global_kkl, only : zeta2
  use global_kkl, only : zeta3
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

  !variables required by sa source term calculation
  use global_vars,only : tv
  use global_sa , only : cb1
  use global_sa , only : cb2
  use global_sa , only : cw1
  use global_sa , only : cw2
  use global_sa , only : cw3
  use global_sa , only : cv1
  use global_sa , only : sigma_sa
  use global_sa , only : kappa_sa
  use global_vars,only : gradtv_x
  use global_vars,only : gradtv_y
  use global_vars,only : gradtv_z

  use global_vars, only : xnx, xny, xnz !face unit normal x
  use global_vars, only : ynx, yny, ynz !face unit normal y
  use global_vars, only : znx, zny, znz !face unit normal z
  use global_vars, only : xA, yA, zA    !face area
  use global_vars  ,only : tkl
  use global_vars  ,only : mu_t
  use global_vars  ,only : KL_residue
  use global_vars  ,only : tv_residue


  implicit none
  private
  public :: add_source_term_residue

  contains

    
    subroutine add_source_term_residue()

      implicit none

      call dmsg(1, 'source', 'add_source_term_residue')

      select case (trim(turbulence))

        case ('none')
          !do nothing
          continue

        case ('sa')
          call add_sa_source()

        case ('saBC')
          call add_saBC_source()

        case ('sst')
          call add_sst_source()

        case ('kkl')
          call add_kkl_source()

        case DEFAULT
          !call turbulence_read_error()
          Fatal_error

      end select

    end subroutine add_source_term_residue


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
            vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                            + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                            + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                             )&
                       )

            CD = 2*density(i,j,k)*sigma_w2*(gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                          + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                          + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                           )/tw(i,j,k)
            !CD = max(CD, 1e-20)
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
            P_k = min(P_k,20.0*D_k)
            P_w = (density(i,j,k)*gama/sst_mu(i,j,k))*P_k

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

      !change for transition modeling
      real :: vort
      real :: Rev
      real :: Rev1
      real :: ReThc
      real :: Tu
      real :: term1


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
            D_k = (cmu**0.75)*density(i,j,k)*(tk(i,j,k)**2.5)/max(tkl(i,j,k),1.e-20)
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

            d2udy2 =(-(gradu_y(i-1,j  ,k  )+gradu_y(i,j,k))*xny(i,j,k)*xA(i,j,k) &
                     -(gradu_y(i  ,j-1,k  )+gradu_y(i,j,k))*yny(i,j,k)*yA(i,j,k) &
                     -(gradu_y(i  ,j  ,k-1)+gradu_y(i,j,k))*zny(i,j,k)*zA(i,j,k) &
                     +(gradu_y(i+1,j  ,k  )+gradu_y(i,j,k))*xny(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradu_y(i  ,j+1,k  )+gradu_y(i,j,k))*yny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradu_y(i  ,j  ,k+1)+gradu_y(i,j,k))*zny(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2udz2 =(-(gradu_z(i-1,j  ,k  )+gradu_z(i,j,k))*xnz(i,j,k)*xA(i,j,k) &
                     -(gradu_z(i  ,j-1,k  )+gradu_z(i,j,k))*ynz(i,j,k)*yA(i,j,k) &
                     -(gradu_z(i  ,j  ,k-1)+gradu_z(i,j,k))*znz(i,j,k)*zA(i,j,k) &
                     +(gradu_z(i+1,j  ,k  )+gradu_z(i,j,k))*xnz(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradu_z(i  ,j+1,k  )+gradu_z(i,j,k))*ynz(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradu_z(i  ,j  ,k+1)+gradu_z(i,j,k))*znz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            ! gradient of v component
            d2vdx2 =(-(gradv_x(i-1,j  ,k  )+gradv_x(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradv_x(i  ,j-1,k  )+gradv_x(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradv_x(i  ,j  ,k-1)+gradv_x(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradv_x(i+1,j  ,k  )+gradv_x(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradv_x(i  ,j+1,k  )+gradv_x(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradv_x(i  ,j  ,k+1)+gradv_x(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2vdy2 =(-(gradv_y(i-1,j  ,k  )+gradv_y(i,j,k))*xny(i,j,k)*xA(i,j,k) &
                     -(gradv_y(i  ,j-1,k  )+gradv_y(i,j,k))*yny(i,j,k)*yA(i,j,k) &
                     -(gradv_y(i  ,j  ,k-1)+gradv_y(i,j,k))*zny(i,j,k)*zA(i,j,k) &
                     +(gradv_y(i+1,j  ,k  )+gradv_y(i,j,k))*xny(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradv_y(i  ,j+1,k  )+gradv_y(i,j,k))*yny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradv_y(i  ,j  ,k+1)+gradv_y(i,j,k))*zny(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2vdz2 =(-(gradv_z(i-1,j  ,k  )+gradv_z(i,j,k))*xnz(i,j,k)*xA(i,j,k) &
                     -(gradv_z(i  ,j-1,k  )+gradv_z(i,j,k))*ynz(i,j,k)*yA(i,j,k) &
                     -(gradv_z(i  ,j  ,k-1)+gradv_z(i,j,k))*znz(i,j,k)*zA(i,j,k) &
                     +(gradv_z(i+1,j  ,k  )+gradv_z(i,j,k))*xnz(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradv_z(i  ,j+1,k  )+gradv_z(i,j,k))*ynz(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradv_z(i  ,j  ,k+1)+gradv_z(i,j,k))*znz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))


            !gradients of w components
            d2wdx2 =(-(gradw_x(i-1,j  ,k  )+gradw_x(i,j,k))*xnx(i,j,k)*xA(i,j,k) &
                     -(gradw_x(i  ,j-1,k  )+gradw_x(i,j,k))*ynx(i,j,k)*yA(i,j,k) &
                     -(gradw_x(i  ,j  ,k-1)+gradw_x(i,j,k))*znx(i,j,k)*zA(i,j,k) &
                     +(gradw_x(i+1,j  ,k  )+gradw_x(i,j,k))*xnx(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradw_x(i  ,j+1,k  )+gradw_x(i,j,k))*ynx(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradw_x(i  ,j  ,k+1)+gradw_x(i,j,k))*znx(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2wdy2 =(-(gradw_y(i-1,j  ,k  )+gradw_y(i,j,k))*xny(i,j,k)*xA(i,j,k) &
                     -(gradw_y(i  ,j-1,k  )+gradw_y(i,j,k))*yny(i,j,k)*yA(i,j,k) &
                     -(gradw_y(i  ,j  ,k-1)+gradw_y(i,j,k))*zny(i,j,k)*zA(i,j,k) &
                     +(gradw_y(i+1,j  ,k  )+gradw_y(i,j,k))*xny(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradw_y(i  ,j+1,k  )+gradw_y(i,j,k))*yny(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradw_y(i  ,j  ,k+1)+gradw_y(i,j,k))*zny(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            d2wdz2 =(-(gradw_z(i-1,j  ,k  )+gradw_z(i,j,k))*xnz(i,j,k)*xA(i,j,k) &
                     -(gradw_z(i  ,j-1,k  )+gradw_z(i,j,k))*ynz(i,j,k)*yA(i,j,k) &
                     -(gradw_z(i  ,j  ,k-1)+gradw_z(i,j,k))*znz(i,j,k)*zA(i,j,k) &
                     +(gradw_z(i+1,j  ,k  )+gradw_z(i,j,k))*xnz(i+1,j  ,k  )*xA(i+1,j  ,k  ) &
                     +(gradw_z(i  ,j+1,k  )+gradw_z(i,j,k))*ynz(i  ,j+1,k  )*yA(i  ,j+1,k  ) &
                     +(gradw_z(i  ,j  ,k+1)+gradw_z(i,j,k))*znz(i  ,j  ,k+1)*zA(i  ,j  ,k+1) &
                    )/(2*volume(i,j,k))

            udd = sqrt( (d2udx2+d2udy2+d2udz2)**2 &
                      + (d2vdx2+d2vdy2+d2vdz2)**2 &
                      + (d2wdx2+d2wdy2+d2wdz2)**2 )

            ud  = sqrt(2*(s11**2 + s12**2 + s13**2 &
                         +s21**2 + s22**2 + s23**2 &
                         +s31**2 + s32**2 + s33**2 ))

            Lvk = kappa*abs(ud/max(udd,1.e-20))

            fp = min(max(P_k/D_k, 0.5),1.0)
            ! Lvk limiter
            Lvk = max(Lvk, tkl(i,j,k)/max((tk(i,j,k)*c11),1.e-20))
            Lvk = min(Lvk, c12*kappa*dist(i,j,k)*fp)

            eta = density(i,j,k)*dist(i,j,k)*sqrt(0.3*tk(i,j,k))/(20*mu(i,j,k))
            fphi = (1 + cd1*eta)/(1 + eta**4)
            cphi2 = zeta3
            cphi1 = (zeta1 - zeta2*((tkl(i,j,k)/max(tk(i,j,k)*Lvk,1.e-20))**2))

            P_kl = cphi1*tkl(i,j,k)*P_k/max(tk(i,j,k),1.e-20)
            D_kl = cphi2*density(i,j,k)*(tk(i,j,k)**1.5)

            !-------------transition moddeling--
!            Tu = min(100.0*sqrt(2.0*tk(i,j,k)/3.0)/(dist(i,j,k)*vort/0.3), 100.0)
!            ReThc = 163.0 + exp(6.91 - Tu)
!            ! __ vorticity __
!            vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
!                            + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
!                            + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
!                             )&
!                       )
!            Rev = density(i,j,k)*dist(i,j,k)*dist(i,j,k)*vort/mu(i,j,k)
!            Rev1 = density(i,j,k)*((Lvk)**2)*vort/mu(i,j,k)
!            term1 = max((Rev/(2.2*ReThc))-1.0, 0.0)
!            !intermittency(i,j,k) = 1.0 - exp(-sqrt(term1))
!            intermittency(i,j,k) = Rev/Rev1
!            !P_k = intermittency(i,j,k)*P_k
            !--- end transition modeling ---

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


    subroutine add_sa_source()
      implicit none
      integer :: i,j,k

      real :: CD1
      real :: CD2
      real :: fv1
      real :: fv2
      real :: fw
      real :: g
      real :: Scap
      real :: r
      real :: vort
      real :: S_v
      real :: D_v
      real :: P_v
      real :: lamda
      real :: kd2
      real :: xi
      real :: nu
      real :: gradrho_x
      real :: gradrho_y
      real :: gradrho_z
      real, dimension(6) :: RhoFace
      real, dimension(6) :: Area
      real, dimension(6,3) :: Normal

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            RhoFace(1) = density(i-1,j  ,k  )+density(i,j,k)
            RhoFace(2) = density(i  ,j-1,k  )+density(i,j,k)
            RhoFace(3) = density(i  ,j  ,k-1)+density(i,j,k)
            RhoFace(4) = density(i+1,j  ,k  )+density(i,j,k)
            RhoFace(5) = density(i  ,j+1,k  )+density(i,j,k)
            RhoFace(6) = density(i  ,j  ,k+1)+density(i,j,k)

            Area(1) = xA(i,j,k)
            Area(2) = yA(i,j,k)
            Area(3) = zA(i,j,k)
            Area(4) = xA(i+1,j  ,k  )
            Area(5) = yA(i  ,j+1,k  )
            Area(6) = zA(i  ,j  ,k+1)

            Normal(1,1:3) = xn(i,j,k,:)
            Normal(2,1:3) = yn(i,j,k,:)
            Normal(3,1:3) = zn(i,j,k,:)
            Normal(4,1:3) = xn(i+1,j  ,k  ,:)
            Normal(5,1:3) = yn(i  ,j+1,k  ,:)
            Normal(6,1:3) = zn(i  ,j  ,k+1,:)

            gradrho_x = (-(RhoFace(1))*Normal(1,1)*Area(1) &
                         -(RhoFace(2))*Normal(2,1)*Area(2) &
                         -(RhoFace(3))*Normal(3,1)*Area(3) &
                         +(RhoFace(4))*Normal(4,1)*Area(4) &
                         +(RhoFace(5))*Normal(5,1)*Area(5) &
                         +(RhoFace(6))*Normal(6,1)*Area(6) &
                        )/(2*volume(i,j,k))

            gradrho_y = (-(RhoFace(1))*Normal(1,2)*Area(1) &
                         -(RhoFace(2))*Normal(2,2)*Area(2) &
                         -(RhoFace(3))*Normal(3,2)*Area(3) &
                         +(RhoFace(4))*Normal(4,2)*Area(4) &
                         +(RhoFace(5))*Normal(5,2)*Area(5) &
                         +(RhoFace(6))*Normal(6,2)*Area(6) &
                        )/(2*volume(i,j,k))

            gradrho_z = (-(RhoFace(1))*Normal(1,3)*Area(1) &
                         -(RhoFace(2))*Normal(2,3)*Area(2) &
                         -(RhoFace(3))*Normal(3,3)*Area(3) &
                         +(RhoFace(4))*Normal(4,3)*Area(4) &
                         +(RhoFace(5))*Normal(5,3)*Area(5) &
                         +(RhoFace(6))*Normal(6,3)*Area(6) &
                        )/(2*volume(i,j,k))


            ! __ vorticity __
            vort = sqrt(  ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                         + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                         + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                          )&
                       )

            ! ___ cross diffusion ___
            CD1 = cb2*((gradtv_x(i,j,k)*gradtv_x(i,j,k))&
                    +  (gradtv_y(i,j,k)*gradtv_y(i,j,k))&
                    +  (gradtv_z(i,j,k)*gradtv_z(i,j,k))&
                     )

            ! ___ addition cross diffusion result conservative form of tv ___
            CD2 =    ((gradrho_x*gradtv_x(i,j,k))&
                    + (gradrho_y*gradtv_y(i,j,k))&
                    + (gradrho_z*gradtv_z(i,j,k))&
                     )

            kd2  = (kappa_sa*dist(i,j,k))**2
            nu   = mu(i,j,k)/density(i,j,k)
            xi   = tv(i,j,k)/nu

            ! ___ functions ___
            fv1  = (xi**3)/((xi**3) + (cv1**3))
            fv2  = 1.0 - xi/(1.0 + (xi*fv1))

            ! ___ Shear stress for production ___
            scap = max(vort + (tv(i,j,k)*fv2/(kd2)), 0.3*vort)

            ! ___ wall function ___
            r    = min(tv(i,j,k)/(Scap*kd2), 10.0)
            g    = r + cw2*((r**6) - r)
            fw   = g*( (1.0+(cw3**6))/((g**6)+(cw3**6)) )**(1.0/6.0)

            ! ____ Dissipation term ___
            D_v = density(i,j,k)*cw1*fw*((tv(i,j,k)/dist(i,j,k))**2)

            ! ____ PRODUCTION term____
            P_v = density(i,j,k)*cb1*Scap*tv(i,j,k)

            ! ____ cross diffusion term ___
            lamda = density(i,j,k)*CD1/sigma_sa - CD2*(nu+tv(i,j,k))/sigma_sa

            S_v = (P_v - D_v  + lamda)*volume(i,j,k)
            tv_residue(i, j, k)   = tv_residue(i, j, k) - S_v

          end do
        end do
      end do

    end subroutine add_sa_source

    subroutine add_saBC_source()
      implicit none
      integer :: i,j,k

      real :: CD1
      real :: CD2
      real :: fv1
      real :: fv2
      real :: fw
      real :: g
      real :: Scap
      real :: r
      real :: vort
      real :: S_v
      real :: D_v
      real :: P_v
      real :: lamda
      real :: kd2
      real :: xi
      real :: nu
      real :: gradrho_x
      real :: gradrho_y
      real :: gradrho_z
      real, dimension(6) :: RhoFace
      real, dimension(6) :: Area
      real, dimension(6,3) :: Normal
      ! transition modeling variables
      real :: chi_1=0.002
      real :: chi_2=5.0
      real :: ReThc
      real :: ReTh
      real :: Term1
      real :: Term2
      real :: nu_BC
      real :: nu_cr
      real :: Local_vel_mag
      real :: nu_t
      real :: rey
      real :: u,v,w

      rey = density_inf*vel_mag*1.0/mu_ref

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            RhoFace(1) = density(i-1,j  ,k  )+density(i,j,k)
            RhoFace(2) = density(i  ,j-1,k  )+density(i,j,k)
            RhoFace(3) = density(i  ,j  ,k-1)+density(i,j,k)
            RhoFace(4) = density(i+1,j  ,k  )+density(i,j,k)
            RhoFace(5) = density(i  ,j+1,k  )+density(i,j,k)
            RhoFace(6) = density(i  ,j  ,k+1)+density(i,j,k)

            Area(1) = xA(i,j,k)
            Area(2) = yA(i,j,k)
            Area(3) = zA(i,j,k)
            Area(4) = xA(i+1,j  ,k  )
            Area(5) = yA(i  ,j+1,k  )
            Area(6) = zA(i  ,j  ,k+1)

            Normal(1,1:3) = xn(i,j,k,:)
            Normal(2,1:3) = yn(i,j,k,:)
            Normal(3,1:3) = zn(i,j,k,:)
            Normal(4,1:3) = xn(i+1,j  ,k  ,:)
            Normal(5,1:3) = yn(i  ,j+1,k  ,:)
            Normal(6,1:3) = zn(i  ,j  ,k+1,:)

            gradrho_x = (-(RhoFace(1))*Normal(1,1)*Area(1) &
                         -(RhoFace(2))*Normal(2,1)*Area(2) &
                         -(RhoFace(3))*Normal(3,1)*Area(3) &
                         +(RhoFace(4))*Normal(4,1)*Area(4) &
                         +(RhoFace(5))*Normal(5,1)*Area(5) &
                         +(RhoFace(6))*Normal(6,1)*Area(6) &
                        )/(2*volume(i,j,k))

            gradrho_y = (-(RhoFace(1))*Normal(1,2)*Area(1) &
                         -(RhoFace(2))*Normal(2,2)*Area(2) &
                         -(RhoFace(3))*Normal(3,2)*Area(3) &
                         +(RhoFace(4))*Normal(4,2)*Area(4) &
                         +(RhoFace(5))*Normal(5,2)*Area(5) &
                         +(RhoFace(6))*Normal(6,2)*Area(6) &
                        )/(2*volume(i,j,k))

            gradrho_z = (-(RhoFace(1))*Normal(1,3)*Area(1) &
                         -(RhoFace(2))*Normal(2,3)*Area(2) &
                         -(RhoFace(3))*Normal(3,3)*Area(3) &
                         +(RhoFace(4))*Normal(4,3)*Area(4) &
                         +(RhoFace(5))*Normal(5,3)*Area(5) &
                         +(RhoFace(6))*Normal(6,3)*Area(6) &
                        )/(2*volume(i,j,k))


            ! __ vorticity __
            vort = sqrt(  ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                         + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                         + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                          )&
                       )

            ! ___ cross diffusion ___
            CD1 = cb2*((gradtv_x(i,j,k)*gradtv_x(i,j,k))&
                    +  (gradtv_y(i,j,k)*gradtv_y(i,j,k))&
                    +  (gradtv_z(i,j,k)*gradtv_z(i,j,k))&
                     )

            ! ___ addition cross diffusion result conservative form of tv ___
            CD2 =    ((gradrho_x*gradtv_x(i,j,k))&
                    + (gradrho_y*gradtv_y(i,j,k))&
                    + (gradrho_z*gradtv_z(i,j,k))&
                     )

            kd2  = (kappa_sa*dist(i,j,k))**2
            nu   = mu(i,j,k)/density(i,j,k)
            xi   = tv(i,j,k)/nu

            ! ___ functions ___
            fv1  = (xi**3)/((xi**3) + (cv1**3))
            fv2  = 1.0 - xi/(1.0 + (xi*fv1))

            ! ___ Shear stress for production ___
            scap = max(vort + (tv(i,j,k)*fv2/(kd2)), 0.3*vort)
            !scap = max(vort + (tv(i,j,k)*fv2/(kd2)), 1.e-10)

            ! ___ wall function ___
            r    = min(tv(i,j,k)/(Scap*kd2), 10.0)
            g    = r + cw2*((r**6) - r)
            fw   = g*( (1.0+(cw3**6))/((g**6)+(cw3**6)) )**(1.0/6.0)

            ! ____ Dissipation term ___
            D_v = density(i,j,k)*cw1*fw*((tv(i,j,k)/dist(i,j,k))**2)

            ! ____ PRODUCTION term____
            P_v = density(i,j,k)*cb1*Scap*tv(i,j,k)

            ! ____ cross diffusion term ___
            lamda = density(i,j,k)*CD1/sigma_sa - CD2*(nu+tv(i,j,k))/sigma_sa

            !--- start of  BC modeling  ----
            chi_1 = 0.002
            chi_2 = 5.0
            !Local_vel_mag = sqrt(sum(qp(i,j,k,2:4)**2))
            u = qp(i,j,k,2)
            v = qp(i,j,k,3)
            w = qp(i,j,k,4)
            Local_vel_mag = sqrt(u*u + v*v + w*w)
            nu_t = tv(i,j,k)*fv1
            nu_cr = chi_2!/rey
            nu_BC = nu_t/(Local_vel_mag*dist(i,j,k))

            ReThc = 803.73*((0.18 + 0.6067)**(-1.027))
            ReTh  = (dist(i,j,k)**2)*vort/(2.193*nu)

            Term1 = sqrt(max(ReTh - ReThc, 0.0)/(chi_1*ReThc))
            Term2 = sqrt(max(nu_BC - nu_cr,0.0)/(nu_cr))

            intermittency(i,j,k) = 1.0 - exp(-(Term1+Term2))

            P_v = intermittency(i,j,k)*P_v
            ! --- end of BC modeling --- !

            S_v = (P_v - D_v  +lamda)*volume(i,j,k)
            tv_residue(i, j, k)   = tv_residue(i, j, k) - S_v

          end do
        end do
      end do

    end subroutine add_saBC_source

end module source


