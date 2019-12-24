  !< Add source's contribution to the residual
module source
  !< Add source's contribution to the residual
  !-----------------------------------------------------------------
#include "debug.h"
#include "error.h"
  use vartypes
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
  use viscosity  ,only :   mu
  use wall_dist  ,only :   dist
  use gradients  ,only :   gradu_x
  use gradients  ,only :   gradu_y
  use gradients  ,only :   gradu_z
  use gradients  ,only :   gradv_x
  use gradients  ,only :   gradv_y
  use gradients  ,only :   gradv_z
  use gradients  ,only :   gradw_x
  use gradients  ,only :   gradw_y
  use gradients  ,only :   gradw_z
  use gradients  ,only :   gradtk_x
  use gradients  ,only :   gradtk_y
  use gradients  ,only :   gradtk_z
  use gradients  ,only :   gradtw_x
  use gradients  ,only :   gradtw_y
  use gradients  ,only :   gradtw_z
  use gradients  ,only :   gradtgm_x
  use gradients  ,only :   gradtgm_y
  use gradients  ,only :   gradtgm_z

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
  use global_sa , only : cb1
  use global_sa , only : cb2
  use global_sa , only : cw1
  use global_sa , only : cw2
  use global_sa , only : cw3
  use global_sa , only : cv1
  use global_sa , only : sigma_sa
  use global_sa , only : kappa_sa
  use global_sa , only : cv1_3
  use global_sa , only : cw3_6
  use gradients ,only : gradtv_x
  use gradients ,only : gradtv_y
  use gradients ,only : gradtv_z
  use viscosity  ,only : mu_t
  use CC, only : DCCVnX
  use CC, only : DCCVnY
  use CC, only : DCCVnZ
  use CC, only : CCnormalX
  use CC, only : CCnormalY
  use CC, only : CCnormalZ

  use CC         , only : find_DCCVn
  use utils,       only: alloc

  implicit none
  private

  public :: add_source_term_residue
!  public :: Setup_source
!  public :: destroy_source

  contains

    
    subroutine add_source_term_residue(qp, residue, cells, Ifaces,Jfaces,Kfaces,scheme,flow, dims)
      !< Call to add different source terms to the residual of different equations.

      implicit none
      type(schemetype), intent(in) :: scheme
      type(flowtype), intent(in) :: flow
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal

      DebugCall('add_source_term_residue')

      select case (trim(scheme%turbulence))

        case ('none')
          !do nothing
          continue

        case ('sa')
          select case(trim(scheme%transition))
            case('none')
              call add_sa_source(qp, residue, cells, Ifaces,Jfaces,Kfaces,dims)
            case('bc')
              call add_saBC_source(qp, residue, cells, Ifaces,Jfaces,Kfaces,flow, dims)
            case DEFAULT
              Fatal_error
          end select

        case ('sst', 'sst2003')
          select case(trim(scheme%transition))
            case('none')
              call add_sst_source(qp, residue, cells,scheme, dims)
            case('lctm2015')
              call add_sst_source_lctm2015(qp, residue, cells,Ifaces,Jfaces,Kfaces,scheme, dims)
            case('bc')
              call add_sst_bc_source(qp, residue, cells,flow,dims)
            case DEFAULT
              Fatal_error
          end select

        case ('kkl')
          call add_kkl_source(qp, residue, cells, Ifaces,Jfaces,Kfaces, dims)

        case DEFAULT
          Fatal_error

      end select

    end subroutine add_source_term_residue


!    subroutine Setup_source()
!      !< Allcoate memory to the required by the variable
!        implicit none
!        !nothing
!    end subroutine Setup_source
!
!
!    subroutine destroy_source()
!      !< deallocate memory before stoping the solver
!        implicit none
!        !nothing
!    end subroutine destroy_source

    subroutine add_sst_source(qp, residue, cells,scheme,dims)
      !< Add residual due to source terms of the SST turbulence model
      implicit none
      type(extent), intent(in) :: dims
      type(schemetype), intent(in) :: scheme
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      !< Store residue at each cell-center
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
      integer :: limiter
      real :: divergence
      real :: density
      real :: tk
      real :: tw
      
      if(trim(scheme%turbulence) == 'sst2003')then
        limiter = 10
        gama1 = 5.0/9.0
        gama2 = 0.44
      else 
        limiter = 20
      end if


      do k = 1,dims%kmx-1
        do j = 1,dims%jmx-1
          do i = 1,dims%imx-1

            density = qp(i,j,k,1)
            tk      = qp(i,j,k,6)
            tw      = qp(i,j,k,7)

            ! __ vorticity __
            vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                            + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                            + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                             )&
                       )

            CD = 2*density*sigma_w2*(gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                   + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                   + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                    )/tw
            CD = max(CD, 10.0**(-limiter))
            F1 = sst_F1(i,j,k)


            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w1*(1. - F1)
            gama        =       gama1*F1  +       gama2*(1. - F1)
            beta        =       beta1*F1  +       beta2*(1. - F1)



            ! ____ Dissipation term ___
            D_k = bstar*density*tw*tk
            D_w = beta*density*tw**2

            ! ____ PRODUCTION term____ 
            divergence = gradu_x(i,j,k) + gradv_y(i,j,k) + gradw_z(i,j,k)
            P_k = mu_t(i,j,k)*(vort**2) -((2.0/3.0)*density*tk*divergence)
            P_k = min(P_k,limiter*D_k)
            P_w = (density*gama/mu_t(i,j,k))*P_k

            ! ____ cross diffusion term ___
            lamda = (1. - F1)*CD

            S_k = P_k - D_k           !Source term TKE
            S_w = P_w - D_w  +lamda   !source term omega

            S_k = S_k * cells(i, j, k)%volume
            S_w = S_w * cells(i, j, k)%volume

            residue(i, j, k, 6) = residue(i, j, k, 6) - S_k
            residue(i, j, k, 7) = residue(i, j, k, 7) - S_w

          end do
        end do
      end do

    end subroutine add_sst_source


    subroutine add_sst_source_lctm2015(qp, residue, cells, Ifaces, Jfaces, Kfaces, scheme, dims)
      !< Add residual due to source terms of the LCTM2015 transition model
      implicit none
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      integer :: i,j,k

      real :: CD
      real :: F1
      real :: vort
      real :: S_K
      real :: S_w
      real :: S_gm
      real :: D_k
      real :: D_w
      real :: D_gm
      real :: P_k
      real :: P_w
      real :: P_gm
      real :: lamda
      real :: Fonset1
      real :: Fonset2
      real :: Fonset3
      real :: Fonset
      real :: Rev
      Real :: RT
      real :: Fturb
      real :: Re_theta
      real :: TuL
      real :: strain
      real :: intermittency
      real :: Pk_lim
      real :: Fon_lim
      real :: lamd
      real :: Fpg
      real :: divergence
      real :: dvdy
      integer :: limiter
      real :: density
      real :: tk
      real :: tw

      if(trim(scheme%turbulence) == 'sst2003')then
        limiter = 10
        gama1 = 5.0/9.0
        gama2 = 0.44
      else 
        limiter = 20
      end if

      !for pressure gradient calculation
      call find_DCCVn(qp, cells, Ifaces, Jfaces, Kfaces, dims)

      do k = 1,dims%kmx-1
        do j = 1,dims%jmx-1
          do i = 1,dims%imx-1

            density       = qp(i,j,k,1)
            tk            = qp(i,j,k,6)
            tw            = qp(i,j,k,7)
            intermittency = qp(i,j,k,8)

            ! __ vorticity __
            vort = sqrt(     ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                            + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                            + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                             )&
                       )

            strain = sqrt(     (((gradw_y(i,j,k) + gradv_z(i,j,k))**2) &
                              + ((gradu_z(i,j,k) + gradw_x(i,j,k))**2) &
                              + ((gradv_x(i,j,k) + gradu_y(i,j,k))**2) &
                              + 2*(gradu_x(i,j,k)**2) &
                              + 2*(gradv_y(i,j,k)**2) &
                              + 2*(gradw_z(i,j,k)**2) &
                               )&
                         )

            CD = 2*density*sigma_w2*(gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                   + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                   + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                    )/tw
            CD = max(CD, 10.0**(-limiter))
            F1 = sst_F1(i,j,k)


            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w1*(1. - F1)
            gama        =       gama1*F1  +       gama2*(1. - F1)
            beta        =       beta1*F1  +       beta2*(1. - F1)



            ! ____ Dissipation term ___
            D_k = bstar*density*tw*tk
            D_w = beta*density*tw**2

            ! ____ PRODUCTION term____ 
            divergence = gradu_x(i,j,k) + gradv_y(i,j,k) + gradw_z(i,j,k)
            P_k = mu_t(i,j,k)*(vort*strain) - ((2.0/3.0)*density*tk*divergence)
            P_k = min(P_k, limiter*D_k)
            P_w = (density*gama/mu_t(i,j,k))*P_k

            ! ____ cross diffusion term ___
            lamda = (1. - F1)*CD

            ! ____Transition modeling  ____
              ! --pressure gradient 
            dvdy = DCCVnX(i,j,k)*CCnormalX(i,j,k) &
                 + DCCVnY(i,j,k)*CCnormalY(i,j,k) &
                 + DCCVnZ(i,j,k)*CCnormalZ(i,j,k)
            lamd =(-7.57e-3)*(dvdy*dist(i,j,k)*dist(i,j,k)*density/mu(i,j,k)) + 0.0128
            lamd = min(max(lamd, -1.0), 1.0)
            if(lamd>=0.0)then
                Fpg = min(1.0 + 14.68*lamd, 1.5)
            else
                Fpg = min(1.0 - 7.34*lamd, 3.0)
            end if
            Fpg = max(Fpg, 0.0)
              ! --gradient
            TuL = min(100.0*sqrt(2.0*tk/3.0)/(tw*dist(i,j,k)),100.0)
            Re_theta = 100.0 + 1000.0*exp(-TuL*Fpg)
            !Re_theta = 100.0 + 1000.0*exp(-TuL)
            Rev = density*dist(i,j,k)*dist(i,j,k)*strain/mu(i,j,k)
            RT = density*tk/(mu(i,j,k)*tw)
            Fturb = exp(-(0.5*Rt)**4)
            Fonset1 = Rev/(2.2*Re_theta)
            Fonset2 = min(Fonset1, 2.0)
            Fonset3 = max(1.0 - (RT/3.5)**3, 0.0)
            Fonset  = max(Fonset2 - Fonset3, 0.0)
            P_gm = 100*density*strain*intermittency*(1.0 - intermittency)*Fonset
            D_gm = 0.06*density*vort*intermittency*Fturb*((50.0*intermittency) - 1.0)

            Fon_lim = min(max((Rev/(2.2*1100.0))-1.0, 0.0), 3.0)
            Pk_lim = 5*max(intermittency - 0.2,0.0)*(1.0 - intermittency)*Fon_lim*max(3*mu(i,j,k) - mu_t(i,j,k), 0.0)*strain*vort
            S_k = intermittency*P_k - max(intermittency,0.1)*D_k  +  Pk_lim      !Source term gm
            S_W = P_w - D_w  + lamda        !Source term gm
            S_gm = P_gm - D_gm           !Source term gm

            S_k = S_k * cells(i, j, k)%volume
            S_w = S_w * cells(i, j, k)%volume
            S_gm= S_gm* cells(i, j, k)%Volume

            residue(i, j, k, 6) = residue(i, j, k, 6) - S_k
            residue(i, j, k, 7) = residue(i, j, k, 7) - S_w
            residue(i, j, k, 8) = residue(i, j, k, 8) - S_gm

          end do
        end do
      end do

    end subroutine add_sst_source_lctm2015


    ! SST-BC model
    subroutine add_sst_bc_source(qp, residue, cells, flow, dims)
      !< Add residual due to source terms of the SST-BC transition model
      implicit none
      type(extent), intent(in) :: dims
      type(flowtype), intent(in) :: flow
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
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
      real :: TuL
      !--------BC model -----
      real :: chi_1=0.002
      real :: chi_2=5.0
      real :: nu_BC
      real :: nu_cr
      real :: nu_t
      real :: re_v
      real :: re_theta
      real :: re_theta_t
      real :: term1
      real :: term2
      real :: term_exponential
      real :: gamma_BC
      real :: vmag
      real :: density
      real :: tk
      real :: tw


      do k = 1,dims%kmx-1
        do j = 1,dims%jmx-1
          do i = 1,dims%imx-1

            density = qp(i,j,k,1)
            tk      = qp(i,j,k,6)
            tw      = qp(i,j,k,7)
            ! __ vorticity __
            vort = sqrt( ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
                        + (gradu_z(i,j,k)- gradw_x(i,j,k))**2 &
                        + (gradv_x(i,j,k)- gradu_y(i,j,k))**2 &
                         )&
                       )

            CD = 2*density*sigma_w2*(gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                   + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                   + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                    )/tw
            !CD = max(CD, 1e-20)
            F1 = sst_F1(i,j,k)


            sigma_k     =    sigma_k1*F1  +    sigma_k2*(1. - F1)
            sigma_w     =    sigma_w1*F1  +    sigma_w1*(1. - F1)
            gama        =       gama1*F1  +       gama2*(1. - F1)
            beta        =       beta1*F1  +       beta2*(1. - F1)



            ! ____ Dissipation term ___
            D_k = bstar*density*tw*tk
            D_w = beta*density*tw**2

            ! ____ PRODUCTION term____ 
            P_k = mu_t(i,j,k)*(vort**2)
            P_k = min(P_k,20.0*D_k)
            P_w = (density*gama/mu_t(i,j,k))*P_k

            ! ____ cross diffusion term ___
            lamda = (1. - F1)*CD

            ! ____Transition modeling  ____
            !------ BC model ---
            vmag = sqrt(SUM(qp(i,j,k,2:4)**2))
            chi_1 = 0.002
            chi_2 = 5.0

            nu_t = mu_t(i,j,k)/density
            nu_cr = chi_2/flow%Reynolds_number
            nu_bc = nu_t/(vmag*dist(i,j,k))

            TuL = flow%tu_inf !local turbulence intensity might not work for BC model
            re_v = density*dist(i,j,k)*dist(i,j,k)*vort/mu(i,j,k)
            re_theta = re_v/2.193
            re_theta_t = (803.73*((TuL + 0.6067)**(-1.027)))

            term1 = sqrt(max(re_theta-re_theta_t,0.)/(chi_1*re_theta_t))
            term2 = sqrt(max(nu_BC-nu_cr,0.0)/nu_cr)
            term_exponential = (term1 + term2)
            gamma_BC = 1.0 - exp(-term_exponential)

            P_k = gamma_BC*P_k

            S_k = P_k - D_k           !Source term TKE
            S_w = P_w - D_w  +lamda   !source term omega

            S_k = S_k * cells(i, j, k)%volume
            S_w = S_w * cells(i, j, k)%volume

            residue(i, j, k, 6) = residue(i, j, k, 6) - S_k
            residue(i, j, k, 7) = residue(i, j, k, 7) - S_w

          end do
        end do
      end do

    end subroutine add_sst_bc_source


    subroutine add_kkl_source(qp, residue, cells, Ifaces,Jfaces,Kfaces, dims)
      !< Add residual due to source terms of the k-kL turbulence model
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
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
      real :: density
      real :: tk
      real :: tkl


      do k = 1,dims%kmx-1
        do j = 1,dims%jmx-1
          do i = 1,dims%imx-1

            density = qp(i,j,k,1)
            tk      = qp(i,j,k,6)
            tkl     = qp(i,j,k,7)

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

            Tau11 = mu_t(i,j,k)*(2*S11 - (2.0/3.0)*delv) - (2.0/3.0)*density*tk
            Tau12 = mu_t(i,j,k)*(2*S12)
            Tau13 = mu_t(i,j,k)*(2*S13)
            Tau21 = mu_t(i,j,k)*(2*S21)
            Tau22 = mu_t(i,j,k)*(2*S22 - (2.0/3.0)*delv) - (2.0/3.0)*density*tk
            Tau23 = mu_t(i,j,k)*(2*S23)
            Tau31 = mu_t(i,j,k)*(2*S31)
            Tau32 = mu_t(i,j,k)*(2*S32)
            Tau33 = mu_t(i,j,k)*(2*S33 - (2.0/3.0)*delv) - (2.0/3.0)*density*tk

            P_k = 0.
            P_k = P_k + Tau11*gradu_x(i,j,k) + Tau12*gradu_y(i,j,k) + Tau13*gradu_z(i,j,k)
            P_k = P_k + Tau21*gradv_x(i,j,k) + Tau22*gradv_y(i,j,k) + Tau23*gradv_z(i,j,k)
            P_k = P_k + Tau31*gradw_x(i,j,k) + Tau32*gradw_y(i,j,k) + Tau33*gradw_z(i,j,k)
            D_k = (cmu**0.75)*density*(tk**2.5)/max(tkl,1.e-20)
            P_k = min(P_k, 20*D_k)

            ! calculation of Lvk
            ! first get second order gradients 
            d2udx2 =(-(gradu_x(i-1,j  ,k  )+gradu_x(i,j,k))*Ifaces(i,j,k)%nx*Ifaces(i,j,k)%A &
                     -(gradu_x(i  ,j-1,k  )+gradu_x(i,j,k))*Jfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                     -(gradu_x(i  ,j  ,k-1)+gradu_x(i,j,k))*Kfaces(i,j,k)%nx*Kfaces(i,j,k)%A &
                     +(gradu_x(i+1,j  ,k  )+gradu_x(i,j,k))*Ifaces(i+1,j  ,k  )%nx*Ifaces(i+1,j  ,k  )%A &
                     +(gradu_x(i  ,j+1,k  )+gradu_x(i,j,k))*Jfaces(i  ,j+1,k  )%nx*Jfaces(i  ,j+1,k  )%A &
                     +(gradu_x(i  ,j  ,k+1)+gradu_x(i,j,k))*Kfaces(i  ,j  ,k+1)%nx*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            d2udy2 =(-(gradu_y(i-1,j  ,k  )+gradu_y(i,j,k))*Ifaces(i,j,k)%ny*Ifaces(i,j,k)%A &
                     -(gradu_y(i  ,j-1,k  )+gradu_y(i,j,k))*Jfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                     -(gradu_y(i  ,j  ,k-1)+gradu_y(i,j,k))*Kfaces(i,j,k)%ny*Kfaces(i,j,k)%A &
                     +(gradu_y(i+1,j  ,k  )+gradu_y(i,j,k))*Ifaces(i+1,j  ,k  )%ny*Ifaces(i+1,j  ,k  )%A &
                     +(gradu_y(i  ,j+1,k  )+gradu_y(i,j,k))*Jfaces(i  ,j+1,k  )%ny*Jfaces(i  ,j+1,k  )%A &
                     +(gradu_y(i  ,j  ,k+1)+gradu_y(i,j,k))*Kfaces(i  ,j  ,k+1)%ny*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            d2udz2 =(-(gradu_z(i-1,j  ,k  )+gradu_z(i,j,k))*Ifaces(i,j,k)%nz*Ifaces(i,j,k)%A &
                     -(gradu_z(i  ,j-1,k  )+gradu_z(i,j,k))*Jfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                     -(gradu_z(i  ,j  ,k-1)+gradu_z(i,j,k))*Kfaces(i,j,k)%nz*Kfaces(i,j,k)%A &
                     +(gradu_z(i+1,j  ,k  )+gradu_z(i,j,k))*Ifaces(i+1,j  ,k  )%nz*Ifaces(i+1,j  ,k  )%A &
                     +(gradu_z(i  ,j+1,k  )+gradu_z(i,j,k))*Jfaces(i  ,j+1,k  )%nz*Jfaces(i  ,j+1,k  )%A &
                     +(gradu_z(i  ,j  ,k+1)+gradu_z(i,j,k))*Kfaces(i  ,j  ,k+1)%nz*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            ! gradient of v component
            d2vdx2 =(-(gradv_x(i-1,j  ,k  )+gradv_x(i,j,k))*Ifaces(i,j,k)%nx*Ifaces(i,j,k)%A &
                     -(gradv_x(i  ,j-1,k  )+gradv_x(i,j,k))*Jfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                     -(gradv_x(i  ,j  ,k-1)+gradv_x(i,j,k))*Kfaces(i,j,k)%nx*Kfaces(i,j,k)%A &
                     +(gradv_x(i+1,j  ,k  )+gradv_x(i,j,k))*Ifaces(i+1,j  ,k  )%nx*Ifaces(i+1,j  ,k  )%A &
                     +(gradv_x(i  ,j+1,k  )+gradv_x(i,j,k))*Jfaces(i  ,j+1,k  )%nx*Jfaces(i  ,j+1,k  )%A &
                     +(gradv_x(i  ,j  ,k+1)+gradv_x(i,j,k))*Kfaces(i  ,j  ,k+1)%nx*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            d2vdy2 =(-(gradv_y(i-1,j  ,k  )+gradv_y(i,j,k))*Ifaces(i,j,k)%ny*Ifaces(i,j,k)%A &
                     -(gradv_y(i  ,j-1,k  )+gradv_y(i,j,k))*Jfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                     -(gradv_y(i  ,j  ,k-1)+gradv_y(i,j,k))*Kfaces(i,j,k)%ny*Kfaces(i,j,k)%A &
                     +(gradv_y(i+1,j  ,k  )+gradv_y(i,j,k))*Ifaces(i+1,j  ,k  )%ny*Ifaces(i+1,j  ,k  )%A &
                     +(gradv_y(i  ,j+1,k  )+gradv_y(i,j,k))*Jfaces(i  ,j+1,k  )%ny*Jfaces(i  ,j+1,k  )%A &
                     +(gradv_y(i  ,j  ,k+1)+gradv_y(i,j,k))*Kfaces(i  ,j  ,k+1)%ny*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            d2vdz2 =(-(gradv_z(i-1,j  ,k  )+gradv_z(i,j,k))*Ifaces(i,j,k)%nz*Ifaces(i,j,k)%A &
                     -(gradv_z(i  ,j-1,k  )+gradv_z(i,j,k))*Jfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                     -(gradv_z(i  ,j  ,k-1)+gradv_z(i,j,k))*Kfaces(i,j,k)%nz*Kfaces(i,j,k)%A &
                     +(gradv_z(i+1,j  ,k  )+gradv_z(i,j,k))*Ifaces(i+1,j  ,k  )%nz*Ifaces(i+1,j  ,k  )%A &
                     +(gradv_z(i  ,j+1,k  )+gradv_z(i,j,k))*Jfaces(i  ,j+1,k  )%nz*Jfaces(i  ,j+1,k  )%A &
                     +(gradv_z(i  ,j  ,k+1)+gradv_z(i,j,k))*Kfaces(i  ,j  ,k+1)%nz*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)


            !gradients of w components
            d2wdx2 =(-(gradw_x(i-1,j  ,k  )+gradw_x(i,j,k))*Ifaces(i,j,k)%nx*Ifaces(i,j,k)%A &
                     -(gradw_x(i  ,j-1,k  )+gradw_x(i,j,k))*Jfaces(i,j,k)%nx*Jfaces(i,j,k)%A &
                     -(gradw_x(i  ,j  ,k-1)+gradw_x(i,j,k))*Kfaces(i,j,k)%nx*Kfaces(i,j,k)%A &
                     +(gradw_x(i+1,j  ,k  )+gradw_x(i,j,k))*Ifaces(i+1,j  ,k  )%nx*Ifaces(i+1,j  ,k  )%A &
                     +(gradw_x(i  ,j+1,k  )+gradw_x(i,j,k))*Jfaces(i  ,j+1,k  )%nx*Jfaces(i  ,j+1,k  )%A &
                     +(gradw_x(i  ,j  ,k+1)+gradw_x(i,j,k))*Kfaces(i  ,j  ,k+1)%nx*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            d2wdy2 =(-(gradw_y(i-1,j  ,k  )+gradw_y(i,j,k))*Ifaces(i,j,k)%ny*Ifaces(i,j,k)%A &
                     -(gradw_y(i  ,j-1,k  )+gradw_y(i,j,k))*Jfaces(i,j,k)%ny*Jfaces(i,j,k)%A &
                     -(gradw_y(i  ,j  ,k-1)+gradw_y(i,j,k))*Kfaces(i,j,k)%ny*Kfaces(i,j,k)%A &
                     +(gradw_y(i+1,j  ,k  )+gradw_y(i,j,k))*Ifaces(i+1,j  ,k  )%ny*Ifaces(i+1,j  ,k  )%A &
                     +(gradw_y(i  ,j+1,k  )+gradw_y(i,j,k))*Jfaces(i  ,j+1,k  )%ny*Jfaces(i  ,j+1,k  )%A &
                     +(gradw_y(i  ,j  ,k+1)+gradw_y(i,j,k))*Kfaces(i  ,j  ,k+1)%ny*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            d2wdz2 =(-(gradw_z(i-1,j  ,k  )+gradw_z(i,j,k))*Ifaces(i,j,k)%nz*Ifaces(i,j,k)%A &
                     -(gradw_z(i  ,j-1,k  )+gradw_z(i,j,k))*Jfaces(i,j,k)%nz*Jfaces(i,j,k)%A &
                     -(gradw_z(i  ,j  ,k-1)+gradw_z(i,j,k))*Kfaces(i,j,k)%nz*Kfaces(i,j,k)%A &
                     +(gradw_z(i+1,j  ,k  )+gradw_z(i,j,k))*Ifaces(i+1,j  ,k  )%nz*Ifaces(i+1,j  ,k  )%A &
                     +(gradw_z(i  ,j+1,k  )+gradw_z(i,j,k))*Jfaces(i  ,j+1,k  )%nz*Jfaces(i  ,j+1,k  )%A &
                     +(gradw_z(i  ,j  ,k+1)+gradw_z(i,j,k))*Kfaces(i  ,j  ,k+1)%nz*Kfaces(i  ,j  ,k+1)%A &
                    )/(2*cells(i,j,k)%volume)

            udd = sqrt( (d2udx2+d2udy2+d2udz2)**2 &
                      + (d2vdx2+d2vdy2+d2vdz2)**2 &
                      + (d2wdx2+d2wdy2+d2wdz2)**2 )

            ud  = sqrt(2*(s11**2 + s12**2 + s13**2 &
                         +s21**2 + s22**2 + s23**2 &
                         +s31**2 + s32**2 + s33**2 ))

            Lvk = kappa*abs(ud/max(udd,1.e-20))

            fp = min(max(P_k/D_k, 0.5),1.0)
            ! Lvk limiter
            Lvk = max(Lvk, tkl/max((tk*c11),1.e-20))
            Lvk = min(Lvk, c12*kappa*dist(i,j,k)*fp)

            eta = density*dist(i,j,k)*sqrt(0.3*tk)/(20*mu(i,j,k))
            fphi = (1 + cd1*eta)/(1 + eta**4)
            cphi2 = zeta3
            cphi1 = (zeta1 - zeta2*((tkl/max(tk*Lvk,1.e-20))**2))

            P_kl = cphi1*tkl*P_k/max(tk,1.e-20)
            D_kl = cphi2*density*(tk**1.5)


            S_k  = P_k  - D_k  - 2*mu(i,j,k)*tk/(dist(i,j,k)**2)       !Source term TKE
            S_kl = P_kl - D_kl - 6*mu(i,j,k)*tkl*fphi/(dist(i,j,k)**2) !source term KL

            S_k  = S_k  * cells(i, j, k)%volume
            S_kl = S_kl * cells(i, j, k)%volume

            residue(i, j, k, 6) = residue(i, j, k, 6) - S_k
            residue(i, j, k, 7) = residue(i, j, k, 7) - S_kl

          end do
        end do
      end do

    end subroutine add_kkl_source


    subroutine add_sa_source(qp, residue, cells, Ifaces,Jfaces,Kfaces, dims)
      !< Add residual due to source terms of SA turbulence model
      implicit none
      type(extent), intent(in) :: dims
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
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
      real :: density
      real :: tv

      do k = 1,dims%kmx-1
        do j = 1,dims%jmx-1
          do i = 1,dims%imx-1

            density = qp(i,j,k,1)
            tv      = qp(i,j,k,6)

            RhoFace(1) = qp(i-1,j  ,k  ,1)+density
            RhoFace(2) = qp(i  ,j-1,k  ,1)+density
            RhoFace(3) = qp(i  ,j  ,k-1,1)+density
            RhoFace(4) = qp(i+1,j  ,k  ,1)+density
            RhoFace(5) = qp(i  ,j+1,k  ,1)+density
            RhoFace(6) = qp(i  ,j  ,k+1,1)+density

            Area(1) = Ifaces(i,j,k)%A
            Area(2) = Jfaces(i,j,k)%A
            Area(3) = Kfaces(i,j,k)%A
            Area(4) = Ifaces(i+1,j  ,k  )%A
            Area(5) = Jfaces(i  ,j+1,k  )%A
            Area(6) = Kfaces(i  ,j  ,k+1)%A

            Normal(1,1:3) = (/Ifaces(i,j,k)%nx,Ifaces(i,j,k)%ny,Ifaces(i,j,k)%nz/)
            Normal(2,1:3) = (/Jfaces(i,j,k)%nx,Jfaces(i,j,k)%ny,Jfaces(i,j,k)%nz/)
            Normal(3,1:3) = (/Kfaces(i,j,k)%nx,Kfaces(i,j,k)%nx,Kfaces(i,j,k)%nx/)
            Normal(4,1:3) = (/Ifaces(i+1,j  ,k)%nx,Ifaces(i+1,j  ,k)%ny,Ifaces(i+1,j  ,k)%nz/)
            Normal(5,1:3) = (/Jfaces(i  ,j+1,k)%nx,Jfaces(i  ,j+1,k)%ny,Jfaces(i  ,j+1,k)%nz/)
            Normal(6,1:3) = (/Kfaces(i  ,j,k+1)%nx,Kfaces(i  ,j,k+1)%ny,Kfaces(i  ,j,k+1)%nz/)

            gradrho_x = (-(RhoFace(1))*Normal(1,1)*Area(1) &
                         -(RhoFace(2))*Normal(2,1)*Area(2) &
                         -(RhoFace(3))*Normal(3,1)*Area(3) &
                         +(RhoFace(4))*Normal(4,1)*Area(4) &
                         +(RhoFace(5))*Normal(5,1)*Area(5) &
                         +(RhoFace(6))*Normal(6,1)*Area(6) &
                        )/(2.0*cells(i,j,k)%volume)

            gradrho_y = (-(RhoFace(1))*Normal(1,2)*Area(1) &
                         -(RhoFace(2))*Normal(2,2)*Area(2) &
                         -(RhoFace(3))*Normal(3,2)*Area(3) &
                         +(RhoFace(4))*Normal(4,2)*Area(4) &
                         +(RhoFace(5))*Normal(5,2)*Area(5) &
                         +(RhoFace(6))*Normal(6,2)*Area(6) &
                        )/(2.0*cells(i,j,k)%volume)

            gradrho_z = (-(RhoFace(1))*Normal(1,3)*Area(1) &
                         -(RhoFace(2))*Normal(2,3)*Area(2) &
                         -(RhoFace(3))*Normal(3,3)*Area(3) &
                         +(RhoFace(4))*Normal(4,3)*Area(4) &
                         +(RhoFace(5))*Normal(5,3)*Area(5) &
                         +(RhoFace(6))*Normal(6,3)*Area(6) &
                        )/(2.0*cells(i,j,k)%volume)


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
            nu   = mu(i,j,k)/density
            xi   = tv/nu

            ! ___ functions ___
            fv1  = (xi**3)/((xi**3) + (cv1**3))
            fv2  = 1.0 - xi/(1.0 + (xi*fv1))

            ! ___ Shear stress for production ___
            scap = max(vort + (tv*fv2/(kd2)), 0.3*vort)

            ! ___ wall function ___
            r    = min(tv/(Scap*kd2), 10.0)
            g    = r + cw2*((r**6) - r)
            fw   = g*( (1.0+(cw3**6))/((g**6)+(cw3**6)) )**(1.0/6.0)

            ! ____ Dissipation term ___
            D_v = density*cw1*fw*((tv/dist(i,j,k))**2)

            ! ____ PRODUCTION term____
            P_v = density*cb1*Scap*tv

            ! ____ cross diffusion term ___
            lamda = density*CD1/sigma_sa - CD2*(nu+tv)/sigma_sa

            S_v = (P_v - D_v  + lamda)*cells(i,j,k)%volume

            residue(i, j, k, 6)   = residue(i, j, k, 6) - S_v

          end do
        end do
      end do

    end subroutine add_sa_source

    subroutine add_saBC_source(qp, residue, cells, Ifaces,Jfaces,Kfaces, flow, dims)
      !< Add residual due to source terms of SABC transition model
      implicit none
      type(extent), intent(in) :: dims
      type(flowtype), intent(in) :: flow
      real, dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real, dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
      !< Input varaible which stores I faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
      !< Input varaible which stores J faces' area and unit normal
      type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
      !< Input varaible which stores K faces' area and unit normal
      integer :: i,j,k

      real :: CD1
      real :: CD2
      real :: fv1
      real :: fv2
      real :: fw
      real :: g
      real :: r
      real :: S_v
      real :: lamda
      real :: dist_i
      real :: dist_i_2
      real :: Ji
      real :: Ji_2
      real :: Ji_3
      real :: S
      real :: Omega
      real :: k2
      real :: inv_k2_d2
      real :: Shat
      real :: inv_Shat
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
      real :: nu_BC
      real :: nu_cr
      real :: nu_t
      real :: u,v,w
      real :: glim
      real :: g_6
      real :: vmag
      real :: Production
      real :: Destruction
      real :: re_v
      real :: re_theta
      real :: re_theta_t
      real :: term1
      real :: term2
      real :: term_exponential
      real :: gamma_BC
      real :: tu
      real :: tv
      real :: density

      tu = flow%tu_inf

      do k = 1,dims%kmx-1
        do j = 1,dims%jmx-1
          do i = 1,dims%imx-1

            !Local_vel_mag
            density= qp(i,j,k,1)
            u      = qp(i,j,k,2)
            v      = qp(i,j,k,3)
            w      = qp(i,j,k,4)
            tv     = qp(i,j,k,6)
            vmag = sqrt(u*u + v*v + w*w)

            RhoFace(1) = qp(i-1,j  ,k  ,1)+density
            RhoFace(2) = qp(i  ,j-1,k  ,1)+density
            RhoFace(3) = qp(i  ,j  ,k-1,1)+density
            RhoFace(4) = qp(i+1,j  ,k  ,1)+density
            RhoFace(5) = qp(i  ,j+1,k  ,1)+density
            RhoFace(6) = qp(i  ,j  ,k+1,1)+density

            Area(1) = Ifaces(i,j,k)%A
            Area(2) = Jfaces(i,j,k)%A
            Area(3) = Kfaces(i,j,k)%A
            Area(4) = Ifaces(i+1,j  ,k  )%A
            Area(5) = Jfaces(i  ,j+1,k  )%A
            Area(6) = Kfaces(i  ,j  ,k+1)%A

            Normal(1,1:3) = (/Ifaces(i,j,k)%nx,Ifaces(i,j,k)%ny,Ifaces(i,j,k)%nz/)
            Normal(2,1:3) = (/Jfaces(i,j,k)%nx,Jfaces(i,j,k)%ny,Jfaces(i,j,k)%nz/)
            Normal(3,1:3) = (/Kfaces(i,j,k)%nx,Kfaces(i,j,k)%nx,Kfaces(i,j,k)%nx/)
            Normal(4,1:3) = (/Ifaces(i+1,j  ,k)%nx,Ifaces(i+1,j  ,k)%ny,Ifaces(i+1,j  ,k)%nz/)
            Normal(5,1:3) = (/Jfaces(i  ,j+1,k)%nx,Jfaces(i  ,j+1,k)%ny,Jfaces(i  ,j+1,k)%nz/)
            Normal(6,1:3) = (/Kfaces(i  ,j,k+1)%nx,Kfaces(i  ,j,k+1)%ny,Kfaces(i  ,j,k+1)%nz/)

            gradrho_x = (-(RhoFace(1))*Normal(1,1)*Area(1) &
                         -(RhoFace(2))*Normal(2,1)*Area(2) &
                         -(RhoFace(3))*Normal(3,1)*Area(3) &
                         +(RhoFace(4))*Normal(4,1)*Area(4) &
                         +(RhoFace(5))*Normal(5,1)*Area(5) &
                         +(RhoFace(6))*Normal(6,1)*Area(6) &
                        )/(2*cells(i,j,k)%volume)

            gradrho_y = (-(RhoFace(1))*Normal(1,2)*Area(1) &
                         -(RhoFace(2))*Normal(2,2)*Area(2) &
                         -(RhoFace(3))*Normal(3,2)*Area(3) &
                         +(RhoFace(4))*Normal(4,2)*Area(4) &
                         +(RhoFace(5))*Normal(5,2)*Area(5) &
                         +(RhoFace(6))*Normal(6,2)*Area(6) &
                        )/(2*cells(i,j,k)%volume)

            gradrho_z = (-(RhoFace(1))*Normal(1,3)*Area(1) &
                         -(RhoFace(2))*Normal(2,3)*Area(2) &
                         -(RhoFace(3))*Normal(3,3)*Area(3) &
                         +(RhoFace(4))*Normal(4,3)*Area(4) &
                         +(RhoFace(5))*Normal(5,3)*Area(5) &
                         +(RhoFace(6))*Normal(6,3)*Area(6) &
                        )/(2*cells(i,j,k)%volume)


            ! __ vorticity __
            Omega = sqrt(  ((gradw_y(i,j,k)- gradv_z(i,j,k))**2 &
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

            dist_i = dist(i,j,k)
            dist_i_2 = dist_i*dist_i
            k2 = kappa_sa*kappa_sa
            nu   = mu(i,j,k)/density
            Ji   = tv/nu
            Ji_2 = Ji*Ji
            Ji_3 = Ji_2*ji


            ! ___ functions ___
            fv1  = (Ji_3)/((Ji_3) + (cv1_3))
            fv2  = 1.0 - Ji/(1.0 + (Ji*fv1))

            ! ___ Shear stress for production ___
            S = Omega
            inv_k2_d2 = 1.0/(k2*dist_i_2)
            Shat      = S + tv*fv2*inv_k2_d2
            Shat      = max(Shat, 1.0e-10)
            inv_Shat  = 1.0/Shat

            ! ____ PRODUCTION term____
            chi_1 = 0.002
            chi_2 = 5.0

            nu_t = tv*fv1
            nu_cr = chi_2/flow%Reynolds_number
            nu_bc = nu_t/(vmag*dist_i)

            re_v = dist_i_2*Omega/nu
            re_theta = re_v/2.193
            re_theta_t = (803.73*((tu+0.6067)**(-1.027)))
            !re_theta_t = 163.0 + exp(6.91-0.18)


            term1 = sqrt(max(re_theta-re_theta_t,0.)/(chi_1*re_theta_t))
            term2 = sqrt(max(nu_BC-nu_cr,0.0)/nu_cr)
            term_exponential = (term1 + term2)
            gamma_BC = 1.0 - exp(-term_exponential)

            Production = gamma_BC*cb1*Shat*tv*cells(i,j,k)%volume

            ! ___ Destruction term___ !
            r    = min(tv*inv_Shat*inv_k2_d2, 10.0)
            g    = r + cw2*((r**6) - r)
            g_6  = g**6
            glim = ((1.0+cw3_6)/(g_6+cw3_6))**(1.0/6.0)
            fw   = g*glim
            Destruction = (cw1*fw*tv*tv/dist_i_2)*(cells(i,j,k)%volume)

            ! ____ cross diffusion term ___
            lamda = (density*CD1/sigma_sa - CD2*(nu+tv)/sigma_sa)*cells(i,j,k)%volume

            S_v = (Production - Destruction  + lamda)
            residue(i, j, k, 6)   = residue(i, j, k, 6) - S_v

          end do
        end do
      end do

    end subroutine add_saBC_source


end module source


