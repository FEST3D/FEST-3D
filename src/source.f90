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

  contains

    
    subroutine add_source_term_residue(qp, residue, cells, Ifaces,Jfaces,Kfaces,scheme,flow, dims)
      !< Call to add different source terms to the residual of different equations.

      implicit none
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
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


    subroutine add_sst_source(qp, residue, cells,scheme,dims)
      !< Add residual due to source terms of the SST turbulence model
      implicit none
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      !< Store residue at each cell-center
      integer :: i,j,k

      real(wp) :: CD
      !< cross diffusion term
      real(wp) :: F1
      !< single cell belding fuction 
      real(wp) :: vort
      !< vorticity magnitude
      real(wp) :: S_k
      !< Total source term of TKE equation
      real(wp) :: S_w
      !< Total source term of omega equation
      real(wp) :: D_k
      !< destruction term of TKE equation
      real(wp) :: D_w
      !< destruction term of omega equation
      real(wp) :: P_k
      !< production term of TKE equation
      real(wp) :: P_w
      !< production term of Omega equation
      real(wp) :: lamda
      !< additional source term in Omega equation
      integer :: limiter
      !< production term limiter
      real(wp) :: divergence
      !< del.V
      real(wp) :: density
      !< single cell density
      real(wp) :: tk
      !< single cell TKE
      real(wp) :: tw
      !< single cell Omega
      
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
      !< finite-volume Schemes
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
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

      real(wp) :: CD
      !< Cross-diffustion term
      real(wp) :: F1
      !< single cell blending function
      real(wp) :: vort
      !< vorticity magnitude
      real(wp) :: S_K
      !< Total source term in TKE equation
      real(wp) :: S_w
      !< Total source term in Omega equation
      real(wp) :: S_gm
      !< Total source term in Gamma equation
      real(wp) :: D_k
      !< Destruction term in TKE equation
      real(wp) :: D_w
      !< Destruction term in Omega equation
      real(wp) :: D_gm
      !< Destruction term in Gamma equation
      real(wp) :: P_k
      !< production term in TKE equation
      real(wp) :: P_w
      !< production term in Omega equation
      real(wp) :: P_gm
      !< production term in Gamma equation
      real(wp) :: lamda
      !< additional source term in Omega equation
      real(wp) :: Fonset1,Fonset2, Fonset3,Fonset
      !< Transition onset term 
      real(wp) :: Rev
      !< Reynodlds number based on vorticity
      real(wp) :: RT
      !< Turbulent reynolds number
      real(wp) :: Fturb
      real(wp) :: Re_theta
      !< Cutt-off reynolds number based on momentum thickness
      real(wp) :: TuL
      !< local turbulence intensity
      real(wp) :: strain
      !< Strain rate magnitude
      real(wp) :: intermittency
      !< intermittency
      real(wp) :: Pk_lim
      !< production lim term
      real(wp) :: Fon_lim
      real(wp) :: lamd
      !< pressure gradient 
      real(wp) :: Fpg
      !< pressure gradient functin
      real(wp) :: divergence
      !< del.V
      real(wp) :: dvdy
      !< pressure gradient sensor
      integer :: limiter
      !< production limiter
      real(wp) :: density
      !< single cell Density
      real(wp) :: tk
      !< single cell TKE
      real(wp) :: tw
      !< single cell Omega

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
      !< Extent of the domain:imx,jmx,kmx
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
      !< Store residue at each cell-center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      integer :: i,j,k

      real(wp) :: CD
      !< cross-diffusion term
      real(wp) :: F1
      !< single cell blending function 
      real(wp) :: vort
      !< vorticity magnitude
      real(wp) :: S_k
      !< Total source term in TKE equation
      real(wp) :: S_w
      !< Total source term in Omega equation
      real(wp) :: D_k
      !< Destruction term in TKE equation
      real(wp) :: D_w
      !< Destruction term in Omega equation
      real(wp) :: P_k
      !< Production term in TKE equation
      real(wp) :: P_w
      !< production term in Omega equation
      real(wp) :: lamda
      !< addtion source term in Omega equation
      real(wp) :: TuL
      !< local turbulence intensity
      !--------BC model -----
      real(wp) :: chi_1=0.002
      real(wp) :: chi_2=5.0
      real(wp) :: nu_BC
      real(wp) :: nu_cr
      real(wp) :: nu_t
      real(wp) :: re_v
      real(wp) :: re_theta
      real(wp) :: re_theta_t
      real(wp) :: term1
      real(wp) :: term2
      real(wp) :: term_exponential
      real(wp) :: gamma_BC
      !< intermittency function
      real(wp) :: vmag
      !< velocity magnitude
      real(wp) :: density
      !< single cell Density
      real(wp) :: tk
      !< single cell TKE
      real(wp) :: tw
      !< single cell omega


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
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
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

      real(wp) :: Tau11
      real(wp) :: Tau12
      real(wp) :: Tau13

      real(wp) :: Tau21
      real(wp) :: Tau22
      real(wp) :: Tau23

      real(wp) :: Tau31
      real(wp) :: Tau32
      real(wp) :: Tau33

      real(wp) :: S11
      real(wp) :: S12
      real(wp) :: S13

      real(wp) :: S21
      real(wp) :: S22
      real(wp) :: S23

      real(wp) :: S31
      real(wp) :: S32
      real(wp) :: S33

      real(wp) :: delv

      real(wp) :: d2udx2
      real(wp) :: d2udy2
      real(wp) :: d2udz2

      real(wp) :: d2vdx2
      real(wp) :: d2vdy2
      real(wp) :: d2vdz2

      real(wp) :: d2wdx2
      real(wp) :: d2wdy2
      real(wp) :: d2wdz2

      real(wp) :: Lvk
      real(wp) :: fp
      real(wp) :: ud
      real(wp) :: udd

      real(wp) :: S_k
      real(wp) :: S_kl
      real(wp) :: D_k
      real(wp) :: D_kl
      real(wp) :: P_k
      real(wp) :: P_kl
      real(wp) :: density
      real(wp) :: tk
      real(wp) :: tkl


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
      !< Extent of the domain:imx,jmx,kmx
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
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

      real(wp) :: CD1
      real(wp) :: CD2
      real(wp) :: fv1
      real(wp) :: fv2
      real(wp) :: fw
      real(wp) :: g
      real(wp) :: Scap
      real(wp) :: r
      real(wp) :: vort
      real(wp) :: S_v
      real(wp) :: D_v
      real(wp) :: P_v
      real(wp) :: lamda
      real(wp) :: kd2
      real(wp) :: xi
      real(wp) :: nu
      real(wp) :: gradrho_x
      real(wp) :: gradrho_y
      real(wp) :: gradrho_z
      real(wp), dimension(6) :: RhoFace
      real(wp), dimension(6) :: Area
      real(wp), dimension(6,3) :: Normal
      real(wp) :: density
      real(wp) :: tv

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
      !< Extent of the domain:imx,jmx,kmx
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      real(wp), dimension(:, :, :, :), intent(inout)  :: residue
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

      real(wp) :: CD1
      real(wp) :: CD2
      real(wp) :: fv1
      real(wp) :: fv2
      real(wp) :: fw
      real(wp) :: g
      real(wp) :: r
      real(wp) :: S_v
      real(wp) :: lamda
      real(wp) :: dist_i
      real(wp) :: dist_i_2
      real(wp) :: Ji
      real(wp) :: Ji_2
      real(wp) :: Ji_3
      real(wp) :: S
      real(wp) :: Omega
      real(wp) :: k2
      real(wp) :: inv_k2_d2
      real(wp) :: Shat
      real(wp) :: inv_Shat
      real(wp) :: nu
      real(wp) :: gradrho_x
      real(wp) :: gradrho_y
      real(wp) :: gradrho_z
      real(wp), dimension(6) :: RhoFace
      real(wp), dimension(6) :: Area
      real(wp), dimension(6,3) :: Normal
      ! transition modeling variables
      real(wp) :: chi_1=0.002
      real(wp) :: chi_2=5.0
      real(wp) :: nu_BC
      real(wp) :: nu_cr
      real(wp) :: nu_t
      real(wp) :: u,v,w
      real(wp) :: glim
      real(wp) :: g_6
      real(wp) :: vmag
      real(wp) :: Production
      real(wp) :: Destruction
      real(wp) :: re_v
      real(wp) :: re_theta
      real(wp) :: re_theta_t
      real(wp) :: term1
      real(wp) :: term2
      real(wp) :: term_exponential
      real(wp) :: gamma_BC
      real(wp) :: tu
      real(wp) :: tv
      real(wp) :: density

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


