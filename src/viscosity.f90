  !< Setup, destroy, calculate molecular and turbulence viscosity
module viscosity
  !< Setup, destroy, calculate molecular and turbulence viscosity
  !-----------------------------------------------------

  use vartypes
  use wall_dist  , only : dist
  use global_kkl   , only : cmu
  use global_sst   , only : bstar
  use global_sst   , only : a1
  use global_sst   , only : sst_F1
  use global_sst   , only : sigma_w2
  use global_sa    , only : cv1

  ! gradients
  use gradients, only : gradu_x
  use gradients, only : gradu_y
  use gradients, only : gradu_z 
  use gradients, only : gradv_x 
  use gradients, only : gradv_y
  use gradients, only : gradv_z
  use gradients, only : gradw_x
  use gradients, only : gradw_y
  use gradients, only : gradw_z
  use gradients, only : gradT_x
  use gradients, only : gradT_y
  use gradients, only : gradT_z
  use gradients, only : gradtk_x
  use gradients, only : gradtk_y
  use gradients, only : gradtk_z
  use gradients, only : gradtw_x
  use gradients, only : gradtw_y
  use gradients, only : gradtw_z
  use copy_bc       , only : copy1
  use utils       , only :   alloc

#include "error.inc"

  implicit none
  private
  real(wp), dimension(:, :, :), allocatable, target     :: mu
   !< Cell-center molecular viscosity
  real(wp), dimension(:, :, :), allocatable, target     :: mu_t
   !< Cell-center turbulent viscosity

  public :: setup_viscosity
  public :: calculate_viscosity
  public :: mu,mu_t

  contains

    subroutine calculate_viscosity(qp, scheme, flow, bc, dims)
      !< Calculate molecular and turbulent viscosity
      implicit none
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(boundarytype), intent(in) :: bc
      !< boundary conditions and fixed values
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      integer :: i,j,k
      real(wp) :: T ! molecular viscosity
      real(wp) :: c ! kkl eddy viscosity
      !- sst varibales -!
      real(wp)    :: F
      real(wp)    :: arg2
      real(wp)    :: vort
      real(wp)    :: NUM
      real(wp)    :: DENOM
      ! for arg2
      real(wp) :: var1
      real(wp) :: var2
      !for vorticity
      real(wp) :: wijwij
      real(wp) :: wx
      real(wp) :: wy
      real(wp) :: wz
      !for strain calculation
      real(wp) :: SijSij
      real(wp) :: Sxx, Syy, Szz
      real(wp) :: Sxy, Szx, Syz
      real(wp) :: strain
      !for arg1
      real(wp) :: arg1
      real(wp) :: CD
      real(wp) :: right
      real(wp) :: left

      ! sa variables
      real(wp) :: fv1
      real(wp) :: xi
      real(wp) :: pressure
      real(wp) :: density
      real(wp) :: tk
      real(wp) :: tkl
      real(wp) :: tw
      real(wp) :: tv
      integer :: imx, jmx, kmx

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      !--- calculate_molecular_viscosity ---!
      if (flow%mu_ref/=0.) then
        select case (trim(flow%mu_variation))
          case ('sutherland_law')
            ! apply_sutherland_law
            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx
                  pressure = qp(i,j,k,5)
                  density  = qp(i,j,k,1)
                  T = pressure/(density*flow%R_gas)
                  mu(i,j,k) = flow%mu_ref * ((T/flow%T_ref)**(1.5)) &
                            *((flow%T_ref + flow%Sutherland_temp)&
                            /(T + flow%Sutherland_temp))
                end do
              end do
            end do

          case ('constant')
            !do nothing
            !mu will be equal to mu_ref
            continue

          case DEFAULT
            print*,"mu_variation not recognized:"
            print*, "   found '",trim(flow%mu_variation),"'"
            print*, "accepted values: 1) sutherland_law"
            print*, "                 2) constant"
            Fatal_error
        end select
      end if
      !--- end molecular viscosity calculation---!

      !--- calculate_turbulent_viscosity  ---!
      if (scheme%turbulence/='none') then
        select case (trim(scheme%turbulence))

          case ('none')
            !do nothing
            continue

          case ('sa', 'saBC')
            !call calculate_sa_mu()
            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx
                  tv = qp(i,j,k,6)
                  density = qp(i,j,k,1)
                  ! xsi 
                   xi = tv*density/mu(i,j,k)
                  !calculation fo fv1 function
                  fv1 = (xi**3)/((xi**3) + (cv1**3))
                  mu_t(i,j,k) = density*tv*fv1
                end do
              end do
            end do

            ! populating ghost cell
            do i = 1,6
              select case(bc%id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-1,-2,-3,-4,-6,-7,-8,-9)
                  !call copy1(sa_mu, "symm", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = mu_t(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = mu_t(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = mu_t(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select

                case(-5)
                  !call copy1(sa_mu, "anti", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = -mu_t(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = -mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = -mu_t(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = -mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = -mu_t(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = -mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
              end select
            end do
            !--- end of sa eddy viscosity  ---!

          case ('sst2003')
            !call calculate_sst_mu()
            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx
                  density = qp(i,j,k,1)
                  tk = qp(i,j,k,6)
                  tw = qp(i,j,k,7)

                  ! calculate_arg2()
                  var1 = sqrt(tk)/(bstar*tw*dist(i,j,k))
                  var2 = 500*(mu(i,j,k)/density)/((dist(i,j,k)**2)*tw)
                  arg2 = max(2*var1, var2)

                  ! calculate_f2()
                  F = tanh(arg2**2)

                  ! calculate_vorticity(
                  sxx = (gradu_x(i,j,k))
                  syy = (gradv_y(i,j,k))
                  szz = (gradw_z(i,j,k))
                  syz = (gradw_y(i,j,k) + gradv_z(i,j,k))
                  szx = (gradu_z(i,j,k) + gradw_x(i,j,k))
                  sxy = (gradv_x(i,j,k) + gradu_y(i,j,k))

                  SijSij = (2.0*(sxx**2)) + (2.0*(syy**2)) + (2.0*(szz**2)) + syz**2 + szx**2 + sxy**2

                  strain = sqrt(SijSij)

                  NUM = density*a1*tk
                  DENOM = max(max((a1*tw), strain*F),1.0e-10)
                  mu_t(i,j,k) = NUM/DENOM
                  !-- end eddy visocisyt calculation --!
                  !-- calculating blending function F1 --!
                  CD = max(2*density*sigma_w2*(                             & 
                                                      gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                                    + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                                    + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                                     )/tw,                  &
                           1.0e-10)

                  right = 4*(density*sigma_w2*tk)/(CD*(dist(i,j,k)**2))
                  left = max(var1, var2)
                  arg1 = min(left, right)
                  sst_F1(i,j,k) = tanh(arg1**4)
                  !-- end of blending function F1 calculation --!
                end do
              end do
            end do

            select case(trim(scheme%transition))
              case('lctm2015')
                do k = 0,kmx
                  do j = 0,jmx
                    do i = 0,imx
                      !modified blending function (Menter 2015)
                      var1 = density*dist(i,j,k)*sqrt(tk)/mu(i,j,k)
                      var2 = exp(-(var1/120)**8)
                      sst_F1(i,j,k) = max(sst_F1(i,j,k),var2)
                    end do
                  end do
                end do
              case DEFAULT
                !do nothing
                continue
            end select

            ! populating ghost cell
            do i = 1,6
              select case(bc%id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-1,-2,-3,-4,-6,-7,-8,-9)
                  !call copy1(sst_mu, "symm", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = mu_t(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) = sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) = sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = mu_t(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) = sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) = sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = mu_t(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) = sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                        sst_F1(1:imx-1, 1:jmx-1,   kmx  ) = sst_F1(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
                case(-5)
                  !call copy1(sst_mu, "anti", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = -mu_t(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) =  sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = -mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) =  sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = -mu_t(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) =  sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = -mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) =  sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = -mu_t(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) =  sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = -mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                        sst_F1(1:imx-1, 1:jmx-1,   kmx  ) =  sst_F1(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
              end select
            end do
            !--- end of sst2003 eddy viscosity  and blending fucntion calculation ---!

          case ('sst')
            !call calculate_sst_mu()
            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx

                  density = qp(i,j,k,1)
                  tk = qp(i,j,k,6)
                  tw = qp(i,j,k,7)
                  ! calculate_arg2()
                  var1 = sqrt(tk)/(bstar*tw*dist(i,j,k))
                  var2 = 500*(mu(i,j,k)/density)/((dist(i,j,k)**2)*tw)
                  arg2 = max(2*var1, var2)

                  ! calculate_f2()
                  F = tanh(arg2**2)

                  ! calculate_vorticity(
                  wx = (gradw_y(i,j,k) - gradv_z(i,j,k))
                  wy = (gradu_z(i,j,k) - gradw_x(i,j,k))
                  wz = (gradv_x(i,j,k) - gradu_y(i,j,k))

                  wijwij = wx**2 + wy**2 + wz**2

                  vort = sqrt(wijwij)

                  NUM = density*a1*tk
                  DENOM = max(max((a1*tw), vort*F),1.e-20)
                  mu_t(i,j,k) = NUM/DENOM
                  !-- end eddy visocisyt calculation --!
                  !-- calculating blending function F1 --!
                  CD = max(2*density*sigma_w2*(                             & 
                                                      gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                                    + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                                    + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                                     )/tw,                  &
                           1e-20)

                  right = 4*(density*sigma_w2*tk)/(CD*(dist(i,j,k)**2))
                  left = max(var1, var2)
                  arg1 = min(left, right)
                  sst_F1(i,j,k) = tanh(arg1**4)
                  !-- end of blending function F1 calculation --!
                end do
              end do
            end do

            select case(trim(scheme%transition))
              case('lctm2015')
                do k = 0,kmx
                  do j = 0,jmx
                    do i = 0,imx
                      !modified blending function (Menter 2015)
                      var1 = density*dist(i,j,k)*sqrt(tk)/mu(i,j,k)
                      var2 = exp(-(var1/120)**8)
                      sst_F1(i,j,k) = max(sst_F1(i,j,k),var2)
                    end do
                  end do
                end do
              case DEFAULT
                !do nothing
                continue
            end select

            ! populating ghost cell
            do i = 1,6
              select case(bc%id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-1,-2,-3,-4,-6,-7,-8,-9)
                  !call copy1(sst_mu, "symm", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = mu_t(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) = sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) = sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = mu_t(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) = sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) = sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = mu_t(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) = sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                        sst_F1(1:imx-1, 1:jmx-1,   kmx  ) = sst_F1(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
                case(-5)
                  !call copy1(sst_mu, "anti", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = -mu_t(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) =  sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = -mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) =  sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = -mu_t(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) =  sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = -mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) =  sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = -mu_t(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) =  sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = -mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                        sst_F1(1:imx-1, 1:jmx-1,   kmx  ) =  sst_F1(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
              end select
            end do
            !--- end of sst eddy viscosity  and blending fucntion calculation ---!


          case ('kkl')
            !--- calculate_kkl_mu()
            c = cmu**0.25

            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx
                  density = qp(i,j,k,1)
                  tk  = qp(i,j,k,6)
                  tkl = qp(i,j,k,7)
                  mu_t(i,j,k) = c*density*tkl/(max(sqrt(tk),1.e-20))
                  if(tkl<1.e-14 .or. tk<1.e-14) &
                    mu_t(i,j,k) =0.0 
                end do
              end do
            end do

            ! populating ghost cell
            do i = 1,6
              select case(bc%id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-4:-1,-6,-8,-9)
                  !call copy1(kkl_mu, "symm", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = mu_t(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = mu_t(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = mu_t(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
                case(-5)
                  !call copy1(kkl_mu, "anti", face_names(i))
                  select case(bc%face_names(i))
                    case("imin")
                        mu_t(      0, 1:jmx-1, 1:kmx-1) = -mu_t(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        mu_t(  imx  , 1:jmx-1, 1:kmx-1) = -mu_t( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        mu_t(1:imx-1,       0, 1:kmx-1) = -mu_t(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        mu_t(1:imx-1,   jmx  , 1:kmx-1) = -mu_t(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        mu_t(1:imx-1, 1:jmx-1,       0) = -mu_t(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        mu_t(1:imx-1, 1:jmx-1,   kmx  ) = -mu_t(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
              end select
            end do
            !--- end of kkl eddy viscosity calculation ---!

          case DEFAULT 
            Fatal_error

        end select
      end if
      !--- end turbulent viscosity calculation---!
      !--- check on viscosity ---!
      if(any(isnan(mu))) then
        Fatal_error
      end if

    end subroutine calculate_viscosity

    subroutine setup_viscosity(scheme,flow, dims)
      !< Allocate and pointer for molecular and turbulent viscosity
      implicit none
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      integer :: imx, jmx, kmx

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx
      !setup_molecular_viscosity()
      if (flow%mu_ref/=0.) then
        call alloc(mu, -2, imx+2, -2, jmx+2, -2, kmx+2)
        mu = flow%mu_ref !intialize
      end if

      !--- setup_turbulent_viscosity ---!
      if (scheme%turbulence/='none') then
        call alloc(mu_t, -2,imx+2, -2,jmx+2, -2,kmx+2)


        select case (trim(scheme%turbulence))

          case ('none', 'sa', 'saBC', 'kkl')
            !do nothing
            continue

          case ('sst', 'sst2003')
            !-- sst blending funciton F1 --!
            call alloc(sst_F1, -2,imx+2, -2,jmx+2, -2,kmx+2)
            sst_F1=0.
            !-- sst blnding function setup compete--!

          case DEFAULT 
            Fatal_error

        end select
      end if
      ! --- end turbulent viscosity setup ---!

    end subroutine setup_viscosity


end module viscosity
