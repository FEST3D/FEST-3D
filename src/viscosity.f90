  !< Setup, destroy, calculate molecular and turbulence viscosity
module viscosity
  !< Setup, destroy, calculate molecular and turbulence viscosity
  !----------------------------------------------------
  ! author   - Jatinder Pal Singh Sandhu
  ! obective - setup     ,
  !            destroy   , and 
  !            calculate 
  ! molecular and turbulence viscosity based on switch
  !-----------------------------------------------------

  use global_vars , only : turbulence
  use global_vars , only : transition
  use global_vars , only : mu
  use global_vars , only : mu_ref
  use global_vars , only : Sutherland_temp
  use global_vars , only : T_ref
  use global_vars , only : R_gas
  use global_vars , only : mu_variation
  use global_vars , only : process_id

  use global_vars , only : pressure
  use global_vars , only : density
  use global_vars  , only : sst_mu
  use global_vars  , only : kkl_mu
  use global_vars  , only : sa_mu
  use global_vars  , only : mu_t
  use global_vars  , only : tw
  use global_vars  , only : tk
  use global_vars  , only : tkl
  use global_vars  , only : tv

  use global_vars , only : imx
  use global_vars , only : jmx
  use global_vars , only : kmx

  use global_vars  , only : id
  use global_vars  , only : face_names
  use global_vars  , only : dist
  use global_kkl   , only : cmu
  use global_sst   , only : bstar
  use global_sst   , only : a1
  use global_sst   , only : sst_F1
  use global_sst   , only : sigma_w2
  use global_sa    , only : cv1

  ! gradients
  use global_vars, only : gradu_x
  use global_vars, only : gradu_y
  use global_vars, only : gradu_z 
  use global_vars, only : gradv_x 
  use global_vars, only : gradv_y
  use global_vars, only : gradv_z
  use global_vars, only : gradw_x
  use global_vars, only : gradw_y
  use global_vars, only : gradw_z
  use global_vars, only : gradT_x
  use global_vars, only : gradT_y
  use global_vars, only : gradT_z
  use global_vars, only : gradtk_x
  use global_vars, only : gradtk_y
  use global_vars, only : gradtk_z
  use global_vars, only : gradtw_x
  use global_vars, only : gradtw_y
  use global_vars, only : gradtw_z

  use copy_bc       , only : copy1

  use utils       , only : dmsg
  use utils       , only :   alloc
  use utils       , only : dealloc

#include "error.inc"

  implicit none
  private

  public :: setup_viscosity
  public :: destroy_viscosity
  public :: calculate_viscosity

  contains

    subroutine calculate_viscosity()
      !< Calculate molecular and turbulent viscosity
      implicit none
      integer :: i,j,k
      real :: T ! molecular viscosity
      real :: c ! kkl eddy viscosity
      !- sst varibales -!
      real    :: F
      real    :: arg2
      real    :: vort
      real    :: NUM
      real    :: DENOM
      ! for arg2
      real :: var1
      real :: var2
      !for vorticity
      real :: wijwij
      real :: wx
      real :: wy
      real :: wz
      !for strain calculation
      real :: SijSij
      real :: Sxx, Syy, Szz
      real :: Sxy, Szx, Syz
      real :: strain
      !for arg1
      real :: arg1
      real :: CD
      real :: right
      real :: left

      ! sa variables
      real :: fv1
      real :: xi

      !--- calculate_molecular_viscosity ---!
      if (mu_ref/=0.) then
        select case (trim(mu_variation))
          case ('sutherland_law')
            ! apply_sutherland_law
            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx
                  T = pressure(i,j,k)/(density(i,j,k)*R_gas)
                  mu(i,j,k) = mu_ref * ((T/T_ref)**(1.5)) &
                            *((T_ref + Sutherland_temp)&
                            /(T + Sutherland_temp))
                end do
              end do
            end do

          case ('constant')
            !do nothing
            !mu will be equal to mu_ref
            continue

          case DEFAULT
            print*,"mu_variation not recognized:"
            print*, "   found '",trim(mu_variation),"'"
            print*, "accepted values: 1) sutherland_law"
            print*, "                 2) constant"
            Fatal_error
        end select
      end if
      !--- end molecular viscosity calculation---!

      !--- calculate_turbulent_viscosity  ---!
      if (turbulence/='none') then
        select case (trim(turbulence))

          case ('none')
            !do nothing
            continue

          case ('sa', 'saBC')
            !call calculate_sa_mu()
            do k = 0,kmx
              do j = 0,jmx
                do i = 0,imx
                  ! xsi 
                   xi = tv(i,j,k)*density(i,j,k)/mu(i,j,k)
                  !calculation fo fv1 function
                  fv1 = (xi**3)/((xi**3) + (cv1**3))
                  sa_mu(i,j,k) = density(i,j,k)*tv(i,j,k)*fv1
                end do
              end do
            end do

            ! populating ghost cell
            do i = 1,6
              select case(id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-1,-2,-3,-4,-6,-7,-8,-9)
                  !call copy1(sa_mu, "symm", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        sa_mu(      0, 1:jmx-1, 1:kmx-1) = sa_mu(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        sa_mu(  imx  , 1:jmx-1, 1:kmx-1) = sa_mu( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        sa_mu(1:imx-1,       0, 1:kmx-1) = sa_mu(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        sa_mu(1:imx-1,   jmx  , 1:kmx-1) = sa_mu(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        sa_mu(1:imx-1, 1:jmx-1,       0) = sa_mu(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        sa_mu(1:imx-1, 1:jmx-1,   kmx  ) = sa_mu(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select

                case(-5)
                  !call copy1(sa_mu, "anti", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        sa_mu(      0, 1:jmx-1, 1:kmx-1) = -sa_mu(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        sa_mu(  imx  , 1:jmx-1, 1:kmx-1) = -sa_mu( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        sa_mu(1:imx-1,       0, 1:kmx-1) = -sa_mu(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        sa_mu(1:imx-1,   jmx  , 1:kmx-1) = -sa_mu(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        sa_mu(1:imx-1, 1:jmx-1,       0) = -sa_mu(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        sa_mu(1:imx-1, 1:jmx-1,   kmx  ) = -sa_mu(1:imx-1, 1:jmx-1,  kmx-1)
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

                  ! calculate_arg2()
                  var1 = sqrt(tk(i,j,k))/(bstar*tw(i,j,k)*dist(i,j,k))
                  var2 = 500*(mu(i,j,k)/density(i,j,k))/((dist(i,j,k)**2)*tw(i,j,k))
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

                  NUM = density(i,j,k)*a1*tk(i,j,k)
                  DENOM = max(max((a1*tw(i,j,k)), strain*F),1.0e-10)
                  sst_mu(i,j,k) = NUM/DENOM
                  !-- end eddy visocisyt calculation --!
                  !-- calculating blending function F1 --!
                  CD = max(2*density(i,j,k)*sigma_w2*(                             & 
                                                      gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                                    + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                                    + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                                     )/tw(i,j,k),                  &
                           1.0e-10)

                  right = 4*(density(i,j,k)*sigma_w2*tk(i,j,k))/(CD*(dist(i,j,k)**2))
                  left = max(var1, var2)
                  arg1 = min(left, right)
                  sst_F1(i,j,k) = tanh(arg1**4)
                  !-- end of blending function F1 calculation --!
                end do
              end do
            end do

            select case(trim(transition))
              case('lctm2015')
                do k = 0,kmx
                  do j = 0,jmx
                    do i = 0,imx
                      !modified blending function (Menter 2015)
                      var1 = density(i,j,k)*dist(i,j,k)*sqrt(tk(i,j,k))/mu(i,j,k)
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
              select case(id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-1,-2,-3,-4,-6,-7,-8,-9)
                  !call copy1(sst_mu, "symm", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        sst_mu(      0, 1:jmx-1, 1:kmx-1) = sst_mu(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) = sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        sst_mu(  imx  , 1:jmx-1, 1:kmx-1) = sst_mu( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) = sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        sst_mu(1:imx-1,       0, 1:kmx-1) = sst_mu(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) = sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        sst_mu(1:imx-1,   jmx  , 1:kmx-1) = sst_mu(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) = sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        sst_mu(1:imx-1, 1:jmx-1,       0) = sst_mu(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) = sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        sst_mu(1:imx-1, 1:jmx-1,   kmx  ) = sst_mu(1:imx-1, 1:jmx-1,  kmx-1)
                        sst_F1(1:imx-1, 1:jmx-1,   kmx  ) = sst_F1(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
                case(-5)
                  !call copy1(sst_mu, "anti", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        sst_mu(      0, 1:jmx-1, 1:kmx-1) = -sst_mu(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) =  sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        sst_mu(  imx  , 1:jmx-1, 1:kmx-1) = -sst_mu( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) =  sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        sst_mu(1:imx-1,       0, 1:kmx-1) = -sst_mu(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) =  sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        sst_mu(1:imx-1,   jmx  , 1:kmx-1) = -sst_mu(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) =  sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        sst_mu(1:imx-1, 1:jmx-1,       0) = -sst_mu(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) =  sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        sst_mu(1:imx-1, 1:jmx-1,   kmx  ) = -sst_mu(1:imx-1, 1:jmx-1,  kmx-1)
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

                  ! calculate_arg2()
                  var1 = sqrt(tk(i,j,k))/(bstar*tw(i,j,k)*dist(i,j,k))
                  var2 = 500*(mu(i,j,k)/density(i,j,k))/((dist(i,j,k)**2)*tw(i,j,k))
                  arg2 = max(2*var1, var2)

                  ! calculate_f2()
                  F = tanh(arg2**2)

                  ! calculate_vorticity(
                  wx = (gradw_y(i,j,k) - gradv_z(i,j,k))
                  wy = (gradu_z(i,j,k) - gradw_x(i,j,k))
                  wz = (gradv_x(i,j,k) - gradu_y(i,j,k))

                  wijwij = wx**2 + wy**2 + wz**2

                  vort = sqrt(wijwij)

                  NUM = density(i,j,k)*a1*tk(i,j,k)
                  DENOM = max(max((a1*tw(i,j,k)), vort*F),1.e-20)
                  sst_mu(i,j,k) = NUM/DENOM
                  !-- end eddy visocisyt calculation --!
                  !-- calculating blending function F1 --!
                  CD = max(2*density(i,j,k)*sigma_w2*(                             & 
                                                      gradtk_x(i,j,k)*gradtw_x(i,j,k)&
                                                    + gradtk_y(i,j,k)*gradtw_y(i,j,k)&
                                                    + gradtk_z(i,j,k)*gradtw_z(i,j,k)&
                                                     )/tw(i,j,k),                  &
                           1e-20)

                  right = 4*(density(i,j,k)*sigma_w2*tk(i,j,k))/(CD*(dist(i,j,k)**2))
                  left = max(var1, var2)
                  arg1 = min(left, right)
                  sst_F1(i,j,k) = tanh(arg1**4)
                  !-- end of blending function F1 calculation --!
                end do
              end do
            end do

            select case(trim(transition))
              case('lctm2015')
                do k = 0,kmx
                  do j = 0,jmx
                    do i = 0,imx
                      !modified blending function (Menter 2015)
                      var1 = density(i,j,k)*dist(i,j,k)*sqrt(tk(i,j,k))/mu(i,j,k)
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
              select case(id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-1,-2,-3,-4,-6,-7,-8,-9)
                  !call copy1(sst_mu, "symm", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        sst_mu(      0, 1:jmx-1, 1:kmx-1) = sst_mu(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) = sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        sst_mu(  imx  , 1:jmx-1, 1:kmx-1) = sst_mu( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) = sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        sst_mu(1:imx-1,       0, 1:kmx-1) = sst_mu(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) = sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        sst_mu(1:imx-1,   jmx  , 1:kmx-1) = sst_mu(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) = sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        sst_mu(1:imx-1, 1:jmx-1,       0) = sst_mu(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) = sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        sst_mu(1:imx-1, 1:jmx-1,   kmx  ) = sst_mu(1:imx-1, 1:jmx-1,  kmx-1)
                        sst_F1(1:imx-1, 1:jmx-1,   kmx  ) = sst_F1(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
                case(-5)
                  !call copy1(sst_mu, "anti", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        sst_mu(      0, 1:jmx-1, 1:kmx-1) = -sst_mu(     1, 1:jmx-1, 1:kmx-1)
                        sst_F1(      0, 1:jmx-1, 1:kmx-1) =  sst_F1(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        sst_mu(  imx  , 1:jmx-1, 1:kmx-1) = -sst_mu( imx-1, 1:jmx-1, 1:kmx-1)
                        sst_F1(  imx  , 1:jmx-1, 1:kmx-1) =  sst_F1( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        sst_mu(1:imx-1,       0, 1:kmx-1) = -sst_mu(1:imx-1,      1, 1:kmx-1)
                        sst_F1(1:imx-1,       0, 1:kmx-1) =  sst_F1(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        sst_mu(1:imx-1,   jmx  , 1:kmx-1) = -sst_mu(1:imx-1,  jmx-1, 1:kmx-1)
                        sst_F1(1:imx-1,   jmx  , 1:kmx-1) =  sst_F1(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        sst_mu(1:imx-1, 1:jmx-1,       0) = -sst_mu(1:imx-1, 1:jmx-1,      1)
                        sst_F1(1:imx-1, 1:jmx-1,       0) =  sst_F1(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        sst_mu(1:imx-1, 1:jmx-1,   kmx  ) = -sst_mu(1:imx-1, 1:jmx-1,  kmx-1)
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
                  kkl_mu(i,j,k) = c*density(i,j,k)*tkl(i,j,k)&
                    /(max(sqrt(tk(i,j,k)),1.e-20))
                  if(tkl(i,j,k)<1.e-14 .or. tk(i,j,k)<1.e-14) &
                    kkl_mu(i,j,k) =0.0 
                end do
              end do
            end do

            ! populating ghost cell
            do i = 1,6
              select case(id(i))
                case(-10,0:)
                  !interface
                  continue

                case(-4:-1,-6,-8,-9)
                  !call copy1(kkl_mu, "symm", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        kkl_mu(      0, 1:jmx-1, 1:kmx-1) = kkl_mu(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        kkl_mu(  imx  , 1:jmx-1, 1:kmx-1) = kkl_mu( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        kkl_mu(1:imx-1,       0, 1:kmx-1) = kkl_mu(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        kkl_mu(1:imx-1,   jmx  , 1:kmx-1) = kkl_mu(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        kkl_mu(1:imx-1, 1:jmx-1,       0) = kkl_mu(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        kkl_mu(1:imx-1, 1:jmx-1,   kmx  ) = kkl_mu(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
                case(-5)
                  !call copy1(kkl_mu, "anti", face_names(i))
                  select case(face_names(i))
                    case("imin")
                        kkl_mu(      0, 1:jmx-1, 1:kmx-1) = -kkl_mu(     1, 1:jmx-1, 1:kmx-1)
                    case("imax")
                        kkl_mu(  imx  , 1:jmx-1, 1:kmx-1) = -kkl_mu( imx-1, 1:jmx-1, 1:kmx-1)
                    case("jmin")
                        kkl_mu(1:imx-1,       0, 1:kmx-1) = -kkl_mu(1:imx-1,      1, 1:kmx-1)
                    case("jmax")
                        kkl_mu(1:imx-1,   jmx  , 1:kmx-1) = -kkl_mu(1:imx-1,  jmx-1, 1:kmx-1)
                    case("kmin")
                        kkl_mu(1:imx-1, 1:jmx-1,       0) = -kkl_mu(1:imx-1, 1:jmx-1,      1)
                    case("kmax")
                        kkl_mu(1:imx-1, 1:jmx-1,   kmx  ) = -kkl_mu(1:imx-1, 1:jmx-1,  kmx-1)
                    case DEFAULT
                      print*, "ERROR: wrong face for boundary condition"
                      Fatal_error
                  end select
              end select
            end do
            !--- end of kkl eddy viscosity calculation ---!

          case DEFAULT 
            !call turbulence_read_error()
            Fatal_error

        end select
      end if
      !--- end turbulent viscosity calculation---!
      !--- check on viscosity ---!
      if(any(isnan(mu))) then
        Fatal_error
      end if

    end subroutine calculate_viscosity

    subroutine setup_viscosity()
      !< Allocate and pointer for molecular and turbulent viscosity

      !setup_molecular_viscosity()
      if (mu_ref/=0.) then
        call alloc(mu, -2, imx+2, -2, jmx+2, -2, kmx+2)
        mu = mu_ref !intialize
      end if

      !--- setup_turbulent_viscosity ---!
      if (turbulence/='none') then
        call alloc(mu_t, -2,imx+2, -2,jmx+2, -2,kmx+2)

        select case (trim(turbulence))

          case ('none')
            !do nothing
            continue

          case ('sa', 'saBC')
            sa_mu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu_t(:,:,:)

          case ('sst', 'sst2003')
            sst_mu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu_t(:,:,:)
            !-- sst blending funciton F1 --!
            call alloc(sst_F1, -2,imx+2, -2,jmx+2, -2,kmx+2)
            sst_F1=0.
            !-- sst blnding function setup compete--!

          case ('kkl')
            kkl_mu(-2:imx+2,-2:jmx+2,-2:kmx+2) => mu_t(:,:,:)

          case DEFAULT 
            !call turbulence_read_error()
            Fatal_error

        end select
      end if
      ! --- end turbulent viscosity setup ---!

    end subroutine setup_viscosity

    subroutine destroy_viscosity()
      !< Deallocate and nullify viscosity (turbulent/molecular)

      ! destroy_molecular_viscosity ---!
      if (mu_ref/=0.) then
        call dealloc(mu)
      end if

      !--- destroy_turbulent_viscosity ---!
      if (turbulence/='none') then
        select case (trim(turbulence))

          case ('none')
            !do nothing
            continue

          case('sa', 'saBC')
            nullify(sa_mu)

          case ('sst', 'sst2003')
            nullify(sst_mu)
            !-- blending funciton F1 --!
            call dealloc(sst_F1)
            !--- sst blending funciton destoryed--!

          case ('kkl')
            nullify(kkl_mu)

          case DEFAULT 
            !call turbulence_read_error()
            Fatal_error

        end select
        call dealloc(mu_t)
      end if
      !--- end of turublent viscosity destruction ---!

    end subroutine destroy_viscosity

end module viscosity
