module sst_viscosity

  !------------------------------------------------------------------
  ! author    - Jatinder Pal Singh Sandhu
  ! objective - calcualte sst viscosity
  ! source    - https://turbmodels.larc.nasa.gov/sst.html#sst-2003
  ! BUg       - double calculation of F1 and vorticity variable 
  !             (in sst_source module)
  !-------------------------------------------------------------------

  use global_vars  , only : imx
  use global_vars  , only : jmx
  use global_vars  , only : kmx
  use global_vars  , only : mu
  use global_vars  , only : sst_mu
  use global_vars  , only : density
  use global_vars  , only : tk
  use global_vars  , only : tw
  use global_vars  , only : dist
  use global_sst   , only : bstar
  use global_sst   , only : a1


  use global_vars   , only : gradu_y
  use global_vars   , only : gradu_z
  use global_vars   , only : gradv_x
  use global_vars   , only : gradv_z
  use global_vars   , only : gradw_x
  use global_vars   , only : gradw_y
  use global_vars   , only : id
  use global_vars   , only : face_names
  use copy_bc       , only : copy1

  private

      
  public :: calculate_sst_mu

  contains

    subroutine calculate_sst_mu()
      implicit none
      integer :: i,j,k
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

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1

            ! calculate_arg2()
            var1 = 2*sqrt(tk(i,j,k))/(bstar*tw(i,j,k)*dist(i,j,k))
            var2 = 500*(mu(i,j,k)/density(i,j,k))/((dist(i,j,k)**2)*tw(i,j,k))
            arg2 = max(var1, var2)

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
          end do
        end do
      end do

      ! populating ghost cell
      do i = 1,6
        select case(id(i))
          case(-4:-1,-6,-8)
            call copy1(sst_mu, "symm", face_names(i))
          case(-5)
            call copy1(sst_mu, "anti", face_names(i))
        end select
      end do


    end subroutine calculate_sst_mu

end module sst_viscosity
