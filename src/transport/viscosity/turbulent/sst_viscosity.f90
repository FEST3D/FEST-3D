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

  private

  integer :: i,j,k
  real    :: F
  real    :: arg2
  real    :: vort
  real    :: NUM
  real    :: DENOM
      
  public :: calculate_sst_mu

  contains

    subroutine calculate_sst_mu()
      implicit none

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            call calculate_arg2()
            call calculate_f2()
            call calculate_vorticity()
            NUM = density(i,j,k)*a1*tk(i,j,k)
            DENOM = max((a1*tw(i,j,k)), vort*F)
            sst_mu(i,j,k) = NUM/DENOM
          end do
        end do
      end do

    end subroutine calculate_sst_mu


    subroutine calculate_f2()
      implicit none

      F = tanh(arg2**2)

    end subroutine calculate_f2


    subroutine calculate_arg2()
      implicit none
      real :: var1
      real :: var2

      var1 = 2*sqrt(tk(i,j,k))/(bstar*tw(i,j,k)*dist(i,j,k))
      var2 = 500*mu(i,j,k)*(dist(i,j,k)**2)*tw(i,j,k)
      arg2 = max(var1, var2)

    end subroutine calculate_arg2


    subroutine calculate_vorticity()
      implicit none
      real :: wijwij
      real :: wx
      real :: wy
      real :: wz

      wx = 0.5*(gradw_y(i,j,k) - gradv_z(i,j,k))
      wy = 0.5*(gradu_z(i,j,k) - gradw_x(i,j,k))
      wz = 0.5*(gradv_x(i,j,k) - gradu_y(i,j,k))

      wijwij = wx**2 + wy**2 + wz**2

      vort = sqrt(2*wijwij)

    end subroutine calculate_vorticity


end module sst_viscosity
