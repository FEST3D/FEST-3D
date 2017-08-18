module kkl_viscosity

  !------------------------------------------------------------------
  ! author    - Jatinder Pal Singh Sandhu
  ! objective - calcualte sst viscosity
  ! source    - https://turbmodels.larc.nasa.gov/sst.html#sst-2003
  ! BUg       - double calculation of F1 and vorticity variable 
  !             (in sst_source module)
  !-------------------------------------------------------------------

  use global_kkl   , only : cmu
  use global_vars  , only : imx
  use global_vars  , only : jmx
  use global_vars  , only : kmx
  use global_vars  , only : mu
  use global_vars  , only : kkl_mu
  use global_vars  , only : density
  use global_vars  , only : tk
  use global_vars  , only : tkl
  use global_vars   , only : id
  use global_vars   , only : face_names
  use copy_bc       , only : copy1

  private

      
  public :: calculate_kkl_mu

  contains

    subroutine calculate_kkl_mu()
      implicit none
      integer :: i,j,k
      real :: c
      c = cmu**0.25

      do k = 1,kmx-1
        do j = 1,jmx-1
          do i = 1,imx-1
            kkl_mu(i,j,k) = c*density(i,j,k)*tkl(i,j,k)/(max(sqrt(tk(i,j,k)),1.e-20))
            if(tkl(i,j,k)<1.e-14 .or. tk(i,j,k)<1.e-14) kkl_mu(i,j,k) =0.0 
          end do
        end do
      end do

      ! populating ghost cell
      do i = 1,6
        select case(id(i))
          case(-4:-1,-6,-8)
            call copy1(kkl_mu, "symm", face_names(i))
          case(-5)
            call copy1(kkl_mu, "anti", face_names(i))
        end select
      end do


    end subroutine calculate_kkl_mu

end module kkl_viscosity
